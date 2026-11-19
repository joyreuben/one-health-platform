"""One Health Intelligence Platform: analysis engine.

Sequence validation, NCBI BLAST, AMR classification, epidemic/credibility
analysis (Claude), confidence scoring, HL7 generation, geocoding cache and
service-status probes. No page layout lives here; see app.py / ui_theme.py.
"""
import os
import subprocess
import time
import threading
import atexit
import traceback
import json
import hashlib
import random
from collections import deque
from typing import List, Dict, Any, Tuple, Optional
from datetime import datetime
from io import StringIO, BytesIO
import concurrent.futures
import re
import html
import shutil
import tempfile

import requests
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

import feeds
from storage import Store

try:
    from anthropic import (Anthropic, APIConnectionError, APIStatusError,
                           AuthenticationError, NotFoundError, PermissionDeniedError)
    _anthropic_available = True
except ImportError:
    _anthropic_available = False

__version__ = "4.0.0"

# =====================================================================
# ANTHROPIC CLIENT
# =====================================================================
ANTHROPIC_MODEL = "claude-sonnet-4-6"  # single source of truth for the model id

def _get_api_key() -> str:
    # st.secrets raises when no secrets.toml exists; that must not prevent
    # falling back to the ANTHROPIC_API_KEY environment variable.
    try:
        k = st.secrets.get("ANTHROPIC_API_KEY", "")
    except Exception:
        k = ""
    return k or os.getenv("ANTHROPIC_API_KEY", "")


@st.cache_resource
def _make_anthropic_client():
    """Create the client once per process instead of on every rerun.
    Returns (client, error_message)."""
    if not _anthropic_available:
        return None, "anthropic package not installed (add to requirements.txt)"
    api_key = _get_api_key()
    if not api_key:
        return None, "ANTHROPIC_API_KEY not set (secrets or env)"
    try:
        return Anthropic(api_key=api_key), ""
    except Exception as exc:
        return None, str(exc)


anthropic_client, _claude_error = _make_anthropic_client()


# =====================================================================
# 0a. LOCAL DEBUG LOGGING (FIX: portable path, size-capped, best-effort)
# =====================================================================
# BUG FIXED: original code hardcoded a Windows-only, user-specific path
# ("C:\Users\JOY\...") which silently failed on any other OS/deployment
# (failures were swallowed by a bare except:pass) and grew unbounded,
# writing raw article/sequence text to disk indefinitely.
LOG_DIR = os.path.join(os.getcwd(), "logs")
CLAUDE_LOG_PATH = os.path.join(LOG_DIR, "claude_raw.log")
LOG_MAX_BYTES = 2_000_000  # cap ~2MB; rotate by dropping the older half


def _append_claude_log(tag: str, attempt: int, raw: str):
    """Best-effort local debug log for Claude API responses.
    Never raises. Uses a portable, app-relative path and caps file size
    so it can't grow forever or leak indefinitely-retained user content."""
    try:
        os.makedirs(LOG_DIR, exist_ok=True)
        entry = f"[{datetime.now().isoformat()}] {tag} attempt {attempt}\n{raw or '<EMPTY>'}\n---\n"
        if os.path.exists(CLAUDE_LOG_PATH) and os.path.getsize(CLAUDE_LOG_PATH) > LOG_MAX_BYTES:
            with open(CLAUDE_LOG_PATH, 'r', encoding='utf-8', errors='replace') as f:
                data = f.read()
            data = data[len(data) // 2:]
            with open(CLAUDE_LOG_PATH, 'w', encoding='utf-8') as f:
                f.write(data)
        with open(CLAUDE_LOG_PATH, 'a', encoding='utf-8') as f:
            f.write(entry)
    except Exception:
        pass


def _claude_json(prompt: str, tag: str, attempts: int = 3, max_tokens: int = 1024) -> Dict[str, Any]:
    """Send a prompt to Claude and return the first JSON object in its reply.
    Retries (with backoff) on empty / non-JSON replies and API errors;
    raises the last error if every attempt fails. The SDK itself already
    retries transient HTTP errors (429/5xx)."""
    last_exc: Exception = ValueError("no attempts made")
    backoff = 1
    for attempt in range(1, attempts + 1):
        try:
            response = anthropic_client.messages.create(
                model=ANTHROPIC_MODEL,
                max_tokens=max_tokens,
                messages=[{"role": "user", "content": prompt}],
            )
            raw = "".join(getattr(b, "text", "") for b in response.content
                          if getattr(b, "type", "") == "text").strip()
            _append_claude_log(tag, attempt, raw)
            if not raw:
                raise ValueError("Empty Claude response content")
            start = raw.find("{")
            if start == -1:
                raise ValueError("No JSON in response")
            # raw_decode handles braces inside JSON strings correctly,
            # unlike naive brace counting.
            obj, _ = json.JSONDecoder().raw_decode(raw[start:])
            if not isinstance(obj, dict):
                raise ValueError("JSON response is not an object")
            return obj
        except Exception as exc:
            last_exc = exc
            if attempt < attempts:
                time.sleep(backoff)
                backoff *= 2
    raise last_exc


def _to_int(v, default=0) -> int:
    try:
        return int(float(str(v).replace(",", "").strip()))
    except Exception:
        m = re.search(r"\d+", str(v))
        return int(m.group()) if m else default


def _to_unit_float(v, default=0.0) -> float:
    try:
        return float(np.clip(float(v), 0.0, 1.0))
    except Exception:
        return default


# =====================================================================
# 0. PERSISTENT STORE (analyses, audit trail, feed items)
# The audit trail used to live in st.session_state and vanished on every
# refresh. It is now a SQLite database shared by all sessions (storage.py).
# =====================================================================
@st.cache_resource
def get_store() -> Store:
    return Store()


def init_session_state():
    """Per-browser-session state. Must be called by the app on every run:
    this module is imported once per process, so module-level code here
    would only initialise the *first* session."""
    st.session_state.setdefault("result_cache", {})
    st.session_state.setdefault("last_request_time", 0)

# ---------------------------------------------------------------------
# Geocoding cache
#
# NOTE ON SCOPE: the location cache, its locks and the geocode executor
# are process-wide (shared by every session on this server), because a
# geocoded lat/lon for "Lagos" is the same fact for every user. They live
# in a st.cache_resource object (see _geo_state) so they survive
# Streamlit's per-rerun re-execution of this script.
# =====================================================================
LOCATION_CACHE_PATH = os.path.join(os.getcwd(), "location_cache.json")


def _migrate_cache_entry(loc, raw_entry):
    """Normalise old/new cache formats into the current schema:
    {'value': {'lon':..,'lat':..} | None, 'attempts': int, 'last_try': float, 'status': str}
    """
    if raw_entry is None:
        return {'value': None, 'attempts': 1, 'last_try': 0.0, 'status': 'failed'}
    if isinstance(raw_entry, dict) and 'value' in raw_entry:
        return raw_entry  # already current schema
    if isinstance(raw_entry, dict) and 'lon' in raw_entry and 'lat' in raw_entry:
        # old schema: bare {'lon':..,'lat':..}
        return {'value': raw_entry, 'attempts': 0, 'last_try': time.time(), 'status': 'ok'}
    return {'value': None, 'attempts': 1, 'last_try': 0.0, 'status': 'failed'}


def _load_location_cache() -> Dict[str, Dict[str, Any]]:
    try:
        if os.path.exists(LOCATION_CACHE_PATH):
            with open(LOCATION_CACHE_PATH, 'r', encoding='utf-8') as f:
                raw = json.load(f)
            return {loc: _migrate_cache_entry(loc, v) for loc, v in raw.items()}
    except Exception:
        pass
    return {}


@st.cache_resource
def _geo_state() -> Dict[str, Any]:
    """Process-wide geocoding state.

    Streamlit re-executes this script on every rerun, so plain module
    globals were rebuilt each time: a new executor (and atexit hook) per
    rerun, the cache reloaded from disk, and background threads writing
    into a stale dict that then overwrote newer entries on disk.
    st.cache_resource keeps one shared instance for the whole process.

    The cache lock is an RLock because _save_location_cache() takes it
    and is called by code that already holds it (a plain Lock deadlocked
    the first time a location was geocoded or the cache was cleared)."""
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)
    atexit.register(lambda: executor.shutdown(wait=False))
    return {
        "cache": _load_location_cache(),
        "cache_lock": threading.RLock(),
        "rate_lock": threading.Lock(),
        "last_time": [0.0],
        "executor": executor,
    }


_GEO = _geo_state()
LOCATION_CACHE: Dict[str, Dict[str, Any]] = _GEO["cache"]
LOCATION_CACHE_LOCK = _GEO["cache_lock"]
GEOCODE_LOCK = _GEO["rate_lock"]


def _save_location_cache():
    try:
        with LOCATION_CACHE_LOCK:
            with open(LOCATION_CACHE_PATH, 'w', encoding='utf-8') as f:
                json.dump(LOCATION_CACHE, f, ensure_ascii=False, indent=2)
    except Exception:
        pass


def _geocode_nominatim(location: str):
    """Simple Nominatim lookup returning {'lon': ..., 'lat': ...} or None.
    Respects a 1s-ish rate limit with a module-level lock for thread safety."""
    if not location or not isinstance(location, str):
        return None
    last_time = _GEO["last_time"]
    with GEOCODE_LOCK:
        elapsed = time.time() - last_time[0]
        if elapsed < 1.1:
            time.sleep(1.1 - elapsed)
        last_time[0] = time.time()

    params = {'q': location, 'format': 'json', 'limit': 1, 'addressdetails': 0}
    headers = {'User-Agent': 'OneHealthPlatform/1.0 (contact: no-reply@example.com)'}
    try:
        r = requests.get('https://nominatim.openstreetmap.org/search', params=params, headers=headers, timeout=10)
        r.raise_for_status()
        data = r.json()
        if data:
            return {'lon': float(data[0]['lon']), 'lat': float(data[0]['lat'])}
    except Exception:
        return None
    return None


def _geocode_and_store(location: str):
    res = _geocode_nominatim(location)
    with LOCATION_CACHE_LOCK:
        prev = LOCATION_CACHE.get(location, {})
        prev_attempts = prev.get('attempts', 0) if isinstance(prev, dict) else 0
        LOCATION_CACHE[location] = {
            'value': res,
            'attempts': prev_attempts if res else prev_attempts + 1,
            'last_try': time.time(),
            'status': 'ok' if res else 'failed',
        }
        _save_location_cache()


def maybe_geocode_in_background(location: str):
    """Start background geocoding for a location if not already cached.

    BUG FIXED: previously a failed lookup was cached as a bare `None` and
    `location in LOCATION_CACHE` returning True meant it would *never* be
    retried again, even after a purely transient network error. Now we
    track attempts/last_try and allow retries with exponential backoff
    (capped) instead of a silent permanent failure.
    """
    if not location:
        return
    if not st.session_state.get('enable_geocoding', True):
        return

    now = time.time()
    with LOCATION_CACHE_LOCK:
        entry = LOCATION_CACHE.get(location)
        if entry is not None:
            if entry.get('value') is not None:
                return  # already resolved successfully
            attempts = entry.get('attempts', 0)
            if attempts >= 5:
                return  # give up after 5 failed attempts
            backoff = min(3600, 30 * (2 ** attempts))  # 30s, 60s, 120s ... capped at 1h
            if now - entry.get('last_try', 0) < backoff:
                return  # too soon to retry
        LOCATION_CACHE[location] = {
            'value': None,
            'attempts': entry.get('attempts', 0) if entry else 0,
            'last_try': now,
            'status': 'pending',
        }
        _save_location_cache()

    try:
        _GEO["executor"].submit(_geocode_and_store, location)
    except Exception:
        try:
            _geocode_and_store(location)
        except Exception:
            pass


# =====================================================================
# 1. SEQUENCE VALIDATION
# =====================================================================
class SequenceValidator:
    VALID_DNA_CHARS     = set('ATGCNRYSWKMBDHV')
    VALID_PROTEIN_CHARS = set('ACDEFGHIKLMNPQRSTVWY')
    PROTEIN_ONLY_CHARS  = VALID_PROTEIN_CHARS - VALID_DNA_CHARS  # {'E','F','I','L','P','Q'}

    MIN_LENGTH = 20
    MAX_LENGTH = 500_000

    @staticmethod
    def validate_and_clean(seq):
        if not seq or not seq.strip():
            return False, "❌ Sequence is empty", ""
        lines     = seq.strip().splitlines()
        seq_lines = lines[1:] if lines[0].startswith('>') else lines
        # Strip ALL whitespace (incl. Windows '\r' and tabs, which were
        # previously reported as "invalid characters") and a trailing
        # protein stop symbol '*'.
        cleaned   = re.sub(r'\s+', '', ''.join(seq_lines)).upper().rstrip('*')
        if not cleaned:
            return False, "❌ Empty after header removal", ""
        if len(cleaned) < SequenceValidator.MIN_LENGTH:
            return False, f"❌ Too short (min {SequenceValidator.MIN_LENGTH} bp)", ""
        if len(cleaned) > SequenceValidator.MAX_LENGTH:
            return False, f"❌ Too long (max {SequenceValidator.MAX_LENGTH:,} bp)", ""
        bad = set(cleaned) - (SequenceValidator.VALID_DNA_CHARS | SequenceValidator.VALID_PROTEIN_CHARS)
        if bad:
            return False, f"❌ Invalid characters: {', '.join(sorted(bad))}", ""
        return True, "", cleaned

    @staticmethod
    def is_protein(seq):
        return bool(set(seq.upper()) & SequenceValidator.PROTEIN_ONLY_CHARS)

    @staticmethod
    def is_likely_coding(seq):
        dna = seq.upper()
        if len(dna) < 90:
            return False
        if dna[:3] in ("ATG", "GTG", "TTG") and dna[-3:] in ("TAA", "TAG", "TGA"):
            return True
        return False


# =====================================================================
# 2. RESFINDER WEB API
# =====================================================================
class ResFinderWebClient:
    API_URL = "https://cge.food.dtu.dk/services/ResFinder/api/v1/run"
    TIMEOUT = 120
    API_DISABLED = True
    DISABLED_REASON = "Disabled in this build (config flag), not a connectivity failure."

    @classmethod
    def is_available(cls):
        return False

    @classmethod
    def analyze(cls, sequence, sample_id):
        return None


# =====================================================================
# 2b. AMRFINDERPLUS (NCBI's curated AMR gene caller, run locally)
# Replaces "regex over BLAST hit titles" as the primary classifier: it
# searches NCBI's Reference Gene Catalog with curated per-gene cutoffs
# and returns in seconds. Web BLAST remains as a fallback.
#
# Where it runs:
#   - Linux/macOS (or a container): `amrfinder` on PATH is used directly.
#   - Windows: through WSL, from the conda env named by OHP_AMRFINDER_ENV
#     (default ~/miniforge3/envs/ohp-amrfinder), optional OHP_WSL_DISTRO.
# =====================================================================
class AMRFinderClient:
    WSL_DISTRO = os.getenv("OHP_WSL_DISTRO", "")
    WSL_ENV = os.getenv("OHP_AMRFINDER_ENV", "$HOME/miniforge3/envs/ohp-amrfinder")
    TIMEOUT = 300
    THREADS = 4

    # method (without the trailing X/P/N search-type letter) -> evidence weight
    METHOD_WEIGHTS = {
        "ALLELE": 1.00, "EXACT": 1.00, "POINT": 0.90, "BLAST": 0.85,
        "PARTIAL_CONTIG_END": 0.65, "PARTIAL": 0.60, "HMM": 0.50, "INTERNAL_STOP": 0.35,
    }
    CRITICAL_SUBCLASSES = ("CARBAPENEM", "COLISTIN", "VANCOMYCIN", "METHICILLIN", "CEPHALOSPORIN")
    HIGH_CLASSES = ("BETA-LACTAM", "TETRACYCLINE", "GLYCOPEPTIDE")

    @classmethod
    def _backend(cls) -> Optional[str]:
        if shutil.which("amrfinder"):
            return "native"
        if os.name == "nt" and shutil.which("wsl"):
            return "wsl"
        return None

    @classmethod
    def _exec(cls, flag_args: str, fasta: Optional[str] = None, timeout: Optional[int] = None):
        """Run amrfinder. `flag_args` is built only from our own constants;
        sequence data goes in via a temp file / stdin, never the command line.
        Returns (returncode, stdout, stderr)."""
        backend = cls._backend()
        if backend is None:
            raise RuntimeError("AMRFinderPlus is not installed")
        timeout = timeout or cls.TIMEOUT
        flags = 0x08000000 if os.name == "nt" else 0  # CREATE_NO_WINDOW
        if backend == "native":
            path = None
            try:
                args = ["amrfinder"] + flag_args.split()
                if fasta is not None:
                    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False, encoding="utf-8") as fh:
                        fh.write(fasta)
                        path = fh.name
                    args = [a.replace("@FASTA@", path) for a in args]
                p = subprocess.run(args, capture_output=True, timeout=timeout)
            finally:
                if path and os.path.exists(path):
                    os.unlink(path)
        else:
            # The FASTA is piped on stdin into a temp file inside WSL, which
            # avoids Windows<->Linux path translation and quoting entirely.
            script = (f"export CONDA_PREFIX={cls.WSL_ENV}; "
                      "export PATH=$CONDA_PREFIX/bin:/usr/local/bin:/usr/bin:/bin; ")
            if fasta is not None:
                script += ("T=$(mktemp --suffix=.fa); cat > $T; "
                           f"amrfinder {flag_args.replace('@FASTA@', '$T')}; rc=$?; rm -f $T; exit $rc")
            else:
                script += f"amrfinder {flag_args}"
            cmd = ["wsl.exe"] + (["-d", cls.WSL_DISTRO] if cls.WSL_DISTRO else []) + ["-e", "bash", "-c", script]
            p = subprocess.run(cmd, input=(fasta or "").encode("utf-8"), capture_output=True,
                               timeout=timeout, creationflags=flags)
        return (p.returncode, p.stdout.decode("utf-8", "replace"), p.stderr.decode("utf-8", "replace"))

    @classmethod
    def info(cls) -> Dict[str, Any]:
        return _amrfinder_info()

    @classmethod
    def is_available(cls) -> bool:
        return bool(_amrfinder_info().get("ok"))

    @staticmethod
    def parse_tsv(text: str) -> List[Dict[str, Any]]:
        """Parse AMRFinderPlus TSV output (handles v3 and v4 column names)."""
        lines = [l for l in text.splitlines() if l.strip()]
        if len(lines) < 2:
            return []
        header = [h.strip() for h in lines[0].split("\t")]

        def pick(row, *names, default=""):
            for n in names:
                v = row.get(n, "")
                if v and v != "NA":
                    return v
            return default

        def num(v):
            try:
                return float(v)
            except Exception:
                return 0.0

        out = []
        for line in lines[1:]:
            cells = line.split("\t")
            row = {h: (cells[i].strip() if i < len(cells) else "") for i, h in enumerate(header)}
            out.append({
                "seq_id":    pick(row, "Contig id", "Protein id", "Protein identifier"),
                "symbol":    pick(row, "Element symbol", "Gene symbol", default="unknown"),
                "name":      pick(row, "Element name", "Sequence name"),
                "scope":     pick(row, "Scope"),
                "type":      pick(row, "Type", "Element type"),
                "subtype":   pick(row, "Subtype", "Element subtype"),
                "class":     pick(row, "Class"),
                "subclass":  pick(row, "Subclass"),
                "method":    pick(row, "Method"),
                "coverage":  num(pick(row, "% Coverage of reference", "% Coverage of reference sequence")),
                "identity":  num(pick(row, "% Identity to reference", "% Identity to reference sequence")),
                "accession": pick(row, "Closest reference accession", "Accession of closest sequence"),
                "ref_name":  pick(row, "Closest reference name", "Name of closest sequence"),
                "start":     pick(row, "Start"),
                "stop":      pick(row, "Stop"),
                "strand":    pick(row, "Strand"),
            })
        return out

    @classmethod
    def run(cls, records: Dict[str, str], protein: bool, organism: Optional[str] = None,
            timeout: Optional[int] = None, threads: Optional[int] = None) -> Dict[str, List[Dict[str, Any]]]:
        """Analyse several sequences in one call. `records` maps a
        whitespace-free id to a cleaned sequence. Returns id -> findings.
        `organism` is an AMRFinderPlus organism name (see `amrfinder -l`,
        e.g. "Klebsiella_pneumoniae"); it enables point-mutation screening
        for that organism."""
        if not records:
            return {}
        fasta = "".join(f">{rid}\n{seq}\n" for rid, seq in records.items())
        flag = "-p" if protein else "-n"
        org = ""
        if organism:
            if not re.fullmatch(r"[A-Za-z_]{3,60}", organism):
                raise ValueError(f"Invalid AMRFinderPlus organism name: {organism!r}")
            org = f" -O {organism}"
        rc, out, err = cls._exec(f"{flag} @FASTA@{org} --threads {threads or cls.THREADS}",
                                 fasta=fasta, timeout=timeout)
        if rc != 0:
            tail = " ".join(err.strip().splitlines()[-3:])[:300]
            raise RuntimeError(f"AMRFinderPlus failed (exit {rc}): {tail}")
        grouped: Dict[str, List[Dict[str, Any]]] = {rid: [] for rid in records}
        for f in cls.parse_tsv(out):
            grouped.setdefault(f["seq_id"], []).append(f)
        return grouped

    @classmethod
    def concern(cls, klass: str, subclass: str) -> str:
        k, s = (klass or "").upper(), (subclass or "").upper()
        if any(c in s for c in cls.CRITICAL_SUBCLASSES) or "COLISTIN" in k:
            return "CRITICAL"
        if any(c in k for c in cls.HIGH_CLASSES):
            return "HIGH"
        return "MODERATE"

    @classmethod
    def method_weight(cls, method: str) -> float:
        m = (method or "").upper()
        if m in cls.METHOD_WEIGHTS:          # e.g. INTERNAL_STOP, HMM: no suffix to strip
            return cls.METHOD_WEIGHTS[m]
        return cls.METHOD_WEIGHTS.get(re.sub(r"[XPN]$", "", m), 0.60)

    @classmethod
    def to_gene_info(cls, f: Dict[str, Any]) -> Dict[str, Any]:
        resistance = (f["subclass"] or f["class"] or "Undetermined").replace("/", " / ").title()
        return {
            "marker_found":     f["symbol"],
            "display_name":     f["symbol"],
            "element_name":     f["name"],
            "class":            (f["class"] or "Unknown").title(),
            "mechanism":        f["subtype"] or f["type"] or "Unknown",
            "resistance":       resistance,
            "clinical_concern": cls.concern(f["class"], f["subclass"]),
            "method":           f["method"],
            "source":           "AMRFINDERPLUS",
        }


_AMR_INFO_TTL_S = 600
_AMR_INFO_CACHE: Dict[str, Any] = {"at": 0.0, "value": None}
_AMR_INFO_LOCK = threading.Lock()


def _amrfinder_info(refresh: bool = False) -> Dict[str, Any]:
    """Is AMRFinderPlus runnable, and which software/database version?

    Cached in-process for 10 minutes. This deliberately does not use
    st.cache_data: that cache is inactive outside a running Streamlit app
    (scripts, tests, the feed thread), which made every call start a new
    `amrfinder -V` process, i.e. one WSL launch per sequence in a batch."""
    with _AMR_INFO_LOCK:
        cached = _AMR_INFO_CACHE["value"]
        if cached is not None and not refresh and time.time() - _AMR_INFO_CACHE["at"] < _AMR_INFO_TTL_S:
            return cached
        value = _probe_amrfinder()
        _AMR_INFO_CACHE.update(at=time.time(), value=value)
        return value


def _probe_amrfinder() -> Dict[str, Any]:
    backend = AMRFinderClient._backend()
    if backend is None:
        return {"ok": False, "error": "not installed (no amrfinder on PATH, no WSL)"}
    try:
        rc, out, err = AMRFinderClient._exec("-V", timeout=150)
    except subprocess.TimeoutExpired:
        return {"ok": False, "backend": backend, "error": "timed out starting AMRFinderPlus"}
    except Exception as exc:
        return {"ok": False, "backend": backend, "error": str(exc)[:160]}
    text = out + "\n" + err
    sw = re.search(r"Software version:\s*(\S+)", text)
    db = re.search(r"Database version:\s*(\S+)", text)
    if rc != 0 or not sw or not db:
        if "No valid AMRFinder database" in text:
            reason = "database not downloaded (run: amrfinder -u)"
        elif "not found" in text or "No such file" in text:
            reason = "not installed in the configured environment"
        else:
            reason = (text.strip().splitlines() or ["unknown error"])[-1][:160]
        return {"ok": False, "backend": backend, "error": reason}
    return {"ok": True, "backend": backend, "software": sw.group(1), "database": db.group(1)}


# =====================================================================
# 3. AMR VALIDATION ENGINE
# =====================================================================
class AMRValidationEngine:

    LOCAL_MARKERS = {
        "ndm":   (r"\bndm-?[0-9]?\b|metallo.beta.lactamase",
                  "Metallo-beta-lactamase", "Antibiotic inactivation", "Carbapenems (last-resort)", "CRITICAL"),
        "kpc":   (r"\bkpc-?[0-9]?\b",
                  "Carbapenemase", "Antibiotic inactivation", "Carbapenems (last-resort)", "CRITICAL"),
        # NOTE: titles are passed through _normalise() first, which turns
        # '-' and '_' into spaces ("mcr-1" -> "mcr 1", "CTX-M" -> "ctx m"),
        # so separators must allow whitespace, not just '-'.
        "mcr":   (r"\bmcr[-\s]?[0-9]\b|colistin.resistance|mobile.colistin.resistance",
                  "Colistin Resistance", "Lipid A modification", "Colistin (last-resort)", "CRITICAL"),
        "ctx-m": (r"ctx[-\s]?m\b|esbl|extended.spectrum.beta.lactamase",
                  "Beta-lactamase (ESBL)", "Antibiotic inactivation", "Cephalosporins", "CRITICAL"),
        "mecA":  (r"\bmeca\b|penicillin.binding.protein.2a|pbp2a|pbp.2a",
                  "PBP Alteration", "Target alteration", "Methicillin", "CRITICAL"),
        "vanA":  (r"\bvan[-_ ]?a\b|\bvan[-_ ]?b\b|\bvana\b|\bvanb\b|\btn1546\b|d.ala.{0,3}d.ala.{0,4}ligase|"
                  r"d.alanine.{0,4}d.alanine.{0,4}ligase|"
                  r"d.ala.{0,3}d.lac.ligase|glycopeptide.resistance|m97297",
                  "Glycopeptide Resistance", "Cell wall alteration", "Vancomycin", "CRITICAL"),
        "tetA":  (r"\btet[a-eg-z]?\b"
                  r"|tetracycline.{0,20}resist"
                  r"|tetracycline.{0,10}efflux"
                  r"|tet\s*[abcdefg]\b",
                  "Tetracycline Resistance", "Efflux pump", "Tetracyclines", "HIGH"),
        # FIX (bug #6): "tetB" appeared as its own key in SEQUENCE_MOTIFS
        # but had no matching LOCAL_MARKERS entry, so a tetB motif hit fell
        # through to the generic "Undetermined" default. Note tetA's regex
        # above (\btet[a-eg-z]?\b) already covers "tetb" for *title*
        # matching (Layer 2); this explicit entry plus MARKER_ALIASES below
        # covers the *sequence motif* path (Layer 4) too.
        "tetB":  (r"\btet[b]?\b|tetracycline.{0,20}resist|tetracycline.{0,10}efflux",
                  "Tetracycline Resistance", "Efflux pump", "Tetracyclines", "HIGH"),
        "ermB":  (r"\berm[a-z]?\b|rrna.methyltransferase|23s.rrna.methylase|macrolide.resistance",
                  "Macrolide Resistance", "Target methylation", "Macrolides", "MODERATE"),
        "aacA":  (r"\baac\([0-9]|aminoglycoside.{0,4}acetyltransferase|aminoglycoside.{0,4}transferase",
                  "Aminoglycoside AT", "Inactivation", "Aminoglycosides", "MODERATE"),
        "qnr":   (r"\bqnr[a-z0-9]?\b|quinolone resistance protein",
                  "Quinolone Resistance", "Target protection", "Fluoroquinolones", "MODERATE"),
        "blaTEM":(r"\bbla.?tem\b|\bbla.?shv\b|beta.lactamase",
                  "Beta-lactamase", "Antibiotic inactivation", "Penicillins", "HIGH"),
    }

    LOCAL_MARKER_DISPLAY = {
        "ndm":   "NDM (Metallo-β-lactamase)",
        "kpc":   "KPC (Carbapenemase)",
        "mcr":   "MCR (Mobile Colistin Resistance)",
        "ctx-m": "CTX-M (ESBL)",
        "mecA":  "mecA (PBP2a, methicillin resistance)",
        "vanA":  "vanA/vanB (Vancomycin resistance)",
        "tetA":  "tetA (Tetracycline efflux pump)",
        "tetB":  "tetB (Tetracycline efflux pump)",
        "ermB":  "ermB (rRNA methyltransferase)",
        "aacA":  "aac (Aminoglycoside acetyltransferase)",
        "qnr":   "qnr (Quinolone resistance protein)",
        "blaTEM":"blaTEM / blaSHV (Beta-lactamase)",
    }

    # Aliases so a marker discovered under one name can borrow another
    # marker's classification metadata when it has none of its own.
    MARKER_ALIASES = {
        "tetB": "tetB",  # now has its own entry; kept explicit for clarity
    }

    @staticmethod
    def _execute_cloud_blast(blast_program, database_target, sequence,
                              extra_params, sequence_hash, timeout_sec=300):
        # FIX (bug #7): identity_threshold used to be baked into the cache
        # key even though it only affects *post-hoc filtering* of hits, not
        # the actual NCBIWWW.qblast request. That meant moving the sidebar
        # slider forced a brand-new, identical 1-5 minute NCBI query instead
        # of reusing the cached raw XML and simply re-filtering it. The raw
        # BLAST result is now cached independent of the threshold; filtering
        # happens afterward in run_blast().
        param_parts = []
        if extra_params:
            param_parts = [f"{k}:{repr(extra_params[k])}" for k in sorted(extra_params)]
        cache_parts = [blast_program, database_target, sequence_hash] + param_parts
        cache_key = '|'.join(cache_parts)
        if cache_key in st.session_state.result_cache:
            st.toast("📦 Using cached result", icon="⚡")
            return st.session_state.result_cache[cache_key]
        elapsed = time.time() - st.session_state.last_request_time
        if elapsed < 2:
            time.sleep(2 - elapsed)

        def _worker():
            handle = NCBIWWW.qblast(program=blast_program, database=database_target,
                                    sequence=sequence, format_type="XML", **extra_params)
            raw = handle.read()
            return raw.decode("utf-8") if isinstance(raw, bytes) else raw

        # Not a `with` block: ThreadPoolExecutor.__exit__ calls
        # shutdown(wait=True), which blocked until qblast finished and made
        # the timeout useless. On timeout we abandon the worker thread.
        ex = concurrent.futures.ThreadPoolExecutor(max_workers=1)
        fut = ex.submit(_worker)
        try:
            xml_data = fut.result(timeout=timeout_sec)
        except concurrent.futures.TimeoutError:
            fut.cancel()
            raise TimeoutError(f"NCBI BLAST timeout (>{timeout_sec}s)")
        finally:
            ex.shutdown(wait=False)

        st.session_state.result_cache[cache_key] = xml_data
        st.session_state.last_request_time = time.time()
        return xml_data

    @classmethod
    def run_blast(cls, sequence, use_cloud=True, identity_threshold=90.0):
        is_valid, err, seq_clean = SequenceValidator.validate_and_clean(sequence)
        if not is_valid:
            raise ValueError(err)

        seq_hash        = hashlib.sha256(seq_clean.encode()).hexdigest()[:8]
        is_protein      = SequenceValidator.is_protein(seq_clean)
        likely_coding   = SequenceValidator.is_likely_coding(seq_clean)
        blast_program   = "blastp" if is_protein else ("blastx" if likely_coding else "blastn")
        database_target = "nr" if blast_program in ("blastp", "blastx") else "nt"
        st.info(f"🧬 {'Protein' if is_protein else ('Protein-coding DNA' if likely_coding else 'DNA')} sequence → {blast_program.upper()} vs {database_target}")

        if use_cloud:
            extra_params = {"expect": 1e-10 if blast_program == 'blastn' else 1e-3}
            if len(seq_clean) < 30:
                extra_params["word_size"] = 7 if blast_program == 'blastn' else 2
            try:
                with st.spinner("⏳ NCBI BLAST (1–5 min)..."):
                    xml_data = cls._execute_cloud_blast(
                        blast_program, database_target, seq_clean,
                        extra_params, seq_hash, timeout_sec=300)
                hits = cls._parse_blast_records(NCBIXML.parse(StringIO(xml_data)), seq_clean)
                score_thresh = identity_threshold if blast_program == 'blastn' else 35
                evalue_thresh = 1e-5 if blast_program == 'blastn' else 1e-2
                hits = [h for h in hits
                        if h["identity_percentage"] >= score_thresh
                        and h["e_value"] <= evalue_thresh]
                if hits:
                    st.success(f"✅ {len(hits)} matches passed filters")
                    return hits, "NCBI_CLOUD"

                # Only blastx queries are nucleotide; a blastp query is a
                # protein sequence and cannot be searched with blastn.
                if blast_program == 'blastx':
                    st.info("ℹ️ No protein-level hits; retrying BLASTN vs nt as fallback")
                    try:
                        xml_data2 = cls._execute_cloud_blast(
                            'blastn', 'nt', seq_clean, extra_params, seq_hash + '_blastn', timeout_sec=300)
                        hits2 = cls._parse_blast_records(NCBIXML.parse(StringIO(xml_data2)), seq_clean)
                        hits2 = [h for h in hits2
                                 if h["identity_percentage"] >= identity_threshold
                                 and h["e_value"] <= 1e-5]
                        if hits2:
                            st.success(f"✅ {len(hits2)} matches passed filters (blastn fallback)")
                            return hits2, "NCBI_CLOUD_BLASTN_FALLBACK"
                    except Exception:
                        st.warning("⚠️ BLASTN fallback failed or timed out")

                if blast_program == 'blastn' and len(seq_clean) >= 90:
                    st.info("ℹ️ No strong nucleotide hits; trying BLASTX vs nr for protein-coding gene detection")
                    extra_params2 = {"expect": 1e-3}
                    if len(seq_clean) < 150:
                        extra_params2["word_size"] = 2
                    try:
                        xml_data3 = cls._execute_cloud_blast(
                            'blastx', 'nr', seq_clean, extra_params2, seq_hash + '_blastx', timeout_sec=300)
                        hits3 = cls._parse_blast_records(NCBIXML.parse(StringIO(xml_data3)), seq_clean)
                        hits3 = [h for h in hits3
                                 if h["identity_percentage"] >= 35
                                 and h["e_value"] <= 1e-2]
                        if hits3:
                            st.success(f"✅ {len(hits3)} matches passed filters (blastx fallback)")
                            return hits3, "NCBI_CLOUD_BLASTX_FALLBACK"
                    except Exception:
                        st.warning("⚠️ BLASTX fallback failed or timed out")

                st.warning("⚠️ No hits passed thresholds")
                return [], "NCBI_CLOUD"
            except TimeoutError as exc:
                st.warning(f"⚠️ {exc}")
                if blast_program == 'blastn' and len(seq_clean) >= 90:
                    st.info("ℹ️ BLASTN timed out; retrying BLASTX vs nr")
                    try:
                        extra_params2 = {"expect": 1e-3}
                        if len(seq_clean) < 150:
                            extra_params2["word_size"] = 2
                        xml_data3 = cls._execute_cloud_blast(
                            'blastx', 'nr', seq_clean, extra_params2, seq_hash + '_blastx', timeout_sec=300)
                        hits3 = cls._parse_blast_records(NCBIXML.parse(StringIO(xml_data3)), seq_clean)
                        hits3 = [h for h in hits3
                                 if h["identity_percentage"] >= 35
                                 and h["e_value"] <= 1e-2]
                        if hits3:
                            st.success(f"✅ {len(hits3)} matches passed filters (blastx timeout fallback)")
                            return hits3, "NCBI_CLOUD_BLASTX_TIMEOUT_FALLBACK"
                    except Exception:
                        st.warning("⚠️ BLASTX timeout fallback failed")
                elif blast_program == 'blastx':
                    st.info("ℹ️ Protein-level search timed out; retrying BLASTN vs nt")
                    try:
                        xml_data2 = cls._execute_cloud_blast(
                            'blastn', 'nt', seq_clean, extra_params, seq_hash + '_blastn', timeout_sec=300)
                        hits2 = cls._parse_blast_records(NCBIXML.parse(StringIO(xml_data2)), seq_clean)
                        hits2 = [h for h in hits2
                                 if h["identity_percentage"] >= identity_threshold
                                 and h["e_value"] <= 1e-5]
                        if hits2:
                            st.success(f"✅ {len(hits2)} matches passed filters (blastn timeout fallback)")
                            return hits2, "NCBI_CLOUD_BLASTN_TIMEOUT_FALLBACK"
                    except Exception:
                        st.warning("⚠️ BLASTN timeout fallback failed")
                st.error("❌ NCBI BLAST timeout (>300s) — try again with a shorter query or wait longer")
                return [], "NCBI_CLOUD_TIMEOUT"
            except Exception as exc:
                st.error(f"❌ NCBI BLAST error: {exc}")
                raise

        st.warning("⚠️ Demo mode — randomly generated")
        return cls._generate_demo_hits(seq_clean), "DEMO_MODE"

    SEQUENCE_MOTIFS = {
        "tetA":  ([
            "ATGAAATCTAACAATGCG",
            "TCGATGCTGTAGGCATAGGC",
            "GCGCGGTGGCGTGGTATGC",
        ], 1),
        "mecA":  (["ATGAAAAAGATAAAAATTGTTCAA", "GCTTTAGTTGCAATACTG"], 1),
        "vanA":  (["ATGGGTTCCGATGATGCGGGAG",  "GCAGCAACGAGTCAGCAAC"],  1),
        "ndm":   (["ATGGAATTGCCCAATATTATGCAC","GGTTTGGCGATCTGGTTTTC"],  1),
        "kpc":   (["ATGTCACTGTATCGCCGTCTAGTT","GCAGTCTAGTTCCGCTAAG"],   1),
        "mcr":   (["ATGCGGTTGATCTTCCTGTTATCT","GCTATGCGTTGGTCGGCG"],    1),
        "ctx-m": (["ATGTTAAGCCGTTTCGCAACGCAG","GAACGTTGCCGCAGCCAG"],    1),
        "blaTEM":(["ATGAGTATTCAACATTTCCGTGTC","TTCGGGGAAATGTGCGCG"],    1),
        "ermB":  (["ATGACAAATTTAAAAATAGATGTT","AAGACCCAAGAAATCGAG"],     1),
        "tetB":  (["ATGAAACTCTTGAAATTCATCTTC","TTCGGGATCGTCAACTTC"],    1),
        "qnr":   (["ATGGAAATAATATTGGCGGGTTTTAACCTGAACT"], 1),
    }

    @classmethod
    def _match_title(cls, text: str):
        for marker, (pattern, cls_, mech, res, concern) in cls.LOCAL_MARKERS.items():
            if re.search(pattern, text, re.IGNORECASE):
                display = cls.LOCAL_MARKER_DISPLAY.get(marker, marker)
                return {"marker_found": marker, "display_name": display, "class": cls_, "mechanism": mech,
                        "resistance": res, "clinical_concern": concern}
        return None

    @classmethod
    def _normalise(cls, text: str) -> str:
        text = text.lower()
        text = re.sub(r'[_\-]', ' ', text)
        return re.sub(r'\s+', ' ', text).strip()

    @classmethod
    def _entrez_fetch_gene_name(cls, accession: str) -> str:
        if not accession:
            return ""
        # blastx/blastp hits are protein accessions; querying only the
        # nucleotide db silently returned nothing for them.
        for db in ("nucleotide", "protein"):
            try:
                time.sleep(0.35)
                r = requests.get(
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                    params={"db": db, "id": accession, "retmode": "json"},
                    timeout=15)
                r.raise_for_status()
                data = r.json()
                uids = list(data.get("result", {}).get("uids", []))
                if not uids:
                    continue
                rec = data["result"][uids[0]]
                parts = [rec.get("title", ""), rec.get("subname", "")]
                text = " ".join(p for p in parts if isinstance(p, str) and p)
                if text:
                    return text
            except Exception:
                continue
        return ""

    @classmethod
    def _sequence_motif_scan(cls, sequence: str):
        """Last-resort classifier: scan the query sequence itself for known
        conserved k-mer signatures when BLAST titles are all uninformative.

        BUG FIXED (#6): a motif marker (e.g. "tetB") that had no matching
        LOCAL_MARKERS key used to silently fall back to a generic
        "Undetermined" tuple, losing real classification info. We now
        resolve through MARKER_ALIASES first."""
        seq_upper = sequence.upper()
        for marker, (motifs, min_hits) in cls.SEQUENCE_MOTIFS.items():
            hits = sum(1 for m in motifs if m.upper() in seq_upper)
            if hits >= min_hits:
                lookup_key = cls.MARKER_ALIASES.get(marker, marker)
                _, cls_, mech, res, concern = cls.LOCAL_MARKERS.get(
                    lookup_key,
                    (None, "Resistance Gene", "Unknown", "Undetermined", "HIGH"))
                display = cls.LOCAL_MARKER_DISPLAY.get(lookup_key, marker)
                return {"marker_found": marker, "display_name": display, "class": cls_, "mechanism": mech,
                        "resistance": res, "clinical_concern": concern,
                        "source": "SEQUENCE_MOTIF"}
        return None

    @classmethod
    def classify_gene(cls, sequence, sample_id, hit_title, all_hits=None):
        """
        Four-layer classification, most reliable → least reliable:
        Layer 1 — ResFinder API (purpose-built AMR DB, disabled for now)
        Layer 2 — Pattern match across ALL BLAST hit titles + sample ID
        Layer 3 — Entrez eSummary fetch for the best-hit accession
        Layer 4 — Sequence-level conserved k-mer motif scan
        Returns NOVEL_VARIANT only if all four layers fail.

        Every branch now tags a "source" key so ConfidenceScorer can weight
        confidence by how reliable the classification path was (bug #10).
        """
        rf = ResFinderWebClient.analyze(sequence, sample_id)
        if rf:
            rf.setdefault("source", "RESFINDER")
            return rf

        titles_to_scan = [hit_title] if not all_hits else (
            [h["title"] for h in all_hits[:20]] + [hit_title])
        titles_to_scan.append(sample_id)

        for title in titles_to_scan:
            result = cls._match_title(cls._normalise(title))
            if result:
                result["source"] = "LOCAL_DATABASE"
                return result

        if all_hits:
            accession = all_hits[0].get("accession", "")
            if accession:
                with st.spinner(f"🔍 Fetching gene name from NCBI for {accession}…"):
                    entrez_text = cls._entrez_fetch_gene_name(accession)
                if entrez_text:
                    result = cls._match_title(cls._normalise(entrez_text))
                    if result:
                        result["source"] = "ENTREZ_FETCH"
                        return result

        motif_result = cls._sequence_motif_scan(sequence)
        if motif_result:
            return motif_result

        return {"marker_found": "NOVEL_VARIANT", "class": "Requires curation",
                "mechanism": "Unknown", "resistance": "Undetermined",
                "clinical_concern": "INVESTIGATE", "source": "MANUAL_REVIEW"}

    @staticmethod
    def _generate_demo_hits(sequence):
        km = {
            "meca": ("Staphylococcus aureus mecA",     "NG_047923.1", (97.5,99.9),(1e-150,1e-100)),
            "van":  ("Enterococcus faecium vanA",      "NG_055612.1", (94.2,98.8),(1e-120,1e-80)),
            "tet":  ("Klebsiella pneumoniae tetA",     "NG_048211.3", (95.0,99.2),(1e-100,1e-50)),
            "ctx":  ("Escherichia coli blaCTX-M-15",  "NG_048931.1", (96.0,100.0),(1e-140,1e-90)),
        }
        sl = sequence.lower()
        title, acc, ir, er = next((v for k,v in km.items() if k in sl), km["ctx"])
        identity = random.uniform(*ir)
        e_value  = 10 ** random.uniform(np.log10(er[0]), np.log10(er[1]))
        return [{"title": title, "accession": acc, "e_value": e_value,
                 "score": float(random.randint(300,500)), "identity_percentage": identity,
                 "alignment_length": random.randint(180,280),
                 "query_coverage": random.uniform(85,100)}]

    @staticmethod
    def _parse_blast_records(blast_records, query_sequence):
        parsed, qlen = [], len(query_sequence)
        try:
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        id_pct = (hsp.identities / hsp.align_length * 100) if hsp.align_length > 0 else 0
                        # Coverage from query coordinates (always in query
                        # units). align_length is in amino acids for blastx,
                        # which understated coverage ~3x against a DNA query.
                        q_span = abs(hsp.query_end - hsp.query_start) + 1

                        raw_title = getattr(alignment, "title", "") or getattr(alignment, "hit_def", "") or alignment.accession
                        parts = re.split(r'\s+>\s+|;|\||,', raw_title)
                        clean_title = next((p.strip() for p in parts if p and p.strip()), alignment.accession)
                        clean_title = re.sub(r'\s+', ' ', re.sub(r'[^\x20-\x7E]', ' ', clean_title)).strip()

                        parsed.append({
                            "title":               clean_title,
                            "accession":           alignment.accession,
                            "e_value":             hsp.expect,
                            "score":               float(hsp.score),
                            "identity_percentage": id_pct,
                            "alignment_length":    hsp.align_length,
                            "query_coverage":      min(100.0, q_span / qlen * 100) if qlen > 0 else 0,
                        })
        except Exception as exc:
            st.error(f"BLAST parse error: {exc}")
            return []

        GENOME_NOISE = re.compile(
            r'chromosome\s+\d|complete\s+(genome|sequence)|whole.genome|'
            r'plasmid\s+\w+,\s+complete|contig\s+\d|scaffold\s+\d|'
            r'NODE_\d|genome\s+assembly',
            re.IGNORECASE)

        GENE_SIGNAL = re.compile(
            r'\b(tet[a-z]|van[a-z]|mec[a-z]|bla\w*|ndm|kpc|mcr|ctx.?m|erm[a-z]|'
            r'qnr|aac|sul|resistance\s+gene|resistance\s+protein|'
            r'tetracycline|aminoglycoside|beta.lactam|carbapenem)\b',
            re.IGNORECASE)

        def _sort_key(h):
            title = h["title"]
            is_noise  = bool(GENOME_NOISE.search(title))
            has_gene  = bool(GENE_SIGNAL.search(title))
            tier = 0 if has_gene else (2 if is_noise else 1)
            return (tier, h["e_value"])

        return sorted(parsed, key=_sort_key)


# =====================================================================
# 3b. FASTA PARSING (FIX: robust decode + now actually wired into the UI)
# =====================================================================
def parse_fasta_files(uploaded_files):
    """Parse one or more uploaded FASTA files into {id: sequence}.

    BUG FIXED (#2): the previous version called
    `file.getvalue().decode("utf-8")` *outside* any try/except, so a
    non-UTF-8 upload (e.g. a stray binary file) crashed the whole app.
    Decoding is now defensive (errors='replace') and this function is
    now actually called from the Genomic Surveillance tab's batch
    uploader below, instead of being dead code."""
    sequences = {}
    for file in uploaded_files:
        try:
            content = file.getvalue().decode("utf-8", errors="replace")
        except Exception:
            continue
        try:
            parsed_any = False
            for record in SeqIO.parse(StringIO(content), "fasta"):
                sequences[record.id] = str(record.seq)
                parsed_any = True
            if not parsed_any:
                raise ValueError("no FASTA records found")
        except Exception:
            name = file.name.rsplit(".", 1)[0]
            raw = "".join(l for l in content.splitlines() if not l.startswith(">"))
            if raw:
                sequences[name] = raw
    return sequences


# =====================================================================
# 4. EPIDEMIC DETECTION ENGINE
# =====================================================================
class EpidemicDetectionEngine:
    DEMO_ARTICLES = {
        "covid_2020": "WHO reports 15 cases of COVID-19 in Beijing hospital. Cases show severe respiratory symptoms. Hospital declared outbreak. Baseline: <1 case/day. Current trend: Rising. Contact tracing initiated.",
        "mpox_2022":  "3 suspected mpox cases reported in Lagos clinic. Patients show characteristic rash and fever. Regional health alert issued. Baseline: 0 cases. Current: 3 cases/week. Spreading to community.",
        "normal":     "Routine hospital report: 2 influenza cases detected in pediatric ward. Expected seasonal variation. Baseline: 5-10 cases/week. Normal levels.",
        "fake_1":     "BREAKING: Scientists confirm 5G towers are secretly spreading a new virus across 47 countries. Millions infected. Government covering it up. Share before deleted!!!",
        "fake_2":     "Local hospital cures 100% of Ebola patients overnight using common household bleach. Big Pharma doesn't want you to know. Source: anonymous doctor.",
    }

    PROMPT = """Analyze this health/disease report for epidemic signals.

Report:
{article}

Respond ONLY with valid JSON (no markdown, no extra text):
{{
  "location": "extracted location or Unknown",
  "disease": "disease name or Unknown",
  "cases": 0,
  "severity": "mild|moderate|severe",
  "trend": "stable|rising|falling|critical",
  "baseline_estimate": "e.g. 1-3 cases/day",
  "risk_score": 0.0,
  "summary": "2-3 sentence summary"
}}"""

    @classmethod
    def analyze_article(cls, article_text):
        # No st.* calls here: this also runs in the feed worker thread. The
        # caller sees which path was taken from the returned "source".
        if not anthropic_client:
            return cls._demo_analysis(article_text)

        try:
            result = _claude_json(cls.PROMPT.format(article=article_text), "EPIDEMIC_ANALYSIS")
        except Exception as exc:
            # Credibility dampening is applied once, in run_epidemic_pipeline,
            # for UNVERIFIABLE articles.
            demo = cls._demo_analysis(article_text)
            demo["source"] = "CLAUDE_API_ERROR+KEYWORD_FALLBACK"
            demo["error"] = str(exc)[:200]
            return demo

        return {
            "location":   str(result.get("location") or "Unknown"),
            "disease":    str(result.get("disease") or "Unknown"),
            "cases":      _to_int(result.get("cases", 0)),
            "severity":   str(result.get("severity") or "unknown"),
            "trend":      str(result.get("trend") or "unknown"),
            "baseline":   str(result.get("baseline_estimate") or "Unknown"),
            "risk_score": _to_unit_float(result.get("risk_score", 0.0)),
            "summary":    str(result.get("summary") or ""),
            "source":     "CLAUDE_API",
        }

    @staticmethod
    def _demo_analysis(text):
        tl = text.lower()
        if "ebola" in tl:
            return {"location":"Unknown","disease":"Ebola","cases":3,"severity":"severe",
                    "trend":"critical","baseline":"0","risk_score":0.95,
                    "summary":"Suspected Ebola cases — critical response required.",
                    "source":"KEYWORD_FALLBACK"}
        if "covid" in tl or "sars-cov" in tl:
            return {"location":"Unknown","disease":"COVID-19","cases":15,"severity":"moderate",
                    "trend":"rising","baseline":"<1/day","risk_score":0.82,
                    "summary":"Rising COVID-19 cases above baseline. Outbreak criteria met.",
                    "source":"KEYWORD_FALLBACK"}
        if "mpox" in tl or "monkeypox" in tl:
            return {"location":"Unknown","disease":"Mpox","cases":3,"severity":"moderate",
                    "trend":"rising","baseline":"0","risk_score":0.70,
                    "summary":"Mpox cluster detected above zero baseline. Monitor closely.",
                    "source":"KEYWORD_FALLBACK"}
        return {"location":"Unknown","disease":"Unknown","cases":2,"severity":"mild",
                "trend":"stable","baseline":"5-10/week","risk_score":0.15,
                "summary":"No significant outbreak signal detected. Routine monitoring recommended.",
                "source":"KEYWORD_FALLBACK"}


# =====================================================================
# 4b. ARTICLE CREDIBILITY CHECKER
# =====================================================================
class ArticleCredibilityChecker:

    SATIRE_MARKERS = [
        "theonion", "babylon bee", "the babylon bee", "waterfordwhispersnews",
        "thespoof", "newsbiscuit", "duffelblog", "clickhole", "the beaverton",
        "reductress", "worldnewsdailyreport", "empirenews", "nationalreport",
    ]

    IMPLAUSIBLE_PATTERNS = [
        r"\b\d{7,}\s*(cases|deaths|infected)\b",
        r"100\s*%\s*(mortality|death rate|fatality)",
        r"(cured?|cures?)\s+(covid|cancer|ebola|hiv|aids)\s+in\s+\d+\s+(hour|minute|day)",
        r"(5g|microchip|chemtrail|government.{0,20}hiding|big pharma.{0,20}secret)",
        r"share before (it.?s )?deleted",
        r"doctors? (don.t|do not|won.t|will not) want you to know",
    ]

    CREDIBILITY_SIGNALS = [
        r"\b(hospital|clinic|ward|facility|laboratory|lab|health (centre|center))\b",
        r"\b(patient|case|cases|individual|person|persons|people)\b",
        r"\b(city|country|region|district|province|state|county|area|reported in)\b",
        r"\b(WHO|CDC|NCDC|ECDC|ministry of health|public health|health authority|official)\b",
        r"\b\d{1,2}[\s/\-]\w+[\s/\-]\d{4}|\b(january|february|march|april|may|june|july|august|september|october|november|december)\b",
    ]

    CREDIBILITY_PROMPT = """You are a medical misinformation analyst. Your job is to determine whether this health article is genuine or fabricated/fake.

Article:
{article}

Evaluate carefully:
1. Does it cite a real, verifiable source (WHO, CDC, named hospital, named official with title)?
2. Are the statistics internally consistent and epidemiologically plausible?
3. Does the language match legitimate public health reporting, or is it sensationalist/clickbait?
4. Are there signs of AI generation, fabrication, or conspiracy theory framing?
5. Does it have specific dates, named locations, and named individuals (real reports do; fake ones usually don't)?
6. Are the claims scientifically possible?

Respond ONLY with valid JSON — no markdown, no explanation outside the JSON:
{{
  "verdict": "REAL" or "LIKELY_REAL" or "UNVERIFIABLE" or "LIKELY_FAKE" or "FAKE",
  "credibility_score": 0.0 to 1.0,
  "flags": ["specific concern 1", "specific concern 2"],
  "reasoning": "2-3 sentence explanation of your verdict",
  "source_cited": true or false,
  "specific_details": true or false
}}

Be strict. Public health decisions depend on this assessment. When in doubt, mark UNVERIFIABLE not REAL."""

    @classmethod
    def check(cls, article_text: str) -> Dict[str, Any]:
        heuristic = cls._heuristic_check(article_text)
        if heuristic["verdict"] == "FAKE" and heuristic["credibility_score"] <= 0.1:
            heuristic["method"] = "HEURISTIC"
            return heuristic
        if anthropic_client:
            api_result = cls._api_check(article_text)
            if api_result:
                combined_flags = list(set(
                    heuristic.get("flags", []) + api_result.get("flags", [])))
                api_result["flags"] = combined_flags
                api_result["method"] = "CLAUDE_API"
                if heuristic["verdict"] == "FAKE":
                    api_result["verdict"] = "FAKE"
                    api_result["credibility_score"] = min(
                        api_result["credibility_score"], 0.1)
                return api_result
        heuristic["method"] = "HEURISTIC"
        return heuristic

    @classmethod
    def _heuristic_check(cls, text: str) -> Dict[str, Any]:
        tl   = text.lower()
        flags = []
        score = 1.0

        for marker in cls.SATIRE_MARKERS:
            if marker in tl:
                flags.append(f"Known satire/misinformation source: '{marker}'")
                score -= 0.9

        for pattern in cls.IMPLAUSIBLE_PATTERNS:
            if re.search(pattern, tl, re.IGNORECASE):
                flags.append(f"Implausible or conspiracy claim detected")
                score -= 0.5

        SENSATIONAL_PATTERNS = [
            r"\b(breaking|shocking|miracle|cure|secret|uncovered|guarantee|share before|must read)\b",
            r"100\s*%",
            r"\b(breaking:)\b"
        ]
        for pat in SENSATIONAL_PATTERNS:
            if re.search(pat, tl, re.IGNORECASE):
                flags.append("Sensational language detected — lowers credibility")
                score -= 0.35

        if re.search(r'100\s*%.*credible|100% credible', tl):
            flags.append("Explicit 100% credibility claim — flagged as suspicious")
            score -= 0.9

        found_signals = sum(
            1 for p in cls.CREDIBILITY_SIGNALS
            if re.search(p, tl, re.IGNORECASE))
        if found_signals < 2:
            flags.append(
                f"Only {found_signals}/5 credibility signals present "
                f"(real reports cite sources, locations, and case counts)")
            score -= 0.25

        if len(text.strip()) < 80:
            flags.append("Text too short to be a credible health report")
            score -= 0.3

        caps_ratio = sum(1 for c in text if c.isupper()) / max(len(text), 1)
        if caps_ratio > 0.35:
            flags.append(f"Excessive capitalisation ({caps_ratio*100:.0f}%) — clickbait signal")
            score -= 0.15

        if not re.search(r'\d+', text):
            flags.append("No numerical data — real outbreak reports always include case counts")
            score -= 0.15

        if text.count('!') >= 3:
            flags.append("Multiple exclamation marks — sensationalist tone")
            score -= 0.2

        score = float(np.clip(score, 0.0, 1.0))

        if   score >= 0.75: verdict = "LIKELY_REAL"
        elif score >= 0.50: verdict = "UNVERIFIABLE"
        elif score >= 0.20: verdict = "LIKELY_FAKE"
        else:               verdict = "FAKE"

        return {
            "verdict":           verdict,
            "credibility_score": score,
            "flags":             flags,
            "reasoning":         "Fast heuristic analysis based on structural and content signals.",
            "source_cited":      bool(re.search(
                r'\b(WHO|CDC|NCDC|ministry|hospital|official|reported by)\b',
                text, re.IGNORECASE)),
            "specific_details":  bool(
                re.search(r'\b\d{4}\b', text) and
                re.search(r'\b[A-Z][a-z]{2,}\b', text)),
        }

    @classmethod
    def _api_check(cls, text: str) -> Optional[Dict[str, Any]]:
        try:
            result = _claude_json(cls.CREDIBILITY_PROMPT.format(article=text), "CREDIBILITY_CHECK")
        except Exception:
            return None

        def _as_bool(v):
            # bool("false") is True; handle string booleans explicitly.
            return v.strip().lower() == "true" if isinstance(v, str) else bool(v)

        verdict = str(result.get("verdict", "UNVERIFIABLE")).upper()
        if verdict not in ("REAL", "LIKELY_REAL", "UNVERIFIABLE", "LIKELY_FAKE", "FAKE"):
            verdict = "UNVERIFIABLE"
        flags = result.get("flags", [])
        if not isinstance(flags, list):
            flags = [str(flags)]
        return {
            "verdict":           verdict,
            "credibility_score": _to_unit_float(result.get("credibility_score", 0.5), 0.5),
            "flags":             [str(f) for f in flags],
            "reasoning":         str(result.get("reasoning", "")),
            "source_cited":      _as_bool(result.get("source_cited", False)),
            "specific_details":  _as_bool(result.get("specific_details", False)),
        }

    @staticmethod
    def render_verdict(result: Dict[str, Any]) -> bool:
        verdict = result.get("verdict", "UNVERIFIABLE")
        score   = result.get("credibility_score", 0.5)
        flags   = result.get("flags", [])
        method  = result.get("method", "UNKNOWN")

        STYLES = {
            "REAL":          ("✅ CREDIBLE SOURCE",       "#0D3320", "#4CAF50"),
            "LIKELY_REAL":   ("✅ LIKELY CREDIBLE",       "#0D2A1A", "#7BC67E"),
            "UNVERIFIABLE":  ("⚠️ UNVERIFIABLE",          "#2A2500", "#FFC107"),
            "LIKELY_FAKE":   ("🚨 LIKELY FAKE",           "#2A1500", "#FF9800"),
            "FAKE":          ("🚫 FAKE / MISINFORMATION", "#2A0000", "#E53935"),
        }
        label, bg, border = STYLES.get(verdict, STYLES["UNVERIFIABLE"])

        st.markdown(
            f'<div style="background:{bg};border-left:5px solid {border};'
            f'padding:14px 18px;border-radius:8px;margin-bottom:14px;">'
            f'<div style="font-size:16px;font-weight:700;color:{border}">{label}</div>'
            f'<div style="font-size:12px;color:#bbb;margin-top:3px">'
            f'Credibility: {score*100:.0f}%'
            f'&nbsp;·&nbsp;Method: {method}'
            f'&nbsp;·&nbsp;Source cited: {"✅" if result.get("source_cited") else "❌"}'
            f'&nbsp;·&nbsp;Specific details: {"✅" if result.get("specific_details") else "❌"}'
            f'</div>'
            f'<div style="font-size:13px;color:#ddd;margin-top:8px">'
            f'{html.escape(str(result.get("reasoning","")))}</div>'
            f'</div>',
            unsafe_allow_html=True)

        if flags:
            with st.expander(f"⚑ {len(flags)} flag(s) detected — click to expand"):
                for f in flags:
                    st.markdown(f"- {f}")

        blocked = verdict in ("FAKE", "LIKELY_FAKE")
        if blocked:
            st.error(
                "🚫 **Analysis blocked.** This article has been flagged as fake or "
                "likely misinformation. Generating a risk score or HL7 output from "
                "fabricated data would be dangerous in a clinical surveillance context. "
                "Please submit a verified health report.")
        return blocked


# =====================================================================
# 5. CONFIDENCE SCORER
# =====================================================================
class ConfidenceScorer:
    # FIX (bug #10): the previous formula gave a boolean "reference_verified"
    # credit — a Layer-2 BLAST-title match and a Layer-4 raw-sequence motif
    # guess both scored the same 0.15 "reference verified" weight, overstating
    # confidence for the weakest classification path. Now weighted by the
    # actual classification source/tier.
    SOURCE_TIER_WEIGHTS = {
        "RESFINDER":       1.00,
        "LOCAL_DATABASE":  0.90,
        "ENTREZ_FETCH":    0.75,
        "SEQUENCE_MOTIF":  0.50,
        "MANUAL_REVIEW":   0.20,
    }

    @classmethod
    def amr_confidence(cls, identity, e_value, coverage, source="MANUAL_REVIEW"):
        i = np.clip(identity / 100.0, 0, 1)
        e = np.clip(-np.log10(np.clip(e_value, 1e-300, 1.0)) / 300.0, 0, 1)
        c = np.clip(coverage / 100.0, 0, 1)
        r = cls.SOURCE_TIER_WEIGHTS.get(source, 0.5)
        return float(np.clip(0.40*i + 0.30*e + 0.15*c + 0.15*r, 0, 1))

    @staticmethod
    def amrfinder_confidence(identity, coverage, method):
        """Confidence for an AMRFinderPlus call. There is no E-value; the
        evidence is identity to the curated reference, how much of the
        reference gene is covered, and the detection method (exact allele
        match > BLAST above curated cutoff > partial > HMM-only)."""
        i = np.clip(identity / 100.0, 0, 1)
        c = np.clip(coverage / 100.0, 0, 1)
        m = AMRFinderClient.method_weight(method)
        return float(np.clip(0.45*i + 0.30*c + 0.25*m, 0, 1))

    @staticmethod
    def risk_category(score):
        if score >= 0.85: return "CRITICAL", "🔴", "#FF4B4B"
        if score >= 0.70: return "HIGH",     "🟠", "#FF9800"
        if score >= 0.50: return "MODERATE", "🟡", "#FFC107"
        if score >= 0.30: return "LOW",      "🟢", "#4CAF50"
        return                   "MINIMAL",  "⚪", "#9E9E9E"


# =====================================================================
# 6. HL7 GENERATOR
# =====================================================================
class HL7Generator:
    @staticmethod
    def _sf(v):
        return str(v).replace("\\","\\E\\").replace("^","\\S\\").replace("|","\\F\\").replace("&","\\T\\")

    @staticmethod
    def _ef(e):
        if e is None:
            return "N/A"   # AMRFinderPlus calls have no E-value
        return "< 1e-300" if (e == 0 or e < 1e-300) else f"{e:.2e}"

    @classmethod
    def amr_report(cls, sample_id, gene, confidence, identity, e_value, concern):
        ts  = datetime.now().strftime("%Y%m%d%H%M%S")
        mid = hashlib.md5(f"{sample_id}{ts}".encode()).hexdigest()[:13]
        af  = "H" if (confidence >= 0.70 or concern in ("CRITICAL","HIGH","INVESTIGATE")) else "N"
        sf, ef = cls._sf, cls._ef
        return "\n".join([
            f"MSH|^~\\&|ONE_HEALTH_PLATFORM|REGIONAL_HUB|LABORATORY_IS|HOSPITAL_SYSTEM|{ts}||ORU^R01|{mid}|P|2.5",
            f"PID|1||{sf(sample_id)}^^^MRN||PATIENT^ANONYMOUS|||||||||||||||||||||||||||",
            f"OBR|1||{sf(sample_id)}|AMR_GENOMICS^Antimicrobial Resistance Genomics|F|||{ts}||||||||||||||||F",
            f"OBX|1|ST|AMR_GENE^Resistance Gene|1|{sf(gene)}|||{af}|||F",
            f"OBX|2|NM|CONFIDENCE^Heuristic Confidence Score|1|{confidence:.4f}|1|0.0^1.0|{af}|||F",
            f"OBX|3|NM|BLAST_IDENTITY^BLAST Identity|1|{identity:.2f}|%|0.0^100.0|N|||F",
            f"OBX|4|ST|EVALUE^NCBI E-value|1|{ef(e_value)}||||N|||F",
            f"OBX|5|ST|CLINICAL_CONCERN^Clinical Risk|1|{sf(concern)}||||N|||F",
        ])

    @classmethod
    def amr_genome_report(cls, sample_id, organism, findings, phenotype=None):
        """ORU^R01 for a whole-genome result: one OBX per resistance
        determinant and, when available, one per predicted drug result."""
        ts  = datetime.now().strftime("%Y%m%d%H%M%S")
        mid = hashlib.md5(f"{sample_id}{ts}".encode()).hexdigest()[:13]
        sf  = cls._sf
        lines = [
            f"MSH|^~\\&|ONE_HEALTH_PLATFORM|REGIONAL_HUB|LABORATORY_IS|HOSPITAL_SYSTEM|{ts}||ORU^R01|{mid}|P|2.5",
            f"PID|1||{sf(sample_id)}^^^MRN||PATIENT^ANONYMOUS|||||||||||||||||||||||||||",
            f"OBR|1||{sf(sample_id)}|AMR_WGS^Antimicrobial Resistance Genome Screen|F|||{ts}||||||||||||||||F",
            f"OBX|1|ST|ORGANISM^Organism (as specified)|1|{sf(organism or 'not specified')}||||N|||F",
        ]
        n = 1
        for f in findings:
            n += 1
            concern = AMRFinderClient.concern(f.get("class", ""), f.get("subclass", ""))
            flag = "H" if concern in ("CRITICAL", "HIGH") else "N"
            # components escaped individually, then joined with the ^ separator
            value = "^".join(sf(c) for c in (f["symbol"], f.get("subclass") or f.get("class") or "",
                                              f.get("method", ""), f"{f.get('identity', 0):.2f}"))
            lines.append(f"OBX|{n}|CWE|AMR_DETERMINANT^Resistance determinant|{n - 1}|{value}|||{flag}|||F")
        for p in phenotype or []:
            n += 1
            lines.append(f"OBX|{n}|ST|PREDICTED_SUSCEPTIBILITY^{sf(p['drug'])} predicted from genotype|1|"
                         f"{'R' if p['predicted_resistant'] else 'S'}||||{'H' if p['predicted_resistant'] else 'N'}|||P")
        return "\n".join(lines)

    @classmethod
    def epidemic_alert(cls, report_id, disease, risk_score, cases, trend, location):
        ts  = datetime.now().strftime("%Y%m%d%H%M%S")
        mid = hashlib.md5(f"{report_id}{ts}".encode()).hexdigest()[:13]
        af  = "H" if risk_score >= 0.70 else "N"
        sf  = cls._sf
        return "\n".join([
            f"MSH|^~\\&|ONE_HEALTH_PLATFORM|REGIONAL_HUB|PH_SYSTEM|NATIONAL_SURVEILLANCE|{ts}||ORU^R01|{mid}|P|2.5",
            f"PID|1||{sf(report_id)}^^^RPT|||||||||||||||||||||||||||",
            f"OBR|1||{sf(report_id)}|EPIDEMIC_ALERT^Epidemic Early Warning|F|||{ts}||||||||||||||||F",
            f"OBX|1|ST|DISEASE^Disease Detected|1|{sf(disease)}|||{af}|||F",
            f"OBX|2|NM|RISK_SCORE^Epidemic Risk Score|1|{risk_score:.4f}|1|0.0^1.0|{af}|||F",
            f"OBX|3|NM|CASES^Case Count|1|{cases}|count|||N|||F",
            f"OBX|4|ST|TREND^Epidemiological Trend|1|{sf(trend)}||||N|||F",
            f"OBX|5|ST|LOCATION^Reported Location|1|{sf(location)}||||N|||F",
        ])


# =====================================================================
# 7. SHARED AMR PIPELINE HELPER (used by both single + batch analysis)
# =====================================================================
AMR_BATCH_SLICE = 150   # sequences per AMRFinderPlus call in batch mode

ENGINE_MODES = {
    "amrfinder":       "AMRFinderPlus (fast, validated)",
    "amrfinder+blast": "AMRFinderPlus, then web BLAST if nothing is found",
    "blast":           "Web BLAST only (slow, 1–5 min)",
}

# AMRFinderPlus --organism values (amrfinder -l, v4.2.7). Setting one adds
# organism-specific point-mutation screening and suppresses genes intrinsic
# to that organism.
AMR_ORGANISMS = {"": "Not specified (acquired genes only)"}
AMR_ORGANISMS.update({o: o.replace("_", " ") + (" / Shigella" if o == "Escherichia" else "") for o in (
    "Acinetobacter_baumannii", "Bordetella_pertussis", "Burkholderia_cepacia", "Burkholderia_mallei",
    "Burkholderia_pseudomallei", "Campylobacter", "Citrobacter_freundii", "Clostridioides_difficile",
    "Corynebacterium_diphtheriae", "Enterobacter_asburiae", "Enterobacter_cloacae", "Enterococcus_faecalis",
    "Enterococcus_faecium", "Escherichia", "Haemophilus_influenzae", "Helicobacter_pylori",
    "Klebsiella_oxytoca", "Klebsiella_pneumoniae", "Neisseria_gonorrhoeae", "Neisseria_meningitidis",
    "Pseudomonas_aeruginosa", "Salmonella", "Serratia_marcescens", "Staphylococcus_aureus",
    "Staphylococcus_epidermidis", "Staphylococcus_pseudintermedius", "Streptococcus_agalactiae",
    "Streptococcus_pneumoniae", "Streptococcus_pyogenes", "Vibrio_cholerae", "Vibrio_parahaemolyticus",
    "Vibrio_vulnificus")})

GENOME_MAX_BP = 20_000_000     # whole-genome upload limit (bacterial genomes are ~1–13 Mb)
GENOME_MIN_BP = 1_000
GENOME_TIMEOUT_S = 3600


# ---------------------------------------------------------------------
# Genotype -> predicted susceptibility (rule-based)
#
# These are the rules evaluated in evaluation/isolate_benchmark.py on 99
# public clinical K. pneumoniae genomes (evaluation/RESULTS.md, section 4).
# They are shown in the app ONLY for PHENOTYPE_VALIDATED_ORGANISMS. Change
# them here and the benchmark uses the change too; re-run it before quoting
# figures. AMRFinderPlus "Subclass" is split on "/" into tokens.
# ---------------------------------------------------------------------
PHENOTYPE_VALIDATED_ORGANISMS = {"Klebsiella_pneumoniae"}
INTRINSIC_GENES = {"oqxA", "oqxB", "fosA"}   # chromosomal in K. pneumoniae; allele variants also matched
_CEPH_TOKENS = {"CEPHALOSPORIN", "CARBAPENEM"}
PHENOTYPE_DRUG_RULES = {
    # 3rd/4th-generation cephalosporins: acquired ESBL, AmpC or carbapenemase gene
    "ceftriaxone":   "ceph", "cefotaxime": "ceph", "ceftazidime": "ceph", "cefepime": "ceph",
    # carbapenems: acquired carbapenemase, or porin disruption + ESBL/AmpC
    "meropenem":     "carb", "imipenem": "carb", "ertapenem": "carb",
    "ciprofloxacin": "QUINOLONE", "levofloxacin": "QUINOLONE",
    "gentamicin":    "GENTAMICIN", "tobramycin": "TOBRAMYCIN", "amikacin": "AMIKACIN",
    "trimethoprim-sulfamethoxazole": "TRIMETHOPRIM",
    "tetracycline":  "TETRACYCLINE",
    "colistin":      "COLISTIN",
}


def _phenotype_usable(findings: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out = []
    for f in findings:
        if (f.get("type") or "").upper() != "AMR":
            continue
        if (f.get("method") or "").upper().startswith("INTERNAL_STOP"):
            continue
        if any(re.fullmatch(re.escape(g) + r"\d*", f.get("symbol") or "") for g in INTRINSIC_GENES):
            continue
        out.append(f)
    return out


def _subclass_tokens(f: Dict[str, Any]) -> set:
    return {t.strip().upper() for t in (f.get("subclass") or "").split("/") if t.strip()}


def _is_point(f: Dict[str, Any]) -> bool:
    return (f.get("subtype") or "").upper().startswith("POINT")


def predict_phenotype(drug: str, findings: List[Dict[str, Any]]):
    """Return (1 = predicted resistant / 0 = susceptible, responsible element symbols)."""
    rule = PHENOTYPE_DRUG_RULES[drug]
    fs = _phenotype_usable(findings)
    if rule == "ceph":
        hits = [f for f in fs if not _is_point(f) and _subclass_tokens(f) & _CEPH_TOKENS]
    elif rule == "carb":
        hits = [f for f in fs if not _is_point(f) and "CARBAPENEM" in _subclass_tokens(f)]
        if not hits:
            porin = [f for f in fs if _is_point(f) and "CARBAPENEM" in _subclass_tokens(f)]
            bl = [f for f in fs if not _is_point(f) and "CEPHALOSPORIN" in _subclass_tokens(f)]
            if porin and bl:
                hits = porin + bl
    else:
        hits = [f for f in fs if rule in _subclass_tokens(f)
                or (rule in ("GENTAMICIN", "TOBRAMYCIN", "AMIKACIN") and _subclass_tokens(f) == {"AMINOGLYCOSIDE"})]
    return (1 if hits else 0), sorted({f["symbol"] for f in hits})


def phenotype_profile(findings: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out = []
    for drug in PHENOTYPE_DRUG_RULES:
        p, why = predict_phenotype(drug, findings)
        out.append({"drug": drug, "predicted_resistant": bool(p), "determinants": why})
    return out


def _amrfinder_result(sample_id: str, seq_clean: str, findings: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Build the pipeline result from AMRFinderPlus findings for one sequence."""
    amr = [f for f in findings if (f["type"] or "AMR").upper() == "AMR"] or findings
    amr = sorted(amr, key=lambda f: f["identity"] * f["coverage"], reverse=True)
    top = amr[0]
    gene_info = AMRFinderClient.to_gene_info(top)
    confidence = ConfidenceScorer.amrfinder_confidence(top["identity"], top["coverage"], top["method"])
    risk, emoji, _ = ConfidenceScorer.risk_category(confidence)
    hits = [{"title": f["name"] or f["symbol"], "accession": f["accession"], "e_value": None,
             "score": None, "identity_percentage": f["identity"], "alignment_length": None,
             "query_coverage": f["coverage"]} for f in amr]
    label = f"{top['symbol']} ({top['name']})" if top["name"] else top["symbol"]
    hl7 = HL7Generator.amr_report(sample_id, label, confidence, top["identity"], None,
                                  gene_info["clinical_concern"])
    info = AMRFinderClient.info()
    return {
        "ok": True, "sample_id": sample_id, "seq_clean": seq_clean,
        "hits": hits, "best": hits[0], "gene_info": gene_info, "findings": amr,
        "confidence": confidence, "risk": risk, "emoji": emoji,
        "display_gene": top["symbol"], "hl7": hl7,
        "engine": "AMRFINDERPLUS", "source": "AMRFINDERPLUS", "is_demo": False,
        "tool_version": info.get("software", "?"), "db_version": info.get("database", "?"),
    }


def run_amr_pipeline(sample_id: str, sequence: str, identity_threshold: float,
                     mode: str = "amrfinder", precomputed: Optional[List[Dict[str, Any]]] = None,
                     organism: Optional[str] = None) -> Dict[str, Any]:
    """Detection + classification + scoring + HL7 for one sequence.

    mode: "amrfinder" | "amrfinder+blast" | "blast" (see ENGINE_MODES).
    If AMRFinderPlus isn't installed, every mode falls back to web BLAST.
    `precomputed` lets the batch runner pass findings from a single
    multi-sequence AMRFinderPlus call. `organism` (an AMR_ORGANISMS key)
    adds point-mutation screening."""
    is_valid, err, seq_clean = SequenceValidator.validate_and_clean(sequence)
    if not is_valid:
        return {"ok": False, "error": err, "sample_id": sample_id}

    if mode != "blast":
        if AMRFinderClient.is_available():
            findings: Optional[List[Dict[str, Any]]] = precomputed
            if findings is None:
                try:
                    with st.spinner("🧬 Running AMRFinderPlus…"):
                        findings = AMRFinderClient.run(
                            {"query": seq_clean}, SequenceValidator.is_protein(seq_clean),
                            organism=organism or None).get("query", [])
                except Exception as exc:
                    if mode == "amrfinder":
                        return {"ok": False, "error": f"❌ {exc}", "sample_id": sample_id}
                    st.warning(f"⚠️ {exc}. Falling back to web BLAST.")
            if findings:
                return _amrfinder_result(sample_id, seq_clean, findings)
            if findings is not None and mode == "amrfinder":
                info = AMRFinderClient.info()
                return {"ok": False, "negative": True, "sample_id": sample_id, "engine": "AMRFINDERPLUS",
                        "error": (f"No known resistance gene detected by AMRFinderPlus "
                                  f"{info.get('software', '?')} (database {info.get('database', '?')}).")}
        else:
            st.info("ℹ️ AMRFinderPlus is not available on this machine; using web BLAST instead.")

    hits, source = AMRValidationEngine.run_blast(
        seq_clean, use_cloud=True, identity_threshold=identity_threshold)
    is_demo = source == "DEMO_MODE"

    if not hits:
        return {"ok": False, "error": "No significant matches found",
                "sample_id": sample_id, "source": source, "is_demo": is_demo}

    best = hits[0]
    gene_info = AMRValidationEngine.classify_gene(
        seq_clean, sample_id, best["title"], all_hits=hits)
    confidence = ConfidenceScorer.amr_confidence(
        best["identity_percentage"], best["e_value"],
        best.get("query_coverage", 50.0),
        source=gene_info.get("source", "MANUAL_REVIEW"))
    risk, emoji, _ = ConfidenceScorer.risk_category(confidence)
    display_gene = gene_info.get("display_name", gene_info.get("marker_found", "UNKNOWN"))
    hl7 = HL7Generator.amr_report(
        sample_id, display_gene, confidence, best["identity_percentage"],
        best["e_value"], gene_info.get("clinical_concern", "UNKNOWN"))

    return {
        "ok": True, "sample_id": sample_id, "seq_clean": seq_clean,
        "hits": hits, "best": best, "gene_info": gene_info,
        "confidence": confidence, "risk": risk, "emoji": emoji,
        "display_gene": display_gene, "hl7": hl7,
        "engine": "BLAST", "source": source, "is_demo": is_demo,
    }


def run_amr_batch(sequences: Dict[str, str], identity_threshold: float, mode: str = "amrfinder",
                  progress=None, organism: Optional[str] = None) -> List[Dict[str, Any]]:
    """Analyse many sequences. With AMRFinderPlus, all nucleotide sequences
    go through one call and all protein sequences through another, instead
    of one process start per sequence."""
    precomputed: Dict[str, List[Dict[str, Any]]] = {}
    if mode != "blast" and AMRFinderClient.is_available():
        cleaned = {}
        for i, (sid, seq) in enumerate(sequences.items()):
            ok, _, sc = SequenceValidator.validate_and_clean(seq)
            if ok:
                cleaned[sid] = (f"s{i}", sc, SequenceValidator.is_protein(sc))
        for protein in (False, True):
            group = {rid: sc for rid, sc, p in cleaned.values() if p == protein}
            if not group:
                continue
            # Sent in slices: genes from large allele families (beta-lactamases)
            # cost ~0.3 s each, so one call with ~800 of them took over 4 minutes
            # in benchmarking and could exceed AMRFinderClient.TIMEOUT.
            rids, out, failed = list(group), {}, None
            for start in range(0, len(rids), AMR_BATCH_SLICE):
                part = {rid: group[rid] for rid in rids[start:start + AMR_BATCH_SLICE]}
                try:
                    out.update(AMRFinderClient.run(part, protein, organism=organism or None))
                except Exception as exc:
                    failed = exc
                    break
            if failed is not None:
                st.warning(f"⚠️ Batch AMRFinderPlus run failed ({failed}); retrying sequences one by one.")
                continue
            for sid, (rid, _, p) in cleaned.items():
                if p == protein:
                    precomputed[sid] = out.get(rid, [])

    results, total = [], len(sequences)
    for i, (sid, seq) in enumerate(sequences.items(), start=1):
        try:
            res = run_amr_pipeline(sid, seq, identity_threshold, mode, precomputed=precomputed.get(sid),
                                   organism=organism)
        except Exception as exc:
            res = {"ok": False, "error": f"❌ Error: {exc}", "sample_id": sid}
        results.append(res)
        if progress:
            progress(i, total)
    return results


def validate_genome(contigs: Dict[str, str]):
    """Check an uploaded assembly. Returns (ok, error, {contig_id: cleaned_seq})."""
    if not contigs:
        return False, "❌ No sequences found in the file", {}
    dna_chars = SequenceValidator.VALID_DNA_CHARS
    cleaned: Dict[str, str] = {}
    for cid, seq in contigs.items():
        s = re.sub(r"\s+", "", seq or "").upper()
        if not s:
            continue
        bad = set(s) - dna_chars
        if bad:
            return False, (f"❌ Contig {cid[:40]} has non-DNA characters ({', '.join(sorted(bad)[:8])}). "
                           "Whole-genome mode takes a nucleotide assembly."), {}
        cleaned[cid] = s
    total = sum(map(len, cleaned.values()))
    if total < GENOME_MIN_BP:
        return False, f"❌ Assembly too small ({total:,} bp; minimum {GENOME_MIN_BP:,})", {}
    if total > GENOME_MAX_BP:
        return False, f"❌ Assembly too large ({total:,} bp; maximum {GENOME_MAX_BP:,})", {}
    return True, "", cleaned


def run_amr_genome(sample_id: str, contigs: Dict[str, str], organism: Optional[str] = None) -> Dict[str, Any]:
    """Screen a whole-genome assembly (all contigs = one isolate) with
    AMRFinderPlus. Reports every determinant, point mutations when an
    organism is given, and rule-based predicted susceptibility for organisms
    in PHENOTYPE_VALIDATED_ORGANISMS. No st.* calls."""
    ok, err, cleaned = validate_genome(contigs)
    if not ok:
        return {"ok": False, "error": err, "sample_id": sample_id}
    if not AMRFinderClient.is_available():
        return {"ok": False, "sample_id": sample_id,
                "error": "❌ Whole-genome screening needs AMRFinderPlus, which is not available on this machine "
                         "(see the Status page). Web BLAST cannot search a whole genome."}
    organism = organism or None
    if organism and organism not in AMR_ORGANISMS:
        return {"ok": False, "error": f"❌ Unknown organism: {organism}", "sample_id": sample_id}

    ids = {f"c{i}": cid for i, cid in enumerate(cleaned, start=1)}
    joined = "".join(cleaned.values())
    t0 = time.time()
    grouped = AMRFinderClient.run({rid: cleaned[cid] for rid, cid in ids.items()}, protein=False,
                                  organism=organism, timeout=GENOME_TIMEOUT_S, threads=os.cpu_count() or 4)
    findings = []
    for rid, fs in grouped.items():
        for f in fs:
            findings.append(dict(f, seq_id=ids.get(rid, rid)))
    findings = [f for f in findings if (f["type"] or "AMR").upper() == "AMR"] or findings

    info = AMRFinderClient.info()
    meta = {
        "input_kind": "GENOME", "organism": organism, "contigs": len(cleaned), "genome_bp": len(joined),
        "seconds": round(time.time() - t0, 1),
        "sequence_length": len(joined), "sequence_sha256": hashlib.sha256(joined.encode()).hexdigest(),
    }
    if not findings:
        return dict(meta, ok=False, negative=True, sample_id=sample_id, engine="AMRFINDERPLUS",
                    error=(f"No known resistance determinant detected by AMRFinderPlus "
                           f"{info.get('software', '?')} (database {info.get('database', '?')})."))

    phenotype = phenotype_profile(findings) if organism in PHENOTYPE_VALIDATED_ORGANISMS else None
    result = _amrfinder_result(sample_id, "", findings)
    result.update(meta)
    result["findings"] = sorted(findings, key=lambda f: (f.get("class") or "", f["symbol"]))
    result["phenotype"] = phenotype
    result["display_gene"] = f"{len(findings)} determinant(s)"
    result["hl7"] = HL7Generator.amr_genome_report(sample_id, AMR_ORGANISMS.get(organism or "", organism),
                                                   result["findings"], phenotype)
    return result


def log_amr_result(result: Dict[str, Any]) -> Optional[int]:
    """Save a successful genomic result to the database. Returns its id."""
    if not result.get("ok"):
        return None
    seq = result.get("seq_clean", "")
    details = {
        "engine": result.get("engine"), "gene_info": result["gene_info"], "best": result["best"],
        "findings": result.get("findings", []), "risk": result["risk"],
        "tool_version": result.get("tool_version"), "db_version": result.get("db_version"),
        "sequence_length": result.get("sequence_length", len(seq)),
        "sequence_sha256": result.get("sequence_sha256") or hashlib.sha256(seq.encode()).hexdigest(),
    }
    for k in ("input_kind", "organism", "contigs", "phenotype"):
        if result.get(k) is not None:
            details[k] = result[k]
    return get_store().log_analysis(
        result["sample_id"], "GENOMIC", result["display_gene"], result["confidence"],
        result["source"], result["is_demo"], details=details, hl7=result["hl7"])


# =====================================================================
# 7b. SHARED EPIDEMIC PIPELINE (Epidemic page + outbreak feed)
# No st.* calls: it also runs in the feed worker thread.
# =====================================================================
def run_epidemic_pipeline(report_id: str, article: str, store: Optional[Store] = None,
                          origin: str = "MANUAL") -> Dict[str, Any]:
    """Credibility screen -> outbreak extraction -> risk score -> HL7.
    Saves non-blocked results to the database."""
    cred_result = ArticleCredibilityChecker.check(article)
    blocked = cred_result.get("verdict") in ("FAKE", "LIKELY_FAKE")
    out: Dict[str, Any] = {"ok": True, "report_id": report_id, "cred_result": cred_result, "blocked": blocked}
    if blocked:
        return out

    analysis = EpidemicDetectionEngine.analyze_article(article)
    # Any non-Claude result is a keyword guess, including the API-error fallback.
    is_demo = analysis["source"] != "CLAUDE_API"
    risk_score = analysis["risk_score"]
    dampened_note = None

    if cred_result["verdict"] == "UNVERIFIABLE":
        credibility_score = _to_unit_float(cred_result.get("credibility_score", 0.5), 0.5)
        floor = 0.5 if risk_score >= 0.70 else 0.35
        # The floor keeps serious signals from being dampened to nothing,
        # but must never *raise* a score above its original value.
        new_risk = float(np.clip(max(risk_score * credibility_score, min(floor, risk_score)), 0.0, 1.0))
        dampened_note = (f"⚠️ Risk score dampened from {risk_score*100:.0f}% → {new_risk*100:.0f}% "
                         "due to low article credibility. Treat this result as provisional until verified.")
        risk_score = new_risk

    risk, emoji, _ = ConfidenceScorer.risk_category(risk_score)
    hl7 = HL7Generator.epidemic_alert(report_id, analysis["disease"], risk_score,
                                      analysis["cases"], analysis["trend"], analysis["location"])
    out.update({"analysis": analysis, "is_demo": is_demo, "risk_score": risk_score, "risk": risk,
                "emoji": emoji, "dampened_note": dampened_note, "hl7": hl7})

    location = str(analysis.get("location", "")).strip()
    out["analysis_id"] = (store or get_store()).log_analysis(
        report_id, "EPIDEMIC", analysis["disease"], risk_score, analysis["source"], is_demo,
        # Don't geocode the literal string "Unknown".
        location="" if location.lower() in ("", "unknown", "n/a", "none") else location,
        details={"analysis": analysis, "risk": risk, "origin": origin,
                 "credibility": {k: cred_result.get(k) for k in
                                 ("verdict", "credibility_score", "method", "flags", "reasoning")}},
        hl7=hl7)
    return out


# =====================================================================
# 7c. OUTBREAK FEED (WHO Disease Outbreak News)
# A daemon thread fetches new items every FEED_INTERVAL_S while the app is
# running. New items are analysed automatically only when Claude is
# connected and the "feed_auto_analyze" setting is on, capped per cycle so
# a backlog can never trigger a burst of API spend.
# =====================================================================
FEED_INTERVAL_S = 6 * 3600
FEED_AUTO_ANALYZE_CAP = 5
_FEED_LOCK = threading.Lock()


def analyze_feed_item(item: Dict[str, Any], store: Store) -> Dict[str, Any]:
    try:
        out = run_epidemic_pipeline(item["id"], item["text"], store=store, origin=item.get("source", "FEED"))
        if out["blocked"]:
            store.set_feed_status(item["id"], "blocked",
                                  note=f"Credibility check: {out['cred_result'].get('verdict')}")
        else:
            store.set_feed_status(item["id"], "analyzed", analysis_id=out.get("analysis_id"))
        return out
    except Exception as exc:
        store.set_feed_status(item["id"], "error", note=str(exc)[:200])
        return {"ok": False, "error": str(exc)}


def sync_feed(store: Store, analyze: Optional[bool] = None) -> Dict[str, Any]:
    """Fetch the newest WHO DON items, store unseen ones, optionally analyse
    new ones. `analyze=None` follows the saved auto-analyze setting."""
    summary = {"fetched": 0, "new": 0, "analyzed": 0, "error": ""}
    with _FEED_LOCK:
        try:
            items = feeds.fetch_who_don(limit=10)
            summary["fetched"] = len(items)
            summary["new"] = store.upsert_feed_items(items)
            store.set_setting("feed_last_fetch", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            store.set_setting("feed_last_error", "")
        except Exception as exc:
            summary["error"] = str(exc)[:200]
            store.set_setting("feed_last_error", summary["error"])
            return summary
        if analyze is None:
            analyze = store.get_setting("feed_auto_analyze", "1") == "1"
        if analyze and anthropic_client is not None:
            for item in store.feed_items(limit=FEED_AUTO_ANALYZE_CAP, status="new"):
                if analyze_feed_item(item, store).get("ok"):
                    summary["analyzed"] += 1
    return summary


@st.cache_resource
def start_feed_worker() -> Optional[threading.Thread]:
    """Start the background feed thread once per process."""
    if os.getenv("OHP_DISABLE_FEED_WORKER") == "1":
        return None
    store = get_store()

    def _loop():
        time.sleep(5)
        while True:
            try:
                sync_feed(store)
            except Exception:
                pass
            time.sleep(FEED_INTERVAL_S)

    t = threading.Thread(target=_loop, name="ohp-feed-worker", daemon=True)
    t.start()
    return t


# =====================================================================
# 8. STARTUP UNIT TESTS
# Three states: not-configured (None) / passed (True) / failed (False).
# Looks for tests/ next to this file, not in the launch directory.
# =====================================================================
APP_DIR = os.path.dirname(os.path.abspath(__file__))


def run_startup_tests():
    tests_dir = os.path.join(APP_DIR, 'tests')
    if not os.path.isdir(tests_dir):
        st.session_state['startup_tests_ok'] = None
        st.session_state['startup_test_output'] = (
            "No 'tests/' directory found next to the app — unit tests are "
            "not configured for this deployment.")
        return
    try:
        import unittest, io
        suite = unittest.TestLoader().discover(tests_dir, top_level_dir=APP_DIR)
        if suite.countTestCases() == 0:
            st.session_state['startup_tests_ok'] = None
            st.session_state['startup_test_output'] = "No test cases discovered in 'tests/'."
            return
        stream = io.StringIO()
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
        st.session_state['startup_tests_ok'] = result.wasSuccessful()
        st.session_state['startup_test_output'] = stream.getvalue()
    except Exception as e:
        st.session_state['startup_tests_ok'] = False
        st.session_state['startup_test_output'] = (
            f"Test runner error: {e}\n" + traceback.format_exc())


# =====================================================================
# 9. SERVICE STATUS PROBES
# Real, cheap requests (no user data sent), run in parallel, cached 5 min.
# Claude is verified with a free Models API lookup, which checks both the
# key and that ANTHROPIC_MODEL is available to it.
# =====================================================================
_STATUS_UA = {'User-Agent': 'OneHealthPlatform/1.0 (service status check)'}


def _probe_http(url: str, params: Optional[dict] = None) -> Tuple[str, str]:
    t0 = time.time()
    try:
        r = requests.get(url, params=params, headers=_STATUS_UA, timeout=6)
        ms = (time.time() - t0) * 1000
        if r.status_code < 400:
            return "ok", f"online · {ms:.0f} ms"
        return "warn", f"responded HTTP {r.status_code}"
    except requests.Timeout:
        return "warn", "timed out (>6 s)"
    except requests.RequestException:
        return "warn", "unreachable (network/DNS)"


def _probe_claude() -> Tuple[str, str]:
    if anthropic_client is None:
        return "warn", _claude_error
    try:
        anthropic_client.with_options(timeout=8, max_retries=0).models.retrieve(ANTHROPIC_MODEL)
        return "ok", f"connected · {ANTHROPIC_MODEL}"
    except AuthenticationError:
        return "warn", "API key rejected (401) — check ANTHROPIC_API_KEY"
    except PermissionDeniedError:
        return "warn", "API key lacks permission (403)"
    except NotFoundError:
        return "warn", f"model '{ANTHROPIC_MODEL}' not available to this key"
    except APIConnectionError:
        return "warn", "unreachable (network/DNS)"
    except APIStatusError as exc:
        return "warn", f"responded HTTP {exc.status_code}"
    except Exception as exc:
        return "warn", str(exc)[:120]


# (name, what the app uses it for, probe)
SERVICE_PROBES = [
    ("Claude API", "Epidemic analysis & article credibility", _probe_claude),
    ("WHO outbreak feed", "Disease Outbreak News, fetched automatically",
     lambda: _probe_http(feeds.WHO_DON_API, {"$top": "1"})),
    ("NCBI BLAST", "Fallback sequence search (nt / nr)",
     lambda: _probe_http("https://blast.ncbi.nlm.nih.gov/Blast.cgi", {"CMD": "Info"})),
    ("NCBI Entrez", "Gene-name lookup for BLAST hits",
     lambda: _probe_http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi", {"retmode": "json"})),
    ("Nominatim", "Geocoding outbreak locations for the map",
     lambda: _probe_http("https://nominatim.openstreetmap.org/status", {"format": "json"})),
]


@st.cache_data(ttl=300, show_spinner=False)
def service_status() -> Dict[str, Any]:
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(SERVICE_PROBES)) as ex:
        futures = [(name, purpose, ex.submit(fn)) for name, purpose, fn in SERVICE_PROBES]
        results = [(name, purpose, *f.result()) for name, purpose, f in futures]

    info = _amrfinder_info(refresh=True)   # a recheck should re-probe the tool too
    if info.get("ok"):
        via = " via WSL" if info.get("backend") == "wsl" else ""
        amr = ("ok", f"v{info['software']} · database {info['database']}{via}")
    else:
        amr = ("warn", f"{info.get('error', 'unavailable')} — using web BLAST instead")
    results.insert(0, ("AMRFinderPlus", "Validated resistance-gene detection (local)", *amr))

    try:
        store = get_store()
        results.append(("Database", "Saved analyses, audit trail & feed items",
                        "ok", f"{store.count_analyses()} analyses saved · SQLite"))
    except Exception as exc:
        results.append(("Database", "Saved analyses, audit trail & feed items", "warn", str(exc)[:120]))
    return {"results": results, "checked_at": datetime.now().strftime("%H:%M:%S")}


