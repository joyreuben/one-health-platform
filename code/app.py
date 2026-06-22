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

import requests
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

try:
    from anthropic import Anthropic
    _anthropic_available = True
except ImportError:
    _anthropic_available = False

__version__ = "3.1.2-ONE-HEALTH"

# =====================================================================
# ANTHROPIC CLIENT
# =====================================================================
anthropic_client = None
_claude_error = ""

if _anthropic_available:
    try:
        key = st.secrets.get("ANTHROPIC_API_KEY", os.getenv("ANTHROPIC_API_KEY", ""))
        if key:
            anthropic_client = Anthropic(api_key=key)
        else:
            _claude_error = "ANTHROPIC_API_KEY not set in secrets"
    except Exception as exc:
        _claude_error = str(exc)
else:
    _claude_error = "anthropic package not installed (add to requirements.txt)"


# =====================================================================
# 0. THREAD-SAFE AUDIT LOG
# =====================================================================
class ThreadSafeAuditLog:
    def __init__(self, max_records=500):
        self.log_queue = deque(maxlen=max_records)
        self.lock = threading.Lock()

    def log_entry(self, sample_id, analysis_type, finding,
                  confidence, data_source, is_demo=False):
        with self.lock:
            self.log_queue.append({
                "Timestamp":     datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "Type":          analysis_type,
                "Sample/Report": sample_id,
                "Finding":       finding,
                "Confidence":    f"{confidence * 100:.1f}%",
                "Source":        data_source,
                "Status":        "⚠️ DEMO" if is_demo else "✅ Real",
            })

    def get_logs(self):
        with self.lock:
            return pd.DataFrame(list(self.log_queue)) if self.log_queue else pd.DataFrame()


for key, default in [
    ("audit_log",         ThreadSafeAuditLog()),
    ("result_cache",      {}),
    ("last_request_time", 0),
]:
    if key not in st.session_state:
        st.session_state[key] = default

GEOCODE_LOCK = threading.Lock()
LAST_GEOCODE_TIME = 0.0
GEOCODE_EXECUTOR = None
LOCATION_CACHE = {}
LOCATION_CACHE_LOCK = threading.Lock()


def _shutdown_geocode_executor():
    global GEOCODE_EXECUTOR
    if GEOCODE_EXECUTOR is not None:
        try:
            GEOCODE_EXECUTOR.shutdown(wait=False)
        except Exception:
            pass
        GEOCODE_EXECUTOR = None


def get_geocode_executor():
    global GEOCODE_EXECUTOR
    if GEOCODE_EXECUTOR is None:
        GEOCODE_EXECUTOR = concurrent.futures.ThreadPoolExecutor(max_workers=2)
        atexit.register(_shutdown_geocode_executor)
    return GEOCODE_EXECUTOR


LOCATION_CACHE_PATH = os.path.join(os.getcwd(), "location_cache.json")
if 'location_cache' not in st.session_state:
    try:
        if os.path.exists(LOCATION_CACHE_PATH):
            with open(LOCATION_CACHE_PATH, 'r', encoding='utf-8') as f:
                LOCATION_CACHE = json.load(f)
        else:
            LOCATION_CACHE = {}
    except Exception:
        LOCATION_CACHE = {}
    st.session_state['location_cache'] = LOCATION_CACHE
else:
    LOCATION_CACHE = st.session_state['location_cache']


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
    global LAST_GEOCODE_TIME
    with GEOCODE_LOCK:
        elapsed = time.time() - LAST_GEOCODE_TIME
        if elapsed < 1.1:
            time.sleep(1.1 - elapsed)
        LAST_GEOCODE_TIME = time.time()

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
        LOCATION_CACHE[location] = res
        _save_location_cache()


def maybe_geocode_in_background(location: str):
    """Start background geocoding for a location if not already cached.
    Marks the location as pending (None) immediately to avoid duplicate work.
    Respects the user-configurable 'enable_geocoding' flag in session state."""
    if not location:
        return
    if not st.session_state.get('enable_geocoding', True):
        return
    with LOCATION_CACHE_LOCK:
        if location in LOCATION_CACHE:
            return
        LOCATION_CACHE[location] = None
        _save_location_cache()

    ex = get_geocode_executor()
    try:
        if ex:
            ex.submit(_geocode_and_store, location)
        else:
            _geocode_and_store(location)
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
        lines     = seq.strip().split('\n')
        seq_lines = lines[1:] if lines[0].startswith('>') else lines
        cleaned   = ''.join(seq_lines).upper().replace(' ', '').replace('\n', '')
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

    @classmethod
    def is_available(cls):
        return False

    @classmethod
    def analyze(cls, sequence, sample_id):
        return None


# =====================================================================
# 3. AMR VALIDATION ENGINE
# =====================================================================
class AMRValidationEngine:

    LOCAL_MARKERS = {
        "ndm":   (r"\bndm-?[0-9]?\b|metallo.beta.lactamase",
                  "Metallo-beta-lactamase", "Antibiotic inactivation", "Carbapenems (last-resort)", "CRITICAL"),
        "kpc":   (r"\bkpc-?[0-9]?\b",
                  "Carbapenemase", "Antibiotic inactivation", "Carbapenems (last-resort)", "CRITICAL"),
        "mcr":   (r"\bmcr-?[0-9]\b|colistin.resistance|mobile.colistin.resistance",
                  "Colistin Resistance", "Lipid A modification", "Colistin (last-resort)", "CRITICAL"),
        "ctx-m": (r"ctx-?m|esbl|extended.spectrum.beta.lactamase",
                  "Beta-lactamase (ESBL)", "Antibiotic inactivation", "Cephalosporins", "CRITICAL"),
        "mecA":  (r"\bmeca\b|penicillin.binding.protein.2a|pbp2a|pbp.2a",
                  "PBP Alteration", "Target alteration", "Methicillin", "CRITICAL"),
        "vanA":  (r"\bvan[-_ ]?a\b|\bvan[-_ ]?b\b|\bvana\b|\bvanb\b|\btn1546\b|d\.ala.{0,3}d\.ala.{0,4}ligase|"
                  r"d\.alanine.{0,4}d\.alanine.{0,4}ligase|"
                  r"d\.ala.{0,3}d\.lac.ligase|glycopeptide\.resistance|m97297",
                  "Glycopeptide Resistance", "Cell wall alteration", "Vancomycin", "CRITICAL"),
        # FIX: Broadened tetA pattern to catch BioPython-reformatted titles
        # - Added "tetracycline resistance protein" phrase match
        # - Added "tet [abcdefg]" with optional space
        # - \btet[a-eg-z]?\b avoids matching "text" while catching tetA-tetZ
        "tetA":  (r"\btet[a-eg-z]?\b"
                  r"|tetracycline.{0,20}resist"
                  r"|tetracycline.{0,10}efflux"
                  r"|tet\s*[abcdefg]\b",
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
        "ermB":  "ermB (rRNA methyltransferase)",
        "aacA":  "aac (Aminoglycoside acetyltransferase)",
        "qnr":   "qnr (Quinolone resistance protein)",
        "blaTEM":"blaTEM / blaSHV (Beta-lactamase)",
    }

    @staticmethod
    def _execute_cloud_blast(blast_program, database_target, sequence,
                              extra_params, sequence_hash, timeout_sec=300,
                              identity_threshold=None):
        param_parts = []
        if extra_params:
            param_parts = [f"{k}:{repr(extra_params[k])}" for k in sorted(extra_params)]
        cache_parts = [blast_program, database_target, sequence_hash, str(identity_threshold)] + param_parts
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

        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as ex:
            fut = ex.submit(_worker)
            try:
                xml_data = fut.result(timeout=timeout_sec)
            except concurrent.futures.TimeoutError:
                fut.cancel()
                raise TimeoutError(f"NCBI BLAST timeout (>{timeout_sec}s)")

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
        # Use 'nr' for protein or translated searches (blastp/blastx); use 'nt' for nucleotide searches (blastn)
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
                        extra_params, seq_hash, timeout_sec=300,
                        identity_threshold=identity_threshold)
                hits = cls._parse_blast_records(NCBIXML.parse(StringIO(xml_data)), seq_clean)
                score_thresh = identity_threshold if blast_program == 'blastn' else 35
                evalue_thresh = 1e-5 if blast_program == 'blastn' else 1e-2
                hits = [h for h in hits
                        if h["identity_percentage"] >= score_thresh
                        and h["e_value"] <= evalue_thresh]
                if hits:
                    st.success(f"✅ {len(hits)} matches passed filters")
                    return hits, "NCBI_CLOUD"

                # If protein-level search returned no hits, try nucleotide BLASTN as a fallback.
                if blast_program in ('blastx', 'blastp'):
                    st.info("ℹ️ No protein-level hits; retrying BLASTN vs nt as fallback")
                    try:
                        xml_data2 = cls._execute_cloud_blast(
                            'blastn', 'nt', seq_clean, extra_params, seq_hash + '_blastn', timeout_sec=300,
                            identity_threshold=identity_threshold)
                        hits2 = cls._parse_blast_records(NCBIXML.parse(StringIO(xml_data2)), seq_clean)
                        hits2 = [h for h in hits2
                                 if h["identity_percentage"] >= identity_threshold
                                 and h["e_value"] <= 1e-5]
                        if hits2:
                            st.success(f"✅ {len(hits2)} matches passed filters (blastn fallback)")
                            return hits2, "NCBI_CLOUD_BLASTN_FALLBACK"
                    except Exception:
                        st.warning("⚠️ BLASTN fallback failed or timed out")
 
                # If a DNA query has no strong blastn hits, try translated BLASTX for protein coding genes.
                if blast_program == 'blastn' and len(seq_clean) >= 90:
                    st.info("ℹ️ No strong nucleotide hits; trying BLASTX vs nr for protein-coding gene detection")
                    extra_params2 = {"expect": 1e-3}
                    if len(seq_clean) < 150:
                        extra_params2["word_size"] = 2
                    try:
                        xml_data3 = cls._execute_cloud_blast(
                            'blastx', 'nr', seq_clean, extra_params2, seq_hash + '_blastx', timeout_sec=300,
                            identity_threshold=identity_threshold)
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
                            'blastx', 'nr', seq_clean, extra_params2, seq_hash + '_blastx', timeout_sec=300,
                            identity_threshold=identity_threshold)
                        hits3 = cls._parse_blast_records(NCBIXML.parse(StringIO(xml_data3)), seq_clean)
                        hits3 = [h for h in hits3
                                 if h["identity_percentage"] >= 35
                                 and h["e_value"] <= 1e-2]
                        if hits3:
                            st.success(f"✅ {len(hits3)} matches passed filters (blastx timeout fallback)")
                            return hits3, "NCBI_CLOUD_BLASTX_TIMEOUT_FALLBACK"
                    except Exception:
                        st.warning("⚠️ BLASTX timeout fallback failed")
                elif blast_program in ('blastx', 'blastp'):
                    st.info("ℹ️ Protein-level search timed out; retrying BLASTN vs nt")
                    try:
                        xml_data2 = cls._execute_cloud_blast(
                            'blastn', 'nt', seq_clean, extra_params, seq_hash + '_blastn', timeout_sec=300,
                            identity_threshold=identity_threshold)
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

    # DNA-level k-mer signatures for genes whose BLAST hits are frequently
    # buried inside whole-genome assemblies with uninformative titles.
    # Each tuple: (unique_substrings_present_in_gene, min_matches_required)
    # Substrings chosen from conserved functional domains — long enough to be
    # specific, short enough to tolerate minor sequencing variants.
    SEQUENCE_MOTIFS = {
        "tetA":  ([
            "ATGAAATCTAACAATGCG",   # tetA start region (E.coli canonical)
            "TCGATGCTGTAGGCATAGGC",  # TM helix 1
            "GCGCGGTGGCGTGGTATGC",   # TM helix region
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
        """Run LOCAL_MARKERS patterns against a pre-normalised text string."""
        for marker, (pattern, cls_, mech, res, concern) in cls.LOCAL_MARKERS.items():
            if re.search(pattern, text, re.IGNORECASE):
                display = getattr(cls, "LOCAL_MARKER_DISPLAY", {}).get(marker, marker)
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
        """
        Fetch the NCBI nucleotide record summary for the accession and
        return the raw title/definition string so classify_gene() can
        pattern-match it.  Returns "" on any failure — never raises.

        Uses the public eutils REST endpoint (no API key required for
        low-volume use; respects NCBI's 3 req/s guideline via a 0.35s
        sleep).  Fetches esummary (JSON) which is much smaller than the
        full GenBank flat file and usually contains the gene product name
        in the 'title' field.
        """
        if not accession:
            return ""
        try:
            time.sleep(0.35)   # respect NCBI rate limit
            url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=nucleotide&id={accession}&retmode=json"
            )
            r = requests.get(url, timeout=15)
            r.raise_for_status()
            data  = r.json()
            uids  = list(data.get("result", {}).get("uids", []))
            if not uids:
                return ""
            rec   = data["result"][uids[0]]
            # 'title' is the GenBank DEFINITION line — most informative
            parts = [
                rec.get("title", ""),
                rec.get("subname", ""),
                " ".join(rec.get("extra", {}).get("gene", [])) if isinstance(rec.get("extra"), dict) else "",
            ]
            return " ".join(p for p in parts if p)
        except Exception:
            return ""

    @classmethod
    def _sequence_motif_scan(cls, sequence: str):
        """
        Last-resort classifier: scan the query sequence itself for
        known conserved k-mer signatures when BLAST titles are all
        uninformative.  Returns a result dict or None.
        """
        seq_upper = sequence.upper()
        for marker, (motifs, min_hits) in cls.SEQUENCE_MOTIFS.items():
            hits = sum(1 for m in motifs if m.upper() in seq_upper)
            if hits >= min_hits:
                _, cls_, mech, res, concern = cls.LOCAL_MARKERS.get(
                    marker,
                    (None, "Resistance Gene", "Unknown", "Undetermined", "HIGH"))
                return {"marker_found": marker, "class": cls_, "mechanism": mech,
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
        """
        # Layer 1: ResFinder
        rf = ResFinderWebClient.analyze(sequence, sample_id)
        if rf:
            return rf

        # Layer 2: Scan ALL hit titles (not just the best one).
        # Whole-genome assemblies dominate BLAST rankings for integrated
        # genes, but hits 5–20 often include named plasmid or gene entries.
        titles_to_scan = [hit_title] if not all_hits else (
            [h["title"] for h in all_hits[:20]] + [hit_title])
        # Also include sample_id so a header like ">ISO_tetA_sample" matches
        titles_to_scan.append(sample_id)

        for title in titles_to_scan:
            result = cls._match_title(cls._normalise(title))
            if result:
                result["source"] = "LOCAL_DATABASE"
                return result

        # Layer 3: Entrez eSummary for the best-hit accession.
        # Extracts the GenBank DEFINITION line which often names the gene
        # even when the BLAST title is truncated/generic.
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

        # Layer 4: Sequence k-mer motif scan.
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

                        # FIX: BioPython concatenates ALL accession titles in alignment.title
                        # separated by " >". Use robust splitting on common separators and
                        # fallback to accession when necessary to avoid noisy titles.
                        raw_title = getattr(alignment, "title", "") or getattr(alignment, "hit_def", "") or alignment.accession
                        parts = re.split(r'\s+>\s+|;|\||,', raw_title)
                        clean_title = next((p.strip() for p in parts if p and p.strip()), alignment.accession)
                        # Remove non-printable characters and collapse whitespace
                        clean_title = re.sub(r'\s+', ' ', re.sub(r'[^\x20-\x7E]', ' ', clean_title)).strip()

                        parsed.append({
                            "title":               clean_title,
                            "accession":           alignment.accession,
                            "e_value":             hsp.expect,
                            "score":               float(hsp.score),
                            "identity_percentage": id_pct,
                            "alignment_length":    hsp.align_length,
                            "query_coverage":      (hsp.align_length / qlen * 100) if qlen > 0 else 0,
                        })
        except Exception as exc:
            st.error(f"BLAST parse error: {exc}")
            return []
        # Sort hits so that entries with a recognisable resistance gene name
        # in the title come before generic whole-genome/chromosome assemblies.
        # Within each tier, rank by E-value (ascending).
        #
        # Whole-genome assemblies like "chromosome 1, complete sequence" or
        # "complete genome" often beat named gene entries on raw E-value
        # because the target sequence is much longer, giving a misleadingly
        # good score for a short alignment.  A named gene entry at e-15 is
        # far more informative than a chromosome hit at e-126.
        GENOME_NOISE = re.compile(
            r'chromosome\s+\d|complete\s+(genome|sequence)|whole.genome|'
            r'plasmid\s+\w+,\s+complete|contig\s+\d|scaffold\s+\d|'
            r'NODE_\d|genome\s+assembly',
            re.IGNORECASE)

        GENE_SIGNAL = re.compile(
            r'\b(tet[a-z]|van[a-z]|mec[a-z]|bla|ndm|kpc|mcr|ctx.?m|erm[a-z]|'
            r'qnr|aac|sul|resistance\s+gene|resistance\s+protein|'
            r'tetracycline|aminoglycoside|beta.lactam|carbapenem)\b',
            re.IGNORECASE)

        def _sort_key(h):
            title = h["title"]
            is_noise  = bool(GENOME_NOISE.search(title))
            has_gene  = bool(GENE_SIGNAL.search(title))
            # tier 0 = named gene, tier 1 = unknown, tier 2 = genome noise
            tier = 0 if has_gene else (2 if is_noise else 1)
            return (tier, h["e_value"])

        return sorted(parsed, key=_sort_key)


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
        if not anthropic_client:
            st.warning(f"⚠️ Claude API unavailable ({_claude_error}). Using keyword analysis.")
            return cls._demo_analysis(article_text)

        attempts = 3
        backoff = 1
        claude_log = "C:\\Users\\JOY\\Documents\\AI Ventures\\claude_raw.log"

        for attempt in range(1, attempts + 1):
            try:
                response = anthropic_client.messages.create(
                    model="claude-sonnet-4-6",
                    max_tokens=512,
                    messages=[{"role": "user",
                               "content": cls.PROMPT.format(article=article_text)}]
                )

                # Robustly extract text from possible response shapes
                text_candidates = []
                if hasattr(response, 'content') and response.content:
                    for item in response.content:
                        if isinstance(item, dict):
                            txt = item.get('text') or item.get('content')
                        else:
                            txt = getattr(item, 'text', None) or getattr(item, 'content', None)
                        if txt:
                            text_candidates.append(txt.strip())
                raw = max(text_candidates, key=len) if text_candidates else ''

                # Log raw response to session and file for debugging
                st.session_state.setdefault('claude_raw_responses', []).append({'attempt': attempt, 'raw': raw})
                try:
                    with open(claude_log, 'a', encoding='utf-8') as lf:
                        lf.write(f"[{datetime.now().isoformat()}] Attempt {attempt}\n")
                        lf.write((raw or '<EMPTY>') + "\n---\n")
                except Exception:
                    pass

                if not raw:
                    raise ValueError("Empty Claude response content")

                # Find the first balanced JSON object in the raw text
                start = raw.find('{')
                if start == -1:
                    raise ValueError("No JSON in response")
                end = None
                stack = 0
                for i in range(start, len(raw)):
                    if raw[i] == '{':
                        stack += 1
                    elif raw[i] == '}':
                        stack -= 1
                        if stack == 0:
                            end = i + 1
                            break
                if end is None:
                    raise ValueError("Unbalanced braces in response")
                json_text = raw[start:end]
                result = json.loads(json_text)
                return {
                    "location":  result.get("location", "Unknown"),
                    "disease":   result.get("disease",  "Unknown"),
                    "cases":     int(result.get("cases", 0)),
                    "severity":  result.get("severity", "unknown"),
                    "trend":     result.get("trend",    "unknown"),
                    "baseline":  result.get("baseline_estimate", "Unknown"),
                    "risk_score":float(np.clip(result.get("risk_score", 0.0), 0.0, 1.0)),
                    "summary":   result.get("summary",  ""),
                    "source":    "CLAUDE_API",
                }

            except Exception as exc:
                st.warning(f"⚠️ Claude API error (attempt {attempt}/{attempts}): {exc}")
                if attempt < attempts:
                    time.sleep(backoff)
                    backoff *= 2
                    continue

                # Final fallback after retries
                st.warning(f"⚠️ Claude API failed after {attempts} attempts. Using keyword fallback and heuristic merge.")
                try:
                    heuristic = ArticleCredibilityChecker._heuristic_check(article_text)
                except Exception:
                    heuristic = {"credibility_score": 0.2}
                demo = cls._demo_analysis(article_text)
                credibility_score = float(np.clip(heuristic.get('credibility_score', 0.0), 0.0, 1.0))
                floor = 0.5 if demo.get('risk_score', 0.0) >= 0.70 else 0.35
                adjusted = max(demo.get('risk_score', 0.0) * credibility_score, floor)
                demo['risk_score'] = float(np.clip(adjusted, 0.0, 1.0))
                demo['source'] = 'CLAUDE_API_ERROR+HEURISTIC'
                return demo

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

        # Strong signals of satire/misinformation
        for marker in cls.SATIRE_MARKERS:
            if marker in tl:
                flags.append(f"Known satire/misinformation source: '{marker}'")
                score -= 0.9

        # Implausible or conspiratorial language must heavily penalise
        for pattern in cls.IMPLAUSIBLE_PATTERNS:
            if re.search(pattern, tl, re.IGNORECASE):
                flags.append(f"Implausible or conspiracy claim detected")
                score -= 0.5

        # Detect sensational language (miracle cures, share-this, guaranteed)
        SENSATIONAL_PATTERNS = [
            r"\b(breaking|shocking|miracle|cure|secret|uncovered|guarantee|share before|must read)\b",
            r"100\s*%",
            r"\b(breaking:)\b"
        ]
        for pat in SENSATIONAL_PATTERNS:
            if re.search(pat, tl, re.IGNORECASE):
                flags.append("Sensational language detected — lowers credibility")
                score -= 0.35

        # If the article explicitly claims "100% credible" or similar, treat as highly suspicious
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

        # Ensure score stays within bounds
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
        attempts = 3
        backoff = 1
        claude_log = "C:\\Users\\JOY\\Documents\\AI Ventures\\claude_raw.log"
        for attempt in range(1, attempts + 1):
            try:
                response = anthropic_client.messages.create(
                    model="claude-sonnet-4-6",
                    max_tokens=512,
                    messages=[{"role": "user",
                               "content": cls.CREDIBILITY_PROMPT.format(article=text)}]
                )
                text_candidates = []
                if hasattr(response, 'content') and response.content:
                    for item in response.content:
                        if isinstance(item, dict):
                            txt = item.get('text') or item.get('content')
                        else:
                            txt = getattr(item, 'text', None) or getattr(item, 'content', None)
                        if txt:
                            text_candidates.append(txt.strip())
                raw = max(text_candidates, key=len) if text_candidates else ''

                # Log raw response
                st.session_state.setdefault('claude_raw_responses', []).append({'attempt': attempt, 'raw': raw})
                try:
                    with open(claude_log, 'a', encoding='utf-8') as lf:
                        lf.write(f"[{datetime.now().isoformat()}] API_CHECK attempt {attempt}\n")
                        lf.write((raw or '<EMPTY>') + "\n---\n")
                except Exception:
                    pass

                if not raw:
                    if attempt < attempts:
                        time.sleep(backoff)
                        backoff *= 2
                        continue
                    return None

                start = raw.find('{')
                if start == -1:
                    if attempt < attempts:
                        time.sleep(backoff)
                        backoff *= 2
                        continue
                    return None
                end = None
                stack = 0
                for i in range(start, len(raw)):
                    if raw[i] == '{':
                        stack += 1
                    elif raw[i] == '}':
                        stack -= 1
                        if stack == 0:
                            end = i + 1
                            break
                if end is None:
                    if attempt < attempts:
                        time.sleep(backoff)
                        backoff *= 2
                        continue
                    return None
                result = json.loads(raw[start:end])
                return {
                    "verdict":           result.get("verdict", "UNVERIFIABLE"),
                    "credibility_score": float(np.clip(
                        result.get("credibility_score", 0.5), 0.0, 1.0)),
                    "flags":             result.get("flags", []),
                    "reasoning":         result.get("reasoning", ""),
                    "source_cited":      bool(result.get("source_cited", False)),
                    "specific_details":  bool(result.get("specific_details", False)),
                }
            except Exception:
                if attempt < attempts:
                    time.sleep(backoff)
                    backoff *= 2
                    continue
                return None

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
            f'{result.get("reasoning","")}</div>'
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
    @staticmethod
    def amr_confidence(identity, e_value, coverage, reference_verified):
        i = np.clip(identity / 100.0, 0, 1)
        e = np.clip(-np.log10(np.clip(e_value, 1e-300, 1.0)) / 300.0, 0, 1)
        c = np.clip(coverage / 100.0, 0, 1)
        r = 1.0 if reference_verified else 0.5
        return float(np.clip(0.40*i + 0.30*e + 0.15*c + 0.15*r, 0, 1))

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
        return "< 1e-300" if (e == 0 or e < 1e-300) else f"{e:.2e}"

    @classmethod
    def amr_report(cls, sample_id, gene, confidence, identity, e_value, concern):
        ts  = datetime.now().strftime("%Y%m%d%H%M%S")
        mid = hashlib.md5(f"{sample_id}{ts}".encode()).hexdigest()[:13]
        # FIX: Lowered threshold to 0.70 and added INVESTIGATE to abnormal flag
        # so that any meaningful finding gets flagged H in HL7
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
# 7. FASTA PARSER
# =====================================================================
def parse_fasta_files(uploaded_files):
    sequences = {}
    for file in uploaded_files:
        content = file.getvalue().decode("utf-8")
        try:
            for record in SeqIO.parse(StringIO(content), "fasta"):
                sequences[record.id] = str(record.seq)
        except Exception:
            name = file.name.rsplit(".", 1)[0]
            raw  = "".join(l for l in content.splitlines() if not l.startswith(">"))
            if raw:
                sequences[name] = raw
    return sequences


# =====================================================================
# 8. STREAMLIT APP
# =====================================================================
st.set_page_config(page_title="🧬 One Health Intelligence Platform",
                   layout="wide", initial_sidebar_state="expanded")

if 'default_sample_id' not in st.session_state:
    st.session_state['default_sample_id'] = f"ISO_{datetime.now().strftime('%H%M%S')}"
if 'default_report_id' not in st.session_state:
    st.session_state['default_report_id'] = f"RPT_{datetime.now().strftime('%H%M%S')}"

st.markdown("""
<style>
.demo-alert{background:#2a2a00;color:#FFD700;padding:10px 16px;border-radius:6px;
            border-left:4px solid #FFD700;margin-bottom:12px;font-weight:500}
.status-ok  {background:#0D2A1A;color:#4CAF50;padding:6px 12px;border-radius:5px;
             border-left:3px solid #4CAF50;font-size:13px;margin:4px 0}
.status-warn{background:#2A2500;color:#FFC107;padding:6px 12px;border-radius:5px;
             border-left:3px solid #FFC107;font-size:13px;margin:4px 0}
</style>""", unsafe_allow_html=True)

for k, v in [("amr_engine",    AMRValidationEngine()),
              ("epi_engine",    EpidemicDetectionEngine()),
              ("credibility",   ArticleCredibilityChecker()),
              ("scorer",        ConfidenceScorer()),
              ("hl7",           HL7Generator())]:
    if k not in st.session_state:
        st.session_state[k] = v

# ---------------------------------------------------------------------
# Startup unit tests runner — run tests once per Streamlit session and
# surface results in the UI so regressions are immediately visible.
# ---------------------------------------------------------------------
def _run_startup_tests():
    try:
        import unittest, io
        loader = unittest.TestLoader()
        suite = loader.discover('tests')
        stream = io.StringIO()
        runner = unittest.TextTestRunner(stream=stream, verbosity=2)
        result = runner.run(suite)
        output = stream.getvalue()
        st.session_state['startup_tests_ok'] = result.wasSuccessful()
        st.session_state['startup_test_output'] = output
    except Exception as e:
        st.session_state['startup_tests_ok'] = False
        st.session_state['startup_test_output'] = (
            f"Test runner error: {e}\n" + traceback.format_exc())

# Run once per session
if 'startup_tests_ok' not in st.session_state:
    _run_startup_tests()

# Header
st.title("🧬 One Health Intelligence Platform")
st.markdown(f"**v{__version__}** | Genomic AMR Surveillance + Epidemic Early Warning + Fake Article Detection")
st.markdown("---")

# Sidebar
st.sidebar.header("⚙️ Settings")
tab_choice = st.sidebar.radio(
    "Mode", ["Genomic Surveillance", "Epidemic Detection",
             "Correlation Dashboard", "Audit Trail"])

st.sidebar.markdown("---")
st.sidebar.subheader("🌐 Service Status")
rf_ok = ResFinderWebClient.is_available()
st.sidebar.markdown(
    f'<div class="{"status-ok" if rf_ok else "status-warn"}">{"✅" if rf_ok else "⚠️"} ResFinder API '
    f'{"online" if rf_ok else "unreachable"}</div>', unsafe_allow_html=True)
st.sidebar.markdown(
    f'<div class="{"status-ok" if anthropic_client else "status-warn"}">{"✅" if anthropic_client else "⚠️"} Claude API '
    f'{"connected" if anthropic_client else _claude_error[:35]}</div>', unsafe_allow_html=True)

# Startup unit tests status
tests_ok = st.session_state.get('startup_tests_ok', True)
st.sidebar.markdown(
    f'<div class="{"status-ok" if tests_ok else "status-warn"}">{"✅" if tests_ok else "⚠️"} Unit tests '
    f'{"passed" if tests_ok else "failed"}</div>', unsafe_allow_html=True)
if not tests_ok:
    with st.sidebar.expander("View failing tests"):
        st.code(st.session_state.get('startup_test_output', ''), language='text')

st.sidebar.markdown("---")
identity_threshold = st.sidebar.slider("Identity Threshold (%)", 80, 100, 90)
# Enable optional background geocoding for unknown locations (Nominatim, cached)
enable_geocoding = st.sidebar.checkbox("Enable Nominatim geocoding (background, cached)", value=True)
if 'enable_geocoding' not in st.session_state:
    st.session_state['enable_geocoding'] = enable_geocoding
else:
    # sync UI checkbox with session state for persistent toggling
    st.session_state['enable_geocoding'] = enable_geocoding

with st.sidebar.expander("ℹ️ Confidence formula"):
    st.markdown("""
| Component | Weight |
|---|---|
| BLAST identity | 40% |
| E-value (log) | 30% |
| Query coverage | 15% |
| Reference verified | 15% |
""")

# =====================================================================
# TAB 1: GENOMIC SURVEILLANCE
# =====================================================================
if tab_choice == "Genomic Surveillance":
    st.subheader("🧬 Antimicrobial Resistance Screening")
    col1, col2 = st.columns([1, 2])

    with col1:
        st.markdown("### Input")
        sample_id = st.text_input("Sample ID", value=st.session_state['default_sample_id'])
        sequence  = st.text_area("DNA / Protein Sequence", height=200,
                                 placeholder=">header\nATGCGAC...\nor raw sequence")
        st.markdown("---")
        cr, cc = st.columns(2)
        trigger = cr.button("▶️ Analyze", type="primary", use_container_width=True)
        if cc.button("🔄 Reset", use_container_width=True):
            st.rerun()

    with col2:
        st.markdown("### Results")
        if trigger:
            if not sequence.strip():
                st.error("❌ Enter a sequence")
            else:
                try:
                    is_valid, err, seq_clean = SequenceValidator.validate_and_clean(sequence)
                    if not is_valid:
                        st.error(err)
                    else:
                        hits, source = st.session_state.amr_engine.run_blast(
                            seq_clean, use_cloud=True,
                            identity_threshold=identity_threshold)
                        is_demo = source == "DEMO_MODE"
                        if is_demo:
                            st.markdown('<div class="demo-alert">⚠️ DEMO MODE — Results randomly generated</div>',
                                        unsafe_allow_html=True)
                        if not hits:
                            st.warning("⚠️ No significant matches found")
                        else:
                            best = hits[0]

                            # Warn if the best hit is still a whole-genome
                            # assembly (can happen when ALL 50 hits are
                            # chromosomal assemblies with no named gene entries)
                            GENOME_NOISE_UI = re.compile(
                                r'chromosome\s+\d|complete\s+(genome|sequence)|'
                                r'whole.genome|plasmid\s+\w+,\s+complete',
                                re.IGNORECASE)
                            if GENOME_NOISE_UI.search(best["title"]):
                                st.warning(
                                    "⚠️ Top BLAST hit is a whole-genome assembly, not a "
                                    "named resistance gene entry. The sequence likely "
                                    "integrated into a host chromosome. Classification "
                                    "will use the sample ID and sequence context only.")

                            gene_info = st.session_state.amr_engine.classify_gene(
                                seq_clean, sample_id, best["title"],
                                all_hits=hits)
                            confidence = st.session_state.scorer.amr_confidence(
                                best["identity_percentage"], best["e_value"],
                                best.get("query_coverage", 50.0),
                                gene_info.get("marker_found","UNKNOWN") not in ("UNKNOWN","NOVEL_VARIANT"))
                            risk, emoji, _ = st.session_state.scorer.risk_category(confidence)
                            st.session_state.audit_log.log_entry(
                                sample_id, "GENOMIC",
                                gene_info.get("display_name", gene_info.get("marker_found","UNKNOWN")),
                                confidence, source, is_demo)

                            m1,m2,m3,m4 = st.columns(4)
                            display_gene = gene_info.get("display_name", gene_info.get("marker_found","UNKNOWN"))
                            m1.metric("Gene",       display_gene)
                            m2.metric("Confidence", f"{confidence*100:.1f}%")
                            m3.metric("Identity",   f"{best['identity_percentage']:.1f}%")
                            m4.metric("Risk",       f"{emoji} {risk}")
                            st.markdown("---")
                            st.dataframe(pd.DataFrame([{
                                "Sample":     sample_id,
                                "Gene Class": gene_info.get("class","Unknown"),
                                "Resistance": gene_info.get("resistance","Unknown"),
                                "Concern":    gene_info.get("clinical_concern","Unknown"),
                                "Identity %": f"{best['identity_percentage']:.2f}",
                                "E-value":    f"{best['e_value']:.2e}",
                                "Coverage %": f"{best.get('query_coverage',0):.1f}",
                                "Source":     gene_info.get("source","Unknown"),
                            }]), use_container_width=True, hide_index=True)

                            # Debug expander
                            with st.expander("🔍 Debug: Raw BLAST hit info"):
                                st.write(f"**Raw title (hit #1):** `{best['title']}`")
                                st.write(f"**Accession:** `{best['accession']}`")
                                st.write(f"**Classification source:** `{gene_info.get('source','?')}`")
                                st.write(f"**Total hits scanned:** `{len(hits)}`")
                                if len(hits) > 1:
                                    st.write("**Top 5 hit titles:**")
                                    for i, h in enumerate(hits[:5], 1):
                                        st.write(f"  {i}. `{h['title'][:100]}`")

                            st.markdown("---")
                            st.subheader("📋 HL7 ORU^R01 (v2.5)")
                            hl7 = st.session_state.hl7.amr_report(
                                sample_id, display_gene,
                                confidence, best["identity_percentage"],
                                best["e_value"], gene_info.get("clinical_concern","UNKNOWN"))
                            st.code(hl7, language="text")
                            st.download_button("📥 Download HL7", data=hl7,
                                               file_name=f"{sample_id}_AMR.hl7",
                                               mime="text/plain")
                except Exception as exc:
                    st.error(f"❌ Error: {exc}")
        else:
            st.info("👈 Enter sequence and click Analyze")


# =====================================================================
# TAB 2: EPIDEMIC DETECTION
# =====================================================================
elif tab_choice == "Epidemic Detection":
    st.subheader("🦠 Outbreak Early Warning System")
    col1, col2 = st.columns([1, 2])

    with col1:
        st.markdown("### Input")
        report_id  = st.text_input("Report ID", value=st.session_state['default_report_id'])
        input_mode = st.radio("Input Type", ["Paste Article", "Use Demo"])

        if input_mode == "Paste Article":
            article = st.text_area("Health article / alert text", height=200,
                                   placeholder="Paste a real health article, WHO alert, or hospital report...")
        else:
            demo_map = {
                "COVID-19 Outbreak (Real)": "covid_2020",
                "Mpox Cluster (Real)":      "mpox_2022",
                "Routine Report (Real)":    "normal",
                "5G Conspiracy (Fake)":     "fake_1",
                "Miracle Cure Hoax (Fake)": "fake_2",
            }
            demo_key = st.selectbox("Select Demo Article", list(demo_map.keys()))
            article  = EpidemicDetectionEngine.DEMO_ARTICLES[demo_map[demo_key]]
            st.text_area("Preview", value=article, height=120, disabled=True)

        st.markdown("---")
        cr, cc = st.columns(2)
        trigger = cr.button("▶️ Analyze", type="primary", use_container_width=True)
        if cc.button("🔄 Reset", use_container_width=True):
            st.rerun()

    with col2:
        st.markdown("### Results")
        if trigger:
            if not article.strip():
                st.error("❌ Enter article text")
            else:
                try:
                    st.markdown("#### 🔍 Article Credibility Check")
                    with st.spinner("Checking article credibility..."):
                        cred_result = ArticleCredibilityChecker.check(article)
                    blocked = ArticleCredibilityChecker.render_verdict(cred_result)

                    if not blocked:
                        st.markdown("---")
                        st.markdown("#### 📊 Epidemic Risk Analysis")

                        if cred_result["verdict"] == "UNVERIFIABLE":
                            st.warning(
                                "⚠️ This article could not be fully verified. "
                                "Treat the risk score below with caution and "
                                "cross-reference with official sources before acting.")

                        analysis  = st.session_state.epi_engine.analyze_article(article)
                        is_demo   = analysis["source"] in ("DEMO_MODE", "KEYWORD_FALLBACK")
                        risk_score = analysis["risk_score"]

                        if cred_result["verdict"] == "UNVERIFIABLE":
                            credibility_score = float(np.clip(
                                cred_result.get("credibility_score", 0.5), 0.0, 1.0))
                            floor = 0.5 if risk_score >= 0.70 else 0.35
                            dampened = risk_score * credibility_score
                            risk_score = float(np.clip(max(dampened, floor), 0.0, 1.0))
                            st.info(
                                f"⚠️ Risk score dampened from {analysis['risk_score']*100:.0f}% → "
                                f"{risk_score*100:.0f}% due to low article credibility. "
                                "Treat this result as provisional until verified.")

                        risk, emoji, _ = st.session_state.scorer.risk_category(risk_score)

                        if is_demo:
                            st.markdown(
                                f'<div class="demo-alert">⚠️ {analysis["source"]} — '
                                f'{"Claude API not connected" if not anthropic_client else "keyword-based fallback"}'
                                f'</div>', unsafe_allow_html=True)

                        st.session_state.audit_log.log_entry(
                            report_id, "EPIDEMIC", analysis["disease"],
                            risk_score, analysis["source"], is_demo)

                        m1,m2,m3,m4 = st.columns(4)
                        m1.metric("Disease",    analysis["disease"])
                        m2.metric("Risk Score", f"{risk_score*100:.1f}%")
                        m3.metric("Trend",      analysis["trend"].upper())
                        m4.metric("Severity",   f"{emoji} {risk}")
                        st.markdown("---")

                        st.dataframe(pd.DataFrame([{
                            "Report":      report_id,
                            "Location":    analysis["location"],
                            "Cases":       analysis["cases"],
                            "Baseline":    analysis["baseline"],
                            "Severity":    analysis["severity"],
                            "Credibility": f"{cred_result['credibility_score']*100:.0f}%",
                            "Summary":     analysis["summary"][:60] + "...",
                        }]), use_container_width=True, hide_index=True)

                        st.markdown("### 📋 Recommended Actions")
                        if risk_score >= 0.85:
                            st.error("🔴 CRITICAL → Activate emergency response, notify public health authorities immediately")
                        elif risk_score >= 0.70:
                            st.warning("🟠 HIGH → Increase surveillance, prepare response team, issue health alert")
                        elif risk_score >= 0.50:
                            st.info("🟡 MODERATE → Monitor closely, prepare communication plan")
                        else:
                            st.success("🟢 LOW → Routine monitoring, standard protocols")

                        st.markdown("---")
                        st.subheader("📋 HL7 Epidemic Alert (v2.5)")
                        hl7 = st.session_state.hl7.epidemic_alert(
                            report_id, analysis["disease"], risk_score,
                            analysis["cases"], analysis["trend"], analysis["location"])
                        st.code(hl7, language="text")
                        st.download_button("📥 Download HL7 Alert", data=hl7,
                                           file_name=f"{report_id}_epidemic_alert.hl7",
                                           mime="text/plain")
                except Exception as exc:
                    st.error(f"❌ Error: {exc}")
        else:
            st.info("👈 Choose article and click Analyze")


# =====================================================================
# TAB 3: CORRELATION DASHBOARD
# =====================================================================
elif tab_choice == "Correlation Dashboard":
    st.subheader("📍 One Health Correlation View")
    audit_df = st.session_state.audit_log.get_logs()

    if audit_df.empty:
        st.info("Run analyses in Genomic Surveillance and Epidemic Detection tabs first.")
    else:
        genomic_df  = audit_df[audit_df["Type"] == "GENOMIC"]
        epidemic_df = audit_df[audit_df["Type"] == "EPIDEMIC"]

        c1,c2,c3,c4 = st.columns(4)
        c1.metric("Total Analyses",   len(audit_df))
        c2.metric("Genomic Findings", len(genomic_df))
        c3.metric("Epidemic Signals", len(epidemic_df))
        c4.metric("Real Results",     len(audit_df[audit_df["Status"] == "✅ Real"]))
        st.markdown("---")

        if not genomic_df.empty:
            st.markdown("#### 🧬 Genomic Findings")
            st.dataframe(genomic_df[["Timestamp","Sample/Report","Finding","Confidence","Source","Status"]],
                         use_container_width=True, hide_index=True)
        if not epidemic_df.empty:
            st.markdown("#### 🦠 Epidemic Signals")
            st.dataframe(epidemic_df[["Timestamp","Sample/Report","Finding","Confidence","Source","Status"]],
                         use_container_width=True, hide_index=True)

        if len(audit_df) >= 2:
            st.markdown("#### 📊 Confidence Timeline")
            chart_df = audit_df.copy()
            chart_df["Confidence_num"] = chart_df["Confidence"].str.replace("%","").astype(float)
            fig = go.Figure(go.Bar(
                x=chart_df["Timestamp"], y=chart_df["Confidence_num"],
                marker_color=["#FF4B4B" if float(c.replace("%","")) >= 85
                              else "#FF9800" if float(c.replace("%","")) >= 70
                              else "#FFC107" if float(c.replace("%","")) >= 50
                              else "#4CAF50"
                              for c in chart_df["Confidence"]],
                text=chart_df["Finding"], textposition="outside"))
            fig.update_layout(yaxis_title="Confidence %", height=350, margin=dict(t=20,b=40))
            st.plotly_chart(fig, use_container_width=True)

    st.markdown("---")
    st.markdown("#### 🌍 Global Reference Map")
    st.markdown("Map of recent findings (dynamic). Hover markers for details. Locations are resolved from the audit log using a cached lookup.")

    # Ensure a persistent small cache for location -> (lon, lat)
    loc_cache = st.session_state.setdefault('location_cache', {})
    # Known quick mapping (fallback, editable)
    _known_locations = {
        'new york': {'lon': -74.0060, 'lat': 40.7128},
        'ny': {'lon': -74.0060, 'lat': 40.7128},
        'lagos': {'lon': 3.3792, 'lat': 6.5244},
        'nairobi': {'lon': 36.8219, 'lat': -1.2921},
        'tokyo': {'lon': 139.6917, 'lat': 35.6895},
        'johannesburg': {'lon': 28.0473, 'lat': -26.2041},
        'london': {'lon': -0.1278, 'lat': 51.5074},
        'beijing': {'lon': 116.4074, 'lat': 39.9042},
        'nigeria': {'lon': 8.6753, 'lat': 9.0820},
    }

    audit_df = st.session_state.audit_log.get_logs()
    points_by_type = {}

    if not audit_df.empty:
        for _, row in audit_df.iterrows():
            loc = row.get('Location', '') if isinstance(row.get('Location', ''), str) else ''
            if not loc:
                continue
            # Try cache first
            if loc in loc_cache:
                latlon = loc_cache[loc]
            else:
                latlon = None
                lk = loc.lower()
                for k, v in _known_locations.items():
                    if k in lk:
                        latlon = v
                        break
                # store even if None to avoid repeated work
                loc_cache[loc] = latlon
                # If still unknown, schedule background geocoding (non-blocking)
                if latlon is None:
                    try:
                        maybe_geocode_in_background(loc)
                    except Exception:
                        # best-effort; silently ignore geocode scheduling failures
                        pass
            if not latlon:
                continue

            t = row.get('Type', 'UNKNOWN')
            conf = row.get('Confidence', '0%')
            try:
                confv = float(str(conf).replace('%', ''))
            except Exception:
                confv = 50.0
            size = 8 + (confv / 100.0) * 24
            hover = f"{row.get('Finding','')} · {row.get('Sample/Report','')}<br>Conf: {confv:.1f}% · Source: {row.get('Source','')}"
            points_by_type.setdefault(t, []).append((latlon['lon'], latlon['lat'], hover, size))

    fig2 = go.Figure()
    color_map = {'GENOMIC': '#1976D2', 'EPIDEMIC': '#D32F2F', 'UNKNOWN': '#9E9E9E'}

    if not points_by_type:
        st.info("No geocoded locations available yet. Run analyses or add entries to the location cache.")
    else:
        for t, pts in points_by_type.items():
            lons = [p[0] for p in pts]
            lats = [p[1] for p in pts]
            hs = [p[2] for p in pts]
            sizes = [p[3] for p in pts]
            fig2.add_trace(go.Scattergeo(
                lon=lons,
                lat=lats,
                mode='markers',
                marker=dict(size=sizes, color=color_map.get(t, '#616161'), line=dict(width=0.5, color='black')),
                name=f"{t} ({len(pts)})",
                hoverinfo='text',
                hovertext=hs,
            ))

        fig2.update_layout(
            geo=dict(projection_type='natural earth', showland=True, landcolor='lightgray', showocean=True, oceancolor='lightblue'),
            height=450, margin=dict(t=10, b=10), legend=dict(title='Type (count)')
        )
        st.plotly_chart(fig2, use_container_width=True)

    with st.expander("Location cache (editable)"):
        st.write(loc_cache)
        if st.button("Clear location cache"):
            with LOCATION_CACHE_LOCK:
                LOCATION_CACHE.clear()
                st.session_state['location_cache'] = LOCATION_CACHE
                _save_location_cache()
            st.success("Location cache cleared (disk cache updated).")


# =====================================================================
# TAB 4: AUDIT TRAIL
# =====================================================================
elif tab_choice == "Audit Trail":
    st.subheader("📋 Analysis History")
    audit_df = st.session_state.audit_log.get_logs()

    if not audit_df.empty:
        c1,c2,c3 = st.columns(3)
        c1.metric("Total Analyses", len(audit_df))
        c2.metric("Real Results",   len(audit_df[audit_df["Status"] == "✅ Real"]))
        c3.metric("Demo Results",   len(audit_df[audit_df["Status"] == "⚠️ DEMO"]))
        st.markdown("---")
        real_only = st.checkbox("Show real results only", value=False)
        display_df = audit_df[audit_df["Status"] == "✅ Real"] if real_only else audit_df
        st.dataframe(display_df.iloc[::-1], use_container_width=True, hide_index=True)
        st.download_button("📥 Export Audit CSV",
            data=display_df.to_csv(index=False).encode(),
            file_name=f"audit_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
            mime="text/csv")
    else:
        st.info("No analyses yet")

# Footer
st.markdown("---")
st.markdown(
    f'<div style="text-align:center;color:#666;font-size:11px;">'
    f'<p>🧬 One Health Intelligence Platform v{__version__}</p>'
    f'<p>NCBI BLAST nt/nr DB · Curated Local Marker Matching · Claude API · HL7 v2.5 · Fake Article Detection</p>'
    f'<p>⚕️ Research &amp; Surveillance Use Only — Not for Clinical Diagnosis</p>'
    f'</div>', unsafe_allow_html=True)