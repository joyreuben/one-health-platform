# 🧬 One Health Intelligence Platform

**Version:** 3.2.0-ONE-HEALTH-FIXED
**Stack:** Streamlit · Biopython (NCBI BLAST) · Anthropic Claude API · Plotly · Pandas/NumPy

A research/surveillance dashboard that fuses three workflows under a "One Health" lens (human, animal, and environmental health signals in one place):

1. **Genomic AMR Surveillance** — BLAST a DNA/protein sequence, classify any antimicrobial resistance (AMR) gene it carries, score confidence, and emit an HL7 v2.5 lab report.
2. **Epidemic Early Warning** — Feed in a health article/report, screen it for credibility (fake news / misinformation detection), then have Claude extract outbreak signals (disease, case counts, trend, risk score) and emit an HL7 v2.5 epidemic alert.
3. **Correlation Dashboard & Audit Trail** — Cross-view of every analysis run, plotted on a world map (geocoded via Nominatim) and logged to a running audit table.

> ⚕️ **Research & Surveillance Use Only — Not for Clinical Diagnosis.** This tool is a decision-support / prototyping aid, not a validated diagnostic or public-health reporting system.

---

## Table of Contents

- [Features](#features)
- [Architecture](#architecture)
- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration (secrets)](#configuration-secrets)
- [Running the App](#running-the-app)
- [How Each Tab Works](#how-each-tab-works)
  - [Genomic Surveillance](#1-genomic-surveillance)
  - [Epidemic Detection](#2-epidemic-detection)
  - [Correlation Dashboard](#3-correlation-dashboard)
  - [Audit Trail](#4-audit-trail)
- [AMR Gene Classification Pipeline](#amr-gene-classification-pipeline)
- [Confidence Scoring](#confidence-scoring)
- [HL7 Output](#hl7-output)
- [Caching & Performance](#caching--performance)
- [Geocoding](#geocoding)
- [Logging](#logging)
- [Unit Tests](#unit-tests)
- [Known Limitations](#known-limitations)
- [Troubleshooting](#troubleshooting)
- [Disclaimer](#disclaimer)

---

## Features

- 🧬 **Live NCBI BLAST** (blastn/blastx/blastp, auto-selected from sequence type) against `nt`/`nr`, with automatic cross-program fallback (nucleotide → protein → nucleotide) if the first search comes up empty or times out.
- 🔬 **Four-layer AMR gene classification**: ResFinder API hook (currently disabled) → curated regex marker library across all BLAST hit titles → live NCBI Entrez gene-name fetch → conserved k-mer sequence-motif scan → falls back to "novel variant, needs manual review" only if all four layers fail.
- 📁 **Batch FASTA upload** — analyze many sequences in one pass with a progress bar and a combined HL7 download.
- 🤖 **Claude-powered epidemic signal extraction** (disease, location, case count, severity, trend, baseline, risk score) with retry/backoff and a keyword-based offline fallback if the API is unavailable.
- 🕵️ **Fake-article / misinformation detector** — a fast rule-based heuristic (satire-site list, conspiracy-language patterns, sensationalism signals, credibility signals) combined with a stricter Claude-based verdict; **blocks** downstream risk scoring for content flagged FAKE/LIKELY_FAKE so no HL7 alert is generated from fabricated data.
- 📊 **Weighted confidence scoring** for AMR calls (BLAST identity, E-value, coverage, and classification-source reliability) and for epidemic risk (dampened by article credibility when a report is unverifiable).
- 🧾 **HL7 v2.5 ORU^R01 message generation** for both AMR findings and epidemic alerts, downloadable as `.hl7` files.
- 🌍 **World map view** of audit-log findings, geocoded via a shared/cached OpenStreetMap Nominatim lookup with retry/backoff for failed lookups.
- 📋 **Persistent audit trail** (in-memory, thread-safe, capped ring buffer) exportable to CSV.
- 🧪 **Startup self-test runner** that discovers and runs any `tests/` directory next to the app and surfaces pass/fail/not-configured status in the sidebar.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│ Streamlit UI (4 tabs)                                                │
│  ├─ Genomic Surveillance   ──> run_amr_pipeline()                    │
│  ├─ Epidemic Detection     ──> ArticleCredibilityChecker + EpidemicDetectionEngine │
│  ├─ Correlation Dashboard  ──> audit log + Plotly map/chart           │
│  └─ Audit Trail            ──> ThreadSafeAuditLog export              │
├─────────────────────────────────────────────────────────────────────┤
│ Core engines                                                         │
│  • SequenceValidator        — cleans/validates DNA or protein input  │
│  • AMRValidationEngine      — BLAST orchestration + 4-layer classify │
│  • ResFinderWebClient       — purpose-built AMR DB hook (disabled)   │
│  • ConfidenceScorer         — weighted confidence formulas           │
│  • HL7Generator             — ORU^R01 message builder                │
│  • EpidemicDetectionEngine  — Claude-based outbreak signal extraction│
│  • ArticleCredibilityChecker— heuristic + Claude fake-news detection │
├─────────────────────────────────────────────────────────────────────┤
│ External services                                                    │
│  • NCBI BLAST (qblast) + Entrez eSummary                              │
│  • Anthropic Claude API (model: claude-sonnet-4-6)                   │
│  • OpenStreetMap Nominatim (geocoding)                                │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Requirements

- Python 3.9+
- Packages: `streamlit`, `pandas`, `numpy`, `plotly`, `biopython`, `requests`, `anthropic`

Example `requirements.txt`:

```
streamlit
pandas
numpy
plotly
biopython
requests
anthropic
```

## Installation

```bash
git clone <your-repo-url>
cd <your-repo>
python -m venv venv
source venv/bin/activate      # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Configuration (secrets)

The app needs an Anthropic API key to enable Claude-powered epidemic analysis and article credibility checks. Without it, both features fall back to local heuristics automatically (the app still runs).

**Option A — Streamlit secrets file** (`.streamlit/secrets.toml`):

```toml
ANTHROPIC_API_KEY = "sk-ant-..."
```

**Option B — environment variable:**

```bash
export ANTHROPIC_API_KEY="sk-ant-..."
```

No API key is required for BLAST (NCBI's public `qblast` endpoint) or Nominatim geocoding, but both are rate-limited by their providers — be considerate of usage.

## Running the App

```bash
streamlit run app.py
```

Then open the local URL Streamlit prints (typically `http://localhost:8501`).

On first load, the app:
- Initializes a thread-safe in-memory audit log and result cache.
- Loads (or creates) a persistent `location_cache.json` geocode cache in the working directory.
- Discovers and runs any `tests/` directory (see [Unit Tests](#unit-tests)) and shows the result in the sidebar.
- Checks Claude API connectivity and shows a status line in the sidebar.

---

## How Each Tab Works

### 1. Genomic Surveillance

**Input:** a Sample ID, a raw or FASTA-formatted DNA/protein sequence (20–500,000 bp), and/or a batch of uploaded FASTA files.

**Pipeline** (`run_amr_pipeline`):
1. Validate/clean the sequence (`SequenceValidator`).
2. Detect protein vs. DNA vs. likely protein-coding DNA, and pick the right BLAST program (`blastp`/`blastx`/`blastn`) and database (`nr`/`nt`).
3. Run NCBI BLAST in a background thread with a 300s timeout; cache raw results by `(program, database, sequence hash, params)` so re-filtering with a different identity threshold doesn't re-query NCBI.
4. Filter hits by identity % and E-value; if nothing passes, automatically retry with a fallback BLAST program (protein↔nucleotide) before giving up.
5. Classify the top hit's gene through the [four-layer pipeline](#amr-gene-classification-pipeline).
6. Score confidence, derive a risk category, and generate an HL7 ORU^R01 report.

**Output:** gene name, confidence %, BLAST identity %, risk badge, a detail table, a raw-hit debug expander, and a downloadable HL7 file. Batch uploads get an aggregated results table and a combined HL7 download. Results persist across unrelated UI reruns (e.g., moving the identity-threshold slider) until you click **Reset** or run a new analysis.

### 2. Epidemic Detection

**Input:** a Report ID and either a pasted article/alert or one of five built-in demo articles (3 real-style, 2 deliberately fake, for testing the detector).

**Pipeline:**
1. **Credibility check first** (`ArticleCredibilityChecker`): a heuristic scan (satire-site list, implausible/conspiracy patterns, sensationalism, capitalization ratio, presence of sources/locations/dates/numbers) is combined with a stricter Claude-based verdict. If flagged `FAKE` or `LIKELY_FAKE`, the pipeline **stops here** — no risk score or HL7 alert is produced, by design.
2. If not blocked, `EpidemicDetectionEngine.analyze_article` sends the article to Claude with a strict JSON-only prompt, extracting location, disease, case count, severity, trend, baseline estimate, risk score, and a summary. On API failure it retries 3x with exponential backoff, then falls back to keyword-based analysis, with the resulting risk score dampened by the heuristic credibility score.
3. If the credibility verdict is `UNVERIFIABLE`, the risk score is further dampened (multiplied by credibility score, with a floor) and a warning banner explains the adjustment.
4. A risk category and recommended action tier (Critical/High/Moderate/Low) are shown, along with an HL7 epidemic alert message.

### 3. Correlation Dashboard

Aggregates everything logged in the audit trail:
- Summary metrics (total analyses, genomic vs. epidemic counts, "real" vs. demo results).
- Per-category tables.
- A confidence-over-time bar chart.
- A world map (Plotly `Scattergeo`) plotting geocoded finding locations, color-coded by type (Genomic/Epidemic), sized by confidence. Locations are resolved from a small built-in city lookup first, then a shared, disk-persisted Nominatim geocode cache (background thread, rate-limited, with retry/backoff for failures). The cache is viewable/editable/clearable from an expander, including a "retry failed lookups" action.

### 4. Audit Trail

A running, timestamped table of every analysis (genomic + epidemic), with a toggle to show only "real" (non-demo) results and a CSV export button.

---

## AMR Gene Classification Pipeline

Applied in order, most → least reliable, returning the first successful match:

| Layer | Method | Notes |
|---|---|---|
| 1. ResFinder | Purpose-built AMR reference DB via CGE web API | **Currently disabled** in this build (config flag, not a connectivity issue) |
| 2. Local pattern match | Regex marker library matched against all BLAST hit titles + sample ID | Covers NDM, KPC, MCR, CTX-M, mecA, vanA/B, tetA/B, ermB, aac, qnr, blaTEM/SHV |
| 3. Entrez fetch | Live NCBI `esummary` lookup of the top hit's accession, then re-matched against the same regex library | Only triggered if layer 2 finds nothing |
| 4. Sequence motif scan | Conserved k-mer signature search directly on the query sequence | Last resort when BLAST titles are uninformative (e.g., top hit is a whole-genome assembly) |
| — | `NOVEL_VARIANT` | Returned only if all four layers fail — flagged for manual curation |

---

## Confidence Scoring

**AMR confidence** = `0.40×identity + 0.30×(−log10(E-value)/300) + 0.15×coverage + 0.15×source_reliability`

Source reliability weights: ResFinder `1.00` · Local database match `0.90` · Entrez fetch `0.75` · Sequence motif scan `0.50` · Manual review `0.20`.

**Risk categories** (shared by AMR and epidemic scores): CRITICAL ≥0.85 · HIGH ≥0.70 · MODERATE ≥0.50 · LOW ≥0.30 · MINIMAL below that.

---

## HL7 Output

Both pipelines emit **HL7 v2.5 ORU^R01** messages (`MSH`/`PID`/`OBR`/`OBX` segments), field-escaped for `|^&\` characters, with an abnormal flag (`H`) set when confidence/risk crosses the high-concern threshold. Files are downloadable individually or (for genomic batches) combined into one `.hl7` file.

---

## Caching & Performance

- **BLAST result cache**: keyed by program + database + sequence hash + extra params (not by identity threshold), stored in `st.session_state.result_cache`, so adjusting the sidebar identity slider re-filters cached hits instantly instead of re-querying NCBI.
- **Geocode cache**: module-level (process-wide, shared across all sessions/users on the same server instance), persisted to `location_cache.json`, with attempt counts and exponential backoff (30s → 1h cap, giving up after 5 failed attempts) for locations that fail to resolve.
- **Rate limiting**: a minimum ~1.1s gap between Nominatim requests and ~2s between NCBI BLAST submissions are enforced via locks.

---

## Geocoding

Uses OpenStreetMap's Nominatim service with a descriptive `User-Agent` header, run in a background thread pool so the UI doesn't block. Because the cache is shared across all users of a single deployed instance, a location looked up once benefits everyone — this is intentional (a lat/lon for "Lagos" is the same fact for every user) and is documented in-code rather than isolated per-session.

---

## Logging

A best-effort local debug log of raw Claude API responses is written to `logs/claude_raw.log` (relative to the app's working directory — portable across OS/deployment, unlike a hardcoded local path). The file is size-capped at ~2MB and self-rotates by discarding its older half rather than growing unbounded. All logging failures are swallowed silently so they never crash the app.

---

## Unit Tests

On startup, the app looks for a `tests/` directory next to `app.py`. Three distinct states are shown in the sidebar:

- **Not configured** (ℹ️) — no `tests/` directory or no test cases found; this is normal if you haven't added tests.
- **Passed** (✅)
- **Failed** (⚠️) — expandable panel shows the failing test output.

To add tests, create a `tests/` package (with `__init__.py` and files matching `test_*.py`) alongside the app; `unittest`'s discovery is used automatically.

---

## Known Limitations

- ResFinder API integration is present in code but **disabled by config flag** in this build — AMR classification currently relies on layers 2–4 (local regex, Entrez, motif scan).
- NCBI BLAST searches typically take 1–5 minutes and are subject to NCBI's own rate limits and availability.
- Epidemic risk scoring and article credibility checks rely on Claude's JSON output; a keyword/heuristic fallback is used if the API is unavailable or misbehaves, which is necessarily less accurate.
- The audit log is in-memory only (capped at 500 records) and is **not persisted** across app restarts — only the geocode cache is written to disk.
- Not validated for clinical or regulatory use; HL7 output is illustrative, not a certified interface message.

---

## Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| Sidebar shows "⚠️ Claude API {error}" | `ANTHROPIC_API_KEY` missing/invalid, or `anthropic` package not installed — app still works via fallbacks |
| BLAST times out (>300s) | NCBI load/network issue — try a shorter sequence or retry later; the app auto-tries a cross-program fallback first |
| "No hits passed thresholds" | Lower the sidebar Identity Threshold, or the sequence may be too divergent/novel |
| Map shows "No geocoded locations available yet" | No audit-log entries have a resolvable `Location` yet, or geocoding is disabled in the sidebar |
| Sidebar shows "⚠️ Unit tests failed" | Expand "View failing tests" in the sidebar for the traceback |

---

## Disclaimer

This platform is provided for **research and surveillance prototyping only**. It is **not** a clinical diagnostic tool, does not replace validated laboratory AMR testing (e.g., ResFinder/CARD with curated databases) or official public-health outbreak reporting channels, and its AI-assisted outputs (Claude analyses, fake-article verdicts, gene classifications) should always be independently verified before informing any real-world health decision.
