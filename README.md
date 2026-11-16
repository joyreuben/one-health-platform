# 🧬 One Health Intelligence Platform

**Version:** 4.0.0 · **Stack:** Streamlit · NCBI AMRFinderPlus · Biopython (NCBI BLAST) · Anthropic Claude API · SQLite · Plotly

A research and surveillance web app that brings genomic and epidemic signals together under a One Health lens:

1. **Genomic AMR Surveillance.** Paste a DNA or protein sequence, or upload FASTA files. NCBI AMRFinderPlus finds the antimicrobial-resistance genes it carries, the app scores the evidence and produces an HL7 v2.5 lab report.
2. **Epidemic Early Warning.** Paste a health article or alert. It is screened for misinformation, then Claude extracts the disease, case count, trend and location, and the app rates the outbreak risk and produces an HL7 v2.5 alert.
3. **Outbreak Feed.** WHO Disease Outbreak News is fetched automatically and run through the same pipeline.
4. **Dashboard and Audit Trail.** Every analysis is saved to a database, plotted over time and on a world map, and can be reopened or exported.

> ⚕️ **Research and surveillance use only. Not for clinical diagnosis.** This is a decision-support prototype, not a validated diagnostic or an official public-health reporting system.

---

## Quick start

```bash
git clone <your-repo-url>
cd one-health-platform
python -m venv .venv
.venv\Scripts\activate            # macOS/Linux: source .venv/bin/activate
pip install -r code/requirements.txt
streamlit run code/app.py
```

Open http://localhost:8501. Always launch from the project root: Streamlit reads the theme from `.streamlit/config.toml` in the folder you start it from.

The app runs with nothing else configured, with two reduced capabilities:

| Not set up | What happens |
|---|---|
| Claude API key | Epidemic analysis uses a keyword fallback, labelled DEMO. Feed items are collected but not analyzed. |
| AMRFinderPlus | Gene detection falls back to web BLAST with keyword classification: 1–5 minutes per sequence and less reliable. |

The **Status** page shows live which services are working.

## Configuration

### Claude API key

Either set an environment variable before launching:

```powershell
$env:ANTHROPIC_API_KEY = "sk-ant-..."        # macOS/Linux: export ANTHROPIC_API_KEY="sk-ant-..."
```

or create `.streamlit/secrets.toml` (git-ignored):

```toml
ANTHROPIC_API_KEY = "sk-ant-..."
```

The model is set in one place, `ANTHROPIC_MODEL` in `code/engine.py`.

### AMRFinderPlus

AMRFinderPlus runs on Linux and macOS only. The app finds it in one of two ways:

- **Linux, macOS, containers:** `amrfinder` on the `PATH`.
- **Windows:** through WSL, from a conda environment inside your Linux distribution.

Install (inside Linux or WSL, with [Miniforge](https://github.com/conda-forge/miniforge)):

```bash
mamba create -y -n ohp-amrfinder -c conda-forge -c bioconda ncbi-amrfinderplus
conda activate ohp-amrfinder
amrfinder -u          # download the reference database; re-run to update it
```

Optional environment variables:

| Variable | Default | Purpose |
|---|---|---|
| `OHP_AMRFINDER_ENV` | `$HOME/miniforge3/envs/ohp-amrfinder` | Conda environment used through WSL |
| `OHP_WSL_DISTRO` | default distribution | WSL distribution name |
| `OHP_DATA_DIR` | `<project>/data` | Where the SQLite database is stored |
| `OHP_DISABLE_FEED_WORKER` | unset | Set to `1` to stop the background feed thread |

## Project layout

```
code/
  app.py          Pages and top navigation (Streamlit entry point)
  engine.py       Analysis logic: validation, AMRFinderPlus, BLAST, Claude, scoring, HL7, status probes
  storage.py      SQLite persistence (analyses, audit trail, feed items, settings)
  feeds.py        WHO Disease Outbreak News fetching and parsing
  ui_theme.py     Design: stylesheet, page themes, CSS/SVG animations
  tests/          Offline unit tests
.streamlit/config.toml   Dark theme
data/             Created at runtime: onehealth.db (git-ignored)
logs/             Created at runtime: Claude raw-response debug log (git-ignored)
```

## Pages

| Page | What it does |
|---|---|
| **Home** | Overview of the platform and links to each tool. |
| **Genomic** | Single sequence or FASTA batch. Shows every resistance gene found, drug class, confidence, and the HL7 report. Engine and BLAST threshold are under *Search settings*. |
| **Epidemic** | Paste an article or pick a demo. Credibility verdict, then risk analysis and HL7 alert. Fake or likely-fake articles are blocked from scoring. |
| **Outbreak Feed** | Latest WHO reports with status, analysis results and a manual *Fetch latest now* button. |
| **Dashboard** | Totals, confidence timeline, and a world map of geocoded reports. |
| **Audit Trail** | Every saved analysis; reopen any report's HL7 and details; CSV export; delete all. |
| **Status** | Live checks of AMRFinderPlus, Claude, WHO feed, NCBI, Nominatim and the database, plus unit-test results. |

## How gene detection works

**Primary: AMRFinderPlus.** The sequence is searched against NCBI's curated Reference Gene Catalog using per-gene cutoffs. Nucleotide input is searched in all six frames, so genes on either strand are found. A batch is sent in one call. Typical time is a few seconds.

Confidence for an AMRFinderPlus call:

| Component | Weight |
|---|---|
| Identity to the reference gene | 45% |
| Coverage of the reference gene | 30% |
| Detection method | 25% |

Method weights: exact or allele match 1.00 · point mutation 0.90 · BLAST above curated cutoff 0.85 · partial 0.60–0.65 · HMM only 0.50 · internal stop 0.35.

Drug classes are mapped to a concern level: carbapenem, colistin, vancomycin, methicillin and cephalosporin resistance are CRITICAL; other beta-lactam, tetracycline and glycopeptide resistance HIGH; everything else MODERATE. This mapping is a simple heuristic defined in `AMRFinderClient`, not a clinical standard.

**Fallback: web BLAST.** Used when AMRFinderPlus is unavailable, or when chosen in *Search settings*. Hits from NCBI's public BLAST service are classified by matching hit titles against a small list of gene-name patterns, then an Entrez lookup, then a sequence motif scan. This path is slow and far less reliable; treat its results as indicative only.

A negative AMRFinderPlus result means no match in the reference database. It does not prove susceptibility: organism-specific point mutations are not searched (that needs the organism to be specified) and novel genes are not covered.

## How epidemic analysis works

1. **Credibility screen.** A fast heuristic (satire sources, conspiracy patterns, sensational language, missing specifics) runs first. With Claude connected, Claude also assesses the article. Verdicts of FAKE or LIKELY_FAKE stop the analysis.
2. **Extraction.** Claude returns disease, location, cases, severity, trend, baseline and a 0–1 risk score as JSON.
3. **Dampening.** For UNVERIFIABLE articles the risk score is multiplied by the credibility score, with a floor that never raises the original score.
4. **Output.** Risk level, recommended action, HL7 alert, a saved audit entry and a map point.

Without Claude, step 2 is a keyword lookup and the result is labelled DEMO.

### Outbreak feed

A background thread fetches the ten newest WHO Disease Outbreak News items a few seconds after start-up and then every six hours while the app is running. Unseen items are stored. If Claude is connected and *Analyze new reports automatically* is on, up to five new items per check are analyzed; each analysis makes two Claude API calls billed to your key. The thread only runs while the app process is running; there is no separate scheduler.

## Data and privacy

- **Saved:** the finding, scores, HL7 message and a details record for each analysis; WHO feed items (public text); settings. For genomic results the sequence length and SHA-256 hash are saved, not the sequence.
- **Not saved:** pasted sequences and pasted article text.
- **Sent to third parties:** sequences to NCBI when web BLAST is used; article text to Anthropic when Claude is connected; place names to Nominatim (OpenStreetMap) when geocoding is on. AMRFinderPlus runs locally and sends nothing.
- **Debug log:** raw Claude responses are written to `logs/claude_raw.log`, capped at about 2 MB.
- There are **no user accounts**. Everyone who can open the app shares one database and can delete its contents.

## Tests

```bash
python -m unittest discover -s code/tests -t code
```

53 offline tests cover sequence validation, AMRFinderPlus output parsing and interpretation, the BLAST-title classifier, scoring, HL7 formatting and escaping, the credibility heuristic, Claude response parsing and retries, the epidemic pipeline, storage and feed parsing. They make no network, Claude or AMRFinderPlus calls. The app also runs them at start-up and shows the result on the Status page.

## Evaluation

Benchmarks and their results are in [`evaluation/`](evaluation/RESULTS.md). In brief: 98.9% of 3,212 ResFinder resistance genes were detected with 99.9% agreement on drug class and no calls on 7,530 negative-control sequences; the offline misinformation rules block 26.5% of fake posts while passing 97.4% of real ones; the Claude-based stages have not been evaluated yet. Read the limitations there before quoting any figure.

## Known limitations

- Detection has been benchmarked on reference gene sequences only, not on isolates with known resistance profiles or laboratory susceptibility results.
- Input is single sequences or FASTA. Raw reads (FASTQ) and assembly are not supported.
- Point-mutation resistance is not reported because the organism is not passed to AMRFinderPlus.
- The risk score is a model estimate from one article, not a statistical outbreak detection against a baseline.
- No authentication, roles or per-organisation data separation.
- HL7 messages are structurally valid v2.5 ORU^R01 but use local codes, not LOINC.
- One feed source (WHO). The feed only updates while the app is running.
- Streamlit is suitable for demos and pilots; a multi-user production deployment needs a proper backend.

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| Status shows AMRFinderPlus "database not downloaded" | Run `amrfinder -u` in the conda environment. |
| Status shows AMRFinderPlus "timed out" on Windows | WSL was starting up. Click *Recheck services*. |
| Claude "API key rejected (401)" | The key is wrong or revoked. |
| Epidemic results labelled DEMO | Claude is not connected or the call failed; the banner states which. |
| Plain light theme, no dark design | The app was started from a folder other than the project root. |
| Map is empty | Only epidemic reports that name a place are mapped; new places take a few seconds to geocode. |
