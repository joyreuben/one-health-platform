# Evaluation results

One Health Intelligence Platform v4.0.0 · sections 1–3 run on 2026-11-15 (seed 20261115) · section 4 run 2026-11-17 to 2026-11-19 (seed 20261117)

This document reports what was measured, how, and what the numbers do and do not show. All intervals are 95% Wilson score intervals. Raw outputs are in `evaluation/results/`.

## Summary

| Component | Test | Result |
|---|---|---|
| Resistance-gene detection | 3,212 ResFinder genes | 98.9% detected (3,177; CI 98.5–99.2) |
| | Drug class agrees with ResFinder | 99.9% of detected genes (3,174 of 3,177; CI 99.7–100.0) |
| | 4,318 *E. coli* K-12 coding sequences | 0 calls (CI 0–0.09%) |
| | 3,212 shuffled sequences | 0 calls (CI 0–0.12%) |
| Previous keyword classifier | Same 3,212 genes, given the true gene name | 30.0% recognised (965) |
| Previous motif scan | Same 3,212 genes | 5.5% recognised (178) |
| Misinformation screen (offline rules) | 10,700 labelled COVID-19 posts | Blocks 26.5% of fake posts; passes 97.4% of real posts |
| | 149 WHO outbreak reports | 0 wrongly blocked (CI 0–2.5%) |
| Outbreak extraction, keyword fallback | 149 WHO outbreak reports | Disease correct in 18.1% (27); location never extracted |
| Outbreak extraction, Claude | | **Not evaluated: no API key was available** |
| Genotype → phenotype, clinical isolates | 99 public *K. pneumoniae* genomes, 1,006 isolate–drug lab results, 15 drugs | Categorical agreement 89.2% (CI 87.1–90.9); sensitivity for resistance 95.0% (CI 93.1–96.4); specificity 77.0% (CI 72.1–81.2) |
| | One isolate per clonal group (95 SNP clusters) | Agreement 89.3%; essentially unchanged |

## 1. Resistance-gene detection

### Method

Every sequence was passed through the platform's own batch pipeline (`engine.run_amr_batch`, AMRFinderPlus mode), so the measurement covers input validation, the AMRFinderPlus call, output parsing and interpretation.

- **Tool:** NCBI AMRFinderPlus 4.2.7, reference database 2026-10-29.1, default ("core") scope, nucleotide input, no organism specified.
- **Positive set:** all 3,212 gene entries in the ResFinder database v2.7.0 (DTU). ResFinder is curated by a different group from the NCBI Reference Gene Catalog that AMRFinderPlus searches.
- **Agreement on drug class:** the ResFinder class for a gene (from `phenotypes.txt`) was mapped to AMRFinderPlus class names. A detected gene agrees if any expected class appears among the classes or subclasses AMRFinderPlus reported. The mapping is in `amr_benchmark.py`.
- **Negative sets:** (a) all 4,318 coding sequences of *E. coli* K-12 MG1655 (NC_000913.3) that have a sequence; (b) a base-shuffled copy of every positive, preserving length and composition.
- **Robustness:** 300 detected genes sampled at random, then altered by random nucleotide substitution (1–20% of positions) or truncated from the 3' end.

### Detection of known genes

| | n | Detected | 95% CI |
|---|---|---|---|
| All ResFinder genes | 3,212 | 3,177 (98.9%) | 98.5–99.2 |
| Detected and drug class agrees | 3,212 | 3,174 (98.8%) | 98.4–99.1 |

How the 3,177 detected genes were matched: 2,901 (91.3%) were exact protein matches to a reference (methods ALLELEX and EXACTX), 265 passed the curated BLAST cutoff (BLASTX), and 11 were partial matches. Platform confidence scores: median 1.00, mean 0.995, minimum 0.769. Assigned concern level: 1,555 critical, 762 high, 860 moderate. In 22 entries more than one gene was reported, as expected for multi-gene entries such as *vanHAX*.

By ResFinder drug-class file:

| Class | n | Detected | Class agrees (of detected) |
|---|---|---|---|
| Beta-lactam | 2,031 | 2,028 (99.9%) | 2,028 (100%) |
| Aminoglycoside | 264 | 258 (97.7%) | 258 (100%) |
| Macrolide / lincosamide / streptogramin | 174 | 164 (94.3%) | 162 (98.8%) |
| Tetracycline | 162 | 158 (97.5%) | 158 (100%) |
| Quinolone | 125 | 125 (100%) | 125 (100%) |
| Colistin | 109 | 106 (97.2%) | 106 (100%) |
| Trimethoprim | 109 | 109 (100%) | 109 (100%) |
| Phenicol | 66 | 61 (92.4%) | 61 (100%) |
| Sulphonamide | 54 | 54 (100%) | 54 (100%) |
| Fosfomycin | 40 | 40 (100%) | 40 (100%) |
| Glycopeptide | 24 | 24 (100%) | 23 (95.8%) |
| Oxazolidinone | 22 | 22 (100%) | 22 (100%) |
| Nitroimidazole | 14 | 12 (85.7%) | 12 (100%) |
| Rifampicin | 10 | 10 (100%) | 10 (100%) |
| Fusidic acid | 3 | 3 (100%) | 3 (100%) |
| Pseudomonic acid | 3 | 3 (100%) | 3 (100%) |
| Misc (ionophores) | 2 | 0 | n/a |

Classes with fewer than about 30 genes have wide intervals; see `amr_benchmark.json`.

### The 35 genes not detected

Re-running these 35 with AMRFinderPlus's wider `--plus` scope finds 6 of them: *mcr-9.1*, *mcr-9.2*, *mcr-9.3*, *aac(3)-Ih*, *narA* and *narB*. NCBI places these outside the core set. The other 29 are not found in either scope. They are mostly genes from antibiotic-producing or environmental organisms (for example *grmA*, *kgmB*, *ole(C)*, *otr(C)*, *lmr(A)*), plus a few others such as *penA*, *blaACC-3* and *erm(41)*. These misses reflect differences in what the two databases include, not sequences that were present in the reference and missed. The full list is in `amr_missed_genes_plus_scope.json`.

### The 3 class disagreements

| Gene | ResFinder | AMRFinderPlus |
|---|---|---|
| *bleO* | Glycopeptide | Bleomycin |
| *vga(B)* | Streptogramin A, B | Lincosamide |
| *vga(D)* | Streptogramin A | Lincosamide |

These are differences in how the two databases label the same gene.

### Negative controls

| Set | n | Sequences with any call | 95% CI |
|---|---|---|---|
| *E. coli* K-12 MG1655 coding sequences | 4,318 | 0 | 0–0.09% |
| Shuffled copies of the positives | 3,212 | 0 | 0–0.12% |

The K-12 result depends on scope. K-12 carries intrinsic chromosomal genes (*ampC*/*blaEC*, *acrF*, *emrD*, *mdtM*) that AMRFinderPlus reports only with `--plus`; this was confirmed directly. In its default configuration the platform reports acquired resistance genes and not these intrinsic ones.

### Robustness to mutation and truncation (n = 300)

| Random nucleotide substitutions | Detected | Mean protein identity of detected | Mean confidence |
|---|---|---|---|
| 1% | 300 (100%) | 97.7% | 0.92 |
| 2% | 298 (99.3%) | 95.5% | 0.88 |
| 5% | 51 (17.0%) | 90.5% | 0.84 |
| 10% | 0 | n/a | n/a |
| 15% | 0 | n/a | n/a |
| 20% | 0 | n/a | n/a |

| Fraction of gene kept | Detected | Mean confidence |
|---|---|---|
| 90% | 300 (100%) | 0.93 |
| 75% | 300 (100%) | 0.84 |
| 50% | 294 (98.0%) | 0.76 |
| 25% | 2 (0.7%) | 0.81 |

Drug class agreed for every detected mutated gene. Detection stops once protein identity falls to about 90%, which is consistent with AMRFinderPlus's curated identity cutoffs, and once less than half the gene is present. Confidence falls as identity and coverage fall. Random substitution is a harsher model than natural divergence: it changes amino acids far more often than evolution does and introduces stop codons. The table shows where the cutoff lies; it does not estimate sensitivity to naturally occurring variants.

### Comparison with the platform's previous classifier

Before AMRFinderPlus was integrated, the platform classified genes by matching keywords in BLAST hit titles, with a short sequence-motif scan as a last resort. Evaluated offline on the same 3,212 genes:

| Previous method | Recognised | 95% CI |
|---|---|---|
| Keyword matcher, given the true gene name | 965 (30.0%) | 28.5–31.7 |
| Sequence-motif scan | 178 (5.5%) | 4.8–6.4 |

Giving the keyword matcher the true gene name is a best case; in use it saw whatever title BLAST returned. Neither method made a call on the K-12 set. The web BLAST step itself was not re-run, because it takes minutes per sequence.

### Speed

10,742 sequences (positives and both negative sets) took 443 seconds of pipeline time, about 24 sequences per second on an 8-thread laptop through WSL. Speed varies widely by gene: beta-lactamase alleles, which each match thousands of near-identical references, take about 0.3 seconds each, while most other sequences take a few milliseconds.

## 2. Misinformation screen (offline rules)

### Method

The platform's rule-based credibility check (`ArticleCredibilityChecker._heuristic_check`) was run on each text. A text counts as blocked when the verdict is FAKE or LIKELY_FAKE, which is what stops an article from being scored in the app. The rules were written before, and not fitted to, either dataset.

- **Labelled posts:** CONSTRAINT-2021 COVID-19 Fake News dataset (Patwa et al.), all three splits: 5,100 fake and 5,600 real posts. Median length 169 characters.
- **Real reports:** 149 WHO Disease Outbreak News reports published 2024-03-13 to 2026-11-13.

### Results on labelled posts

| | All 10,700 posts | Official test split (2,140) |
|---|---|---|
| Fake posts blocked (sensitivity) | 26.5% (1,352 of 5,100; CI 25.3–27.7) | 27.5% (CI 24.8–30.3) |
| Real posts passed (specificity) | 97.4% (5,455 of 5,600; CI 97.0–97.8) | 96.7% (CI 95.5–97.6) |
| Blocks that were correct (precision) | 90.3% (1,352 of 1,497; CI 88.7–91.7) | 88.3% (CI 84.3–91.4) |
| Accuracy | 63.6% | 63.7% |
| F1 for the fake class | 0.41 | 0.42 |

Verdicts given to the 5,100 fake posts: FAKE 105, LIKELY_FAKE 1,247, UNVERIFIABLE 1,088, LIKELY_REAL 2,660.

### Results on WHO reports

None of the 149 reports was blocked (CI 0–2.5%). Verdicts: LIKELY_REAL 139, UNVERIFIABLE 10.

### Reading

The rules are conservative: they rarely block genuine content, and when they block they are usually right, but they let about three in four fake posts through. On their own they are not an adequate misinformation filter. In the app they are the first stage; the second stage is an assessment by Claude, which was not evaluated here.

## 3. Outbreak information extraction

### Method

WHO report titles follow the pattern "Disease - Country". The title was parsed into an answer key and withheld; the extractor received only the report body. A disease is correct when the significant words of the shorter name are contained in the longer, after normalising common synonyms. Location was scored only for the 120 reports whose title names a specific place, not "Global situation", "Multi-country" or a WHO region. The answer key covers 58 distinct diseases. Matching is automatic and lenient; per-report outputs are saved for manual review.

### Results

| Extractor | Disease correct | Location correct |
|---|---|---|
| Keyword fallback (used when Claude is not connected) | 27 of 149 (18.1%; CI 12.8–25.1) | 0 of 120 (CI 0–3.1%) |
| Claude | not evaluated | not evaluated |

The keyword fallback recognises only Ebola, COVID-19 and mpox, and never extracts a location. This is the baseline against which the Claude extractor should be compared.

## 4. Genotype–phenotype validation on clinical isolates

### Method

- **Data.** NCBI Pathogen Detection, *Klebsiella* group (PDG000000012.2582, downloaded 2026-11-17). Eligible: *K. pneumoniae*, a genome assembly, a human host or clinical source, and laboratory resistant/susceptible results for at least 3 of the 15 study drugs. 1,002 isolates were eligible.
- **Sample.** At most 8 isolates per BioProject (to spread the sample across laboratories), then 100 drawn at random (seed 20261117). 99 were analysed; one (GCA_004322995.1) did not finish within the 60-minute limit. The 99 come from 43 BioProjects and at least 20 countries (13 without a recorded country), collected 2003–2025 (26 without a date).
- **Genotype.** Each assembly was downloaded from NCBI and run through the platform's AMRFinderPlus module (`AMRFinderClient.run`, nucleotide mode) with the organism set to *Klebsiella_pneumoniae*, which adds point-mutation screening. AMRFinderPlus 4.2.7, database 2026-10-29.1. NCBI's own precomputed genotypes were not used.
- **Phenotype.** The resistant (R) or susceptible (S) result as submitted to NCBI. Submitters used CLSI or EUCAST breakpoints of different years. Intermediate and other categories were excluded. 1,006 isolate–drug results: 680 R, 326 S.
- **Prediction rules.** Fixed in `isolate_benchmark.py` before results were seen. Genes and mutations are mapped to drugs through the AMRFinderPlus subclass: cephalosporins are predicted resistant from an acquired ESBL, AmpC or carbapenemase gene; carbapenems from a carbapenemase gene, or a porin (OmpK35/36) disruption together with an ESBL/AmpC gene; quinolones, aminoglycosides, trimethoprim-sulfamethoxazole (from *dfr*), tetracycline and colistin from any element of the matching subclass. Calls with an internal stop codon and the chromosomal genes intrinsic to *K. pneumoniae* (*oqxA*, *oqxB*, *fosA*) are ignored.
- **Clonal structure.** Isolates were grouped by NCBI SNP cluster; isolates outside any cluster form their own group.

### Correction made after the first run

The first version of the script compared intrinsic gene names exactly, but AMRFinderPlus reports allele variants (*oqxB5*, *oqxA11*, *fosA5*). The pre-specified exclusion therefore never applied, and *oqxB* alone produced false resistant calls for quinolones. The comparison was corrected to match variants and the cached AMRFinderPlus results were re-analysed; nothing else changed. Both outputs are kept:

| | Sensitivity | Specificity | Categorical agreement |
|---|---|---|---|
| As first run (`results/isolate_*_as_run_v1.*`) | 95.2% | 71.8% | 87.6% |
| Corrected (reported below) | 95.0% | 77.0% | 89.2% |

Only ciprofloxacin and levofloxacin results changed.

### Results by drug

Sensitivity is for detecting resistance. Very major error (VME) = resistant isolate predicted susceptible; major error (ME) = susceptible isolate predicted resistant.

| Drug | n | R | Sensitivity (CI) | Specificity (CI) | Agreement | VME | ME |
|---|---|---|---|---|---|---|---|
| Ceftriaxone | 63 | 55 | 96.4% (87.7–99.0) | 100% (67.6–100) | 96.8% | 3.6% | 0% |
| Cefotaxime | 56 | 49 | 100% (92.7–100) | 85.7% (48.7–97.4) | 98.2% | 0% | 14.3% |
| Ceftazidime | 78 | 65 | 98.5% (91.8–99.7) | 84.6% (57.8–95.7) | 96.2% | 1.5% | 15.4% |
| Cefepime | 74 | 53 | 96.2% (87.3–99.0) | 47.6% (28.3–67.6) | 82.4% | 3.8% | 52.4% |
| Meropenem | 90 | 55 | 96.4% (87.7–99.0) | 85.7% (70.6–93.7) | 92.2% | 3.6% | 14.3% |
| Imipenem | 75 | 48 | 95.8% (86.0–98.9) | 70.4% (51.5–84.2) | 86.7% | 4.2% | 29.6% |
| Ertapenem | 81 | 52 | 98.1% (89.9–99.7) | 93.1% (78.0–98.1) | 96.3% | 1.9% | 6.9% |
| Ciprofloxacin | 86 | 62 | 96.8% (89.0–99.1) | 79.2% (59.5–90.8) | 91.9% | 3.2% | 20.8% |
| Levofloxacin | 54 | 37 | 100% (90.6–100) | 58.8% (36.0–78.4) | 87.0% | 0% | 41.2% |
| Gentamicin | 82 | 48 | 91.7% (80.5–96.7) | 94.1% (80.9–98.4) | 92.7% | 8.3% | 5.9% |
| Tobramycin | 55 | 41 | 90.2% (77.5–96.1) | 100% (78.5–100) | 92.7% | 9.8% | 0% |
| Amikacin | 82 | 26 | 100% (87.1–100) | 53.6% (40.7–66.0) | 68.3% | 0% | 46.4% |
| Trimethoprim-sulfamethoxazole | 71 | 56 | 96.4% (87.9–99.0) | 73.3% (48.1–89.1) | 91.5% | 3.6% | 26.7% |
| Tetracycline | 38 | 26 | 73.1% (53.9–86.3) | 91.7% (64.6–98.5) | 79.0% | 26.9% | 8.3% |
| Colistin | 21 | 7 | 28.6% (8.2–64.1) | 92.9% (68.5–98.7) | 71.4% | 71.4% | 7.1% |
| **Pooled** | **1,006** | **680** | **95.0% (93.1–96.4)** | **77.0% (72.1–81.2)** | **89.2% (87.1–90.9)** | **5.0%** | **23.0%** |
| One isolate per SNP cluster | 968 | 647 | 95.1% (93.1–96.5) | 77.6% (72.7–81.8) | 89.3% (87.2–91.1) | 5.0% | 22.4% |

### Where the errors come from

- **False resistant calls** are dominated by *aac(6')-Ib-cr*, which AMRFinderPlus assigns to amikacin: of the 26 amikacin-susceptible isolates wrongly called resistant, 21 had *aac(6')-Ib-cr5* as the only determinant and 2 had *aac(6')-Ib*. This enzyme usually raises amikacin MICs without crossing the breakpoint. Cefepime false positives come from ESBLs in isolates still testing susceptible to cefepime. A few carbapenem false positives come from porin disruption plus ESBL, and from blaOXA-48, which hydrolyses imipenem weakly.
- **Missed resistance** concentrates in colistin (5 of 7 resistant isolates had no known determinant; colistin resistance in *K. pneumoniae* often comes from *mgrB* disruption by insertion sequences, which AMRFinderPlus does not report) and tetracycline (7 of 26).
- A per-drug rule refinement (for example, not counting *aac(6')-Ib* for amikacin) would reduce major errors, but was not applied: changing rules after seeing these results would overstate performance. It should be tested on a separate set of isolates.

### Clonal structure and cross-validation

The 99 isolates fell into 95 SNP clusters; only 7 shared a cluster with another sampled isolate, because sampling was capped per BioProject. Keeping one isolate per cluster changed pooled agreement by 0.1 percentage points.

A logistic regression on gene and mutation presence (12 drugs with ≥10 resistant and ≥10 susceptible isolates; 5-fold, 10 repeats) was evaluated with random and with cluster-aware folds:

| Drug | Random CV balanced accuracy | Cluster-aware | Difference |
|---|---|---|---|
| Ceftazidime | 85.0% | 84.9% | +0.1 |
| Cefepime | 73.8% | 72.8% | +1.0 |
| Meropenem | 78.5% | 79.1% | −0.6 |
| Imipenem | 73.9% | 73.8% | +0.1 |
| Ertapenem | 82.5% | 82.9% | −0.4 |
| Ciprofloxacin | 85.5% | 84.8% | +0.7 |
| Levofloxacin | 81.5% | 78.6% | +2.9 |
| Gentamicin | 80.1% | 80.5% | −0.4 |
| Tobramycin | 86.3% | 85.2% | +1.1 |
| Amikacin | 87.2% | 88.4% | −1.3 |
| Trimethoprim-sulfamethoxazole | 76.5% | 77.5% | −1.0 |
| Tetracycline | 55.5% | 51.9% | +3.6 |

The differences are small and in both directions. **This sample does not show clonal leakage**, because it contains almost no clonal pairs; that is a consequence of the sampling design, not evidence that leakage is unimportant. Testing leakage needs a sample with real clonal structure, such as whole outbreak collections. The learned model did not outperform the fixed rules and is reported only for this comparison.

### Speed

Median 212 s per genome with 2 CPU threads, 4 genomes in parallel, through WSL on an 8-core laptop.

### Reading

On public clinical *K. pneumoniae*, the platform's AMRFinderPlus module found a known determinant for 95% of laboratory-confirmed resistance, with the best results for third-generation cephalosporins, carbapenems (meropenem, ertapenem) and gentamicin/tobramycin. Specificity is lower, mainly because the presence of a gene does not always mean resistance above the clinical breakpoint. Colistin and tetracycline are poorly predicted.

## What has not been evaluated

- **Every Claude-based stage:** credibility assessment, outbreak extraction, and the risk score. The script supports this (`--with-claude`) but needs an API key and is billed to it.
- **The risk score** has no reference standard in this evaluation, with or without Claude.
- **The app's upload screen itself.** Section 4 called the platform's AMRFinderPlus client directly. The app's whole-genome upload (added afterwards) uses the same client, organism setting and prediction rules (`engine.predict_phenotype`; the benchmark imports it, and re-running it reproduced every prediction exactly). One assembly was also checked end to end through the upload path: identical AMR findings.
- **Species other than *K. pneumoniae*,** and isolates from other sources (animal, food, environment).
- **Raw sequencing reads (FASTQ).**
- **The web BLAST fallback path.**
- **Use by real users,** including usability and time saved.

## Limitations of what was evaluated

- ResFinder and the NCBI catalog are curated separately but draw on the same literature, and 91% of detected genes were exact matches. The detection figure measures agreement between two reference databases through this pipeline. It is not clinical sensitivity.
- The negative sets show the pipeline does not call resistance genes in one susceptible laboratory strain or in random sequence. They do not cover close non-resistance homologues from other species.
- CONSTRAINT posts are short social-media texts about one disease. The rules were designed for article-length reports, so performance on news articles may differ in either direction.
- The WHO answer key was produced by parsing titles automatically, and string matching can both over- and under-credit an extractor.
- WHO reports are well-written official documents; results on them say nothing about noisy or local-language sources.
- Section 4 phenotypes were submitted by many laboratories using different breakpoint standards and years, and were not re-interpreted against one standard. Public isolates with AST data over-represent resistant and outbreak isolates (68% of results were resistant). 99 isolates give wide intervals for individual drugs, especially specificity and colistin.
- These are public genomes that NCBI already screens with AMRFinderPlus. The reference catalogue may have been informed by some of them, so performance on new genomes could differ.

## Reproducing

```bash
pip install -r code/requirements.txt          # plus AMRFinderPlus, see README
python evaluation/download_data.py            # public datasets into evaluation/data/
python evaluation/amr_benchmark.py            # about 15 minutes; resumable
python evaluation/epidemic_benchmark.py       # under a minute
python evaluation/epidemic_benchmark.py --with-claude   # optional; needs ANTHROPIC_API_KEY
pip install scikit-learn                      # for section 4 only
python evaluation/isolate_benchmark.py --n 100 --workers 4 --threads 2   # several hours; resumable
```

Section 4 needs two files from NCBI Pathogen Detection in `evaluation/data/isolates/`: the *Klebsiella* metadata TSV saved as `kleb.metadata.tsv` (`/pathogen/Results/Klebsiella/latest_snps/Metadata/`) and the SNP cluster list saved as `kleb.clusters.tsv` (`.../Clusters/*.cluster_list.tsv`). NCBI replaces these files regularly; this run used release PDG000000012.2582, and the selected isolates are listed in `results/isolate_selection.csv`.

The 35-gene `--plus` re-run and the K-12 intrinsic-gene check were done separately; the first is saved as `results/amr_missed_genes_plus_scope.json`.

## Data sources

| Dataset | Version / date | Source | Licence note |
|---|---|---|---|
| ResFinder database | 2.7.0 | bitbucket.org/genomicepidemiology/resfinder_db | Apache-2.0 |
| AMRFinderPlus and Reference Gene Catalog | 4.2.7 / 2026-10-29.1 | NCBI | Public domain |
| *E. coli* K-12 MG1655 | NC_000913.3 | NCBI RefSeq | Public domain |
| CONSTRAINT-2021 COVID-19 Fake News | train, val, test | github.com/diptamath/covid_fake_news | Check the repository's terms before redistributing |
| WHO Disease Outbreak News | 149 reports, 2024–2026 | who.int public API | WHO terms of use apply |
| NCBI Pathogen Detection, *Klebsiella* | PDG000000012.2582 (metadata, SNP clusters, 99 assemblies) | ncbi.nlm.nih.gov/pathogens | Public domain; cite the submitting BioProjects where required |

Cite the original publications for ResFinder, AMRFinderPlus and CONSTRAINT in any manuscript.
