# Evaluation results

One Health Intelligence Platform v4.0.0 · run on 2026-11-15 · random seed 20261115

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

## What has not been evaluated

- **Every Claude-based stage:** credibility assessment, outbreak extraction, and the risk score. The script supports this (`--with-claude`) but needs an API key and is billed to it.
- **The risk score** has no reference standard in this evaluation, with or without Claude.
- **Whole genomes and real isolates.** All genomic tests used single gene sequences. The platform does not accept sequences over 500,000 bases, raw reads or assemblies.
- **Agreement with laboratory susceptibility results.** Nothing here measures whether a detected gene predicts phenotypic resistance.
- **Point-mutation resistance,** which requires specifying the organism.
- **The web BLAST fallback path.**
- **Use by real users,** including usability and time saved.

## Limitations of what was evaluated

- ResFinder and the NCBI catalog are curated separately but draw on the same literature, and 91% of detected genes were exact matches. The detection figure measures agreement between two reference databases through this pipeline. It is not clinical sensitivity.
- The negative sets show the pipeline does not call resistance genes in one susceptible laboratory strain or in random sequence. They do not cover close non-resistance homologues from other species.
- CONSTRAINT posts are short social-media texts about one disease. The rules were designed for article-length reports, so performance on news articles may differ in either direction.
- The WHO answer key was produced by parsing titles automatically, and string matching can both over- and under-credit an extractor.
- WHO reports are well-written official documents; results on them say nothing about noisy or local-language sources.

## Reproducing

```bash
pip install -r code/requirements.txt          # plus AMRFinderPlus, see README
python evaluation/download_data.py            # public datasets into evaluation/data/
python evaluation/amr_benchmark.py            # about 15 minutes; resumable
python evaluation/epidemic_benchmark.py       # under a minute
python evaluation/epidemic_benchmark.py --with-claude   # optional; needs ANTHROPIC_API_KEY
```

The 35-gene `--plus` re-run and the K-12 intrinsic-gene check were done separately; the first is saved as `results/amr_missed_genes_plus_scope.json`.

## Data sources

| Dataset | Version / date | Source | Licence note |
|---|---|---|---|
| ResFinder database | 2.7.0 | bitbucket.org/genomicepidemiology/resfinder_db | Apache-2.0 |
| AMRFinderPlus and Reference Gene Catalog | 4.2.7 / 2026-10-29.1 | NCBI | Public domain |
| *E. coli* K-12 MG1655 | NC_000913.3 | NCBI RefSeq | Public domain |
| CONSTRAINT-2021 COVID-19 Fake News | train, val, test | github.com/diptamath/covid_fake_news | Check the repository's terms before redistributing |
| WHO Disease Outbreak News | 149 reports, 2024–2026 | who.int public API | WHO terms of use apply |

Cite the original publications for ResFinder, AMRFinderPlus and CONSTRAINT in any manuscript.
