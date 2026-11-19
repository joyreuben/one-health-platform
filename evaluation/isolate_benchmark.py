"""
Genotype-phenotype validation on public clinical Klebsiella pneumoniae isolates.

Data: NCBI Pathogen Detection (Klebsiella group). Isolates with a genome
assembly AND laboratory antimicrobial susceptibility results (AST_phenotypes),
from human/clinical sources. NCBI's own AMR_genotypes column is NOT used: the
platform runs AMRFinderPlus itself on each downloaded assembly.

Two analyses:
  1. Rule-based prediction. Fixed, pre-specified rules (DRUG_RULES below,
     written before any result was seen) turn the platform's AMRFinderPlus
     findings into a resistant/susceptible prediction per drug, compared with
     the laboratory result. Reported per drug and pooled, on all isolates and
     on one isolate per clonal group (NCBI SNP cluster) to limit clonal bias.
  2. Learned prediction. A logistic regression on gene/mutation presence,
     evaluated with ordinary random cross-validation versus cluster-aware
     cross-validation (isolates of one SNP cluster never split between
     training and test), to measure clonal leakage.

Usage (from the project root):
  python evaluation/isolate_benchmark.py --n 100
  python evaluation/isolate_benchmark.py --analyze-only    # reuse cached runs

Inputs (downloaded once, git-ignored) in evaluation/data/isolates/:
  kleb.metadata.tsv   https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/latest_snps/Metadata/
  kleb.clusters.tsv   https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/latest_snps/Clusters/ (cluster_list)
Assemblies are fetched from the NCBI Datasets API and cached in asm/;
AMRFinderPlus results are cached in amr/ so an interrupted run resumes.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import math
import os
import random
import sys
import threading
import time
import urllib.request
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

logging.disable(logging.CRITICAL)
os.environ.setdefault("OHP_DISABLE_FEED_WORKER", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "code"))
from engine import (AMRFinderClient, _amrfinder_info, PHENOTYPE_DRUG_RULES,  # noqa: E402
                    INTRINSIC_GENES, predict_phenotype)

DATA = os.path.join(HERE, "data", "isolates")
ASM_DIR = os.path.join(DATA, "asm")
AMR_DIR = os.path.join(DATA, "amr")
RESULTS = os.path.join(HERE, "results")
SEED = 20261117
ORGANISM = "Klebsiella_pneumoniae"
DATASETS_URL = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
                "?include_annotation_type=GENOME_FASTA")
csv.field_size_limit(2**31 - 1)

# ---------------------------------------------------------------------------
# Prediction rules: imported from the platform (engine.PHENOTYPE_DRUG_RULES,
# engine.predict_phenotype), so the app shows exactly the rules measured here.
# They were fixed before results were seen. Ignored: non-AMR elements,
# INTERNAL_STOP calls, and genes intrinsic to K. pneumoniae (oqxA/oqxB, fosA,
# allele variants included). v1 of this script matched intrinsic names
# exactly, so that exclusion silently never applied; v1 output is kept as
# results/isolate_*_as_run_v1.* and both are reported in RESULTS.md.
# ---------------------------------------------------------------------------
DRUG_RULES = PHENOTYPE_DRUG_RULES
INTRINSIC = INTRINSIC_GENES
predict = predict_phenotype


# ---------------------------------------------------------------------------
# Data selection
# ---------------------------------------------------------------------------
def parse_ast(text):
    out = {}
    for item in (text or "").split(","):
        if "=" in item:
            d, c = item.split("=", 1)
            d, c = d.strip(), c.strip().upper()
            if d in DRUG_RULES and c in ("R", "S"):
                out[d] = 1 if c == "R" else 0
    return out


def load_candidates():
    clusters = {}
    with open(os.path.join(DATA, "kleb.clusters.tsv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("gencoll_acc") and row["gencoll_acc"] != "NULL":
                clusters[row["gencoll_acc"]] = row["PDS_acc"].split(".")[0]
    cands, seen = [], set()
    with open(os.path.join(DATA, "kleb.metadata.tsv"), encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = row.get("asm_acc", "NULL")
            if acc in ("NULL", "") or acc in seen:
                continue
            if not (row.get("scientific_name") or "").startswith("Klebsiella pneumoniae"):
                continue
            human = "homo sapiens" in (row.get("host") or "").lower()
            if not (human or (row.get("epi_type") or "").lower() == "clinical"):
                continue
            ast = parse_ast(row.get("AST_phenotypes"))
            if len(ast) < 3:
                continue
            seen.add(acc)
            cands.append({
                "asm_acc": acc, "biosample": row.get("biosample_acc", ""),
                "bioproject": row.get("bioproject_acc", ""),
                "cluster": clusters.get(acc, "singleton:" + acc),
                "country": (row.get("geo_loc_name") or "").split(":")[0],
                "year": (row.get("collection_date") or "")[:4],
                "source": row.get("isolation_source", ""),
                "asm_level": row.get("asm_level", ""),
                "ast": ast,
            })
    return cands


def select(cands, n, per_project):
    rng = random.Random(SEED)
    by_proj = defaultdict(list)
    for c in sorted(cands, key=lambda c: c["asm_acc"]):
        by_proj[c["bioproject"]].append(c)
    pool = []
    for proj in sorted(by_proj):
        items = by_proj[proj][:]
        rng.shuffle(items)
        pool.extend(items[:per_project])
    rng.shuffle(pool)
    return sorted(pool[:n], key=lambda c: c["asm_acc"])


# ---------------------------------------------------------------------------
# Download + AMRFinderPlus (cached, resumable)
# ---------------------------------------------------------------------------
_dl_lock = threading.Lock()


def download(acc):
    path = os.path.join(ASM_DIR, acc + ".zip")
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    for attempt in range(5):
        try:
            with _dl_lock:  # stay well under the API's 5 requests/second
                time.sleep(0.4)
            req = urllib.request.Request(DATASETS_URL.format(acc=acc),
                                         headers={"User-Agent": "ohp-evaluation/1.0"})
            with urllib.request.urlopen(req, timeout=300) as r:
                data = r.read()
            zipfile.ZipFile(io.BytesIO(data)).testzip()
            with open(path + ".part", "wb") as fh:
                fh.write(data)
            os.replace(path + ".part", path)
            return path
        except Exception:
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"download failed: {acc}")


def read_assembly(path):
    z = zipfile.ZipFile(path)
    names = [n for n in z.namelist() if n.endswith((".fna", ".fa", ".fasta"))]
    if not names:
        raise RuntimeError("no FASTA in archive")
    recs, cur, buf = {}, None, []
    for line in z.read(names[0]).decode("utf-8", "replace").splitlines():
        if line.startswith(">"):
            if cur:
                recs[cur] = "".join(buf)
            cur, buf = "c%d" % (len(recs) + 1), []
        else:
            buf.append(line.strip().upper())
    if cur:
        recs[cur] = "".join(buf)
    return recs


def analyse_isolate(acc, threads):
    out = os.path.join(AMR_DIR, acc + ".json")
    if os.path.exists(out):
        with open(out, encoding="utf-8") as fh:
            return json.load(fh)
    recs = read_assembly(download(acc))
    t = time.time()
    res = AMRFinderClient.run(recs, protein=False, organism=ORGANISM, timeout=3600, threads=threads)
    rec = {"asm_acc": acc, "seconds": round(time.time() - t, 1), "contigs": len(recs),
           "length_bp": sum(map(len, recs.values())),
           "findings": [f for v in res.values() for f in v]}
    with open(out + ".part", "w", encoding="utf-8") as fh:
        json.dump(rec, fh)
    os.replace(out + ".part", out)
    return rec


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(100 * max(0, c - h), 2), round(100 * min(1, c + h), 2)]


def pct(k, n):
    return round(100 * k / n, 2) if n else None


def confusion(pairs):
    tp = sum(1 for p, y in pairs if p == 1 and y == 1)
    fn = sum(1 for p, y in pairs if p == 0 and y == 1)
    fp = sum(1 for p, y in pairs if p == 1 and y == 0)
    tn = sum(1 for p, y in pairs if p == 0 and y == 0)
    nr, ns = tp + fn, fp + tn
    return {
        "n": nr + ns, "resistant": nr, "susceptible": ns, "TP": tp, "FN": fn, "FP": fp, "TN": tn,
        "sensitivity": pct(tp, nr), "sensitivity_ci": wilson(tp, nr),
        "specificity": pct(tn, ns), "specificity_ci": wilson(tn, ns),
        "categorical_agreement": pct(tp + tn, nr + ns), "agreement_ci": wilson(tp + tn, nr + ns),
        "very_major_error_rate": pct(fn, nr), "major_error_rate": pct(fp, ns),
    }


def rule_analysis(isolates, amr):
    rows, per_drug = [], defaultdict(list)
    for iso in isolates:
        fs = amr[iso["asm_acc"]]["findings"]
        for drug, y in sorted(iso["ast"].items()):
            p, why = predict(drug, fs)
            per_drug[drug].append((p, y, iso["cluster"]))
            rows.append({"asm_acc": iso["asm_acc"], "cluster": iso["cluster"], "drug": drug,
                         "phenotype": "R" if y else "S", "predicted": "R" if p else "S",
                         "agree": int(p == y), "determinants": ";".join(why)})
    rng = random.Random(SEED)
    drugs, pooled, pooled_dedup = {}, [], []
    for drug in sorted(per_drug):
        trip = per_drug[drug]
        pooled += [(p, y) for p, y, _ in trip]
        by_c = defaultdict(list)
        for p, y, c in trip:
            by_c[c].append((p, y))
        dedup = [rng.choice(v) for _, v in sorted(by_c.items())]
        pooled_dedup += dedup
        drugs[drug] = {"all": confusion([(p, y) for p, y, _ in trip]),
                       "one_per_cluster": confusion(dedup)}
    return rows, {"per_drug": drugs, "pooled_all": confusion(pooled),
                  "pooled_one_per_cluster": confusion(pooled_dedup)}


def ml_analysis(isolates, amr, min_class=10, repeats=10):
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import balanced_accuracy_score, roc_auc_score
    from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

    feats = sorted({f["symbol"] for a in amr.values() for f in a["findings"]
                    if (f.get("type") or "").upper() == "AMR"})
    idx = {s: i for i, s in enumerate(feats)}
    out = {"features": len(feats), "model": "L2 logistic regression, class_weight=balanced, C=1.0",
           "folds": 5, "repeats": repeats, "drugs": {}}
    for drug in DRUG_RULES:
        sub = [i for i in isolates if drug in i["ast"]]
        y = np.array([i["ast"][drug] for i in sub])
        if y.sum() < min_class or (len(y) - y.sum()) < min_class:
            continue
        X = np.zeros((len(sub), len(feats)))
        for r, iso in enumerate(sub):
            for f in amr[iso["asm_acc"]]["findings"]:
                if f["symbol"] in idx:
                    X[r, idx[f["symbol"]]] = 1
        groups = np.array([i["cluster"] for i in sub])
        res = {}
        for name in ("random", "cluster_aware"):
            bas, aucs = [], []
            for rep in range(repeats):
                cv = (StratifiedKFold(5, shuffle=True, random_state=SEED + rep) if name == "random"
                      else StratifiedGroupKFold(5, shuffle=True, random_state=SEED + rep))
                prob = np.zeros(len(y))
                splits = cv.split(X, y) if name == "random" else cv.split(X, y, groups)
                for tr, te in splits:
                    m = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
                    m.fit(X[tr], y[tr])
                    prob[te] = m.predict_proba(X[te])[:, 1]
                bas.append(balanced_accuracy_score(y, prob >= 0.5))
                aucs.append(roc_auc_score(y, prob))
            res[name] = {"balanced_accuracy_mean": round(100 * float(np.mean(bas)), 2),
                         "balanced_accuracy_sd": round(100 * float(np.std(bas)), 2),
                         "auc_mean": round(float(np.mean(aucs)), 3)}
        res["n"], res["resistant"] = len(y), int(y.sum())
        res["clusters"] = len(set(groups))
        res["inflation_pp"] = round(res["random"]["balanced_accuracy_mean"]
                                    - res["cluster_aware"]["balanced_accuracy_mean"], 2)
        out["drugs"][drug] = res
    return out


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--per-project", type=int, default=8,
                    help="max isolates per BioProject, to spread the sample across laboratories")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--analyze-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(ASM_DIR, exist_ok=True)
    os.makedirs(AMR_DIR, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)

    cands = load_candidates()
    chosen = select(cands, a.n, a.per_project)
    print(f"candidates: {len(cands)} clinical K. pneumoniae with assembly + AST; selected {len(chosen)}"
          f" from {len({c['bioproject'] for c in chosen})} BioProjects", flush=True)
    with open(os.path.join(RESULTS, "isolate_selection.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["asm_acc", "biosample", "bioproject", "snp_cluster", "country", "year",
                    "isolation_source", "asm_level", "drugs_tested"])
        for c in chosen:
            w.writerow([c["asm_acc"], c["biosample"], c["bioproject"], c["cluster"], c["country"],
                        c["year"], c["source"], c["asm_level"], len(c["ast"])])

    amr, failed, t0 = {}, {}, time.time()
    if a.analyze_only:
        for c in chosen:
            p = os.path.join(AMR_DIR, c["asm_acc"] + ".json")
            if os.path.exists(p):
                with open(p, encoding="utf-8") as fh:
                    amr[c["asm_acc"]] = json.load(fh)
            else:
                failed[c["asm_acc"]] = "no cached AMRFinderPlus result (failed or not yet run)"
    else:
        with ThreadPoolExecutor(a.workers) as ex:
            futs = {ex.submit(analyse_isolate, c["asm_acc"], a.threads): c["asm_acc"] for c in chosen}
            for k, fut in enumerate(as_completed(futs), 1):
                acc = futs[fut]
                try:
                    amr[acc] = fut.result()
                    msg = f"{amr[acc]['seconds']}s, {len(amr[acc]['findings'])} findings"
                except Exception as e:
                    failed[acc] = str(e)[:200]
                    msg = "FAILED " + failed[acc]
                el = time.time() - t0
                print(f"[{k}/{len(chosen)}] {acc}: {msg} | elapsed {el/60:.1f} min", flush=True)

    done = [c for c in chosen if c["asm_acc"] in amr]
    rows, rules = rule_analysis(done, amr)
    ml = ml_analysis(done, amr)
    with open(os.path.join(RESULTS, "isolate_predictions.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(RESULTS, "isolate_findings.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["asm_acc", "symbol", "type", "subtype", "class", "subclass", "method", "identity", "coverage"])
        for acc in sorted(amr):
            for f in amr[acc]["findings"]:
                w.writerow([acc, f["symbol"], f["type"], f["subtype"], f["class"], f["subclass"],
                            f["method"], f["identity"], f["coverage"]])

    info = _amrfinder_info(refresh=True)
    secs = [amr[k]["seconds"] for k in amr]
    summary = {
        "dataset": "NCBI Pathogen Detection, Klebsiella (PDG000000012.2582), clinical K. pneumoniae "
                   "isolates with assembly and AST phenotypes",
        "seed": SEED, "organism_option": ORGANISM,
        "amrfinder": {"version": info.get("software"), "database": info.get("database")},
        "candidates": len(cands), "selected": len(chosen), "analysed": len(done),
        "failed": failed, "bioprojects": len({c["bioproject"] for c in done}),
        "countries": dict(Counter(c["country"] for c in done).most_common()),
        "snp_clusters": len({c["cluster"] for c in done}),
        "isolates_in_shared_clusters": sum(1 for c in done
                                           if sum(1 for d in done if d["cluster"] == c["cluster"]) > 1),
        "drug_results": len(rows),
        "median_seconds_per_isolate": sorted(secs)[len(secs) // 2] if secs else None,
        "phenotype_note": "Laboratory R/S as submitted to NCBI (CLSI or EUCAST per submitter); "
                          "intermediate and other categories excluded.",
        "rules": {"intrinsic_ignored": sorted(INTRINSIC), "drug_rules": DRUG_RULES},
        "rule_based": rules, "machine_learning": ml,
    }
    with open(os.path.join(RESULTS, "isolate_benchmark.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print("\n=== Rule-based genotype -> phenotype ===")
    print(f"{'drug':32} {'n':>4} {'R':>4} {'sens':>7} {'spec':>7} {'agree':>7} {'VME':>6} {'ME':>6}")
    for drug, r in rules["per_drug"].items():
        x = r["all"]
        print(f"{drug:32} {x['n']:>4} {x['resistant']:>4} {x['sensitivity']!s:>7} {x['specificity']!s:>7} "
              f"{x['categorical_agreement']!s:>7} {x['very_major_error_rate']!s:>6} {x['major_error_rate']!s:>6}")
    for k in ("pooled_all", "pooled_one_per_cluster"):
        x = rules[k]
        print(f"{k:32} {x['n']:>4} {x['resistant']:>4} {x['sensitivity']!s:>7} {x['specificity']!s:>7} "
              f"{x['categorical_agreement']!s:>7} {x['very_major_error_rate']!s:>6} {x['major_error_rate']!s:>6}")
    print(f"\nisolates {len(done)}, SNP clusters {summary['snp_clusters']}, "
          f"in shared clusters {summary['isolates_in_shared_clusters']}")
    print("\n=== Logistic regression: random vs cluster-aware CV (balanced accuracy %) ===")
    for drug, r in ml["drugs"].items():
        print(f"{drug:32} n={r['n']:>3} R={r['resistant']:>3} clusters={r['clusters']:>3} "
              f"random {r['random']['balanced_accuracy_mean']:>6} | cluster-aware "
              f"{r['cluster_aware']['balanced_accuracy_mean']:>6} | inflation {r['inflation_pp']}")


if __name__ == "__main__":
    main()
