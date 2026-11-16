"""Benchmark of the platform's resistance-gene detection.

    python evaluation/download_data.py
    python evaluation/amr_benchmark.py

Every sequence goes through the app's own batch pipeline
(engine.run_amr_batch, AMRFinderPlus mode), so this measures the platform
end to end: validation, the AMRFinderPlus call, parsing and interpretation.

Test sets
  Positives   ResFinder database genes (DTU): acquired resistance genes curated
              independently of the NCBI Reference Gene Catalog that AMRFinderPlus uses.
  Negatives   (a) every coding sequence of E. coli K-12 MG1655 (NC_000913.3)
              (b) a base-shuffled copy of each positive (same length and composition)
  Robustness  a random sample of positives with random substitutions (1-20 %)
              or truncated to 90-25 % of their length

Baseline
  The platform's previous, pre-AMRFinderPlus classifier, evaluated offline:
  its keyword matcher given the true gene name (a best case: it normally sees a
  BLAST hit title), and its sequence-motif scan.

Outputs: evaluation/results/amr_benchmark.json and CSV files alongside it.
"""
import csv
import json
import logging
import math
import os
import random
import sys
import tempfile
import time

logging.disable(logging.CRITICAL)
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "results")
os.environ["OHP_DISABLE_FEED_WORKER"] = "1"
os.environ.setdefault("OHP_DATA_DIR", tempfile.mkdtemp(prefix="ohp_eval_"))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "code"))

import engine  # noqa: E402
from engine import AMRFinderClient, AMRValidationEngine  # noqa: E402

SEED = 20261115
# Beta-lactamase alleles each match thousands of near-identical references, so
# AMRFinderPlus needs ~0.3 s per such gene; small chunks keep every call well
# inside the client timeout, and finished chunks are cached so a run can resume.
CHUNK = 150
ROBUSTNESS_SAMPLE = 300
CACHE = os.path.join(DATA, "cache")
SUBSTITUTION_RATES = [0.01, 0.02, 0.05, 0.10, 0.15, 0.20]
TRUNCATION_KEEP = [0.90, 0.75, 0.50, 0.25]

# ResFinder "Class" vocabulary -> AMRFinderPlus class/subclass tokens that count as agreement.
RESFINDER_TO_AMRFINDER = {
    "aminoglycoside": {"AMINOGLYCOSIDE"},
    "aminocyclitol": {"AMINOGLYCOSIDE", "SPECTINOMYCIN"},
    "beta-lactam": {"BETA-LACTAM"},
    "polymyxin": {"COLISTIN"},
    "fosfomycin": {"FOSFOMYCIN"},
    "steroid antibacterial": {"FUSIDIC ACID"},
    "glycopeptide": {"GLYCOPEPTIDE"},
    "macrolide": {"MACROLIDE"},
    "lincosamide": {"LINCOSAMIDE"},
    "streptogramin a": {"STREPTOGRAMIN"},
    "streptogramin b": {"STREPTOGRAMIN"},
    "pleuromutilin": {"PLEUROMUTILIN"},
    "nitroimidazole": {"NITROIMIDAZOLE"},
    "oxazolidinone": {"OXAZOLIDINONE"},
    "amphenicol": {"PHENICOL"},
    "phenicol": {"PHENICOL"},
    "tetracyclines": {"TETRACYCLINE"},
    "pseudomonic acid": {"MUPIROCIN"},
    "quinolone": {"QUINOLONE"},
    "rifamycin": {"RIFAMYCIN"},
    "folate pathway antagonist": {"TRIMETHOPRIM", "SULFONAMIDE"},
    "tetracycline": {"TETRACYCLINE"},
}
# Fallback when a gene is missing from phenotypes.txt: the ResFinder file it came from.
FILE_TO_AMRFINDER = {
    "aminoglycoside": {"AMINOGLYCOSIDE"}, "beta-lactam": {"BETA-LACTAM"}, "colistin": {"COLISTIN"},
    "fosfomycin": {"FOSFOMYCIN"}, "fusidicacid": {"FUSIDIC ACID"}, "glycopeptide": {"GLYCOPEPTIDE"},
    "macrolide": {"MACROLIDE", "LINCOSAMIDE", "STREPTOGRAMIN", "PLEUROMUTILIN"},
    "nitroimidazole": {"NITROIMIDAZOLE"}, "oxazolidinone": {"OXAZOLIDINONE"}, "phenicol": {"PHENICOL"},
    "pseudomonicacid": {"MUPIROCIN"}, "quinolone": {"QUINOLONE"}, "rifampicin": {"RIFAMYCIN"},
    "sulphonamide": {"SULFONAMIDE"}, "tetracycline": {"TETRACYCLINE"}, "trimethoprim": {"TRIMETHOPRIM"},
}
# Drug classes the legacy keyword markers stand for.
LEGACY_MARKER_CLASS = {
    "ndm": "BETA-LACTAM", "kpc": "BETA-LACTAM", "ctx-m": "BETA-LACTAM", "blaTEM": "BETA-LACTAM",
    "mecA": "BETA-LACTAM", "mcr": "COLISTIN", "vanA": "GLYCOPEPTIDE", "tetA": "TETRACYCLINE",
    "tetB": "TETRACYCLINE", "ermB": "MACROLIDE", "aacA": "AMINOGLYCOSIDE", "qnr": "QUINOLONE",
}


VERSION_KEY = ""   # platform + tool + database versions; part of every cache key


def wilson(k: int, n: int, z: float = 1.96):
    """Wilson 95 % score interval for a proportion."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def prop(k: int, n: int):
    lo, hi = wilson(k, n)
    return {"k": k, "n": n, "pct": round(100 * k / n, 2) if n else None,
            "ci95": [round(100 * lo, 2), round(100 * hi, 2)]}


def read_fasta(path):
    name, seq = None, []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq).upper()
                name, seq = line[1:], []
            elif line:
                seq.append(line)
    if name:
        yield name, "".join(seq).upper()


def load_positives():
    """ResFinder genes: header -> {seq, file_classes, expected tokens}."""
    rf = os.path.join(DATA, "resfinder_db")
    pheno = {}
    with open(os.path.join(rf, "phenotypes.txt"), encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pheno[row["Gene_accession no."].strip()] = row.get("Class", "")
    genes, unknown_classes = {}, set()
    for fn in sorted(os.listdir(rf)):
        if not fn.endswith(".fsa") or fn == "all.fsa":
            continue
        file_class = fn[:-4]
        for header, seq in read_fasta(os.path.join(rf, fn)):
            g = genes.setdefault(header, {"seq": seq, "files": set(), "expected": set(), "resfinder_class": ""})
            g["files"].add(file_class)
    for header, g in genes.items():
        cls = pheno.get(header, "")
        g["resfinder_class"] = cls
        for part in [p.strip().lower() for p in cls.split(",") if p.strip()]:
            if part in RESFINDER_TO_AMRFINDER:
                g["expected"] |= RESFINDER_TO_AMRFINDER[part]
            else:
                unknown_classes.add(part)
        if not g["expected"]:
            for fc in g["files"]:
                g["expected"] |= FILE_TO_AMRFINDER.get(fc, set())
    return genes, sorted(unknown_classes)


def run_pipeline(records: dict, label: str):
    """Run {id: sequence} through the app's batch pipeline in chunks."""
    import hashlib
    os.makedirs(CACHE, exist_ok=True)
    ids = list(records)
    out, elapsed = {}, 0.0
    for i in range(0, len(ids), CHUNK):
        chunk = {k: records[k] for k in ids[i:i + CHUNK]}
        key = hashlib.sha256(json.dumps([VERSION_KEY, sorted(chunk.items())]).encode()).hexdigest()[:20]
        path = os.path.join(CACHE, f"{key}.json")
        if os.path.exists(path):
            with open(path, encoding="utf-8") as f:
                saved = json.load(f)
            results, secs = saved["results"], saved["seconds"]
        else:
            t0 = time.time()
            results = engine.run_amr_batch(chunk, 90, "amrfinder")
            secs = time.time() - t0
            # Input-validation rejections are legitimate results; tool failures are not.
            failed = [r for r in results if not r.get("ok") and not r.get("negative")
                      and ("AMRFinderPlus" in r.get("error", "") or "timed out" in r.get("error", "").lower())]
            if failed:
                sys.exit(f"AMRFinderPlus call failed in '{label}' chunk {i}: {failed[0]['error'][:200]} (re-run to resume)")
            for r in results:
                r.pop("seq_clean", None)
                r.pop("hl7", None)
            with open(path, "w", encoding="utf-8") as f:
                json.dump({"results": results, "seconds": secs}, f)
        elapsed += secs
        for res in results:
            out[res["sample_id"]] = res
        print(f"  {label}: {min(i + CHUNK, len(ids))}/{len(ids)}  ({secs:.0f}s)", flush=True)
    return out, elapsed


def tokens(result):
    t = set()
    for f in result.get("findings", []):
        for field in (f.get("class", ""), f.get("subclass", "")):
            t |= {x.strip().upper() for x in field.split("/") if x.strip()}
    return t


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    s = list(seq)
    for pos in rng.sample(range(len(s)), max(1, round(len(s) * rate))):
        s[pos] = rng.choice([b for b in "ACGT" if b != s[pos]])
    return "".join(s)


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = random.Random(SEED)
    global VERSION_KEY
    AMRFinderClient.THREADS = os.cpu_count() or 4
    AMRFinderClient.TIMEOUT = 900
    info = AMRFinderClient.info()
    if not info.get("ok"):
        sys.exit(f"AMRFinderPlus is not available: {info.get('error')}")
    VERSION_KEY = f"{engine.__version__}|{info['software']}|{info['database']}"
    print(f"AMRFinderPlus {info['software']}, database {info['database']}")

    genes, unknown_classes = load_positives()
    rf_version = open(os.path.join(DATA, "resfinder_db", "VERSION")).read().strip()
    print(f"ResFinder {rf_version}: {len(genes)} unique gene entries")
    if unknown_classes:
        print("ResFinder classes without a mapping (fall back to file name):", unknown_classes)

    # ---------------------------------------------------------- positives
    pos_res, pos_time = run_pipeline({h: g["seq"] for h, g in genes.items()}, "positives")
    rows, per_class = [], {}
    n_detected = n_class_eval = n_class_ok = n_invalid = 0
    for h, g in genes.items():
        r = pos_res[h]
        invalid = not r.get("ok") and not r.get("negative")
        detected = bool(r.get("ok"))
        tk = tokens(r) if detected else set()
        class_ok = bool(g["expected"] & tk) if (detected and g["expected"]) else None
        n_invalid += invalid
        n_detected += detected
        if detected and g["expected"]:
            n_class_eval += 1
            n_class_ok += bool(class_ok)
        for fc in g["files"]:
            pc = per_class.setdefault(fc, {"n": 0, "detected": 0, "class_eval": 0, "class_ok": 0})
            pc["n"] += 1
            pc["detected"] += detected
            pc["class_eval"] += class_ok is not None
            pc["class_ok"] += bool(class_ok)
        top = r["findings"][0] if detected else {}
        rows.append({
            "resfinder_gene": h, "resfinder_file": ";".join(sorted(g["files"])), "resfinder_class": g["resfinder_class"],
            "length_bp": len(g["seq"]), "detected": int(detected), "invalid_input": int(invalid),
            "amrfinder_gene": top.get("symbol", ""), "amrfinder_class": top.get("class", ""),
            "amrfinder_subclass": top.get("subclass", ""), "method": top.get("method", ""),
            "identity": top.get("identity", ""), "coverage": top.get("coverage", ""),
            "n_findings": len(r.get("findings", [])), "class_concordant": "" if class_ok is None else int(class_ok),
            "confidence": round(r["confidence"], 4) if detected else "", "concern": r["gene_info"]["clinical_concern"] if detected else "",
            "note": "" if detected else r.get("error", "")[:120],
        })
    with open(os.path.join(OUT, "amr_positives_per_gene.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    methods = {}
    for r in rows:
        if r["detected"]:
            methods[r["method"]] = methods.get(r["method"], 0) + 1

    # ---------------------------------------------------------- negatives
    k12 = {}
    for i, (header, seq) in enumerate(read_fasta(os.path.join(DATA, "ecoli_k12_mg1655_cds.fna"))):
        gene = header.split("[gene=")[1].split("]")[0] if "[gene=" in header else f"cds{i}"
        k12[f"k12_{i}_{gene}"] = seq
    k12_res, k12_time = run_pipeline(k12, "E. coli K-12 CDS")
    shuffled = {}
    for h, g in genes.items():
        s = list(g["seq"])
        rng.shuffle(s)
        shuffled[f"shuf_{h}"] = "".join(s)
    shuf_res, shuf_time = run_pipeline(shuffled, "shuffled controls")

    fp_rows = []
    for name, res in (("ecoli_k12", k12_res), ("shuffled", shuf_res)):
        for sid, r in res.items():
            if r.get("ok"):
                for f in r["findings"]:
                    fp_rows.append({"set": name, "sequence": sid, "amrfinder_gene": f["symbol"], "name": f["name"],
                                    "class": f["class"], "subclass": f["subclass"], "method": f["method"],
                                    "identity": f["identity"], "coverage": f["coverage"]})
    with open(os.path.join(OUT, "amr_negative_set_calls.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["set", "sequence", "amrfinder_gene", "name", "class", "subclass",
                                          "method", "identity", "coverage"])
        w.writeheader()
        w.writerows(fp_rows)
    k12_valid = sum(1 for r in k12_res.values() if r.get("ok") or r.get("negative"))
    k12_called = sum(1 for r in k12_res.values() if r.get("ok"))
    shuf_valid = sum(1 for r in shuf_res.values() if r.get("ok") or r.get("negative"))
    shuf_called = sum(1 for r in shuf_res.values() if r.get("ok"))

    # --------------------------------------------------------- robustness
    detected_ids = [h for h in genes if pos_res[h].get("ok")]
    sample = rng.sample(detected_ids, min(ROBUSTNESS_SAMPLE, len(detected_ids)))
    robustness = {"sample_size": len(sample), "substitution": [], "truncation": []}
    for rate in SUBSTITUTION_RATES:
        recs = {f"m{i}": mutate(genes[h]["seq"], rate, rng) for i, h in enumerate(sample)}
        res, _ = run_pipeline(recs, f"substitution {rate:.0%}")
        det = [res[f"m{i}"] for i in range(len(sample)) if res[f"m{i}"].get("ok")]
        cls_ok = sum(1 for i, h in enumerate(sample)
                     if res[f"m{i}"].get("ok") and genes[h]["expected"] & tokens(res[f"m{i}"]))
        cls_n = sum(1 for i, h in enumerate(sample) if res[f"m{i}"].get("ok") and genes[h]["expected"])
        robustness["substitution"].append({
            "nucleotide_substitution_rate": rate, "detected": prop(len(det), len(sample)),
            "class_concordant_among_detected": prop(cls_ok, cls_n),
            "mean_confidence": round(sum(r["confidence"] for r in det) / len(det), 4) if det else None,
            "mean_protein_identity": round(sum(r["best"]["identity_percentage"] for r in det) / len(det), 2) if det else None})
    for keep in TRUNCATION_KEEP:
        recs = {f"t{i}": genes[h]["seq"][:max(60, int(len(genes[h]["seq"]) * keep))] for i, h in enumerate(sample)}
        res, _ = run_pipeline(recs, f"truncation keep {keep:.0%}")
        det = [res[f"t{i}"] for i in range(len(sample)) if res[f"t{i}"].get("ok")]
        robustness["truncation"].append({
            "fraction_of_gene_kept": keep, "detected": prop(len(det), len(sample)),
            "mean_confidence": round(sum(r["confidence"] for r in det) / len(det), 4) if det else None})

    # ------------------------------------------------- legacy baseline
    def legacy_title(name):
        hit = AMRValidationEngine._match_title(AMRValidationEngine._normalise(name))
        return hit["marker_found"] if hit else None

    leg_title_hit = leg_title_ok = leg_motif_hit = leg_motif_ok = 0
    for h, g in genes.items():
        gene_name = h.split("_")[0]
        m = legacy_title(gene_name)
        if m:
            leg_title_hit += 1
            leg_title_ok += LEGACY_MARKER_CLASS.get(m) in g["expected"]
        mm = AMRValidationEngine._sequence_motif_scan(g["seq"])
        if mm:
            leg_motif_hit += 1
            leg_motif_ok += LEGACY_MARKER_CLASS.get(mm["marker_found"]) in g["expected"]
    k12_motif_fp = sum(1 for s in k12.values() if AMRValidationEngine._sequence_motif_scan(s))
    k12_title_fp = sum(1 for k in k12 if legacy_title(k.split("_", 2)[2]))

    n = len(genes)
    n_valid = n - n_invalid
    total_seqs = len(genes) + len(k12) + len(shuffled)
    total_time = pos_time + k12_time + shuf_time
    summary = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "platform_version": engine.__version__,
        "amrfinderplus": {"software": info["software"], "database": info["database"]},
        "resfinder_db_version": rf_version,
        "random_seed": SEED,
        "positives": {
            "source": "ResFinder database, per-class FASTA files, unique headers",
            "n": n, "rejected_by_input_validation": n_invalid,
            "detected": prop(n_detected, n_valid),
            "class_concordant_among_detected": prop(n_class_ok, n_class_eval),
            "detected_and_class_concordant": prop(n_class_ok, n_valid),
            "detection_method_counts": dict(sorted(methods.items(), key=lambda kv: -kv[1])),
            "per_resfinder_file": {k: {"n": v["n"], "detected": prop(v["detected"], v["n"]),
                                       "class_concordant_among_detected": prop(v["class_ok"], v["class_eval"])}
                                   for k, v in sorted(per_class.items())},
        },
        "negatives": {
            "ecoli_k12_cds": {"n": len(k12), "valid_input": k12_valid, "with_any_amr_call": prop(k12_called, k12_valid),
                              "no_call": prop(k12_valid - k12_called, k12_valid)},
            "shuffled_positives": {"n": len(shuffled), "valid_input": shuf_valid,
                                   "with_any_amr_call": prop(shuf_called, shuf_valid),
                                   "no_call": prop(shuf_valid - shuf_called, shuf_valid)},
        },
        "robustness": robustness,
        "legacy_baseline": {
            "keyword_matcher_given_true_gene_name": {"recognised": prop(leg_title_hit, n),
                                                     "recognised_and_class_correct": prop(leg_title_ok, n)},
            "sequence_motif_scan": {"recognised": prop(leg_motif_hit, n),
                                    "recognised_and_class_correct": prop(leg_motif_ok, n)},
            "ecoli_k12_calls": {"keyword_matcher_on_gene_name": prop(k12_title_fp, len(k12)),
                                "sequence_motif_scan": prop(k12_motif_fp, len(k12))},
        },
        "throughput": {"sequences": total_seqs, "seconds": round(total_time, 1),
                       "sequences_per_second": round(total_seqs / total_time, 1)},
    }
    with open(os.path.join(OUT, "amr_benchmark.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: summary[k] for k in ("positives", "negatives", "legacy_baseline", "throughput")}, indent=1)[:6000])


if __name__ == "__main__":
    main()
