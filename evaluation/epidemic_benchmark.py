"""Benchmark of the platform's epidemic-intelligence components.

    python evaluation/download_data.py
    python evaluation/epidemic_benchmark.py                 # offline components only
    python evaluation/epidemic_benchmark.py --with-claude   # also the Claude stages (needs ANTHROPIC_API_KEY, costs money)

Part A  Misinformation screen (offline heuristic) on CONSTRAINT-2021: 10,700 COVID-19
        posts labelled real / fake. "Predicted fake" = the platform would block the
        article (verdict FAKE or LIKELY_FAKE).
Part B  The same screen on real WHO Disease Outbreak News reports: how often would a
        genuine official report be wrongly blocked?
Part C  Outbreak extraction on WHO reports. The WHO title ("Disease - Country") is the
        answer key; the title is removed from the text given to the extractor.
        Offline: the keyword fallback. With --with-claude: the Claude extractor too.

With --with-claude, Parts A and B are also run through the full credibility check
(heuristic + Claude) on a random sample, and Part C through Claude on a sample.

Outputs: evaluation/results/epidemic_benchmark.json and CSV files alongside it.
"""
import argparse
import csv
import json
import logging
import math
import os
import random
import re
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
import feeds  # noqa: E402
from engine import ArticleCredibilityChecker as Cred, EpidemicDetectionEngine as Epi  # noqa: E402

SEED = 20261115
BLOCKING = ("FAKE", "LIKELY_FAKE")


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def prop(k, n):
    lo, hi = wilson(k, n)
    return {"k": k, "n": n, "pct": round(100 * k / n, 2) if n else None,
            "ci95": [round(100 * lo, 2), round(100 * hi, 2)]}


def binary_metrics(rows):
    """rows: (is_fake: bool, predicted_fake: bool)."""
    tp = sum(1 for t, p in rows if t and p)
    fn = sum(1 for t, p in rows if t and not p)
    fp = sum(1 for t, p in rows if not t and p)
    tn = sum(1 for t, p in rows if not t and not p)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    return {
        "n": len(rows), "confusion": {"fake_blocked": tp, "fake_passed": fn, "real_blocked": fp, "real_passed": tn},
        "sensitivity_fake_blocked": prop(tp, tp + fn),
        "specificity_real_passed": prop(tn, tn + fp),
        "precision_of_block": prop(tp, tp + fp),
        "accuracy": prop(tp + tn, len(rows)),
        "f1_fake": round(2 * prec * rec / (prec + rec), 4) if prec + rec else 0.0,
    }


# --------------------------------------------------------------- data
def load_constraint():
    rows = []
    for split, name in (("train", "constraint_constraint_train.csv"), ("val", "constraint_constraint_val.csv"),
                        ("test", "constraint_english_test_with_labels.csv")):
        with open(os.path.join(DATA, name), encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("tweet") and r.get("label") in ("real", "fake"):
                    rows.append({"split": split, "id": f"{split}-{r['id']}", "text": r["tweet"],
                                 "fake": r["label"] == "fake"})
    return rows


_STOP = {"disease", "diseases", "virus", "viral", "caused", "by", "infection", "infections", "human", "humans",
         "in", "the", "of", "and", "a", "an", "outbreak", "cases", "case", "update", "situation", "due", "to",
         "with", "associated", "novel", "new", "suspected", "confirmed", "type"}
_SYNONYMS = [  # (pattern, canonical)
    (r"middle east respiratory syndrome( coronavirus)?|mers cov|\bmers\b", "mers"),
    (r"monkeypox", "mpox"),
    (r"coronavirus disease( 2019)?|covid 19|sars cov 2", "covid"),
    (r"avian influenza|bird flu", "avian influenza"),
    (r"marburg virus disease", "marburg"),
    (r"ebola virus disease|ebola disease", "ebola"),
    (r"democratic republic of the congo|\bdrc\b|dr congo", "drc"),
    (r"united states of america|\busa\b|united states|\bu s\b", "usa"),
    (r"united kingdom of great britain and northern ireland|united kingdom|\buk\b", "uk"),
    (r"kingdom of saudi arabia", "saudi arabia"),
    (r"russian federation", "russia"),
    (r"united republic of tanzania", "tanzania"),
]


def norm(s: str) -> str:
    s = re.sub(r"[^a-z0-9 ]+", " ", (s or "").lower())
    s = re.sub(r"\s+", " ", s).strip()
    for pat, canon in _SYNONYMS:
        s = re.sub(pat, canon, s)
    return s


def sig_tokens(s: str):
    return {t for t in norm(s).split() if t not in _STOP}


def disease_match(pred: str, truth: str) -> bool:
    p, t = sig_tokens(pred), sig_tokens(truth)
    if not p or not t or norm(pred) in ("unknown", ""):
        return False
    small, big = (p, t) if len(p) <= len(t) else (t, p)
    return small <= big


def location_match(pred: str, truth: str) -> bool:
    p = norm(pred)
    if not p or p in ("unknown", "n a", "none"):
        return False
    places = [norm(x) for x in re.split(r",| and |;|/", truth) if norm(x)]
    return any(pl in p or p in pl for pl in places)


def split_title(title: str):
    """'Disease - Country' -> (disease, location or None)."""
    # WHO writes the separator as " - ", " – ", "- ", "–" or even "x)-Country".
    # A plain hyphen inside a word ("vaccine-derived", "COVID-19") is not a separator.
    parts = re.split(r"\s*[–—]\s*|\s+-\s*|-\s+|(?<=[a-z\)])-(?=[A-Z])", title, maxsplit=1)
    if len(parts) == 2 and parts[0].strip() and parts[1].strip():
        return parts[0].strip(), parts[1].strip()
    if "," in title:
        a, b = title.split(",", 1)
        return a.strip(), b.strip()
    return title.strip(), None


_NON_SPECIFIC = re.compile(r"global|multi|region|hemisphere|situation|update|worldwide|countries", re.IGNORECASE)


def specific_location(loc) -> bool:
    """True when the title names a country or place, not 'Global situation' / 'Multi-country' / a WHO region."""
    return bool(loc) and not _NON_SPECIFIC.search(loc)


def load_don():
    with open(os.path.join(DATA, "who_don.json"), encoding="utf-8") as f:
        raw = json.load(f)["value"]
    items, seen = [], set()
    for it in raw:
        key = it.get("DonId") or it.get("UrlName") or it.get("Id")
        if key in seen:
            continue
        seen.add(key)
        title = (it.get("OverrideTitle") if it.get("UseOverrideTitle") and it.get("OverrideTitle") else it.get("Title")) or ""
        body = "\n\n".join(x for x in (feeds.html_to_text(it.get("Summary", "")),
                                       feeds.html_to_text(it.get("Overview", "")),
                                       feeds.html_to_text(it.get("Epidemiology", ""))) if x)
        disease, location = split_title(title)
        items.append({"id": it.get("DonId") or it.get("UrlName") or it.get("Id"), "title": title.strip(),
                      "published": (it.get("PublicationDate") or "")[:10], "truth_disease": disease,
                      "truth_location": location, "body": body[:feeds.MAX_TEXT_CHARS],
                      "full_text": (title + "\n\n" + body)[:feeds.MAX_TEXT_CHARS]})
    return [i for i in items if len(i["body"]) > 200]


def verdict_counts(verdicts):
    out = {}
    for v in verdicts:
        out[v] = out.get(v, 0) + 1
    return dict(sorted(out.items()))


def write_csv(name, rows):
    if not rows:
        return
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


# --------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--with-claude", action="store_true", help="also evaluate the Claude stages (billed to your API key)")
    ap.add_argument("--claude-posts", type=int, default=300, help="CONSTRAINT posts for the Claude credibility check")
    ap.add_argument("--claude-reports", type=int, default=60, help="WHO reports for the Claude stages")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = random.Random(SEED)
    summary = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "platform_version": engine.__version__,
               "random_seed": SEED, "claude_evaluated": False}

    # ---- Part A: heuristic screen on CONSTRAINT
    posts = load_constraint()
    a_rows = []
    for p in posts:
        h = Cred._heuristic_check(p["text"])
        a_rows.append({"id": p["id"], "split": p["split"], "label": "fake" if p["fake"] else "real",
                       "verdict": h["verdict"], "credibility_score": round(h["credibility_score"], 3),
                       "blocked": int(h["verdict"] in BLOCKING), "n_flags": len(h["flags"]), "chars": len(p["text"])})
    write_csv("credibility_constraint_per_post.csv", a_rows)

    def metrics_for(rows):
        return binary_metrics([(r["label"] == "fake", bool(r["blocked"])) for r in rows])

    summary["A_heuristic_screen_constraint"] = {
        "dataset": "CONSTRAINT-2021 COVID-19 Fake News (Patwa et al.), train+val+test; the screen is rule-based "
                   "and was not fitted to this data",
        "all": metrics_for(a_rows),
        "official_test_split": metrics_for([r for r in a_rows if r["split"] == "test"]),
        "verdicts_for_fake_posts": verdict_counts(r["verdict"] for r in a_rows if r["label"] == "fake"),
        "verdicts_for_real_posts": verdict_counts(r["verdict"] for r in a_rows if r["label"] == "real"),
        "median_post_length_chars": sorted(r["chars"] for r in a_rows)[len(a_rows) // 2],
    }

    # ---- Part B: heuristic screen on real WHO reports
    don = load_don()
    b_rows = []
    for d in don:
        h = Cred._heuristic_check(d["full_text"])
        b_rows.append({"id": d["id"], "published": d["published"], "title": d["title"], "verdict": h["verdict"],
                       "credibility_score": round(h["credibility_score"], 3), "blocked": int(h["verdict"] in BLOCKING),
                       "flags": " | ".join(h["flags"])})
    write_csv("credibility_who_reports.csv", b_rows)
    summary["B_heuristic_screen_who_reports"] = {
        "dataset": f"{len(don)} WHO Disease Outbreak News reports published "
                   f"{min(d['published'] for d in don)} to {max(d['published'] for d in don)}",
        "reports_per_publication_year": verdict_counts(d["published"][:4] for d in don),
        "wrongly_blocked": prop(sum(r["blocked"] for r in b_rows), len(b_rows)),
        "verdicts": verdict_counts(r["verdict"] for r in b_rows),
    }

    # ---- Part C: extraction, keyword fallback
    # Location is scored only where the title names a specific place.
    for d in don:
        if not specific_location(d["truth_location"]):
            d["truth_location"] = None
    with_loc = [d for d in don if d["truth_location"]]
    c_rows = []
    for d in don:
        k = Epi._demo_analysis(d["body"])
        c_rows.append({"id": d["id"], "title": d["title"], "truth_disease": d["truth_disease"],
                       "truth_location": d["truth_location"] or "",
                       "keyword_disease": k["disease"], "keyword_disease_ok": int(disease_match(k["disease"], d["truth_disease"])),
                       "keyword_location": k["location"],
                       "keyword_location_ok": int(bool(d["truth_location"]) and location_match(k["location"], d["truth_location"]))})
    summary["C_extraction_who_reports"] = {
        "answer_key": "disease and location parsed from the WHO report title; title withheld from the extractor",
        "n_reports": len(don), "n_with_specific_location_in_title": len(with_loc),
        "distinct_diseases_in_answer_key": len({norm(d["truth_disease"]) for d in don}),
        "keyword_fallback": {
            "disease_correct": prop(sum(r["keyword_disease_ok"] for r in c_rows), len(c_rows)),
            "location_correct": prop(sum(r["keyword_location_ok"] for r in c_rows), len(with_loc)),
            "note": "the fallback recognises only Ebola, COVID-19 and mpox and never extracts a location",
        },
    }

    # ---- Claude stages (opt-in)
    if args.with_claude:
        if engine.anthropic_client is None:
            print("! --with-claude given but Claude is not connected:", engine._claude_error)
        else:
            summary["claude_evaluated"] = True
            summary["claude_model"] = engine.ANTHROPIC_MODEL
            fake = [p for p in posts if p["split"] == "test" and p["fake"]]
            real = [p for p in posts if p["split"] == "test" and not p["fake"]]
            half = args.claude_posts // 2
            sample = rng.sample(fake, min(half, len(fake))) + rng.sample(real, min(half, len(real)))
            full_rows = []
            for i, p in enumerate(sample, 1):
                r = Cred.check(p["text"])
                full_rows.append({"id": p["id"], "label": "fake" if p["fake"] else "real", "verdict": r["verdict"],
                                  "method": r.get("method", ""), "credibility_score": round(r["credibility_score"], 3),
                                  "blocked": int(r["verdict"] in BLOCKING)})
                if i % 25 == 0:
                    print(f"  Claude credibility: {i}/{len(sample)}", flush=True)
            write_csv("credibility_constraint_claude_sample.csv", full_rows)
            same_heur = [r for r in a_rows if r["id"] in {x["id"] for x in sample}]
            summary["A_full_check_constraint_sample"] = {
                "sample": f"{len(sample)} posts from the official test split, balanced, seed {SEED}",
                "heuristic_plus_claude": metrics_for(full_rows),
                "heuristic_only_same_sample": metrics_for(same_heur),
                "decided_by": verdict_counts(r["method"] for r in full_rows),
            }
            reports = rng.sample(don, min(args.claude_reports, len(don)))
            ext_rows, blocked = [], 0
            for i, d in enumerate(reports, 1):
                cred = Cred.check(d["full_text"])
                blocked += cred["verdict"] in BLOCKING
                a = Epi.analyze_article(d["body"])
                ext_rows.append({"id": d["id"], "title": d["title"], "truth_disease": d["truth_disease"],
                                 "truth_location": d["truth_location"] or "", "source": a["source"],
                                 "claude_disease": a["disease"], "disease_ok": int(disease_match(a["disease"], d["truth_disease"])),
                                 "claude_location": a["location"],
                                 "location_ok": int(bool(d["truth_location"]) and location_match(a["location"], d["truth_location"])),
                                 "cases": a["cases"], "trend": a["trend"], "severity": a["severity"],
                                 "risk_score": a["risk_score"], "credibility_verdict": cred["verdict"]})
                if i % 10 == 0:
                    print(f"  Claude extraction: {i}/{len(reports)}", flush=True)
            write_csv("extraction_who_reports_claude_sample.csv", ext_rows)
            ok_src = [r for r in ext_rows if r["source"] == "CLAUDE_API"]
            loc_rows = [r for r in ok_src if r["truth_location"]]
            summary["C_extraction_who_reports"]["claude"] = {
                "sample": f"{len(reports)} reports, seed {SEED}", "api_failures": len(ext_rows) - len(ok_src),
                "disease_correct": prop(sum(r["disease_ok"] for r in ok_src), len(ok_src)),
                "location_correct": prop(sum(r["location_ok"] for r in loc_rows), len(loc_rows)),
                "matching": "automatic lenient string match; review extraction_who_reports_claude_sample.csv by hand",
            }
            summary["B_full_check_who_reports_sample"] = {"wrongly_blocked": prop(blocked, len(reports))}

    write_csv("extraction_who_reports_keyword.csv", c_rows)
    with open(os.path.join(OUT, "epidemic_benchmark.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
