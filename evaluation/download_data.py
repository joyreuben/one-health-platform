"""Download the public datasets used by the evaluations into evaluation/data/.

    python evaluation/download_data.py

Sources (all public, fetched without credentials):
  - ResFinder database (DTU, Apache-2.0): acquired resistance gene sequences by drug class.
  - Escherichia coli K-12 MG1655 coding sequences (NCBI RefSeq NC_000913.3): negative controls.
  - CONSTRAINT-2021 COVID-19 Fake News dataset (Patwa et al.): labelled real/fake health posts.
  - WHO Disease Outbreak News (WHO public API): real outbreak reports.
"""
import io
import json
import os
import sys
import zipfile

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
UA = {"User-Agent": "OneHealthPlatform-evaluation/1.0"}


def get(url, **kw):
    r = requests.get(url, headers=UA, timeout=120, **kw)
    r.raise_for_status()
    return r


def resfinder():
    out = os.path.join(DATA, "resfinder_db")
    if os.path.isdir(out) and any(f.endswith(".fsa") for f in os.listdir(out)):
        return "already present"
    r = get("https://bitbucket.org/genomicepidemiology/resfinder_db/get/master.zip")
    os.makedirs(out, exist_ok=True)
    n = 0
    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
        for name in z.namelist():
            base = os.path.basename(name)
            if base.endswith(".fsa") or base in ("phenotypes.txt", "VERSION", "README.md", "LICENSE"):
                with open(os.path.join(out, base), "wb") as f:
                    f.write(z.read(name))
                n += 1
        top = z.namelist()[0].split("/")[0]
    with open(os.path.join(out, "SOURCE.txt"), "w", encoding="utf-8") as f:
        f.write(f"https://bitbucket.org/genomicepidemiology/resfinder_db (archive folder: {top})\n")
    return f"{n} files ({top})"


def ecoli_k12():
    out = os.path.join(DATA, "ecoli_k12_mg1655_cds.fna")
    if os.path.exists(out) and os.path.getsize(out) > 1_000_000:
        return "already present"
    r = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params={"db": "nuccore", "id": "NC_000913.3", "rettype": "fasta_cds_na", "retmode": "text"})
    with open(out, "w", encoding="utf-8") as f:
        f.write(r.text)
    return f"{r.text.count('>')} CDS"


def constraint():
    base = "https://raw.githubusercontent.com/diptamath/covid_fake_news/main/data/"
    done = []
    for name in ("Constraint_Train.csv", "Constraint_Val.csv", "english_test_with_labels.csv"):
        out = os.path.join(DATA, "constraint_" + name.lower())
        if not (os.path.exists(out) and os.path.getsize(out) > 10_000):
            with open(out, "wb") as f:
                f.write(get(base + name).content)
        done.append(f"{name}: {os.path.getsize(out) // 1024} KB")
    return "; ".join(done)


def who_don(n=150):
    out = os.path.join(DATA, "who_don.json")
    if os.path.exists(out) and os.path.getsize(out) > 100_000:
        return "already present"
    items, skip = [], 0
    while len(items) < n:
        page = get("https://www.who.int/api/news/diseaseoutbreaknews",
                   params={"$orderby": "PublicationDate desc", "$top": "50", "$skip": str(skip)}).json().get("value", [])
        if not page:
            break
        items += page
        skip += len(page)
    with open(out, "w", encoding="utf-8") as f:
        json.dump({"value": items[:n]}, f, ensure_ascii=False)
    return f"{len(items[:n])} reports, {items[-1]['PublicationDate'][:10]} to {items[0]['PublicationDate'][:10]}"


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    failed = False
    for label, fn in [("ResFinder database", resfinder), ("E. coli K-12 CDS", ecoli_k12),
                      ("CONSTRAINT fake-news dataset", constraint), ("WHO Disease Outbreak News", who_don)]:
        try:
            print(f"{label}: {fn()}")
        except Exception as exc:
            failed = True
            print(f"{label}: FAILED - {exc}")
    sys.exit(1 if failed else 0)
