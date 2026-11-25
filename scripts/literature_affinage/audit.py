"""Fixed rule checks on the reading output, no model (after the Affinage R1-R10 audit).

A1 the finding gene is one of that paper's listed genes;
A2 every partner gene id is a listed gene of that paper;
A3 kind, status, method tier and compartment take allowed values;
A4 the pmid is in the batch;
A5 the finding sentence names the gene (its ToxoDB name, id number, or a product word);
A7 the subject of the finding sentence (its first word or two) appears in the title or abstract,
   or the gene ToxoDB name or id number does (added after the ATP synthase case, 2026-11-24);
A6 per batch: findings or doubts reach the last quarter of the batch (a reader that stopped
   halfway left batch 018 with nothing after paper 29 of 60, 2026-11-24).
Prints the counts per rule and writes `audit.tsv`; a finding that fails A1-A4 is not loaded.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

KINDS = {"localisation", "partner", "complex", "activity", "phenotype", "host_facing", "structure"}
COMPARTMENTS = {"", "19S proteasome", "20S proteasome", "40S ribosome", "60S ribosome", "ER", "Golgi", "IMC",
                "PM - integral", "PM - peripheral", "apical", "apicoplast", "cytosol", "dense granules",
                "endomembrane vesicles", "micronemes", "mitochondrion - membranes", "mitochondrion - soluble",
                "nucleolus", "nucleus - chromatin", "nucleus - non-chromatin", "rhoptries", "tubulin cytoskeleton",
                "PV lumen", "PVM", "host cytosol", "host nucleus", "conoid", "basal complex", "apical annuli",
                "residual body", "cyst wall"}


def findings_files(work: Path) -> list:
    """The findings of the first pass (batch_NNN) and of the supplement pass (supp_NNN)."""
    return sorted(f for pat in ("batch_[0-9][0-9][0-9].jsonl", "supp_[0-9][0-9][0-9].jsonl")
                  for f in (work / "findings").glob(pat))


def audit(work: Path) -> pd.DataFrame:
    """One row per finding with the rules it fails."""
    rows = []
    for f in findings_files(work):
        batch = {p["pmid"]: p for p in json.load(open(work / "batches" / f.name.replace(".jsonl", ".json")))}
        for i, line in enumerate(open(f)):
            if not line.strip():
                continue
            r = json.loads(line)
            paper = batch.get(str(r.get("pmid")))
            listed = {g["gene_id"]: g for g in paper["genes"]} if paper else {}
            fails = []
            if paper is None:
                fails.append("A4")
            if r.get("gene_id") not in listed:
                fails.append("A1")
            # A supplement paper lists only the genes the first pass missed; its partners may be among
            # the genes read before.
            known = set(listed) | set((paper or {}).get("already_read", []))
            if any(p not in known for p in r.get("partner_gene_ids", []) or []):
                fails.append("A2")
            if (r.get("kind") not in KINDS or r.get("status") not in {"established", "refuted"}
                    or r.get("method_tier") not in (1, 2, 3) or r.get("compartment", "") not in COMPARTMENTS):
                fails.append("A3")
            g = listed.get(r.get("gene_id"), {})
            words = {w.lower() for w in re.findall(r"[A-Za-z0-9]{3,}", g.get("product", ""))} - {"protein", "domain", "containing", "family", "putative"}
            text = r.get("finding", "").lower()
            name = (g.get("name") or "").lower()
            if not ((name and name != "n/a" and name in text) or r.get("gene_id", "x").split("_")[-1] in text
                    or any(w in text for w in words)):
                fails.append("A5")
            ab = ((paper or {}).get("title", "") + " " + (paper or {}).get("abstract", "")).lower()
            subj = re.sub(r"^(the|a|an|deletion of|loss of|depletion of|disruption of|knockout of)\s+", "",
                          text.strip()).split(" ")
            heads = [w.strip(".,;:()") for w in subj[:2] if len(w.strip(".,;:()")) >= 3]
            flat = re.sub(r"[^a-z0-9]", "", ab)                   # AP2XII-2 and AP2XII2 are one name
            squash = lambda w: re.sub(r"[^a-z0-9]", "", w)
            if not ((name and name != "n/a" and squash(name) in flat) or r.get("gene_id", "x").split("_")[-1] in ab
                    or any(squash(h) and squash(h) in flat for h in heads)):
                fails.append("A7")
            rows.append({"batch": f.stem, "line": i, "pmid": r.get("pmid"), "gene_id": r.get("gene_id"),
                         "kind": r.get("kind"), "fails": ",".join(fails), "finding": r.get("finding", "")})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    work = Path(sys.argv[1])
    a = audit(work)
    a.to_csv(work / "audit.tsv", sep="\t", index=False)
    print(f"findings {len(a)} in {a.batch.nunique()} batches")
    for rule in ("A1", "A2", "A3", "A4", "A5", "A7"):
        print(rule, int(a.fails.str.contains(rule).sum()))
    print(a[a.fails != ""].head(12).to_string(max_colwidth=90))
    for f in findings_files(work):
        papers = [p["pmid"] for p in json.load(open(work / "batches" / f.name.replace(".jsonl", ".json")))]
        seen = {json.loads(l)["pmid"] for g in (f, f.with_suffix(".doubts.jsonl")) if g.exists()
                for l in open(g) if l.strip()}
        last = max([i for i, p in enumerate(papers) if p in seen], default=-1)
        if last < 0.75 * len(papers):
            print(f"A6 {f.stem}: nothing after paper {last} of {len(papers)} -- reread the batch")
