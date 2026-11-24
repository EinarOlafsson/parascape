"""Stage 1 batches for the free Toxoplasma Affinage run: papers grouped for reading.

Each batch is one JSON file of papers. Each paper carries its abstract and the ParaScape
genes it names (gene id, ToxoDB name, product). A reader follows READING_RULES.md and
writes `findings/<batch>.jsonl`. A batch whose findings file exists is done; a rerun
skips it.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def main(argv=None) -> int:
    """Write the batch files."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--work", required=True, help="the Stage 0 output folder")
    ap.add_argument("--identity", required=True, help="toxodb_identity.tsv of the release")
    ap.add_argument("--per-batch", type=int, default=60)
    args = ap.parse_args(argv)
    work = Path(args.work)
    gp = pd.read_parquet(work / "gene_papers.parquet")
    ident = (pd.read_csv(args.identity, sep="\t", dtype=str, keep_default_na=False)
             .drop_duplicates("gene_id").set_index("gene_id"))
    papers = {}
    with open(work / "abstracts.jsonl") as f:
        for line in f:
            r = json.loads(line)
            if r["abstract"].strip():
                papers[r["pmid"]] = r
    genes = gp[gp.pmid.isin(papers)].groupby("pmid").gene_id.apply(lambda s: sorted(set(s)))
    # In order of year, then PubMed id
    order = sorted(genes.index, key=lambda p: (papers[p]["year"], p))
    out = work / "batches"; out.mkdir(exist_ok=True); (work / "findings").mkdir(exist_ok=True)
    n = 0
    for b in range(0, len(order), args.per_batch):
        items = []
        for p in order[b:b + args.per_batch]:
            r = papers[p]
            items.append({"pmid": p, "year": r["year"], "journal": r["journal"], "title": r["title"],
                          "abstract": r["abstract"],
                          "genes": [{"gene_id": g,
                                     "name": ident.gene_name.get(g, "") if g in ident.index else "",
                                     "product": ident["product"].get(g, "") if g in ident.index else ""}
                                    for g in genes[p]]})
        (out / f"batch_{n:03d}.json").write_text(json.dumps(items, indent=1))
        n += 1
    print(json.dumps({"papers": len(order), "batches": n, "per_batch": args.per_batch}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
