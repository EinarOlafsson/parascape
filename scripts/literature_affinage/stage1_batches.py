"""Stage 1 batches for the free Toxoplasma Affinage run: papers grouped for reading.

Each batch is one JSON file of papers. Each paper carries its abstract and the ParaScape
genes it names (gene id, ToxoDB name, product). A reader follows READING_RULES.md and
writes `findings/<batch>.jsonl`. A batch whose findings file exists is done; a rerun
skips it.

A paper lists EVERY gene that its title or abstract names (from the release `mentions.parquet`),
not only the genes for which Stage 0 selected it. The first run listed only the latter, so a paper
picked for one gene was read without the well-studied genes it also names: ROP18, MIC2, AMA1 and
GRA2 had more than 20 papers each and were missing from 976 gene-paper pairs in 564 papers that
had been read (the readers reported 214 of these as "unlisted protein" doubts, 2026-11-24).

`--supplement` writes `supp_<n>.json` batches for papers that already have findings, listing
only the genes the first pass left out; `already_read` names the genes read before, for context.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

MAX_GENES_PER_PAPER = 50     # the same catalogue cut as Stage 0


def named_in_abstract(mentions: pd.DataFrame, pmids) -> pd.DataFrame:
    """(gene_id, pmid) pairs where the title or abstract names the gene, for papers of at most 50 genes."""
    m = mentions.dropna(subset=["pmid"]).copy()
    m["pmid"] = m["pmid"].astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    m = m[m.pmid.str.fullmatch(r"\d+")]
    fan = m.groupby("pmid").gene_id.nunique()
    m = m[m.pmid.map(fan) <= MAX_GENES_PER_PAPER]
    m = m[m.pmid.isin(set(pmids)) & m.section.isin(["title", "abstract"])]
    return m[["gene_id", "pmid"]].drop_duplicates().reset_index(drop=True)


def gene_lists(gp: pd.DataFrame, named: pd.DataFrame, pmids) -> pd.Series:
    """pmid -> sorted gene ids: the Stage 0 selection plus every gene the abstract names."""
    both = pd.concat([gp[["gene_id", "pmid"]].astype(str), named.astype(str)])
    both = both[both.pmid.isin(set(pmids))]
    return both.groupby("pmid").gene_id.apply(lambda s: sorted(set(s)))


def _gene(ident, g):
    return {"gene_id": g,
            "name": ident.gene_name.get(g, "") if g in ident.index else "",
            "product": ident["product"].get(g, "") if g in ident.index else ""}


def main(argv=None) -> int:
    """Write the batch files."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--work", required=True, help="the Stage 0 output folder")
    ap.add_argument("--identity", required=True, help="toxodb_identity.tsv of the release")
    ap.add_argument("--mentions", required=True, help="mentions.parquet of the release")
    ap.add_argument("--per-batch", type=int, default=60)
    ap.add_argument("--supplement", action="store_true",
                    help="write supp_<n>.json for genes the first pass did not list")
    args = ap.parse_args(argv)
    work = Path(args.work)
    gp = pd.read_parquet(work / "gene_papers.parquet")
    gp["pmid"] = gp.pmid.astype(str)
    ident = (pd.read_csv(args.identity, sep="\t", dtype=str, keep_default_na=False)
             .drop_duplicates("gene_id").set_index("gene_id"))
    papers = {}
    with open(work / "abstracts.jsonl") as f:
        for line in f:
            r = json.loads(line)
            if r["abstract"].strip():
                papers[r["pmid"]] = r
    named = named_in_abstract(pd.read_parquet(args.mentions), papers)
    genes = gene_lists(gp, named, papers)
    out = work / "batches"; out.mkdir(exist_ok=True); (work / "findings").mkdir(exist_ok=True)

    if args.supplement:
        listed = {}
        for b in sorted(out.glob("batch_*.json")):
            for it in json.loads(b.read_text()):
                listed[it["pmid"]] = {g["gene_id"] for g in it["genes"]}
        todo = {p: [g for g in genes[p] if g not in listed[p]] for p in genes.index if p in listed}
        todo = {p: gs for p, gs in todo.items() if gs}
        prefix, extra = "supp", lambda p: {"already_read": sorted(listed[p])}
    else:
        # The audit reads each findings file against ITS batch file; rewriting a batch that has
        # findings would make every one of them fail A1.
        if any((work / "findings").glob("batch_*.jsonl")):
            sys.exit("findings exist for these batches; use --supplement, or start in a new --work folder")
        todo = dict(genes.items())
        prefix, extra = "batch", lambda p: {}

    order = sorted(todo, key=lambda p: (papers[p]["year"], p))
    n = 0
    for b in range(0, len(order), args.per_batch):
        items = []
        for p in order[b:b + args.per_batch]:
            r = papers[p]
            items.append({"pmid": p, "year": r["year"], "journal": r["journal"], "title": r["title"],
                          "abstract": r["abstract"], "genes": [_gene(ident, g) for g in todo[p]],
                          **extra(p)})
        (out / f"{prefix}_{n:03d}.json").write_text(json.dumps(items, indent=1))
        n += 1
    print(json.dumps({"papers": len(order), "pairs": sum(map(len, todo.values())),
                      "batches": n, "prefix": prefix, "per_batch": args.per_batch}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
