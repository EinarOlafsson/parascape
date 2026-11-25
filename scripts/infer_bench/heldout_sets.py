"""Held-out gene sets of the Infer benchmark, after Affinage (DESIGN.md, "Benchmark changes", 2026-11-24).

A score averaged over all held-out genes is carried by the easy ones. Affinage scores its system on a
set enriched for where it fails, beside controls matched for study attention. Here:

- **paralog**: a gene whose orthogroup holds another Toxoplasma gene (shares data and often the label
  with its copy, so a random fold flatters it);
- **sparse**: the tenth of genes with the fewest measured values over every view;
- **attention stratum**: papers naming the gene (focal + substantive + incidental): none, 1-2, 3-10,
  more than 10. Labels come mostly from well-studied genes, so a score is reported within each
  stratum, and a control set is drawn to match the strata of the failure-enriched set.

`failure_enriched` = paralog or sparse. "No papers" is a stratum, not part of the set: 5,574 genes
have none, so a set that included them would be most of the genome (a first version, counting only
focal and substantive papers, put 7,588 of 8,140 genes in it). Writes `heldout_sets.tsv` (one row per gene).

    python scripts/infer_bench/heldout_sets.py OUT_DIR
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

NODES = "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots/nodes.parquet"
GENES = "/media/carruthers/mnt2/parascape/maxproj_inputs/genes.parquet"
STRATA = ((0, 0, "none"), (1, 2, "1-2"), (3, 10, "3-10"), (11, np.inf, "more than 10"))
SPARSE_QUANTILE = 0.10


def attention_stratum(papers: pd.Series) -> pd.Series:
    """The attention stratum of each gene from its paper count."""
    out = pd.Series("none", index=papers.index, dtype=object)
    for lo, hi, name in STRATA:
        out[(papers >= lo) & (papers <= hi)] = name
    return out


def sets(nodes: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """One row per gene with its paralog, sparse and attention flags."""
    t = nodes[["gene_id", "orthogroup"]].copy()
    og = t["orthogroup"].astype(str)
    t["paralog"] = og.notna() & ~og.isin(["", "nan", "None"]) & (t.groupby("orthogroup")["gene_id"].transform("size") > 1)
    share = genes.set_index("gene_id").notna().mean(axis=1)
    t["measured_share"] = t["gene_id"].map(share).fillna(0.0)
    t["sparse"] = t["measured_share"] <= t["measured_share"].quantile(SPARSE_QUANTILE)
    papers = (nodes["n_papers_focal"].fillna(0) + nodes["n_papers_substantive"].fillna(0)
              + nodes["n_papers_incidental"].fillna(0))
    t["papers"] = papers.astype(int)
    t["attention"] = attention_stratum(papers)
    t["failure_enriched"] = t["paralog"] | t["sparse"]
    return t


def matched_control(t: pd.DataFrame, chosen: pd.Series, pool: pd.Series, seed: int = 0) -> pd.Series:
    """Genes from `pool` (not chosen) drawn to match the attention strata of the `chosen` genes."""
    rng = np.random.default_rng(seed)
    want = t.loc[chosen, "attention"].value_counts()
    picked = []
    for stratum, n in want.items():
        cand = t.index[pool & ~chosen & (t["attention"] == stratum)]
        picked += list(rng.choice(cand, min(n, len(cand)), replace=False))
    out = pd.Series(False, index=t.index)
    out[picked] = True
    return out


def main(argv=None) -> int:
    """Write the gene sets and print their sizes."""
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        raise SystemExit(__doc__)
    os.makedirs(args[0], exist_ok=True)
    t = sets(pd.read_parquet(NODES), pd.read_parquet(GENES))
    t.to_csv(os.path.join(args[0], "heldout_sets.tsv"), sep="\t", index=False)
    print(t[["paralog", "sparse", "failure_enriched"]].sum().to_string())
    print(t["attention"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
