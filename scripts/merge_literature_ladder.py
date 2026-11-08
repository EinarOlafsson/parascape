"""Merge the literature-mined rungs into the compartment ladder (instruction 82).

Inputs: the graded review (results/host_facing_truth_2026_11_04) and the mapped literature rows
(results/literature_set_candidates_2026_11_08/candidates_mapped.tsv, each a verbatim abstract quote with a PMID).
Rules: a gene absent from the review joins the rung of its quote; a review gene on the PVM takes the side its quote
states, because the review recorded only host-facing yes/no; excluded rows are left out. Every change is listed.

    python scripts/merge_literature_ladder.py
"""
from __future__ import annotations

import os

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIT = os.path.join(REPO, "results", "literature_set_candidates_2026_11_08", "candidates_mapped.tsv")
OUT = os.path.join(REPO, "results", "literature_set_candidates_2026_11_08")
RUNG = {"moving_junction": "1 injected at the moving junction",
        "pvm_vacuole_facing": "2 PVM, vacuole-facing",
        "pvm_host_facing": "3 PVM, host-facing",
        "host_cytosol": "5 exported to the host cytosol",
        "host_nucleus": "6 host nucleus",
        "host_organelle_or_pm": "7 host organelle or plasma membrane"}


def main() -> None:
    lit = pd.read_csv(LIT, sep="\t")
    lit = lit[lit.mapped_gene_id.notna() & lit.excluded.isna()]
    lit = lit.groupby("mapped_gene_id").agg(protein=("protein", "first"), category=("category", "first"),
                                            pmids=("pmids", lambda s: ";".join(map(str, s)))).reset_index()
    lit["rung"] = lit.category.map(RUNG)
    lit = lit.rename(columns={"mapped_gene_id": "gene_id"})
    lit.to_csv(os.path.join(OUT, "literature_rungs.tsv"), sep="\t", index=False)
    print(lit.rung.value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
