"""Merge the literature-mined rungs into the compartment ladder (instruction 82).

Inputs: the mapped literature rows of both mining passes -- abstracts (`candidates_mapped.tsv`) and open-access
full text (`fulltext_candidates_mapped.tsv`) -- each a verbatim quote with its PMID or PMCID.
Rules: a gene absent from the review joins the rung of its quote; a review gene on the PVM takes the side its quote
states, because the review recorded only host-facing yes/no; excluded rows are left out. Every change is listed.

    python scripts/merge_literature_ladder.py
"""
from __future__ import annotations

import os

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIT = os.path.join(REPO, "results", "literature_set_candidates_2026_11_08", "candidates_mapped.tsv")
FULLTEXT = os.path.join(REPO, "results", "literature_set_candidates_2026_11_08", "fulltext_candidates_mapped.tsv")
OUT = os.path.join(REPO, "results", "literature_set_candidates_2026_11_08")
RUNG = {"moving_junction": "1 injected at the moving junction",
        "pvm_vacuole_facing": "2 PVM, vacuole-facing",
        "pvm_host_facing": "3 PVM, host-facing",
        "host_nucleus": "6 host nucleus",
        "host_organelle_or_pm": "7 host organelle or plasma membrane",
        "vacuole_resident": "4 vacuole-resident (hard negative)",
        "host_cytosol": "5 exported to the host cytosol"}


def _mapped(path: str, column: str = "pmids") -> pd.DataFrame:
    """The rows of one mined file that have a gene and no exclusion, one row per gene."""
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    frame = frame[(frame.mapped_gene_id != "") & (frame.excluded == "")]
    agg = frame.groupby("mapped_gene_id").agg(protein=("protein", "first"), category=("category", "first"),
                                              pmids=(column, lambda s: ";".join(map(str, s))))
    return agg.reset_index().rename(columns={"mapped_gene_id": "gene_id"})


def main() -> None:
    lit = pd.concat([_mapped(LIT), _mapped(FULLTEXT)], ignore_index=True)
    lit["rung"] = lit.category.map(RUNG)
    # A gene quoted on the PVM in one paper and in the lumen or on the IVN in another cannot take one rung:
    # the two statements disagree about where it is, and the ladder needs one place per gene.
    disagree = lit.groupby("gene_id").rung.nunique()
    conflicted = set(disagree[disagree > 1].index)
    lit["conflict"] = lit.gene_id.isin(conflicted)
    lit[lit.conflict].to_csv(os.path.join(OUT, "literature_conflicts.tsv"), sep="\t", index=False)
    lit = lit[~lit.conflict].drop(columns="conflict").drop_duplicates("gene_id")
    lit.to_csv(os.path.join(OUT, "literature_rungs.tsv"), sep="\t", index=False)
    print(lit.rung.value_counts().sort_index().to_string())
    print(len(conflicted), "genes left out: two papers place them differently")


if __name__ == "__main__":
    main()
