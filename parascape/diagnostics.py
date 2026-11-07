"""Tell a finding from an artefact: missingness and dominance warnings (instruction 60, change 7).

Two measured checks, each returning the numbers and a plain sentence only when the numbers call for it:

* **Missingness** -- are the genes a result names much emptier in the evidence it was built from than
  genes in general? A map places a gene with little measured on it from whatever filled the gaps, so a
  named set that is mostly gaps rests on the filling, not on the measurements.
* **Dominance** -- does one evidence family carry most of a model's importance? Then the answer is that
  one dataset's, however many families were offered, and a reader should know which.

No window here; every rule is testable.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

#: Warn when the named genes miss this much more of the evidence than genes in general (absolute).
MISSING_GAP = 0.15
#: Warn when one family carries at least this share of the positive importance.
DOMINANT_SHARE = 0.5


def missingness(nodes: pd.DataFrame, columns, genes, id_column: str = "gene_id") -> dict:
    """Share of evidence cells missing for the named genes, against all genes, on `columns`."""
    columns = [c for c in columns if c in nodes.columns]
    if not columns or not len(nodes):
        return {}
    table = nodes.set_index(id_column)[columns] if id_column in nodes.columns else nodes[columns]
    named = table.reindex([g for g in map(str, genes) if g in set(map(str, table.index))])
    if not len(named):
        return {}
    overall = float(table.isna().to_numpy().mean())
    mine = float(named.isna().to_numpy().mean())
    gap = mine - overall
    out = {"columns": len(columns), "named": int(len(named)), "named_missing": mine,
           "all_missing": overall, "gap": gap, "warn": gap >= MISSING_GAP}
    if out["warn"]:
        out["sentence"] = (f"The named genes are missing {mine:.0%} of the evidence the map was built "
                           f"from, against {overall:.0%} for genes in general: their places rest on "
                           "filled-in values more than on measurements.")
    return out


def dominance(importance: pd.DataFrame, family: str = "family",
              value: str = "importance_mean") -> dict:
    """The share of positive importance the strongest family carries, from a harness importance table."""
    if importance is None or not len(importance) or value not in importance:
        return {}
    positive = importance[importance[value] > 0]
    total = float(positive[value].sum())
    if total <= 0:
        return {}
    top = positive.sort_values(value, ascending=False).iloc[0]
    share = float(top[value]) / total
    out = {"family": str(top[family]), "share": share, "families": int(len(positive)),
           "warn": share >= DOMINANT_SHARE and len(positive) > 1}
    if out["warn"]:
        out["sentence"] = (f"One evidence family, {top[family]}, carries {share:.0%} of what the model "
                           f"relies on among {len(positive)}: the answer is largely that dataset's.")
    return out
