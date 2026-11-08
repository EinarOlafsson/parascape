"""Measured redundancy between columns, reported -- never used to exclude (instruction 73; Einar 2026-11-07).

`data/measured_groups_tg.tsv` comes from `scripts/build_redundancy_report_notebook.py`: Toxoplasma numeric columns
grouped so that every pair in a group has |Spearman| >= 0.8. Einar chose "report only": an input in the same
measured group as the label but of another information category is named as a warning on the question; the
category rule alone decides what is held out.
"""
from __future__ import annotations

import os
from functools import lru_cache

from . import organisms
import pandas as pd

_FILES = {organisms.TOXOPLASMA: "measured_groups_tg.tsv"}


@lru_cache(maxsize=None)
def groups(organism: str = organisms.TOXOPLASMA) -> pd.DataFrame:
    """column, group, category for every column the report measured (empty when there is no report)."""
    name = _FILES.get(organism)
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", name) if name else ""
    return pd.read_csv(path, sep="\t") if path and os.path.exists(path) else pd.DataFrame(
        columns=["column", "group", "category"])


def overlap(label_column: str, inputs, organism: str = organisms.TOXOPLASMA, table: pd.DataFrame | None = None) -> dict:
    """{input: reason} for inputs in the label column's measured group but of a different category."""
    g = groups(organism) if table is None else table
    if not len(g) or label_column not in set(g["column"]):
        return {}
    row = g[g["column"] == label_column].iloc[0]
    mates = g[(g["group"] == row["group"]) & (g["column"] != label_column)]
    mates = mates[mates["category"].astype(str) != str(row["category"])]
    wanted = set(map(str, inputs))
    return {c: (f"moves with the label {label_column} (|Spearman| >= 0.8, measured) although its category is "
                f"{cat}, not {row['category']}: reported, not held out")
            for c, cat in zip(mates["column"], mates["category"]) if c in wanted}
