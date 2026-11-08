#!/usr/bin/env python3
"""Measured comprehensiveness of every registered dataset (instruction 65.04): genes it covers, from the shipped tables.

Einar's selection rule prefers the more comprehensive of comparable datasets. The review queue of 2026-10-07 held
`measured_comprehensiveness = unknown` for all 162 sources. Here it is measured: for each dataset, the genes of its
organism table with a value in any of its registered columns, and that share of the table. A dataset whose columns
are not in the shipped tables (source-only, protein-keyed host tables, unavailable) stays "unknown" with the reason.
Writes a new queue version beside the old one; nothing else changes.

    python scripts/build_dataset_comprehensiveness.py
"""
from __future__ import annotations

import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
QUEUE = os.path.join(ROOT, "results", "dataset_selection_2026_10_07", "source_review_queue.csv")
OUT = os.path.join(ROOT, "results", "dataset_selection_2026_11_08")


def main() -> int:
    """Measure coverage per dataset and write the updated queue."""
    from parascape import datasets as D
    tables = {"Tg": pd.read_parquet(os.path.join(ROOT, "parascape/data/nodes.parquet")),
              "Pf": pd.read_parquet(os.path.join(ROOT, "parascape/data/pf_nodes.parquet"))}
    by_key = {d.key: d for d in D.REGISTRY}
    q = pd.read_csv(QUEUE)
    genes, share, note = [], [], []
    for row in q.itertuples(index=False):
        d = by_key.get(row.dataset_id)
        table = tables.get(row.organism)
        if d is None or table is None:
            genes.append(None); share.append(None)
            note.append("unknown: not a gene table of a shipped organism" if d is not None else "unknown: not in the registry")
            continue
        cols = [c for c in d.columns if c in table.columns]
        if not cols:
            genes.append(None); share.append(None); note.append("unknown: no registered column in the shipped table")
            continue
        n = int(table[cols].notna().any(axis=1).sum())
        genes.append(n); share.append(round(n / len(table), 4)); note.append(f"measured: {len(cols)} column(s)")
    q["measured_genes"] = genes
    q["measured_share"] = share
    q["measured_comprehensiveness"] = note
    os.makedirs(OUT, exist_ok=True)
    q.to_csv(os.path.join(OUT, "source_review_queue.csv"), index=False)
    done = q["measured_share"].notna()
    print(f"{int(done.sum())} of {len(q)} datasets measured; median share {q.loc[done, 'measured_share'].median():.3f}")
    print(q.loc[done].groupby("organism")["measured_share"].describe().round(3).to_string())
    print(q.loc[~done, "measured_comprehensiveness"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
