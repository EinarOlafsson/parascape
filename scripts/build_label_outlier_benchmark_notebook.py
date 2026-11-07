"""Corrupted-label detection, a SYNTHETIC task (instruction 64.12, outlier ranking), as an executed notebook.

A share of one label's values is deliberately swapped to another class in a copy of the table; label-outlier
ranking runs on the corrupted column, and the corrupted genes are looked for at the top. It measures whether the
method can find planted mistakes -- not whether any real annotation is wrong, which nothing here adjudicates.

    python scripts/build_label_outlier_benchmark_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"label_outlier_benchmark_{DATE}")
NB = os.path.join(REPO, "notebooks", f"label_outlier_benchmark_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Can outlier ranking find planted label mistakes?", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**A synthetic task, stated first.** 5 % of the genes with a value in each label have that value swapped "
          "to another class of the same label (seed 0), in a copy of the table. Label-outlier ranking runs on "
          "the corrupted label; the planted genes are counted among its top calls at depths fixed here (50, "
          "100, 200). This says whether the method finds mistakes we made; it says nothing about which real "
          "annotations are wrong (instruction 64.12).",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, numpy as np, pandas as pd",
            "from parascape import strategies as S",
            "base = S.Context.shipped('Tg')",
            "LABELS = ['compartment', 'dtm_class', 'cellcycle_phase']",
            "rows = []",
            "for label in LABELS:",
            "    nodes = base.nodes.copy()",
            "    values = nodes[label]",
            "    known = values.dropna().index",
            "    rng = np.random.default_rng(0)",
            "    planted = rng.choice(known, size=max(1, int(0.05 * len(known))), replace=False)",
            "    classes = values.dropna().unique()",
            "    for i in planted:",
            "        nodes.at[i, label] = rng.choice([c for c in classes if c != values.at[i]])",
            "    ctx = S.Context(nodes, organism='Tg')",
            "    r = S.get('label_outliers').run(ctx, target=label, top=200)",
            "    table = next(t for t in r.tables.values() if 'gene_id' in getattr(t, 'columns', []))",
            "    ranked = table.gene_id.astype(str).tolist()",
            "    planted_ids = set(nodes.loc[planted, 'gene_id'].astype(str))",
            "    hits = np.isin(ranked, list(planted_ids))",
            "    prevalence = len(planted_ids) / len(known)",
            "    row = dict(label=label, labelled=len(known), planted=len(planted_ids), prevalence=round(prevalence, 3))",
            "    for k in (50, 100, 200):",
            "        row[f'precision_at_{k}'] = round(float(hits[:k].mean()), 3) if len(hits) >= k else np.nan",
            "        row[f'recall_at_{k}'] = round(float(hits[:k].sum() / len(planted_ids)), 3)",
            "    rows.append(row)",
            "table = pd.DataFrame(rows)",
            "table.to_csv(os.path.join(OUT, 'planted_recovery.tsv'), sep='\\t', index=False)",
            "table")
    nb.write(NB)


if __name__ == "__main__":
    main()
