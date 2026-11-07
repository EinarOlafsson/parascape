"""Error dependence between all recorded methods, on matched held-out cohorts (instruction 64.32), as a notebook.

Uses the strictly nested outer-test calls of the 14 recorded methods (instruction 64.34): every method called the
same genes with those genes' answers hidden. For each pair of methods: how correlated their errors are, on the genes
both answered, with an orthogroup bootstrap interval; pairs answering fewer than 100 genes in common are marked
unresolved rather than estimated.

    python scripts/build_method_dependence_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"method_dependence_{DATE}")
NB = os.path.join(REPO, "notebooks", f"method_dependence_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Which methods make the same mistakes?", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**What is measured (instruction 64.32).** On three labels, all 14 recorded methods called the same "
          "genes with those genes' answers hidden (the strictly nested outer tests of 64.34). For every pair: the "
          "correlation of their errors on the genes both answered, with a 95 % interval from resampling whole "
          "orthogroups. Two methods with highly correlated errors are one piece of evidence, not two; a pair with "
          "under 100 genes in common is unresolved.", f"Built {_dt.date.today().isoformat()}.")
    nb.code("import itertools, os, numpy as np, pandas as pd",
            "FOLDERS = {'dtm_class': 'meta_inference_nested_2026_11_07', 'compartment': 'meta_inference_nested_2026_11_07_compartment',",
            "           'cellcycle_phase': 'meta_inference_nested_2026_11_07_cellcycle_phase'}",
            "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet')).set_index('gene_id')",
            "rng = np.random.default_rng(0)",
            "rows = []",
            "for label, folder in FOLDERS.items():",
            "    path = os.path.join(REPO, 'results', folder)",
            "    calls = pd.concat([pd.read_parquet(os.path.join(path, f)).assign(strategy=f[6:-8]) for f in os.listdir(path)",
            "                       if f.startswith('calls_') and f.endswith('.parquet')])",
            "    calls = calls[calls.role == 'outer_test']",
            "    wide = calls.pivot_table(index='gene_id', columns='strategy', values='call', aggfunc='first')",
            "    truth = calls.groupby('gene_id').truth.first().reindex(wide.index)",
            "    groups = nodes['orthogroup'].reindex(wide.index).fillna(pd.Series(wide.index, index=wide.index)).to_numpy()",
            "    for a, b in itertools.combinations(wide.columns, 2):",
            "        both = (wide[a] != 'ABSTAIN') & (wide[b] != 'ABSTAIN')",
            "        n = int(both.sum())",
            "        if n < 100:",
            "            rows.append(dict(label=label, a=a, b=b, genes=n, error_correlation=np.nan, low=np.nan, high=np.nan, status='unresolved'))",
            "            continue",
            "        ea = (wide.loc[both, a] != truth[both]).to_numpy(float); eb = (wide.loc[both, b] != truth[both]).to_numpy(float)",
            "        if ea.std() == 0 or eb.std() == 0:",
            "            rows.append(dict(label=label, a=a, b=b, genes=n, error_correlation=np.nan, low=np.nan, high=np.nan, status='no errors to compare'))",
            "            continue",
            "        phi = float(np.corrcoef(ea, eb)[0, 1])",
            "        g = groups[both.to_numpy()]; ug = np.unique(g); idx_of = {x: np.flatnonzero(g == x) for x in ug}",
            "        boots = []",
            "        for _ in range(200):",
            "            idx = np.concatenate([idx_of[x] for x in rng.choice(ug, size=len(ug))])",
            "            if ea[idx].std() and eb[idx].std():",
            "                boots.append(np.corrcoef(ea[idx], eb[idx])[0, 1])",
            "        rows.append(dict(label=label, a=a, b=b, genes=n, error_correlation=round(phi, 3),",
            "                         low=round(float(np.percentile(boots, 2.5)), 3), high=round(float(np.percentile(boots, 97.5)), 3), status='measured'))",
            "pairs = pd.DataFrame(rows)",
            "pairs.to_csv(os.path.join(OUT, 'pairs.tsv'), sep='\\t', index=False)",
            "summary = pairs.groupby('label').agg(pairs=('a', 'size'), unresolved=('status', lambda s: int((s != 'measured').sum())),",
            "    median_error_correlation=('error_correlation', 'median'), highest=('error_correlation', 'max')).round(3)",
            "summary.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\\t')",
            "summary")
    nb.code("pairs[pairs.status == 'measured'].sort_values('error_correlation', ascending=False).head(12)")
    nb.write(NB)


if __name__ == "__main__":
    main()
