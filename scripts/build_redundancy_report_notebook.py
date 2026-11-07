"""Which columns carry the same information (instruction 73, phase 3: the redundancy report), as a notebook.

Measured, not declared: for every numeric column of the Toxoplasma node table measured on at least 1,000 genes,
the rank correlation with every other column on the genes both measure (pairwise, at least 300 shared genes).
Columns are grouped when |Spearman| >= 0.8 (complete linkage, so every pair inside a group clears the bar), and
the number of independent axes is estimated from the eigenvalues of the correlation matrix (participation ratio).
Groups are then compared with the declared information categories: a group that spans categories is a place where
the provenance-based categories miss shared information.

    python scripts/build_redundancy_report_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"redundancy_report_{DATE}")
NB = os.path.join(REPO, "notebooks", f"redundancy_report_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Which columns carry the same information", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**The question (Einar, 2026-11-07):** does the information map use feature correlation to decide which "
          "features hold the same information? Until now, no -- categories come from where a column was measured "
          "(provenance). This report measures it: rank correlation between every pair of numeric columns, groups "
          "of columns that move together (|Spearman| >= 0.8 for every pair in a group), the effective number of "
          "independent axes, and where the measured groups cross the declared categories.",
          "Nothing is changed by this report; it is the evidence for deciding whether to add a measured redundancy "
          "rule next to the category rule.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, numpy as np, pandas as pd",
            "from scipy.cluster.hierarchy import linkage, fcluster",
            "from scipy.spatial.distance import squareform",
            "from parascape import information_categories as IC",
            "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet'))",
            "num = nodes.select_dtypes('number')",
            "num = num.loc[:, (num.notna().sum() >= 1000) & (num.nunique() > 2)]",
            "print(num.shape[1], 'numeric columns measured on >= 1,000 genes (binary columns left out)')")
    nb.md("## 1 · Pairwise rank correlation", "On the genes both columns measure; pairs sharing fewer than 300 genes "
          "are left undecided (not counted as independent).")
    nb.code("ranks = num.rank()",
            "rho = ranks.corr(method='pearson', min_periods=300)   # Pearson on ranks = Spearman, pairwise complete",
            "undecided = int(rho.isna().values[np.triu_indices(len(rho), 1)].sum())",
            "pairs = rho.where(np.triu(np.ones(rho.shape, bool), 1)).stack()",
            "strong = pairs[pairs.abs() >= 0.8]",
            "print(len(pairs), 'pairs;', len(strong), 'with |Spearman| >= 0.8;', undecided, 'undecided (too few shared genes)')")
    nb.md("## 2 · Groups of columns that move together")
    nb.code("d = (1 - rho.abs()).fillna(1.0).clip(lower=0).to_numpy(copy=True)",
            "np.fill_diagonal(d, 0)",
            "d = (d + d.T) / 2",
            "Z = linkage(squareform(d, checks=False), method='complete')",
            "group = pd.Series(fcluster(Z, t=0.2, criterion='distance'), index=rho.index)",
            "cat = pd.Series({c: '+'.join(IC.categories_of(c)) or 'unannotated' for c in rho.index})",
            "sizes = group.value_counts()",
            "multi = sizes[sizes > 1].index",
            "groups = pd.DataFrame([dict(group=int(g), columns=int(sizes[g]), categories=', '.join(sorted(set(cat[group == g].astype(str)))),",
            "                            example=', '.join(list(group[group == g].index[:5]))) for g in multi]).sort_values('columns', ascending=False)",
            "groups['crosses_categories'] = groups.categories.str.contains(',')",
            "groups.to_csv(os.path.join(OUT, 'groups.tsv'), sep='\\t', index=False)",
            "pd.DataFrame({'column': group.index, 'group': group.values, 'category': cat.values}).to_csv(os.path.join(OUT, 'column_groups.tsv'), sep='\\t', index=False)",
            "print(len(rho), 'columns fall into', group.nunique(), 'groups;', len(multi), 'groups hold more than one column;',",
            "      int(groups.crosses_categories.sum()), 'of those cross a declared category')",
            "groups.head(25)")
    nb.md("## 3 · How many independent axes",
          "Participation ratio of the eigenvalues of the correlation matrix (missing correlations set to 0): "
          "(sum of eigenvalues)^2 / sum of squared eigenvalues. It equals the number of columns if all are "
          "independent and 1 if all carry one signal.")
    nb.code("ev = np.clip(np.linalg.eigvalsh(rho.fillna(0).values), 0, None)",
            "pr = float(ev.sum() ** 2 / (ev ** 2).sum())",
            "axes = dict(columns=len(rho), groups_at_0_8=int(group.nunique()), participation_ratio=round(pr, 1))",
            "pd.Series(axes).to_csv(os.path.join(OUT, 'axes.tsv'), sep='\\t', header=False)",
            "axes")
    nb.md("## 4 · Where measured redundancy crosses the declared categories",
          "Each such group is a place where an input from one category carries the label information of another, "
          "so the category rule alone would not catch it.")
    nb.code("cross = groups[groups.crosses_categories]",
            "cross.to_csv(os.path.join(OUT, 'crossing_groups.tsv'), sep='\\t', index=False)",
            "cross")
    nb.write(NB)


if __name__ == "__main__":
    main()
