"""Per-class label records from per-gene rows (instruction 64.10), as an executed notebook.

For every method and label in the shipped track record (both organisms): coverage, accuracy on answered and on all
genes, chance, macro F1, and per class recall / precision / F1, with the confusion table -- all derived from the
per-gene held-out rows by `label_metrics`. Rare classes the record left out (fewer than the minimum class size)
are listed as evaluation gaps.

    python scripts/build_label_metrics_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"label_metrics_{DATE}")
NB = os.path.join(REPO, "notebooks", f"label_metrics_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Per-class label records", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**What this reads (instruction 64.10).** The shipped track record holds one row per held-out gene for "
          "each method and label (the gene and its orthogroup hidden, five folds). Every number below is computed "
          "from those rows by `parascape.label_metrics`: coverage (share of genes answered), accuracy on answered "
          "and on all genes, chance (the commonest class), macro F1, and per class recall, precision and F1. An "
          "abstention counts against recall; a class never predicted has no precision (unknown, not zero).",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, pandas as pd, numpy as np",
            "from parascape import label_metrics as LM, track_record as T, strategies as S",
            "ledger = T.shipped()",
            "together = ledger[ledger['mode'].astype(str) == 'together']",
            "pairs = together[['organism', 'strategy', 'target']].astype(str).drop_duplicates().sort_values(['organism', 'target', 'strategy'])",
            "print(len(ledger), 'rows;', len(pairs), 'method x label records')")
    nb.md("## 1 · Every method on every label")
    nb.code("overall, classes, confusions = [], [], {}",
            "for org, key, target in pairs.itertuples(index=False):",
            "    rows = LM.track_record_rows(together[together.organism.astype(str) == org], key, target)",
            "    o = LM.overall(rows); o.update(organism=org, strategy=key, target=target); overall.append(o)",
            "    per = LM.by_class(rows).assign(organism=org, strategy=key, target=target); classes.append(per)",
            "    confusions[(org, key, target)] = LM.confusion(rows)",
            "overall = pd.DataFrame(overall)[['organism', 'target', 'strategy', 'genes', 'coverage', 'accuracy_answered', 'accuracy_all', 'chance', 'macro_f1', 'classes']]",
            "classes = pd.concat(classes, ignore_index=True)",
            "overall.to_csv(os.path.join(OUT, 'overall.tsv'), sep='\\t', index=False)",
            "classes.to_csv(os.path.join(OUT, 'by_class.tsv'), sep='\\t', index=False)",
            "with open(os.path.join(OUT, 'confusion.tsv'), 'w') as fh:",
            "    for (org, key, target), table in confusions.items():",
            "        fh.write(f'# {org} {target} {key}\\n'); table.to_csv(fh, sep='\\t'); fh.write('\\n')",
            "overall.assign(beats_chance=overall.accuracy_all > overall.chance).groupby(['organism', 'target']).agg(",
            "    methods=('strategy', 'size'), best_macro_f1=('macro_f1', 'max'), beat_chance=('beats_chance', 'sum')).round(3)")
    nb.md("## 2 · The best method per label, and its weakest classes")
    nb.code("best = overall.loc[overall.groupby(['organism', 'target']).macro_f1.idxmax()]",
            "weak = []",
            "for r in best.itertuples(index=False):",
            "    c = classes[(classes.organism == r.organism) & (classes.strategy == r.strategy) & (classes.target == r.target) & (classes.genes > 0)]",
            "    w = c.nsmallest(2, 'recall')",
            "    weak.append(dict(organism=r.organism, target=r.target, method=r.strategy, macro_f1=round(r.macro_f1, 3),",
            "                     weakest=', '.join(f'{x.label} (recall {x.recall:.0%}, n={x.genes})' for x in w.itertuples(index=False))))",
            "pd.DataFrame(weak)")
    nb.md("## 3 · Classes the record could not evaluate",
          "Classes with fewer genes than the minimum class size (`strategies.MIN_CLASS`) are left out of every hold-out run; they "
          "are evaluation gaps, listed here per label.")
    nb.code("gaps = []",
            "for org in together.organism.astype(str).unique():",
            "    ctx = S.Context.shipped(org)",
            "    for target in together[together.organism.astype(str) == org].target.astype(str).unique():",
            "        counts = ctx.truth(target).value_counts()",
            "        small = counts[counts < S.MIN_CLASS]",
            "        gaps.append(dict(organism=org, target=target, min_class=S.MIN_CLASS, classes_left_out=len(small),",
            "                         genes_left_out=int(small.sum()), which=', '.join(f'{k} ({v})' for k, v in small.head(6).items())))",
            "gaps = pd.DataFrame(gaps); gaps.to_csv(os.path.join(OUT, 'evaluation_gaps.tsv'), sep='\\t', index=False)",
            "gaps")
    nb.write(NB)


if __name__ == "__main__":
    main()
