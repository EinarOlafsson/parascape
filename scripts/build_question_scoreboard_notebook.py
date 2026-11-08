"""The question scoreboard (instruction 73, phase 4), as an executed notebook.

Reads the per-question harness results written by `scripts/score_question_catalogue.py` and keeps those whose label
column is still admitted in the catalogue (bookkeeping columns were excluded after scoring began). Per question:
frozen-test AUC and average precision against prevalence, the cross-validated AUC, whether it may fill (frozen-test
lower bound >= 0.90). Per category: how many questions, median frozen AUC, how many clearly above chance (frozen AUC
>= 0.70), how many may fill.

    python scripts/build_question_scoreboard_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"question_scoreboard_{DATE}")
NB = os.path.join(REPO, "notebooks", f"question_scoreboard_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    nb = ExecutedNotebook("The question scoreboard", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**What this reads (instruction 73, phase 4).** Every admitted catalogue question scored through the "
          "harness (Einar, 2026-11-08: keep all admitted): its own dataset family and information category held "
          "out, orthogroup-grouped nested cross-validation, the frozen test touched once. A question is one yes/no "
          "label: the extreme 20 % of a value column, or one class against the rest.",
          "**How to read a row.** Frozen AUC 0.5 is chance; average precision is compared with the prevalence of "
          "yes; *may fill* means the frozen-test lower bound on precision reached 0.90 at some call threshold.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import glob, json, os, numpy as np, pandas as pd",
            "cat = pd.read_csv(os.path.join(REPO, 'results/question_catalogue_2026_11_08/catalogue.tsv'), sep='\\t')",
            "admitted = set(cat[cat.admitted].column)",
            "rows = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(REPO, 'results/question_scores_2026_11_08/*.json')))]",
            "S = pd.DataFrame(rows)",
            "errors = S[S.get('error').notna()] if 'error' in S else S.iloc[:0]",
            "S = S[S.column.isin(admitted) & (S.get('error').isna() if 'error' in S else True)]",
            "print(len(rows), 'scored;', len(errors), 'errors;', len(S), 'on admitted columns')",
            "S['lift_ap'] = (S.frozen_ap / S.prevalence).round(2)",
            "S.to_csv(os.path.join(OUT, 'scores.tsv'), sep='\\t', index=False)",
            "errors.to_csv(os.path.join(OUT, 'errors.tsv'), sep='\\t', index=False)")
    nb.md("## 1 · Per category")
    nb.code("board = S.groupby('category').agg(questions=('question', 'size'), median_frozen_auc=('frozen_auc', 'median'),",
            "                                  above_070=('frozen_auc', lambda a: int((a >= 0.70).sum())),",
            "                                  near_chance=('frozen_auc', lambda a: int((a < 0.60).sum())),",
            "                                  may_fill=('may_fill', 'sum')).round(3).sort_values('median_frozen_auc', ascending=False)",
            "board.to_csv(os.path.join(OUT, 'by_category.tsv'), sep='\\t')",
            "summary = dict(questions=len(S), median_frozen_auc=round(float(S.frozen_auc.median()), 3),",
            "               above_070=int((S.frozen_auc >= 0.70).sum()), near_chance=int((S.frozen_auc < 0.60).sum()),",
            "               may_fill=int(S.may_fill.sum()), errors=len(errors))",
            "json.dump(summary, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "print(summary)",
            "board")
    nb.md("## 2 · The best-answered and the worst-answered questions")
    nb.code("cols = ['question', 'category', 'labelled', 'positives', 'frozen_auc', 'frozen_ap', 'prevalence', 'may_fill', 'top_families']",
            "S.sort_values('frozen_auc', ascending=False)[cols].head(20)")
    nb.code("S.sort_values('frozen_auc')[cols].head(15)")
    nb.md("## 3 · Questions that may fill",
          "Only these may write a cell into the map as a fact; every other question ranks genes and says how "
          "often such a ranking is right.")
    nb.code("S[S.may_fill][cols]")
    nb.write(NB)


if __name__ == "__main__":
    main()
