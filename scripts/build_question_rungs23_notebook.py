"""Rungs 2-3 of the strategy ladder on the curated-label questions (instruction 73, phase 4), as an executed notebook.

The admitted catalogue questions of the curated-label categories (place, function, phenotype, modification,
interaction) were re-scored on a fresh partition (seed 2) with three strategies the inner folds choose between:
boosted trees (rung 1), one model per evidence family combined by a meta-learner (rung 2), and network propagation
over the evidence graph (rung 3). Both comparisons are reported: on the SAME folds (the inner-fold score of each
strategy, the fair one) and the frozen test of the chosen strategy against the rung-1 frozen test on seed 0 (a
different partition, so differences there are partly the split).

    python scripts/score_question_catalogue.py --models hist_gbt,family_stacking,graph_propagation --seed 2 \
        --categories place,function,phenotype,modification,interaction --out results/question_scores_rungs23_2026_11_08
    python scripts/build_question_rungs23_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"question_rungs23_{DATE}")
NB = os.path.join(REPO, "notebooks", f"question_rungs23_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    nb = ExecutedNotebook("Strategy ladder rungs 2-3 on curated-label questions", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**What this tests (instruction 73, phase 4).** Whether a model per evidence family combined (rung 2) or "
          "network propagation (rung 3) answers the curated-label questions better than boosted trees on all the "
          "evidence at once (rung 1), and whether more cells become fillable at the 90 % bar.",
          "**Fair comparison:** the inner-fold score of every strategy on the same folds (average precision, the "
          "harness selection score). The frozen test is touched once, for the strategy the inner folds chose.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import glob, json, os, numpy as np, pandas as pd",
            "recs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(REPO, 'results/question_scores_rungs23_2026_11_08/*.json')))]",
            "R = pd.DataFrame([r for r in recs if 'error' not in r])",
            "errors = [r for r in recs if 'error' in r]",
            "rung1 = pd.read_csv(os.path.join(REPO, 'results/question_scoreboard_2026_11_08/scores.tsv'), sep='\\t').set_index('question')",
            "print(len(R), 'questions scored;', len(errors), 'errors')",
            "R.groupby('category').size().to_frame('questions')")
    nb.md("## 1 · Same folds: which strategy wins the inner folds")
    nb.code("rows = []",
            "for r in R.itertuples(index=False):",
            "    folds = pd.DataFrame(r.inner_scores).T",
            "    rows.append(dict(question=r.question, category=r.category, **folds.mean().round(4).to_dict()))",
            "inner = pd.DataFrame(rows)",
            "strats = [c for c in inner.columns if c not in ('question', 'category')]",
            "inner['best'] = inner[strats].fillna(-1).idxmax(axis=1).where(inner[strats].notna().any(axis=1), '')",
            "inner['gain_over_rung1'] = (inner[strats].max(axis=1) - inner['boosted trees']).round(4)",
            "inner.to_csv(os.path.join(OUT, 'inner_scores.tsv'), sep='\\t', index=False)",
            "summary = inner.groupby('category')[strats].median().round(3)",
            "summary['questions'] = inner.groupby('category').size()",
            "summary['rung 2 or 3 best'] = inner[inner.best != 'boosted trees'].groupby('category').size().reindex(summary.index).fillna(0).astype(int)",
            "summary")
    nb.md("## 2 · Frozen test of the chosen strategy, and cells fillable at 90 %",
          "Rung 1 numbers are from seed 0, the ladder from seed 2: different partitions.")
    nb.code("R['rung1_frozen_auc'] = R.question.map(rung1.frozen_auc)",
            "R['rung1_may_fill'] = R.question.map(rung1.may_fill)",
            "R['chosen_mostly'] = R.chosen.map(lambda c: pd.Series(c).mode().iloc[0] if len(c) else '')",
            "cols = ['question', 'category', 'chosen_mostly', 'frozen_auc', 'rung1_frozen_auc', 'precision_lower90', 'may_fill', 'rung1_may_fill']",
            "F = R[cols].sort_values('frozen_auc', ascending=False)",
            "F.to_csv(os.path.join(OUT, 'frozen.tsv'), sep='\\t', index=False)",
            "out = dict(questions=len(R), errors=len(errors),",
            "           median_frozen_auc=round(float(R.frozen_auc.median()), 3),",
            "           median_rung1_frozen_auc=round(float(R.rung1_frozen_auc.median()), 3),",
            "           may_fill=int(R.may_fill.sum()), rung1_may_fill=int(R.rung1_may_fill.fillna(False).sum()),",
            "           inner_best=inner.best.value_counts().to_dict())",
            "json.dump(out, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "print(out)",
            "F.head(20).round(3)")
    nb.md("## 3 · Reading",
          "A strategy that wins the inner folds on a category is the one to use there; a category where boosted "
          "trees still win shows that splitting evidence by family or walking the evidence graph adds nothing on "
          "these labels. Cells become fillable only where the frozen lower bound reaches 0.90.")
    nb.write(NB)


if __name__ == "__main__":
    main()
