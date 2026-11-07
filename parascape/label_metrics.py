"""Per-class numbers derived from per-gene held-out rows (instruction 64.10).

Every summary here is computed from rows of (truth, prediction), one per held-out gene, where a missing
prediction is an abstention. The same function reads the 14-strategy track record (`track_record.shipped`) and
the harness rows (`harness.HarnessResult.rows`), so no number exists that the rows do not support.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _calls(rows: pd.DataFrame) -> tuple:
    truth = rows["truth"].astype("object").astype(str)
    prediction = rows["prediction"].astype("object")
    answered = prediction.notna() & (prediction.astype(str) != "None")
    return truth, prediction.where(answered).astype("object"), answered


def by_class(rows: pd.DataFrame) -> pd.DataFrame:
    """Per true class: genes, answered, right, recall, precision of calls naming it, F1.

    Recall counts an abstention as not recovered; precision counts only calls that name the class. Classes
    that were never predicted keep a missing precision rather than zero.
    """
    truth, prediction, answered = _calls(rows)
    out = []
    for label in sorted(set(truth) | set(prediction.dropna().astype(str))):
        is_true = truth == label
        named = answered & (prediction.astype(str) == label)
        right = int((named & is_true).sum())
        recall = right / int(is_true.sum()) if is_true.any() else float("nan")
        precision = right / int(named.sum()) if named.any() else float("nan")
        f1 = 2 * precision * recall / (precision + recall) \
            if np.isfinite(precision) and np.isfinite(recall) and precision + recall > 0 else float("nan")
        out.append(dict(label=label, genes=int(is_true.sum()), answered=int((is_true & answered).sum()),
                        named=int(named.sum()), right=right, recall=recall, precision=precision, f1=f1))
    return pd.DataFrame(out)


def overall(rows: pd.DataFrame) -> dict:
    """Coverage, accuracy on answered genes, accuracy on all, chance (commonest class) and macro F1."""
    truth, prediction, answered = _calls(rows)
    n = len(rows)
    right = int((answered & (prediction.astype(str) == truth)).sum())
    per = by_class(rows)
    chance = float(truth.value_counts(normalize=True).iloc[0]) if n else float("nan")
    return dict(genes=n, answered=int(answered.sum()), coverage=float(answered.mean()) if n else float("nan"),
                accuracy_answered=right / int(answered.sum()) if answered.any() else float("nan"),
                accuracy_all=right / n if n else float("nan"), chance=chance,
                macro_f1=float(per.loc[per["genes"] > 0, "f1"].fillna(0).mean()) if len(per) else float("nan"),
                classes=int((per["genes"] > 0).sum()) if len(per) else 0)


def confusion(rows: pd.DataFrame) -> pd.DataFrame:
    """True class (rows) by call (columns), with an ``abstained`` column."""
    truth, prediction, _answered = _calls(rows)
    call = prediction.fillna("abstained").astype(str)
    return pd.crosstab(truth.rename("truth"), call.rename("call"))


def track_record_rows(ledger: pd.DataFrame, strategy: str, target: str) -> pd.DataFrame:
    """The per-gene fold hold-out rows of one strategy and label, at its first recorded setting."""
    rows = ledger[(ledger["strategy"].astype(str) == strategy) & (ledger["target"].astype(str) == target)
                  & (ledger["mode"].astype(str) == "together")]
    if not len(rows):
        return rows
    first = rows["setting_key"].astype(str).iloc[0]
    return rows[rows["setting_key"].astype(str) == first]
