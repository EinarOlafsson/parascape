"""Per-class numbers from per-gene rows (instruction 64.10)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from parascape import label_metrics as LM


def _rows():
    return pd.DataFrame({"truth": ["a", "a", "a", "b", "b", "c"],
                         "prediction": ["a", "a", None, "a", "b", None]})


def test_recall_counts_abstentions_and_precision_only_named_calls():
    per = LM.by_class(_rows()).set_index("label")
    assert per.loc["a", "recall"] == 2 / 3 and per.loc["a", "precision"] == 2 / 3
    assert per.loc["b", "recall"] == 0.5 and per.loc["b", "precision"] == 1.0
    assert per.loc["c", "recall"] == 0 and np.isnan(per.loc["c", "precision"])   # never predicted: unknown, not 0


def test_overall_coverage_accuracy_and_chance():
    o = LM.overall(_rows())
    assert o["coverage"] == 4 / 6 and o["accuracy_answered"] == 3 / 4 and o["accuracy_all"] == 3 / 6
    assert o["chance"] == 0.5 and o["classes"] == 3


def test_confusion_keeps_abstentions():
    c = LM.confusion(_rows())
    assert c.loc["a", "abstained"] == 1 and c.loc["b", "a"] == 1 and c.loc["c", "abstained"] == 1
