"""Missingness and dominance warnings (instruction 60, change 7)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from parascape import diagnostics as D


def _nodes():
    rng = np.random.default_rng(0)
    ids = [f"g{i}" for i in range(100)]
    a = rng.normal(size=100)
    b = rng.normal(size=100)
    a[:10] = np.nan                 # the first ten genes are mostly empty
    b[:10] = np.nan
    b[50:60] = np.nan
    return pd.DataFrame({"gene_id": ids, "a": a, "b": b})


def test_named_genes_that_are_mostly_gaps_are_flagged():
    out = D.missingness(_nodes(), ["a", "b"], [f"g{i}" for i in range(10)])
    assert out["named_missing"] == 1.0 and out["warn"]
    assert "filled-in values" in out["sentence"]


def test_named_genes_as_measured_as_the_rest_are_not_flagged():
    out = D.missingness(_nodes(), ["a", "b"], [f"g{i}" for i in range(20, 40)])
    assert not out["warn"] and "sentence" not in out


def test_missingness_with_nothing_to_measure_says_nothing():
    assert D.missingness(_nodes(), ["no_such"], ["g1"]) == {}
    assert D.missingness(_nodes(), ["a"], ["not_a_gene"]) == {}


def test_one_family_carrying_most_of_the_importance_is_flagged():
    imp = pd.DataFrame({"family": ["x", "y", "z"], "importance_mean": [0.3, 0.05, -0.01]})
    out = D.dominance(imp)
    assert out["family"] == "x" and out["share"] > 0.8 and out["warn"]


def test_spread_importance_is_not_flagged():
    imp = pd.DataFrame({"family": ["x", "y", "z"], "importance_mean": [0.1, 0.09, 0.08]})
    assert not D.dominance(imp)["warn"]
    assert D.dominance(pd.DataFrame()) == {}


def test_a_harness_result_reports_its_dominance():
    from parascape.harness import HarnessResult
    r = HarnessResult(question_id="q", answer_type="boolean", question_class="B", n_labels=10, seed=0,
                      split_identity="s", bar=0.9, strategy="x", strategy_hash="h", threshold=0.5,
                      cross_validated={},
                      importance=pd.DataFrame({"family": ["a", "b"], "importance_mean": [0.2, 0.01]}))
    assert r.dominance["family"] == "a" and r.dominance["warn"]
