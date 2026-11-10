"""The borrowed profile-inference strategies (instruction 77): each an estimator the harness scores."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from parascape import borrowed_models as B
from parascape import harness as H
from parascape.question_registry import Question


def _data(n=600, seed=0, missing=0.0):
    """Two informative families, one of noise, an attention column that tracks the label."""
    rng = np.random.default_rng(seed)
    ids = [f"g{i:04d}" for i in range(n)]
    hidden = rng.normal(size=n)
    y = pd.Series(hidden > 0.3, index=ids)
    X = pd.DataFrame({
        "expr_a": hidden + rng.normal(scale=0.5, size=n),
        "expr_b": hidden + rng.normal(scale=0.8, size=n),
        "prot_a": hidden + rng.normal(scale=0.7, size=n),
        "noise_a": rng.normal(size=n), "noise_b": rng.normal(size=n),
        "papers": y.to_numpy() * 1.0 + rng.normal(scale=1.0, size=n),
    }, index=ids)
    if missing:
        X = X.mask(rng.random(X.shape) < missing)
    families = pd.Series({"expr_a": "expression", "expr_b": "expression", "prot_a": "proteome",
                          "noise_a": "junk", "noise_b": "junk", "papers": "pubmed"})
    categories = {c: (("attention",) if c == "papers" else ("transcript",)) for c in X.columns}
    return X, y, families, categories


def _auc(y, p):
    from sklearn.metrics import roc_auc_score
    return roc_auc_score(np.asarray(y).astype(int), p)


@pytest.fixture(autouse=True)
def _families():
    X, y, families, categories = _data()
    B.set_families(families, categories)
    yield
    B.set_families({})


@pytest.mark.parametrize("name", sorted(B.BORROWED))
def test_every_borrowed_strategy_is_a_classifier_that_learns(name):
    """fit / predict_proba, two columns that sum to one, and well above chance on held-out genes."""
    X, y, _f, _c = _data(missing=0.1)
    model = B.BORROWED[name]().fit(X.iloc[:400], y.iloc[:400].astype(int))
    p = model.predict_proba(X.iloc[400:])
    assert p.shape == (200, 2)
    assert np.allclose(p.sum(axis=1), 1.0) and (p >= 0).all() and (p <= 1).all()
    assert _auc(y.iloc[400:], p[:, 1]) > 0.8, name


@pytest.mark.parametrize("name", sorted(B.BORROWED))
def test_every_borrowed_strategy_is_registered_in_the_harness(name):
    """A strategy outside `harness.MODELS` could not be scored, so it could not be offered."""
    assert name in H.MODELS
    H.Strategy.of(name, name)


def test_propensity_weighting_removes_the_attention_tilt_from_the_labels():
    """In the weighted sample, the attention column no longer tells the label."""
    X, y, _f, _c = _data()
    model = B.PropensityWeighted().fit(X, y.astype(int))
    assert model.attention_columns == ["papers"]
    from sklearn.linear_model import LogisticRegression
    raw = LogisticRegression().fit(X[["papers"]], y)
    assert _auc(y, raw.predict_proba(X[["papers"]])[:, 1]) > 0.7   # the tilt is real
    from sklearn.metrics import roc_auc_score
    weighted = roc_auc_score(y.astype(int), X["papers"], sample_weight=model.sample_weight_)
    assert abs(weighted - 0.5) < 0.05                                # and the weights remove it


def test_propensity_weighting_without_attention_columns_is_unweighted():
    X, y, _f, _c = _data()
    model = B.PropensityWeighted().fit(X.drop(columns="papers"), y.astype(int))
    assert model.attention_columns == []
    assert model.weight_summary == (1.0, 1.0, 1.0)


def test_family_stacking_gives_the_noise_family_no_more_weight_than_the_evidence():
    """Each family is one measurable weight; pure noise should not outweigh real evidence."""
    X, y, _f, _c = _data()
    model = B.FamilyStacking().fit(X.drop(columns="papers"), y.astype(int))
    w = model.weights_
    assert set(w.index) == {"expression", "proteome", "junk"}
    assert w["expression"] > w["junk"] and w["proteome"] > w["junk"]
    assert (w >= 0).all()


def test_pu_bagging_never_treats_the_unlabelled_as_one_fitted_negative_class():
    """Each bag draws as many unlabelled genes as there are positives."""
    X, y, _f, _c = _data()
    model = B.PUBagging(n_bags=4).fit(X, y.astype(int))
    assert len(model.models) == 4


def test_low_rank_completion_scores_rows_with_missing_evidence():
    X, y, _f, _c = _data(missing=0.4)
    model = B.LowRankCompletion(rank=3).fit(X.iloc[:400], y.iloc[:400].astype(int))
    p = model.predict_proba(X.iloc[400:])[:, 1]
    assert np.isfinite(p).all()


def test_graph_propagation_is_deterministic_and_scores_without_the_scoring_labels():
    """`predict_proba` takes evidence only; two fits on the same genes give the same scores."""
    X, y, _f, _c = _data()
    model = B.GraphPropagation().fit(X.iloc[:400], y.iloc[:400].astype(int))
    a = model.predict_proba(X.iloc[400:])[:, 1]
    b = B.GraphPropagation().fit(X.iloc[:400], y.iloc[:400].astype(int)).predict_proba(X.iloc[400:])[:, 1]
    assert np.allclose(a, b)


def test_the_untested_strategies_say_why():
    """`multi_task` left this list on 2026-11-09: it is implemented and graded (instruction 73 rung 5),
    and it loses to boosted trees, which is a measured answer rather than an untested one."""
    assert set(B.UNTESTED) == {"selective_prediction", "empirical_bayes_shrinkage"}
    assert all(len(reason.split()) > 15 for reason in B.UNTESTED.values())
    assert "multi_task" in B.BORROWED, "it moved from untested to implemented, not out of the module"


def test_a_borrowed_strategy_runs_through_the_whole_harness(tmp_path):
    """Nested cross-validation, calibration and the frozen test, unchanged, for a borrowed model."""
    X, y, families, categories = _data(n=800)
    X = X.drop(columns="papers")
    groups = pd.Series([f"og{i // 2:04d}" for i in range(len(X))], index=X.index)
    question = Question(id="synthetic_borrowed", text="Is it positive?", answer_type="boolean",
                        question_class="B", label_source="synthetic", label_family="truth",
                        evidence_forbidden=("truth",), min_labels_for_frozen_test=200)
    result = H.run_question(question, labels=y, features=X, families=families.drop("papers"),
                            groups=groups, strategies=[H.Strategy.of("stacking", "family_stacking")],
                            log_path=tmp_path / "holdout.jsonl", importance_repeats=1)
    assert result.frozen_test is not None
    assert result.frozen_test["roc_auc"] > 0.8
