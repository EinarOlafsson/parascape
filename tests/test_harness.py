"""The harness enforces what three hand-written scoring runs forgot.

Each test here corresponds to one of the mechanisms instruction 73 asks for, and they are written to
fail if the mechanism is removed rather than to exercise the code path: the stage guard must *raise*
when a stage reaches the frozen partition, the holdout budget must *refuse* a second look by a
changed strategy, a question below its own label threshold must be unable to fill, and the decision
to fill must read the frozen-test bound and not the cross-validated one -- which was precisely the
optimism in the three published scoreboards.
"""
import numpy as np
import pandas as pd
import pytest

from parascape import harness as H
from parascape import splits as S
from parascape.question_registry import Question


def _question(**kwargs):
    """A boolean synthetic question whose own label family is held out, as the registry insists."""
    fields = dict(id="synthetic", text="Is the hidden variable positive?", answer_type="boolean",
                  question_class="B", label_source="synthetic", label_family="truth",
                  evidence_forbidden=("truth",), min_labels_for_frozen_test=200)
    fields.update(kwargs)
    return Question(**fields)


def _data(n=800, noise=0.15, seed=0, informative=True):
    """Synthetic genes in orthogroups of two, with one informative column and two of pure noise."""
    rng = np.random.default_rng(seed)
    ids = [f"g{i:04d}" for i in range(n)]
    groups = pd.Series([f"og{i // 2:04d}" for i in range(n)], index=ids)
    hidden = rng.normal(size=n)
    labels = pd.Series(hidden > 0.0, index=ids)
    column = hidden + rng.normal(scale=noise, size=n) if informative else rng.normal(size=n)
    features = pd.DataFrame({"evidence": column, "noise_a": rng.normal(size=n),
                             "noise_b": rng.normal(size=n)}, index=ids)
    families = pd.Series({"evidence": "assay", "noise_a": "junk", "noise_b": "junk"})
    return labels, features, families, groups


def _strategies():
    return [H.Strategy.of("every family", "hist_gbt", max_iter=60),
            H.Strategy.of("assay only", "hist_gbt", families=("assay",), max_iter=60)]


def _run(tmp_path, question=None, data=None, **kwargs):
    labels, features, families, groups = data or _data()
    return H.run_question(question or _question(), labels=labels, features=features,
                          families=families, groups=groups,
                          strategies=kwargs.pop("strategies", _strategies()),
                          log_path=tmp_path / "holdout.jsonl", **kwargs)


# --------------------------------------------------------------------------- the bound itself
def test_the_clopper_pearson_lower_bound_is_the_exact_one_sided_bound():
    assert H.clopper_pearson_lower(10, 10) == pytest.approx(0.1 ** (1 / 10))
    assert H.clopper_pearson_lower(60, 60) == pytest.approx(0.1 ** (1 / 60))
    assert H.clopper_pearson_lower(0, 20) == 0.0
    assert np.isnan(H.clopper_pearson_lower(0, 0))                 # no calls is not "wrong"
    assert H.clopper_pearson_lower(19, 20) < 19 / 20               # a bound, not the point estimate
    bounds = [H.clopper_pearson_lower(n, n) for n in (10, 20, 40, 80)]
    assert bounds == sorted(bounds)                                # more calls, tighter bound
    assert H.clopper_pearson_lower(21, 21) < 0.90 <= H.clopper_pearson_lower(22, 22)


# --------------------------------------------------------------------------- the stage guard
def test_the_guard_raises_when_a_stage_touches_the_test_partition():
    labels, _features, _families, groups = _data(n=200)
    split = S.make_split(list(labels.index), list(groups), organism="Tg",
                         benchmark_id="guard_probe", seed=0)
    guard = H.StageGuard(split)
    test = list(split.entities("test"))
    train = list(split.entities("train"))
    for stage in ("model_fit", "feature_selection", "setting_selection", "refit", "calibration"):
        with pytest.raises(ValueError, match="Forbidden"):
            guard.guard(stage, test)
    with pytest.raises(ValueError, match="Forbidden"):
        guard.guard("calibration", train)          # calibration may see only the calibration partition
    with pytest.raises(ValueError, match="Forbidden"):
        guard.guard("model_fit", list(split.entities("tune")))
    with pytest.raises(ValueError, match="Unknown fit stage"):
        guard.guard("have_a_look", train)
    assert guard.calls == ()                        # a refused call is not a recorded one
    guard.guard("refit", train + list(split.entities("tune")))
    assert guard.stages == frozenset({"refit"}) and guard.calls[0][0] == "refit"


def test_every_stage_the_harness_runs_is_guarded_and_recorded(tmp_path):
    result = _run(tmp_path)
    stages = {stage for stage, _n in result.guard_calls}
    assert set(H.STAGES_REQUIRED) <= stages
    assert len(result.guard_calls) > len(H.STAGES_REQUIRED)        # once per fold, not once per run
    assert not any(flag.startswith("guard_calls_missing") for flag in result.flags)


def test_a_result_with_no_guard_calls_may_not_fill_however_good_its_bound():
    unaudited = H.HarnessResult(question_id="x", answer_type="boolean", question_class="B",
                                n_labels=800, seed=0, split_identity="0" * 64, bar=0.90,
                                strategy="s", strategy_hash="h", threshold=0.5,
                                cross_validated={"precision_lower90": 0.99},
                                frozen_test={"precision_lower90": 0.99})
    assert not H.may_fill(unaudited)
    audited = H.HarnessResult(**{**unaudited.__dict__,
                                 "guard_calls": tuple((stage, 1) for stage in H.STAGES_REQUIRED)})
    assert H.may_fill(audited)


# --------------------------------------------------------------------------- the holdout budget
def test_the_holdout_log_refuses_a_second_use_by_a_changed_strategy(tmp_path):
    data = _data()
    first = _run(tmp_path, data=data)
    assert first.frozen_test is not None and "frozen_test_refused" not in first.flags
    assert len(H.holdout_log(tmp_path / "holdout.jsonl")) == 1

    changed = _run(tmp_path, data=data,
                   strategies=[H.Strategy.of("every family", "hist_gbt", max_iter=61)])
    assert changed.frozen_test is None
    assert "frozen_test_refused" in changed.flags and "needs_fresh_partition" in changed.flags
    assert "refused" in changed.refusal and not H.may_fill(changed)
    assert len(H.holdout_log(tmp_path / "holdout.jsonl")) == 1     # a refused look is not logged

    repeat = _run(tmp_path, data=data)
    assert repeat.frozen_test is not None and "frozen_test_repeat" in repeat.flags
    assert repeat.strategy_hash == first.strategy_hash
    assert repeat.frozen_test["precision_lower90"] == pytest.approx(
        first.frozen_test["precision_lower90"])
    assert len(H.holdout_log(tmp_path / "holdout.jsonl")) == 2

    fresh = _run(tmp_path, data=data, seed=7,
                 strategies=[H.Strategy.of("every family", "hist_gbt", max_iter=61)])
    assert fresh.frozen_test is not None and "frozen_test_refused" not in fresh.flags


def test_changing_the_evidence_also_spends_the_partition(tmp_path):
    """A narrower forbidden list is a changed strategy even when the model is untouched.

    Found by the harness's own permutation importance on the first real run: widening
    `evidence_forbidden` changes what the model sees while leaving the model, its settings and its
    version identical, so the strategy hash alone would have passed it as a repeat.
    """
    labels, features, families, groups = _data()
    one = _strategies()[:1]
    first = _run(tmp_path, data=(labels, features, families, groups), strategies=one)
    assert first.frozen_test is not None
    narrower = H.run_question(_question(evidence_forbidden=("truth", "junk")), labels=labels,
                              features=features[["evidence"]], families=families, groups=groups,
                              strategies=one, log_path=tmp_path / "holdout.jsonl")
    assert narrower.strategy_hash == first.strategy_hash      # the model is untouched
    assert "frozen_test_refused" in narrower.flags and "changed evidence" in narrower.refusal
    assert narrower.evidence_hash != first.evidence_hash


def test_the_budget_is_keyed_by_question_and_a_hash_of_the_strategy_definition():
    one = H.Strategy.of("a", "hist_gbt", max_iter=60)
    same = H.Strategy.of("a different name", "hist_gbt", max_iter=60)
    other = H.Strategy.of("a", "hist_gbt", max_iter=60, learning_rate=0.1)
    version = H.Strategy.of("a", "hist_gbt", version="2", max_iter=60)
    families = H.Strategy.of("a", "hist_gbt", families=("assay",), max_iter=60)
    assert one.identity == same.identity                 # renaming is not a changed strategy
    assert len({one.identity, other.identity, version.identity, families.identity}) == 4
    with pytest.raises(ValueError, match="Unknown model"):
        H.Strategy.of("a", "telepathy")


# --------------------------------------------------------------------------- the label threshold
def test_a_question_below_the_label_threshold_is_flagged_and_cannot_fill(tmp_path):
    data = _data(n=318)
    result = _run(tmp_path, question=_question(min_labels_for_frozen_test=400), data=data)
    assert result.n_labels == 318
    assert "no_frozen_test" in result.flags and result.frozen_test is None
    assert "rank genes" in result.refusal and "put nothing into the map" in result.refusal
    assert not H.may_fill(result)
    # The cross-validated estimate is still produced -- ranking is what these questions are for.
    assert np.isfinite(result.cross_validated["roc_auc"])
    assert not H.holdout_log(tmp_path / "holdout.jsonl").shape[0]   # the partition was never touched


# --------------------------------------------------------------------------- which number decides
def test_the_frozen_test_bound_not_the_cross_validated_one_decides_filling():
    common = dict(question_id="x", answer_type="boolean", question_class="B", n_labels=800, seed=0,
                  split_identity="0" * 64, bar=0.90, strategy="s", strategy_hash="h", threshold=0.9,
                  guard_calls=tuple((stage, 1) for stage in H.STAGES_REQUIRED))
    optimistic = H.HarnessResult(cross_validated={"precision_lower90": 0.95},
                                 frozen_test={"precision_lower90": 0.88}, **common)
    assert not H.may_fill(optimistic)
    pessimistic = H.HarnessResult(cross_validated={"precision_lower90": 0.70},
                                  frozen_test={"precision_lower90": 0.93}, **common)
    assert H.may_fill(pessimistic)
    assert not H.may_fill(H.HarnessResult(cross_validated={"precision_lower90": 0.99},
                                          frozen_test=None, flags=("no_frozen_test",), **common))


def test_both_numbers_are_returned_side_by_side(tmp_path):
    result = _run(tmp_path)
    summary = result.summary()
    assert len(summary) == 2
    assert list(summary["partition"]) == ["cross_validated (train+tune, nested)",
                                          "frozen_test (touched once)"]
    for column in ("precision", "precision_lower90", "coverage_at_bar", "calibration_error",
                   "roc_auc", "reaches_bar"):
        assert column in summary.columns
    assert np.isfinite(result.calibration["calibration_error"])
    assert not result.baselines.empty and not result.importance.empty
    assert set(result.importance["family"]) == {"assay", "junk"}
    assert result.as_dict()["may_fill"] == result.may_fill


# --------------------------------------------------------------------------- leaks
def test_a_deliberately_leaky_feature_is_refused_rather_than_quietly_dropped(tmp_path):
    labels, features, families, groups = _data()
    features["the_label_itself"] = labels.astype(float)
    families["the_label_itself"] = "truth"              # the family the labels came from
    with pytest.raises(ValueError, match="forbidden evidence families"):
        H.run_question(_question(), labels=labels, features=features, families=families,
                       groups=groups, strategies=_strategies(),
                       log_path=tmp_path / "holdout.jsonl")
    filtered = H.allowed_features(_question(), features, families)
    assert "the_label_itself" not in filtered.columns
    with pytest.raises(ValueError, match="forbidden families"):
        H.run_question(_question(), labels=labels, features=filtered, families=families,
                       groups=groups, strategies=[H.Strategy.of("leaky", "hist_gbt",
                                                                families=("truth",))],
                       log_path=tmp_path / "holdout.jsonl")


def test_nested_cross_validation_never_fits_outside_train_and_tune(tmp_path):
    result = _run(tmp_path)
    labels, _f, _fam, groups = _data()
    split = S.make_split(list(labels.index), list(groups), organism="Tg",
                         benchmark_id="synthetic", seed=0)
    assert result.split_identity == split.identity
    allowed = {"setting_selection": {"train", "tune"}, "refit": {"train", "tune"},
               "calibration": {"calibration"}}
    for stage, _n in result.guard_calls:
        assert S.FIT_ACCESS[stage] == allowed[stage]
    # The frozen partition's size is what it should be and nothing fitted on it.
    assert len(split.entities("test")) > 0
    assert result.frozen_test["n"] == len(split.entities("test"))


# --------------------------------------------------------------------------- the known answer
def test_synthetic_data_whose_answer_is_known_gives_the_expected_bound(tmp_path):
    """One informative column with a 0.10 noise level: a near-perfect question must clear the bar."""
    result = _run(tmp_path, data=_data(n=1600, noise=0.10))
    assert result.cross_validated["roc_auc"] > 0.97
    assert result.frozen_test["precision_lower90"] >= 0.90
    assert result.frozen_test["precision_lower90"] <= result.frozen_test["precision"]
    assert H.may_fill(result)
    assert result.importance.iloc[0]["family"] == "assay"
    best = result.baselines.query("partition == 'frozen_test'")["precision_lower90"].max()
    assert result.frozen_test["precision_lower90"] >= best - 0.05
    assert result.baselines.query("baseline == 'ortholog copy (abstains)'")["coverage"].max() == 0.0


def test_the_cross_validated_estimate_is_the_optimistic_one(tmp_path):
    """The reason the frozen partition exists, shown on data whose answer is known.

    At a 0.15 noise level the nested cross-validated bound clears 0.90 and the frozen-test bound
    does not: the folds that produced the first number also chose the strategy and the threshold.
    This is the same relationship the three published scoreboards stand in to a frozen test, which
    is why they are relabelled as cross-validated estimates rather than results.
    """
    result = _run(tmp_path, data=_data(n=800, noise=0.15))
    assert result.cross_validated["precision_lower90"] > result.frozen_test["precision_lower90"]
    assert result.cross_validated["reaches_bar"] and not result.frozen_test["reaches_bar"]
    assert not H.may_fill(result)


def test_a_question_with_no_signal_cannot_fill(tmp_path):
    """Pure noise: ranking is chance, no threshold clears the bar, and nothing may be filled."""
    result = _run(tmp_path, data=_data(informative=False, seed=3))
    assert 0.35 < result.cross_validated["roc_auc"] < 0.65
    assert result.frozen_test["precision_lower90"] < 0.90
    assert not H.may_fill(result)


def test_the_frozen_test_is_reported_within_each_stratum(tmp_path):
    """Instruction 64.16: the same threshold and calibration, read again inside each stratum."""
    labels, features, families, groups = _data()
    strata = pd.Series(["well-studied" if i % 2 else "understudied" for i in range(len(labels))],
                       index=labels.index)
    result = _run(tmp_path, data=(labels, features, families, groups), strata=strata)
    by = result.frozen_test["by_stratum"]
    assert set(by) == {"understudied", "well-studied"}
    assert sum(v["n"] for v in by.values()) == result.frozen_test["n"]


def test_per_gene_rows_reproduce_every_partition_number(tmp_path):
    """Instruction 64.10: the harness keeps one row per scored gene, and its summaries derive from them."""
    from parascape import label_metrics as LM
    r = _run(tmp_path)
    assert set(r.rows["partition"]) == {"cross_validated", "frozen_test"}
    for partition, block in (("cross_validated", r.cross_validated), ("frozen_test", r.frozen_test)):
        rows = r.rows[r.rows["partition"] == partition]
        assert len(rows) == block["n"]
        assert int(rows["called"].sum()) == block["calls"]
        assert int((rows["called"] & rows["correct"]).sum()) == block["correct"]
        assert LM.overall(rows)["accuracy_all"] == pytest.approx(block["accuracy"])
