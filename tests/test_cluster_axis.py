"""The cluster-against-axis engine: the three rules it exists to enforce."""
import numpy as np
import pandas as pd
import pytest

from parascape import cluster_axis as CA


@pytest.fixture
def world():
    """600 genes, a moved cluster, an untouched cluster, and a cluster nobody measured."""
    rng = np.random.default_rng(0)
    genes = [f"G{i:04d}" for i in range(600)]
    expression = pd.Series(rng.lognormal(size=600), index=genes)
    axis = pd.Series(rng.normal(scale=1.0, size=600), index=genes)
    moved = genes[:40]
    axis.loc[moved] = axis.loc[moved] + 6.0          # a real, large effect
    markers = {
        "moved": moved,
        "quiet": genes[100:140],
        "barely_measured": [f"UNSEEN{i}" for i in range(38)] + genes[200:202],
    }
    return markers, axis, expression


def test_a_moved_cluster_is_found(world):
    markers, axis, expression = world
    res = CA.test_clusters_against_axis(markers, axis, expression, n_null=300, axis_name="drug")
    assert "moved" in res.reached
    assert res.table.loc["moved", "q"] < 0.05


def test_an_untouched_cluster_is_not_reported_as_reached(world):
    markers, axis, expression = world
    res = CA.test_clusters_against_axis(markers, axis, expression, n_null=300)
    assert "quiet" not in res.reached
    assert bool(res.table.loc["quiet", "testable"]) is True


def test_a_cluster_the_axis_did_not_measure_is_unknown_not_zero(world):
    """The rule that matters most: an empty cell must not mean two different things."""
    markers, axis, expression = world
    res = CA.test_clusters_against_axis(markers, axis, expression, n_null=200)
    row = res.table.loc["barely_measured"]
    assert bool(row["testable"]) is False
    assert bool(row["reaches"]) is False
    assert np.isnan(row["observed"]) and np.isnan(row["q"])
    assert "barely_measured" in res.unknown
    assert "markers are in the axis" in row["reason"]


def test_every_tested_cluster_carries_a_detection_limit(world):
    markers, axis, expression = world
    res = CA.test_clusters_against_axis(markers, axis, expression, n_null=300)
    tested = res.table[res.table["testable"].to_numpy(dtype=bool)]
    assert len(tested) >= 2
    assert (tested["detection_limit"] > 0).all()
    assert tested["over_detection_limit"].notna().all()


def test_the_null_is_matched_on_size_and_expression_not_on_category_size():
    """A small set of highly expressed genes must not win by being small and bright."""
    rng = np.random.default_rng(1)
    genes = [f"G{i:04d}" for i in range(500)]
    expression = pd.Series(np.linspace(0.01, 100.0, 500), index=genes)
    # the axis is a pure function of expression: brighter genes move more, for everyone
    axis = pd.Series(np.log1p(expression.to_numpy()), index=genes)
    bright = genes[-20:]                       # small AND the brightest
    res = CA.test_clusters_against_axis({"bright": bright}, axis, expression, n_null=400)
    # matched on expression, the brightest genes are compared with other bright genes, so this
    # is NOT an effect. An unmatched null would call it one.
    assert "bright" not in res.reached


def test_a_cluster_with_no_shared_genes_raises_rather_than_scoring_nothing():
    axis = pd.Series([1.0, 2.0], index=["A", "B"])
    expression = pd.Series([1.0, 2.0], index=["X", "Y"])
    with pytest.raises(ValueError, match="share no genes"):
        CA.test_clusters_against_axis({"c": ["A"]}, axis, expression)


def test_direction_can_be_asked_for_instead_of_size(world):
    markers, axis, expression = world
    signed = CA.test_clusters_against_axis(markers, axis, expression, n_null=300, absolute=False)
    assert signed.table.loc["moved", "observed"] > signed.table.loc["moved", "null_mean"]


def test_the_summary_names_all_three_outcomes(world):
    markers, axis, expression = world
    res = CA.test_clusters_against_axis(markers, axis, expression, n_null=200, axis_name="BKI-1708")
    text = res.summary()
    assert "BKI-1708" in text
    assert "not measured well enough" in text


def test_correction_is_applied_across_clusters(world):
    markers, axis, expression = world
    res = CA.test_clusters_against_axis(markers, axis, expression, n_null=300)
    tested = res.table[res.table["testable"].to_numpy(dtype=bool)]
    assert (tested["q"] >= tested["p"] - 1e-9).all(), "a corrected value is never below its raw one"


def test_an_unmatchable_set_says_why_rather_than_claiming_an_effect():
    """A set at the edge of the expression range cannot be matched, and must say so."""
    rng = np.random.default_rng(2)
    genes = [f"G{i:04d}" for i in range(500)]
    expression = pd.Series(np.linspace(0.01, 100.0, 500), index=genes)
    axis = pd.Series(np.log1p(expression.to_numpy()), index=genes)
    res = CA.test_clusters_against_axis({"bright": genes[-20:]}, axis, expression, n_null=300)
    row = res.table.loc["bright"]
    assert bool(row["testable"]) is True, "it WAS measured; the problem is the matching"
    assert bool(row["expression_matched"]) is False
    assert bool(row["reaches"]) is False
    assert "expression" in row["reason"]
    assert any("could not be centred" in n for n in res.notes)


def test_a_set_in_the_middle_of_the_range_is_matched_fine():
    rng = np.random.default_rng(3)
    genes = [f"G{i:04d}" for i in range(500)]
    expression = pd.Series(np.linspace(0.01, 100.0, 500), index=genes)
    axis = pd.Series(rng.normal(size=500), index=genes)
    res = CA.test_clusters_against_axis({"middle": genes[240:280]}, axis, expression, n_null=300)
    assert bool(res.table.loc["middle", "expression_matched"]) is True
