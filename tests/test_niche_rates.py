"""Two rates on every location call (Einar 2026-11-26): the stated rate and the rate on other genes."""
import numpy as np
import pandas as pd

from parascape import niche_rates as NR


def _nodes(n=200):
    rng = np.random.default_rng(0)
    niches = np.array(["cytosol", "apicoplast", "dense granules", "unassigned"])
    return pd.DataFrame({"gene_id": [f"g{i}" for i in range(n)], "lopit_map": niches[rng.integers(0, 4, n)]})


def test_a_location_answer_gets_both_rates_and_the_other_gene_rate_is_measured():
    nodes = _nodes()
    calls = pd.DataFrame({"gene_id": nodes.gene_id, "predicted": "cytosol", "cluster_precision": 0.95})
    out, table = NR.add_two_rates(calls, nodes, "predicted", "cluster_precision")
    assert {"rate on markers", "rate on other genes", "hyperLOPIT classifier says"} <= set(out.columns)
    checkable = nodes.lopit_map != "unassigned"
    expected = float((nodes.lopit_map[checkable] == "cytosol").mean())
    assert abs(out["rate on other genes"].iloc[0] - expected) < 1e-9
    assert out["rate on markers"].iloc[0] == 0.95
    assert table["checkable"].sum() == int(checkable.sum())


def test_harness_wording_names_the_class():
    nodes = _nodes()
    calls = pd.DataFrame({"gene_id": nodes.gene_id, "answer": "lopit_map = apicoplast", "probability": 0.6})
    out, table = NR.add_two_rates(calls, nodes, "answer", "probability")
    assert table is not None and "rate on other genes" in out


def test_an_answer_that_is_not_about_niches_is_left_alone():
    nodes = _nodes()
    calls = pd.DataFrame({"gene_id": nodes.gene_id, "predicted": "EC 2", "cluster_precision": 0.9})
    out, table = NR.add_two_rates(calls, nodes, "predicted", "cluster_precision")
    assert table is None and list(out.columns) == list(calls.columns)


def test_an_answer_in_coarse_compartments_is_checked_against_the_coarsened_classifier():
    nodes = _nodes()
    coarse = {"cytosol": "cytosol", "apicoplast": "apicoplast", "dense granules": "apical_secretory",
              "unassigned": "unassigned"}
    nodes["compartment_best"] = nodes.lopit_map.map(coarse)
    calls = pd.DataFrame({"gene_id": nodes.gene_id, "predicted": "apical_secretory", "cluster_precision": 0.8})
    out, table = NR.add_two_rates(calls, nodes, "predicted", "cluster_precision")
    assert table is not None
    checkable = nodes.lopit_map != "unassigned"
    assert abs(out["rate on other genes"].iloc[0]
               - float((nodes.lopit_map[checkable] == "dense granules").mean())) < 1e-9


def test_a_gene_hyperlopit_never_measured_is_not_counted_as_a_disagreement():
    nodes = _nodes()
    nodes.loc[:99, "lopit_map"] = np.nan
    calls = pd.DataFrame({"gene_id": nodes.gene_id[:100], "predicted": "cytosol", "cluster_precision": 0.95})
    out, table = NR.add_two_rates(calls, nodes, "predicted", "cluster_precision")
    assert table["checkable"].sum() == 0
    assert out["rate on other genes"].isna().all()


def test_the_report_says_when_nothing_can_be_checked():
    nodes = _nodes()
    nodes["lopit_map"] = np.where(np.arange(len(nodes)) < 150, np.nan, nodes.lopit_map)
    calls = pd.DataFrame({"gene_id": nodes.gene_id[:100], "predicted": "cytosol", "cluster_precision": 0.95})
    _out, table = NR.add_two_rates(calls, nodes, "predicted", "cluster_precision")
    assert "upper bound" in NR.report(table)
