"""THE REACH MODULE: a single-cell cluster block met with a bulk treatment, in a real window.

Einar's words, which are the acceptance criteria: *"we should try to build a mechanism in
parascape that would allow the user to, in parascape, do the analysis we did, using the scRNAseq
data together with the bulk drug RNA-seq."*

So the tests here are about the mechanism being USABLE by somebody who did not write the script:
the bradyzoite subtype block is on offer without being typed in, the stored contrasts are on
offer, and the three refusals the analysis depends on survive the trip through the interface.
A unit test of the engine is in `test_cluster_axis.py`; this file is about the module.
"""
import numpy as np
import pandas as pd
import pytest
from PySide6 import QtWidgets

from parascape import cluster_axis as CA
from parascape import modules as MOD
from parascape import reach_module as RM


@pytest.fixture(scope="module")
def qapp():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def window(qapp):
    from parascape import app as A
    try:
        win = A.Window(species="Toxoplasma gondii")
    except Exception as error:                      # pragma: no cover - needs a GL context
        pytest.skip(f"cannot create the window offscreen: {error}")
    yield win
    from parascape.map_module import wait_for_projections
    wait_for_projections()
    win.close()
    qapp.processEvents()


@pytest.fixture(scope="module")
def module(window):
    shell = getattr(window, "reach_module", None)
    if shell is None or getattr(shell, "module", None) is None:
        pytest.skip("the Reach module did not build in this window")
    return shell.module


# ----------------------------------------------------------------- it is there, and reachable

def test_the_module_is_in_the_catalogue_and_built_by_its_own_shell():
    assert MOD.get("reach").band == "Ask"
    assert MOD.shell_module("reach") == "reach_module"


def test_typing_what_it_is_called_finds_it():
    kind, found, view = MOD.search_all("Reach")[0]
    assert (kind, found.key, view) == ("module", "reach", None)


def test_the_bradyzoite_subtypes_are_on_offer_without_anybody_typing_a_column(module):
    """The analysis this exists to repeat starts from a published block of five columns."""
    assert "brady_subtypes" in module.blocks
    found = module.blocks["brady_subtypes"]
    assert sorted(found["columns"]) == ["A", "B", "C", "D", "E"]
    offered = [module.block.itemData(i) for i in range(module.block.count())]
    assert "brady_subtypes" in offered
    assert "" in offered, "a table of your own is always the last offer"


def test_a_block_of_two_quantities_is_not_offered_as_a_block_of_clusters(module):
    """`stage_enriched_derived` beside `stage_margin_derived` is which stage and how far ahead,
    not two stages, and offering it would invite an answer about nothing."""
    assert "stage_enriched" not in module.blocks


def test_the_stored_contrasts_are_on_offer(module):
    if module.catalogue.empty:
        pytest.skip("no expression traits release in this checkout")
    assert module.contrast.count() > 0
    column = module.contrast.currentData()
    assert str(column) in set(module.catalogue.index)


# ------------------------------------------------------------------------ it runs end to end

def test_it_tests_every_cluster_and_names_what_it_could_not_see(module):
    if module.catalogue.empty:
        pytest.skip("no expression traits release in this checkout")
    module.block.setCurrentIndex([module.block.itemData(i) for i in
                                  range(module.block.count())].index("brady_subtypes"))
    module.n_null.setValue(200)
    result = module.run()
    assert result is not None, "the five subtypes against a stored contrast must produce a result"
    assert len(result.table) >= 2
    shown = RM.as_shown(result)
    assert list(shown.columns)[0] == "cluster"
    assert "verdict" in shown.columns
    assert set(shown["verdict"]) <= {RM.REACHES, RM.NO_EFFECT, RM.UNKNOWN, RM.UNMATCHED}
    assert module.table.rowCount() == len(shown)
    assert module.export_button.isEnabled()


def test_the_prose_says_what_the_baseline_was_and_never_the_size_of_the_cluster(module):
    if module.catalogue.empty:
        pytest.skip("no expression traits release in this checkout")
    result = module.run()
    text = RM.reading(result)
    assert "matched" in text and "expression" in text
    assert "floor" not in text.lower(), "Einar's standing instruction on the word"


def test_the_markers_of_the_chosen_cluster_are_what_the_next_steps_act_on(module):
    if module.catalogue.empty:
        pytest.skip("no expression traits release in this checkout")
    module.run()
    genes = module.genes()
    assert genes, "a tested cluster carries its markers on"
    assert set(genes) <= set(module.nodes.gene_id.astype(str))


# ------------------------------------------------- the three refusals survive the interface

def test_a_cluster_the_contrast_barely_measured_is_unknown_and_never_no_effect():
    """The refusal that matters most, read through the words the table prints."""
    genes = [f"TGME49_{i:06d}" for i in range(400)]
    rng = np.random.default_rng(0)
    markers = {"measured": genes[:60], "barely": [f"ABSENT_{i}" for i in range(50)]}
    axis = pd.Series(rng.normal(size=400), index=genes)
    expression = pd.Series(rng.lognormal(size=400), index=genes)
    result = CA.test_clusters_against_axis(markers, axis, expression, n_null=200)
    shown = RM.as_shown(result).set_index("cluster")
    assert shown.loc["barely", "verdict"] == RM.UNKNOWN
    assert pd.isna(shown.loc["barely", "observed"]), "no number is better than a wrong zero"
    assert RM.UNKNOWN == "n.m."


def test_a_cluster_that_is_measured_and_unmoved_carries_what_could_have_been_seen():
    genes = [f"TGME49_{i:06d}" for i in range(400)]
    rng = np.random.default_rng(1)
    result = CA.test_clusters_against_axis(
        {"quiet": genes[:60], "other": genes[100:160]},
        pd.Series(rng.normal(size=400), index=genes),
        pd.Series(rng.lognormal(size=400), index=genes), n_null=200)
    shown = RM.as_shown(result).set_index("cluster")
    assert set(shown["verdict"]) == {RM.NO_EFFECT}
    assert (shown["detection_limit"] > 0).all()
    assert "could have seen" in RM.reading(result)


def test_a_marker_rule_that_marks_nothing_says_so_rather_than_returning_an_answer(module, tmp_path):
    """At a margin no gene clears, the module refuses instead of testing empty sets.

    Written first against the bradyzoite block at a margin of ten, which still leaves 72
    markers in subtype A -- the block is more specific than the test assumed. A table whose
    clusters differ by a fifth is the case the refusal is actually for.
    """
    if module.catalogue.empty:
        pytest.skip("no expression traits release in this checkout")
    genes = [f"TGME49_{i:06d}" for i in range(80)]
    flat = tmp_path / "flat.csv"
    pd.DataFrame({"gene_id": genes,
                  "cluster_1": np.linspace(1.0, 2.0, 80),
                  "cluster_2": np.linspace(1.2, 2.4, 80)}).to_csv(flat, index=False)
    module.load(str(flat))
    module.margin.setValue(10.0)
    assert module.run() is None, "nothing marks either cluster, so there is nothing to test"
    module.margin.setValue(CA.MIN_MARGIN)
    module.block.setCurrentIndex([module.block.itemData(i) for i in
                                  range(module.block.count())].index("brady_subtypes"))


# ----------------------------------------------------------------------- a table of your own

def test_a_published_block_of_your_own_is_read_and_refused_honestly(module, tmp_path):
    good = tmp_path / "clusters.csv"
    pd.DataFrame({"gene_id": ["TGME49_000001", "TGME49_000002"],
                  "cluster_1": [10.0, 1.0], "cluster_2": [1.0, 10.0]}).to_csv(good, index=False)
    frame = module.load(str(good))
    assert list(frame.columns)[0] == "gene_id"
    assert module.source["clusters"] == ("cluster_1", "cluster_2")

    one = tmp_path / "one.csv"
    pd.DataFrame({"gene_id": ["TGME49_000001"], "cluster_1": [1.0]}).to_csv(one, index=False)
    with pytest.raises(ValueError, match="runner-up"):
        module.load(str(one))

    unnamed = tmp_path / "unnamed.csv"
    pd.DataFrame({"x": ["a"], "c1": [1.0], "c2": [2.0]}).to_csv(unnamed, index=False)
    with pytest.raises(ValueError, match="gene_id"):
        module.load(str(unnamed))
