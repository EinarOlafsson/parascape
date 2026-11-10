"""Leakage attacks, one per space (instruction 53, WP7).

The category rule and the measured-redundancy rule are written once in the harness and are meant to hold in
ANY space. That is a claim, and a claim about leakage is worth only the attacks it survives. Each test below
tries to get something past the closure in a way that has either happened here or is one rename away from
happening, and asserts it is caught.

The attacks are deliberately the cheap ones, because those are what actually leak: a column nobody
categorised, a coordinate standing in for gene neighbourhood, and a value copied along an ortholog link and
then used to predict the thing it was copied from.
"""
from __future__ import annotations

import glob
import os

import pandas as pd
import pytest

from parascape import information_categories as IC

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPACES = os.path.join(ROOT, "results", "host_gene_spaces_2026_11_09")


def _tables():
    """(code, frame) for every further space built in this checkout."""
    for path in sorted(glob.glob(os.path.join(SPACES, "*_gene_nodes.parquet"))):
        yield os.path.basename(path)[:2].capitalize(), pd.read_parquet(path)


def test_every_further_space_is_fully_categorised():
    """The first attack is doing nothing: an uncategorised column is INVISIBLE to the closure.

    The closure excludes evidence by category. A column with no category is not excluded by anything, so a
    space nobody annotated has no leakage protection at all -- which is exactly what the five new spaces
    looked like on 2026-11-09 before this rule was written.
    """
    tables = dict(_tables())
    assert tables, "no further space is built in this checkout"
    for code, frame in tables.items():
        missing = [c for c in frame.columns if c != "gene_id" and not IC.categories_of(c, code)]
        assert not missing, f"{code}: {missing} carry no information category, so nothing can exclude them"


def test_coordinates_and_annotation_are_never_a_feature():
    """Gene order is not biology: a model reading a start coordinate is reading the neighbourhood."""
    for code, frame in _tables():
        for column in ("seqid", "start", "end", "strand", "biotype", "gene_name"):
            if column not in frame.columns:
                continue
            assert IC.categories_of(column, code) == ("identity",), f"{code}.{column}"
    # `identity` is excluded from every question, which is the property that makes the above matter.
    assert "identity" in IC.NEVER_A_FEATURE


def test_a_value_transferred_along_an_ortholog_cannot_predict_what_it_came_from():
    """The attack this codebase is one rename away from.

    The mouse space carries locations copied from human orthologs, named `*_from_human_ortholog` and graded
    inferred. If such a column were categorised as anything but `place`, a question about location could use
    it -- and it IS that location, copied. The link itself is `conservation`, the same kind of evidence as an
    orthogroup, so a question about conservation excludes it too.
    """
    assert IC.categories_of("location_from_human_ortholog", "Mm") == ("place",)
    assert IC.categories_of("hpa_location_from_human_ortholog", "Mm") == ("place",)
    assert IC.categories_of("human_ortholog", "Rn") == ("conservation",)
    assert IC.categories_of("human_ortholog_high_confidence", "Rn") == ("conservation",)
    excluded = IC.excluded_for("place", ["location_from_human_ortholog", "bgee_top_score"], "Mm")
    assert "location_from_human_ortholog" in excluded
    assert "bgee_top_score" not in excluded, "expression is not location, and must stay usable"
    assert "human_ortholog" in IC.excluded_for("conservation", ["human_ortholog", "bgee_top_score"], "Rn")


def test_expression_in_a_vector_or_parasite_space_is_transcript_evidence():
    """A new space's own columns must land in the vocabulary, not beside it."""
    assert IC.categories_of("expr_midgut_36h", "Ag") == ("transcript",)
    assert IC.categories_of("expr_pv_liver", "Pv") == ("transcript",)
    for column in ("expr_midgut_36h", "expr_pv_liver"):
        assert column in IC.excluded_for("transcript", [column, "human_ortholog"], "Ag")


def test_an_unknown_category_name_raises_rather_than_excluding_nothing():
    """Silently excluding nothing is the failure this vocabulary replaces."""
    with pytest.raises(Exception):
        IC.excluded_for("not_a_category", ["bgee_top_score"], "Rn")
