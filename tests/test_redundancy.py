"""Measured redundancy is reported, never used to exclude (instruction 73, Einar 2026-11-07: report only)."""
import pandas as pd

from parascape import redundancy as R


def test_inputs_moving_with_the_label_across_categories_are_named():
    table = pd.DataFrame({"column": ["label", "rpf_a", "rna_b", "other"], "group": [1, 1, 1, 2],
                          "category": ["transcript", "translation", "transcript", "translation"]})
    out = R.overlap("label", ["rpf_a", "rna_b", "other"], table=table)
    assert set(out) == {"rpf_a"}                     # same category is the category rule's job; group 2 unrelated
    assert "reported, not held out" in out["rpf_a"]
    assert R.overlap("unmeasured", ["rpf_a"], table=table) == {}


def test_the_shipped_report_finds_the_translation_transcript_crossing():
    g = R.groups("Tg")
    assert len(g) > 300
    crossing = g.groupby("group")["category"].nunique()
    assert (crossing > 1).any()
