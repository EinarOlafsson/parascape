"""The label catalogue: every label is described, and a name says its evidence when it was copied."""
import importlib.util
import os

import numpy as np
import pandas as pd
import pytest

from parascape import label_catalogue as LC

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(rel):
    spec = importlib.util.spec_from_file_location(os.path.basename(rel)[:-3], os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_every_entry_is_complete_and_uses_a_known_evidence_type():
    for name, e in LC.ENTRIES.items():
        assert e["evidence"] in LC.EVIDENCE, name
        assert e["sources"] and e["made"] and e["column"] and e["left_out"], name


def test_copied_labels_say_so_in_their_name_and_only_they_do():
    for name, e in LC.ENTRIES.items():
        assert ("copied by orthology" in name) == (e["evidence"] == LC.COPIED), name


def test_the_interaction_build_names_have_entries():
    B = _load("scripts/labels/build_interaction_labels.py")
    for name in B.LABELS.values():
        assert name in LC.ENTRIES


def test_catalogue_counts_classes_and_refuses_an_undescribed_label():
    y = np.array(["a", "a", "b", None], dtype=object)
    L = {"hyperLOPIT niche": (y, "categorical", ["__place__"])}
    t = LC.catalogue(L, {"hyperLOPIT niche": ["lopit_x"]})
    row = t.iloc[0]
    assert row["genes_labelled"] == 3 and row["classes"] == "a (2); b (1)"
    assert row["left_out_columns"] == 1 and row["evidence"] == LC.OBSERVED
    with pytest.raises(KeyError):
        LC.catalogue({"made-up label": (y, "categorical", [])})


def test_sources_resolve_to_a_name():
    for e in LC.ENTRIES.values():
        for k in e["sources"]:
            assert LC.source_record(k)["name"] != "UNKNOWN SOURCE", k


def test_modules_drop_small_communities():
    B = _load("scripts/labels/build_interaction_labels.py")
    big = [(f"g{i}", f"g{j}") for i in range(12) for j in range(i + 1, 12)]
    small = [("s1", "s2"), ("s2", "s3")]
    pairs = pd.DataFrame(big + small, columns=["gene_a", "gene_b"])
    lab = B.modules(pairs, "module", {})
    assert set(lab.index) == {f"g{i}" for i in range(12)}
    assert lab.nunique() == 1
