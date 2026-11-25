"""The literature findings (free Affinage run, 2026-11-24): what ships and how it is used."""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parascape import datasets, literature_findings as L  # noqa: E402


def test_the_shipped_tables_load_and_match_the_registry():
    t = L.labels(log=lambda *_: None)
    assert len(t) > 500 and t.index.is_unique
    assert tuple(t.columns) == L.COLUMNS == datasets.get("literature_affinage").columns
    assert t.index.str.startswith("TGME49_").all()


def test_labels_rest_only_on_established_findings_whose_protein_the_abstract_names():
    f = L.findings()
    used = f[(f.status == "established") & f.named_in_abstract]
    t = L.labels(log=lambda *_: None)
    with_place = t.index[t.lit_compartment != ""]
    placed = set(used[(used.compartment != "") & (used.method_tier <= 2)].gene_id)
    assert set(with_place) <= placed
    host = set(t.index[t.lit_host_facing])
    assert host <= set(used[used.kind == "host_facing"].gene_id)


def test_every_finding_cites_a_paper_and_partner_ids_are_genes():
    f = L.findings()
    assert f.pmid.str.fullmatch(r"\d+").all()
    assert f.source_group.eq(L.SOURCE_GROUP).all()
    p = L.pairs()
    assert (p.a < p.b).all() and p.a.str.startswith("TGME49_").all() and p.b.str.startswith("TGME49_").all()


def test_resolve_maps_old_ids_and_keeps_one_row_per_gene():
    t = L.labels(log=lambda *_: None, resolve=lambda g: "TGME49_000000" if g.endswith("0") else g)
    assert t.index.is_unique


def test_a_missing_table_leaves_the_columns_out(tmp_path):
    assert L.labels(log=lambda *_: None, root=str(tmp_path)).empty
    assert L.findings(root=str(tmp_path)).empty


def test_each_label_carries_both_axes_of_confidence():
    """Method tier (how strong the experiment) and papers (how much evidence), per label (Affinage)."""
    t = pd.read_parquet(os.path.join(L.HERE, "literature_gene_labels.parquet"))
    placed = t[t["lit_compartment"] != ""]
    assert placed["lit_compartment_tier"].isin([1, 2]).all(), "a place rests on imaging or stronger"
    assert (placed["lit_compartment_papers"] >= 1).all()
    hf = t[t["lit_host_facing"]]
    assert hf["lit_host_facing_tier"].notna().all() and (hf["lit_host_facing_papers"] >= 1).all()
