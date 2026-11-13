"""The host PROTEIN spaces (instruction 53, WP9; Einar, 2026-11-13): UniProt-keyed, opened through lenses.

The tables come from the data directory (the downloaded release when PARASCAPE_CACHE is set); the layouts are
read from there or, until the next release ships them, from the source tree. Every test that needs a file
skips when it is absent, so the suite says what it could not check rather than failing for a missing download.
"""
from __future__ import annotations

import hashlib
import os
import shutil

import numpy as np
import pandas as pd
import pytest

from parascape import organisms as O
from parascape.spaces import host_proteins as HP

ROWS = {O.HUMAN_PROTEINS: 21_020, O.MOUSE_PROTEINS: 16_473}


def _need(code):
    if not os.path.exists(HP.nodes_path(code)):
        pytest.skip(f"{HP.get(code).nodes} is not in this data directory")


@pytest.mark.parametrize("code", list(HP.SPACES))
def test_each_protein_table_is_a_valid_uniprot_space(code):
    _need(code)
    nodes = HP.read_nodes(code)
    assert len(nodes) == ROWS[code]
    assert nodes.gene_id.is_unique and nodes.gene_id.map(HP.get(code).matches).all()
    # Nothing but the identifier column is renamed: every measured column survives.
    raw = pd.read_parquet(HP.nodes_path(code))
    assert list(nodes.columns) == ["gene_id"] + list(raw.columns[1:])
    for column in raw.columns[1:]:
        assert raw[column].notna().sum() == nodes[column].notna().sum(), column


def test_an_ensembl_gene_id_is_not_a_protein_space_id():
    """The point of two kinds of host space: a gene id must not pass as a protein id, nor the reverse."""
    from parascape.spaces.human import HUMAN_SPACE
    assert not HP.HUMAN_PROTEIN_SPACE.matches("ENSG00000075624")
    assert not HUMAN_SPACE.matches("P60709")
    assert HP.HUMAN_PROTEIN_SPACE.matches("P60709") and HP.HUMAN_PROTEIN_SPACE.matches("A0A024R1R8")


@pytest.mark.parametrize("code,lens", [(c, l) for c, ls in HP.LENSES.items() for l in ls])
def test_every_built_lens_loads_beside_the_table_it_was_built_from(code, lens):
    _need(code)
    nodes, xyz, meta = HP.load(code, lens)
    assert xyz.shape == (len(nodes), 3) and np.isfinite(xyz).all()
    assert meta["space"] == code[:2] and meta["lens"] == lens
    assert meta["genes"] == len(nodes)
    assert meta["node_file_sha256"] == hashlib.sha256(open(HP.nodes_path(code), "rb").read()).hexdigest()


def test_a_layout_built_from_another_table_is_refused(tmp_path, monkeypatch):
    code = O.HUMAN_PROTEINS
    _need(code)
    changed = pd.read_parquet(HP.nodes_path(code))
    changed.loc[0, "fibroblast_tpm"] = 12345.0
    path = tmp_path / HP.get(code).nodes
    changed.to_parquet(path)
    layout = HP.layout_path(code, "tissue")
    shutil.copy(layout, tmp_path / os.path.basename(layout))
    monkeypatch.setattr(HP, "nodes_path", lambda c: str(path))
    with pytest.raises(ValueError, match="different"):
        HP.load(code, "tissue")


def test_the_refused_mouse_lens_stays_refused_with_its_reason():
    assert (O.MOUSE_PROTEINS, "parasite_facing") in HP.REFUSED
    assert "parasite_facing" not in HP.LENSES[O.MOUSE_PROTEINS]
    assert "0.9 per cent" in HP.REFUSED[(O.MOUSE_PROTEINS, "parasite_facing")]
    with pytest.raises(ValueError, match="refused"):
        HP.load(O.MOUSE_PROTEINS, "parasite_facing")


def test_protein_spaces_are_declared_beside_the_registry_not_in_it():
    """Registering them would put a UniProt pattern into `detect`, which cannot tell human from mouse."""
    assert not set(HP.SPACES) & set(O.SPACES)
    assert set(O.codes(available=True)) <= set(O.SPACES)
    assert O.detect(["P60709", "Q9Y6K9", "A0A024R1R8"]) is None
    assert O.by_species("Homo sapiens host proteins") is HP.HUMAN_PROTEIN_SPACE
    assert O.by_species("Mus musculus host proteins") is HP.MOUSE_PROTEIN_SPACE


def test_available_lists_a_space_only_with_its_table_and_a_layout(monkeypatch, tmp_path):
    monkeypatch.setattr(HP, "nodes_path", lambda c: str(tmp_path / "absent.parquet"))
    assert HP.available() == [] and HP.lenses(O.HUMAN_PROTEINS) == []


@pytest.fixture(scope="module")
def qapp():
    pytest.importorskip("PySide6")
    from PySide6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def test_a_host_protein_space_opens_in_the_window_through_each_lens(qapp):
    from parascape.app import Window, HOST_PROTEIN_SPECIES
    name = HP.HUMAN_PROTEIN_SPACE.species
    if O.HUMAN_PROTEINS not in HP.available():
        pytest.skip("the human protein table or its layouts are absent")
    assert name in HOST_PROTEIN_SPECIES
    w = Window(species=name, lens="tissue")
    try:
        assert w.n == ROWS[O.HUMAN_PROTEINS] and w.lens == "tissue"
        assert w.edges == {} and "by tissue" in w.windowTitle()
        # Picking a protein shows the protein panel: its name, its UniProt link, and "unknown" for a
        # missing value -- never the parasite panel calling it uncharacterised.
        row = int(np.flatnonzero(w.nodes.fibroblast_tpm.isna().to_numpy())[0])
        w.show_detail(row)
        html = w.detail.toHtml()
        assert "uniprot.org" in html and "unknown" in w.detail.toPlainText()
        assert "uncharacterised" not in w.detail.toPlainText()
        # The menu offers each lens; the one open is checked; the refused mouse lens is disabled.
        actions = w.host_protein_actions
        assert actions[(O.HUMAN_PROTEINS, "tissue")].isChecked()
        assert not actions[(O.MOUSE_PROTEINS, "parasite_facing")].isEnabled()
        other = w.open_species(name, lens="parasite_facing")
        try:
            assert other is not w and other.lens == "parasite_facing"
            assert np.array_equal(other.nodes.gene_id, w.nodes.gene_id)
            assert not np.allclose(other.xyz, w.xyz)
            assert other.open_species(name, lens="parasite_facing") is other
        finally:
            other.close()
    finally:
        w.close()


def test_a_host_protein_window_is_never_reopened_at_start_up(qapp):
    """Choosing one must not move the next start-up (or later tests) onto a protein table."""
    from parascape.app import available_species
    assert HP.HUMAN_PROTEIN_SPACE.species not in available_species()
    assert HP.MOUSE_PROTEIN_SPACE.species not in available_species()


def test_an_installed_host_gene_pack_opens_as_its_own_map(qapp, tmp_path, monkeypatch):
    """Route 1 as far as the data on this machine goes: the Ensembl-keyed gene space opens from its pack."""
    from parascape import space_catalogue as SC
    from parascape.app import Window, HOST_GENE_SPECIES, species_code
    monkeypatch.setenv("PARASCAPE_STATE", str(tmp_path))
    if not SC.local_archives("Mm"):
        pytest.skip("no mouse gene pack in this checkout")
    SC.install("Mm")
    name = "Mus musculus genes"
    assert name in HOST_GENE_SPECIES and species_code(name) == O.MOUSE
    w = Window(species=name)
    try:
        assert w.n == 18_428 and w.nodes.gene_id.str.fullmatch(r"ENSMUSG\d{11}").all()
        assert w.edges == {}
        w.show_detail(0)
        assert "ensembl.org" in w.detail.toHtml() and "host gene" in w.detail.toPlainText()
    finally:
        w.close()
