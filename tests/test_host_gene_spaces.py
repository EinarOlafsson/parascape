"""The host GENE spaces: two lenses, the protein-gene link rule, and the three Map gaps.

Instruction 53, WP9 route 1 (Einar decided "both" on 2026-11-13). Everything that needs an installed
pack skips when there is none, so the suite says what it could not check rather than failing for a
missing download.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import zipfile

import numpy as np
import pandas as pd
import pytest

from parascape import host_links, packs
from parascape.spaces import host_genes as HG

GENES = {"Hs": 58_988, "Mm": 18_428}


@pytest.fixture(scope="module")
def installed(tmp_path_factory):
    """Install whichever host gene packs this checkout holds, into a state directory of their own."""
    from parascape import space_catalogue as SC
    state = tmp_path_factory.mktemp("packs")
    os.environ["PARASCAPE_STATE"] = str(state)
    have = []
    for code in ("Hs", "Mm"):
        if SC.local_archives(code):
            SC.install(code)
            have.append(code)
    yield have
    os.environ.pop("PARASCAPE_STATE", None)


def _need(code, installed):
    if code not in installed:
        pytest.skip(f"no {code} gene pack in this checkout")
    if not HG.lenses(code):
        pytest.skip(f"the installed {code} pack carries no lens layout")


# --------------------------------------------------------------------------- the licence class
def _space():
    from parascape.spaces.human import HUMAN_SPACE
    return HUMAN_SPACE


def _pack(tmp_path, licences, sources=None):
    nodes = tmp_path / "hs_gene_nodes.parquet"
    pd.DataFrame({"gene_id": [f"ENSG{i:011d}" for i in range(3)],
                  "tg_lung_log2fc": [1.0, 2.0, 3.0]}).to_parquet(nodes, index=False)
    graph = tmp_path / "hs_gene_graph_tissue.npz"
    np.savez(graph, xyz=np.zeros((3, 3), np.float32),
             gene_ids=np.array([f"ENSG{i:011d}" for i in range(3)], dtype=str))
    out = tmp_path / "pack.zip"
    return packs.build(_space(), "2026.11.13",
                       {"hs_gene_nodes.parquet": nodes, "hs_gene_graph_tissue.npz": graph},
                       licences, out, sources), out


CITED = {"gene_id": "O", "tg_lung_log2fc": "G"}
SOURCE = {"tg_lung_log2fc": {"accession": "GSE119855",
                             "citation": "Singhania A et al. Nature Communications 2019"}}


def test_a_cited_public_deposit_column_ships_when_its_accession_and_paper_are_recorded(tmp_path):
    """Einar, 2026-11-13: GEO columns may be redistributed, with citation. G is that licence class."""
    sha, out = _pack(tmp_path, CITED, SOURCE)
    manifest = packs.verify(out, sha)
    assert manifest["column_licenses"]["tg_lung_log2fc"] == "G"
    assert manifest["column_sources"]["tg_lung_log2fc"]["accession"] == "GSE119855"


@pytest.mark.parametrize("sources,missing", [
    (None, "no recorded source"),
    ({"tg_lung_log2fc": {"citation": "Singhania A et al. Nature Communications 2019"}},
     "needs an accession"),
    ({"tg_lung_log2fc": {"accession": "GSE119855", "citation": "x"}}, "needs a citation"),
    ({"tg_lung_log2fc": {"accession": "not an accession", "citation": "a long enough citation"}},
     "needs an accession"),
])
def test_a_cited_column_without_both_an_accession_and_a_citation_is_refused(tmp_path, sources,
                                                                           missing):
    """The permission is to redistribute a deposit WITH its attribution, so G without it is refused."""
    with pytest.raises(ValueError, match=missing):
        _pack(tmp_path, CITED, sources)


def test_the_unreviewed_and_local_only_codes_are_still_refused(tmp_path):
    for code in ("V", "X", "CC-BY"):
        with pytest.raises(ValueError, match="O, SA or G"):
            _pack(tmp_path, {"gene_id": "O", "tg_lung_log2fc": code}, SOURCE)


def test_an_installed_pack_is_rechecked_for_its_citations(tmp_path):
    """The manifest is what an installed pack is read back from, so the check lives there too."""
    sha, out = _pack(tmp_path, CITED, SOURCE)
    manifest = json.loads(zipfile.ZipFile(out).read(packs.MANIFEST))
    manifest["column_sources"] = {}
    with pytest.raises(ValueError, match="no recorded source"):
        packs._manifest(json.dumps(manifest))


# --------------------------------------------------------------------------- the link rule
def test_the_link_tables_are_many_to_many_in_both_directions(installed):
    _need("Hs", installed)
    pairs = host_links.read_links("Hs").dropna(subset=["gene_id"])
    assert (pairs.groupby("host_id").gene_id.nunique() > 1).sum() > 0
    assert (pairs.groupby("gene_id").host_id.nunique() > 1).sum() > 0


def test_a_selection_carries_as_the_union_of_everything_it_names(installed):
    _need("Hs", installed)
    pairs = host_links.read_links("Hs").dropna(subset=["gene_id"])
    counts = pairs.groupby("host_id").gene_id.nunique()
    fans = str(counts.idxmax())
    expected = set(pairs.loc[pairs.host_id.astype(str) == fans, "gene_id"].astype(str))
    carried = host_links.genes_of("Hs", [fans])
    assert set(carried.ids) == expected and len(expected) > 1
    assert carried.fanned_out == 1 and carried.one_to_one == 0 and carried.unresolved == []
    assert carried.grew and "named several" in carried.summary()


def test_an_endpoint_that_names_nothing_carries_nothing_and_is_counted(installed):
    _need("Hs", installed)
    table = host_links.read_links("Hs")
    orphan = str(table.loc[table.gene_id.isna(), "host_id"].iloc[0])
    carried = host_links.genes_of("Hs", [orphan])
    assert carried.ids == [] and carried.unresolved == [orphan]
    assert "named none and carried nothing" in carried.summary()


def test_carrying_works_the_other_way_round_too(installed):
    _need("Hs", installed)
    pairs = host_links.read_links("Hs").dropna(subset=["gene_id"])
    counts = pairs.groupby("gene_id").host_id.nunique()
    shared = str(counts.idxmax())
    carried = host_links.proteins_of("Hs", [shared])
    expected = set(pairs.loc[pairs.gene_id.astype(str) == shared, "host_id"].astype(str))
    assert set(carried.ids) == expected and len(expected) > 1
    assert carried.direction == "genes to proteins"


def test_carrying_out_and_back_is_not_a_round_trip_and_says_so(installed):
    """One protein, several genes, several proteins back. The result reports the growth, not identity."""
    _need("Hs", installed)
    pairs = host_links.read_links("Hs").dropna(subset=["gene_id"])
    counts = pairs.groupby("gene_id").host_id.nunique()
    shared = str(counts.idxmax())
    one = str(pairs.loc[pairs.gene_id.astype(str) == shared, "host_id"].iloc[0])
    there = host_links.genes_of("Hs", [one])
    back = host_links.proteins_of("Hs", there.ids)
    # Against the ORIGINAL selection, not against the hop. `grew` compares one carry with what was
    # put into that carry, so the second hop (5 genes in, 5 proteins out) reports no growth while
    # the round trip has turned one protein into five. The test says what it means.
    assert one in back.ids and len(back.ids) > len(there.asked) == 1
    assert back.shared and back.fanned_out


def test_the_mouse_link_is_one_to_one_where_it_exists(installed):
    """A rule written only against this table would look right and be wrong in human."""
    _need("Mm", installed)
    pairs = host_links.read_links("Mm").dropna(subset=["gene_id"])
    assert (pairs.groupby("host_id").gene_id.nunique() == 1).all()
    assert (pairs.groupby("gene_id").host_id.nunique() == 1).all()


def test_carrying_nothing_carries_nothing(installed):
    _need("Hs", installed)
    carried = host_links.genes_of("Hs", [])
    assert carried.ids == [] and not carried.grew and "0 carried" in carried.summary()


# --------------------------------------------------------------------------- the lenses
@pytest.mark.parametrize("code", list(HG.LENSES))
def test_both_lenses_of_a_gene_space_load_beside_the_table_they_were_built_from(code, installed):
    _need(code, installed)
    from parascape import paths
    for lens in HG.LENSES[code]:
        nodes, xyz, meta = HG.load(code, lens)
        assert len(nodes) == GENES[code]
        assert xyz.shape == (len(nodes), 3) and np.isfinite(xyz).all()
        assert meta["space"] == code and meta["lens"] == lens and meta["features"]
        node_file = paths.space_file(HG.space(code), HG.space(code).nodes)
        assert meta["node_file_sha256"] == hashlib.sha256(open(node_file, "rb").read()).hexdigest()


def test_the_two_lenses_place_the_genes_differently(installed):
    _need("Hs", installed)
    _, tissue, tmeta = HG.load("Hs", "tissue")
    _, facing, fmeta = HG.load("Hs", "parasite_facing")
    assert not np.allclose(tissue, facing)
    assert set(tmeta["features"]) != set(fmeta["features"])


def test_a_layout_built_from_another_table_is_refused(installed, tmp_path, monkeypatch):
    _need("Mm", installed)
    from parascape import paths
    space = HG.space("Mm")
    real_nodes = paths.space_file(space, space.nodes)
    changed = pd.read_parquet(real_nodes)
    changed.loc[0, "bgee_top_score"] = 12345.0
    changed.to_parquet(tmp_path / space.nodes, index=False)
    shutil.copyfile(paths.space_file(space, space.graph), tmp_path / space.graph)
    monkeypatch.setattr(paths, "space_file", lambda s, name: str(tmp_path / name))
    with pytest.raises(ValueError, match="different"):
        HG.load("Mm", "tissue")


def test_a_space_with_no_installed_pack_offers_no_lens(monkeypatch, tmp_path):
    from parascape import paths
    def absent(space, name):
        raise FileNotFoundError(f"{space.species} has no installed data pack")
    monkeypatch.setattr(paths, "space_file", absent)
    assert HG.lenses("Hs") == [] and HG.available() == []


def test_the_gene_lenses_make_the_same_two_claims_as_the_protein_lenses():
    """Same keys, same wording: two kinds of host space must not read as two different questions."""
    from parascape.spaces import host_proteins as HP
    assert set(HG.LENSES["Hs"]) == {"tissue", "parasite_facing"}
    assert HG.LENS_NAMES == HP.LENS_NAMES


def test_the_mouse_parasite_facing_lens_is_built_for_genes_though_it_was_refused_for_proteins(
        installed):
    """The protein table carries one such column; the gene table carries ten, so the lens is real."""
    from parascape import organisms as O
    from parascape.spaces import host_proteins as HP
    _need("Mm", installed)
    assert (O.MOUSE_PROTEINS, "parasite_facing") in HP.REFUSED
    _nodes, _xyz, meta = HG.load("Mm", "parasite_facing")
    assert len(meta["features"]) >= 5
    assert meta["genes_with_any_value"] / meta["genes"] > 0.5


# --------------------------------------------------------------------------- the window
@pytest.fixture(scope="module")
def qapp():
    pytest.importorskip("PySide6")
    from PySide6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def test_a_host_gene_space_opens_through_each_lens(qapp, installed):
    from parascape.app import Window, HOST_GENE_SPECIES
    _need("Mm", installed)
    name = next(n for n, w in HOST_GENE_SPECIES.items() if w["code"] == "Mm")
    w = Window(species=name, lens="tissue")
    try:
        assert w.n == GENES["Mm"] and w.lens == "tissue" and w.edges == {}
        assert "by tissue" in w.windowTitle() and w.lens_features
        assert w.host_lens_actions[("Mm", "tissue")].isChecked()
        other = w.open_species(name, lens="parasite_facing")
        try:
            assert other is not w and other.lens == "parasite_facing"
            assert np.array_equal(other.nodes.gene_id, w.nodes.gene_id)
            assert not np.allclose(other.xyz, w.xyz)
        finally:
            other.close()
    finally:
        w.close()


def test_the_default_colour_of_a_host_map_is_a_measurement_the_lens_used(qapp, installed):
    """The third host gap: it used to be the first category column, a flag true for a few hundred rows."""
    from parascape.app import Window, HOST_GENE_SPECIES, BIN_PREFIX
    _need("Mm", installed)
    name = next(n for n, w in HOST_GENE_SPECIES.items() if w["code"] == "Mm")
    w = Window(species=name, lens="parasite_facing")
    try:
        assert w.category.startswith(BIN_PREFIX)
        assert w.category[len(BIN_PREFIX):] in w.lens_features
        assert w.category_values().notna().any()
    finally:
        w.close()


def test_the_human_protein_default_colour_is_no_longer_the_senegal_flag(qapp):
    from parascape import organisms as O
    from parascape.app import Window, BIN_PREFIX
    from parascape.spaces import host_proteins as HP
    if O.HUMAN_PROTEINS not in HP.available():
        pytest.skip("the human protein table or its layouts are absent")
    w = Window(species=HP.HUMAN_PROTEIN_SPACE.species, lens="tissue")
    try:
        assert w.category != "rbc_surface_found_senegal"
        assert w.category.startswith(BIN_PREFIX)
    finally:
        w.close()


def test_the_map_module_names_the_lens_features_and_refuses_to_rebuild(qapp, installed):
    """The first two host gaps: the similarity box lists the lens, and Build this map is off."""
    from parascape.app import Window, HOST_GENE_SPECIES
    from parascape.map_module import MapModule
    _need("Mm", installed)
    name = next(n for n, w in HOST_GENE_SPECIES.items() if w["code"] == "Mm")
    w = Window(species=name, lens="parasite_facing")
    try:
        module = MapModule(w)
        assert module.host_lens == "parasite_facing"
        assert module.similarity == {}
        shown = [module.lens_list.item(i).text() for i in range(module.lens_list.count())]
        assert shown == w.lens_features and shown
        build = module.shell.actions["Build this map"]
        assert not build.isEnabled() and "parasite measurement blocks" in build.toolTip()
        assert not module.shell.actions["Find the best map"].isEnabled()
        assert not module.live.isChecked() and not module.live.isEnabled()
    finally:
        w.close()


def test_a_parasite_map_still_offers_its_measurement_blocks(qapp):
    from parascape.app import Window
    from parascape.map_module import MapModule
    w = Window()
    try:
        module = MapModule(w)
        assert module.host_lens == "" and module.similarity
        assert module.shell.actions["Build this map"].isEnabled()
        assert module.live.isChecked()
    finally:
        w.close()


def test_a_selection_carries_from_a_protein_window_to_its_genes_and_back(qapp, installed):
    from parascape import organisms as O
    from parascape.app import Window, HOST_GENE_SPECIES
    from parascape.spaces import host_proteins as HP
    _need("Mm", installed)
    if O.MOUSE_PROTEINS not in HP.available():
        pytest.skip("the mouse protein table or its layouts are absent")
    pairs = host_links.read_links("Mm").dropna(subset=["gene_id"])
    w = Window(species=HP.MOUSE_PROTEIN_SPACE.species, lens="tissue")
    try:
        wanted = set(pairs.host_id.astype(str).head(40))
        rows = np.flatnonzero(w.nodes.gene_id.astype(str).isin(wanted).to_numpy())[:5]
        w.on_gated(rows)
        picked = set(w.nodes.gene_id.astype(str).iloc[rows])
        genes = Window.carry_selection(w)
        try:
            assert genes is not w
            assert genes.species == next(n for n, s in HOST_GENE_SPECIES.items()
                                         if s["code"] == "Mm")
            expected = set(pairs.loc[pairs.host_id.astype(str).isin(picked), "gene_id"].astype(str))
            assert set(genes.nodes.gene_id.astype(str).iloc[genes.gated]) == expected
            # And back the other way, from the gene window to its proteins.
            proteins = genes.carry_selection()
            try:
                assert proteins is not genes
                assert proteins.species == HP.MOUSE_PROTEIN_SPACE.species
                assert picked <= set(proteins.nodes.gene_id.astype(str).iloc[proteins.gated])
            finally:
                proteins.close()
        finally:
            genes.close()
    finally:
        w.close()


def test_carrying_nothing_leaves_the_window_where_it_is(qapp, installed):
    from parascape.app import Window, HOST_GENE_SPECIES
    _need("Mm", installed)
    name = next(n for n, w in HOST_GENE_SPECIES.items() if w["code"] == "Mm")
    w = Window(species=name)
    try:
        w.gated = None
        w.sel = None
        assert w.carry_selection() is w
    finally:
        w.close()


def test_a_parasite_window_has_no_space_to_carry_a_selection_to(qapp):
    from parascape.app import Window
    w = Window()
    try:
        w.on_gated(np.array([0, 1, 2]))
        assert w.carry_selection() is w
    finally:
        w.close()
