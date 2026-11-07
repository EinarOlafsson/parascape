"""Cross-space evidence links (instruction 64.24): organism-qualified, relation types kept apart."""
from __future__ import annotations

from parascape import bridges as B


def test_every_link_names_both_organisms_and_one_relation_type():
    links = B.all_links()
    assert list(links.columns) == B.COLUMNS
    assert links["source_organism"].notna().all() and links["target_organism"].notna().all()
    kinds = set(links["relation"])
    assert kinds <= {"measured interaction (curated literature)", "measured interaction (crosslinking MS)",
                     "orthology (shared OrthoMCL group)"}


def test_measured_and_inferred_links_are_never_merged():
    cov = B.coverage().set_index("relation")
    assert "orthology (shared OrthoMCL group)" in cov.index
    assert cov.loc["orthology (shared OrthoMCL group)", "links"] > 1000
    host = B.host_links()
    assert (host["relation"].str.startswith("measured interaction")).all()


def test_a_host_species_is_never_guessed():
    host = B.host_links()
    lit = host[host["relation"] == "measured interaction (curated literature)"]
    assert (lit["target_organism"] == B.HOST_UNRECORDED).all()


def _human_pack(tmp_path, monkeypatch):
    """A two-gene human pack installed into a scratch cache, with one protein linked to two genes."""
    import numpy as np
    import pandas as pd
    from parascape import packs, paths
    from parascape.spaces import build_pack
    from parascape.spaces.human import HUMAN_SPACE
    monkeypatch.setenv(paths.ENV_STATE, str(tmp_path / "state"))
    data = tmp_path / "data"
    data.mkdir()
    nodes = pd.DataFrame({"gene_id": ["ENSG00000000001", "ENSG00000000002"], "gene_name": ["AAA", "BBB"]})
    nodes.to_parquet(data / HUMAN_SPACE.nodes, index=False)
    np.savez(data / HUMAN_SPACE.graph, xyz=np.zeros((2, 3)), gene_ids=nodes.gene_id.to_numpy(dtype=str))
    links = pd.DataFrame({"host_id": ["P11111", "P11111", "P22222"],
                          "gene_id": ["ENSG00000000001", "ENSG00000000002", None]})
    links.to_parquet(data / "hs_protein_gene_links.parquet", index=False)
    archive = tmp_path / "hs.zip"
    sha = build_pack(HUMAN_SPACE, data, "1", {c: "O" for c in nodes}, archive,
                     extra={"hs_protein_gene_links.parquet": data / "hs_protein_gene_links.parquet"})
    packs.install(archive, sha, "Hs")


def _link(target_organism, target_id, source_id="PF3D7_0000001", source_organism="Pf"):
    return dict(source_organism=source_organism, source_id=source_id, target_organism=target_organism,
                target_id=target_id)


def test_a_host_protein_opens_as_host_genes_in_its_own_space(tmp_path, monkeypatch):
    _human_pack(tmp_path, monkeypatch)
    c = B.counterpart(_link("human", "AAA (P11111)"))
    assert c["space"] == "Hs" and c["genes"] == ["ENSG00000000001", "ENSG00000000002"]
    assert len(c["rows"]) == 2 and c["space_genes"] == 2      # the human denominator, not the parasite one
    lost = B.counterpart(_link("human", "XXX (P22222)"))
    assert lost["space"] == "Hs" and lost["genes"] == [] and "no reviewed gene" in lost["status"]


def test_without_a_host_pack_or_a_recorded_species_no_space_is_chosen(tmp_path, monkeypatch):
    from parascape import paths
    monkeypatch.setenv(paths.ENV_STATE, str(tmp_path / "empty"))
    assert B.counterpart(_link("human", "AAA (P11111)"))["space"] is None
    assert "not installed" in B.counterpart(_link("human", "AAA (P11111)"))["status"]
    unrecorded = B.counterpart(_link(B.HOST_UNRECORDED, "IRGB6", "TGME49_000001", "Tg"))
    assert unrecorded["space"] is None and "does not say" in unrecorded["status"]


def test_a_parasite_end_opens_in_its_own_map_from_either_side():
    o = B.orthology_links()
    if not len(o):
        return
    row = o.iloc[0]
    c = B.counterpart(row)
    assert c["space"] == "Pf" and c["genes"] == [row["target_id"]] and len(c["rows"]) == 1
    back = B.counterpart(row, from_gene=row["target_id"])
    assert back["space"] == "Tg" and back["genes"] == [row["source_id"]]
