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
