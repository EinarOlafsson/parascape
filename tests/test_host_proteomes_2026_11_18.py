"""Instruction 41, 2026-11-18: six host slots filled from public proteomes, and what they are not.

Each column ships only with its accession and its paper. The things worth asserting are the ones that
would be wrong to forget: a fibroblast that is not the foreskin line, a bone marrow that is one autopsy
study, a surface list with no negatives, a salivary-gland contrast that no protein passes, and a
Plasmodium hold that none of this touched.
"""
from __future__ import annotations

import json
import os

import pandas as pd
import pytest

from parascape import datasets, deposits, organisms

#: key -> (accession, PMID, organism, columns)
ENTRIES = {
    "host_fibroblast_proteome": ("PXD032893", "35743299", organisms.HUMAN, ("fibroblast_protein_ppm",)),
    "host_monocyte_proteome": ("PXD000561", "24870542", organisms.HUMAN, ("monocyte_protein_ppm",)),
    "host_bone_marrow_proteome": ("PXD020192", "33107741", organisms.HUMAN, ("bone_marrow_protein_ppm",)),
    "host_dermis_proteome": ("PXD027093", "35817779", organisms.HUMAN,
                             ("dermis_superficial_protein_ppm", "dermis_deep_protein_ppm")),
    "host_monocyte_surface": ("PXD013832", "32165698", organisms.HUMAN,
                              ("monocyte_surface_classical_intensity",)),
    "host_anopheles_salivary_gland_proteome": (
        "PXD057628", "41266344", organisms.ANOPHELES,
        ("anopheles_salivary_gland_protein_ppm", "anopheles_salivary_gland_infection_log2fc",
         "anopheles_salivary_gland_infection_padj")),
}
#: slot -> the columns that fill it
SLOTS = {
    "host proteome · human fibroblast": ("fibroblast_protein_ppm",),
    "host proteome · human monocyte": ("monocyte_protein_ppm",),
    "host surface / receptor repertoire · human monocyte": ("monocyte_surface_classical_intensity",),
    "host proteome · human bone marrow": ("bone_marrow_protein_ppm",),
    "host proteome · human dermis": ("dermis_superficial_protein_ppm", "dermis_deep_protein_ppm"),
    "host proteome · Anopheles salivary gland": ("anopheles_salivary_gland_protein_ppm",),
    "host response to infection · Anopheles salivary gland": (
        "anopheles_salivary_gland_infection_log2fc", "anopheles_salivary_gland_infection_padj"),
}
RAW = "/media/carruthers/mnt3/claude/toxoplasma_projects/datasets"
#: These columns are in the repository's own data directory. They are not in the published 0.0.13.0 release,
#: so a run pinned to that release through PARASCAPE_CACHE must still read them from the checkout.
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "parascape", "data")


def _table(code):
    return pd.read_parquet(os.path.join(DATA, organisms.HOST_TABLES[code]))


@pytest.mark.parametrize("key", sorted(ENTRIES))
def test_every_column_carries_its_accession_and_its_paper(key):
    """No accession or no citation means the column does not ship (Einar, 2026-11-13)."""
    accession, pmid, organism, columns = ENTRIES[key]
    entry = datasets.get(key)
    assert entry.accession == accession and entry.pmid == pmid
    assert entry.organism == organism
    assert entry.citation and accession in entry.citation and f"PMID {pmid}" in entry.citation
    assert tuple(entry.columns) == columns
    assert "Redistributed with its accession and paper cited" in entry.note
    for column in columns:
        assert datasets.provenance(column, organism).key == key


def _deposit(key):
    return pd.read_csv(os.path.join(DATA, f"deposit_{key}.tsv"), sep="\t").set_index("host_id")


@pytest.mark.parametrize("key", sorted(ENTRIES))
def test_the_shipped_host_table_holds_each_column_and_the_deposit_matches_it(key):
    accession, pmid, organism, columns = ENTRIES[key]
    table = _table(organism).set_index("host_id")
    deposit = _deposit(key)
    for column in columns:
        assert column in table.columns
        assert table[column].notna().sum() == deposit[column].notna().sum() > 0
        assert set(deposit.index[deposit[column].notna()]) == set(table.index[table[column].notna()])


def test_no_host_table_grew_and_the_plasmodium_hold_is_untouched():
    """The host tables ARE the gene space of each species; a measurement outside it is dropped, not added."""
    assert len(_table("Hs")) == 21_020 and len(_table("Mm")) == 16_473
    assert len(_table(organisms.ANOPHELES)) == 12_574
    assert set(deposits.HELD_OUT_OF_NODES) == {"pf_liver_transcriptome", "pf_translation_efficiency",
                                               "pf_dbsp_occupancy", "pf_midgut_transcriptome"}
    assert not {d.key for d in deposits.DEPOSITS if d.organism == "Pf"} & set(ENTRIES)


def test_each_slot_is_answered_by_its_own_columns():
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
    import generate_slot_table as G
    for slot, columns in SLOTS.items():
        assert tuple(G.HOST_COLUMNS[slot]) == columns


def test_the_fibroblast_is_labelled_as_skin_vehicle_cells_and_not_the_foreskin_line():
    note = datasets.get("host_fibroblast_proteome").note
    assert "PRIMARY SKIN" in note and "VEHICLE" in note and "not the foreskin line" in note
    assert "one deposit and not a consensus" in note


def test_bone_marrow_says_it_is_one_shallow_study_and_the_join_loses_mass():
    note = datasets.get("host_bone_marrow_proteome").note
    assert "ONE autopsy study" in note and "76.7%" in note
    assert _deposit("host_bone_marrow_proteome").bone_marrow_protein_ppm.notna().sum() < 3000


def test_the_dermis_keeps_its_two_layers_apart():
    d = _deposit("host_dermis_proteome")
    assert {"dermis_superficial_protein_ppm", "dermis_deep_protein_ppm"} <= set(d.columns)
    assert "NOT used" in datasets.get("host_dermis_proteome").note


def test_the_monocyte_surface_list_has_no_zeros_because_it_has_no_negatives():
    """A protein the experiment did not quantify is unknown, so no value is zero and none is invented."""
    values = _deposit("host_monocyte_surface").monocyte_surface_classical_intensity.dropna()
    assert (values > 0).all() and 400 < len(values) < 500
    note = datasets.get("host_monocyte_surface").note
    assert "UNKNOWN and not absent" in note and "not proven exposed" in note


def test_the_salivary_gland_is_a_berghei_infection_at_21_days_and_no_protein_passes_005():
    d = _deposit("host_anopheles_salivary_gland_proteome")
    assert d.index.str.fullmatch(r"AGAP\d{6}").all() and d.index.is_unique
    padj = d.anopheles_salivary_gland_infection_padj.dropna()
    assert len(padj) > 2000 and padj.min() > 0.05, "the registry note says nothing passes 0.05"
    note = datasets.get("host_anopheles_salivary_gland_proteome").note
    assert "21 days" in note and "BERGHEI" in note and "NO protein passes 0.05" in note
    assert d.anopheles_salivary_gland_protein_ppm.dropna().gt(0).all()


@pytest.mark.skipif(not os.path.isdir(RAW), reason="the raw dataset archive is not on this machine")
def test_the_symbol_join_keeps_at_least_93_per_cent_of_each_source_file():
    """Reported per column because the join is by symbol; a floor so a regression is visible."""
    from parascape import host
    index = host.uniprot_index(RAW, "human", log=lambda *a: None)["by_symbol"]
    for name in ("9606-SKIN_FIBROBLAST-integrated.txt", "9606-Adult_Monocytes_Kim_2014_SEQUEST.txt",
                 "9606-PXD020192_bone_marrow.txt", "9606-PXD027093_deep_dermis.txt",
                 "9606-PXD027093_superficial_dermis.txt"):
        t = pd.read_csv(os.path.join(RAW, "host", "paxdb", name), sep="\t", comment="#", usecols=[0, 1, 2],
                        names=["g", "e", "a"]).dropna(subset=["g", "a"])
        assert t.g.astype(str).str.strip().map(index).notna().mean() >= 0.93, name


def test_a_machine_without_the_archive_gets_empty_frames_not_half_built_columns():
    for fn in (deposits.human_fibroblast_proteome, deposits.human_monocyte_proteome,
               deposits.human_bone_marrow_proteome, deposits.human_dermis_proteome,
               deposits.human_monocyte_surface, deposits.anopheles_salivary_gland_proteome):
        assert fn("/nonexistent").empty
