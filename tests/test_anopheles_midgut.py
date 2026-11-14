"""Instruction 41: the Anopheles midgut host transcriptome, and what it is not.

`Pf_host transcriptome · Anopheles midgut` is filled from GSE284537 (PMID 41125888), the only
genuinely midgut series in the archive and a single-cell one. Einar admitted it on 2026-11-13 as a
per-sample pseudobulk. These tests assert the things that would be wrong to forget: that a pseudobulk
is labelled as a pseudobulk, that the knockdown arm did not get into the slot, that the table is
keyed on AGAP and not on UniProt, and that the mapping is INCOMPLETE -- a reader who assumed the
column covered the genome would be reading 40 absent genes as zero.
"""
from __future__ import annotations

import pandas as pd
import pytest

from parascape import datasets, deposits, host, organisms, paths, slots

SLOT = "host transcriptome · Anopheles midgut"
COLUMNS = ("anopheles_midgut_36h_cpm", "anopheles_midgut_48h_cpm", "anopheles_midgut_4d_cpm",
           "anopheles_midgut_7d_cpm")


@pytest.fixture(scope="module")
def table():
    return pd.read_parquet(paths.cache_file(organisms.HOST_TABLES[organisms.ANOPHELES]))


def test_the_anopheles_host_table_keys_on_agap_and_not_on_uniprot(table):
    """The identifier contract this source changed, asserted rather than described.

    A. gambiae has about 300 reviewed UniProt entries against 14,559, which is why the midgut slots
    read `unreachable` from 2026-08-22. A Swiss-Prot key here would silently cover two per cent of
    the organism, so this table keys on VectorBase gene ids and the test says so.
    """
    assert organisms.HOST_TABLES[organisms.ANOPHELES] == "ag_host_genes.parquet"
    assert table.host_id.str.fullmatch(r"AGAP\d{6}").all()
    assert table.host_id.is_unique
    for code in ("Hs", "Mm"):
        other = pd.read_parquet(paths.cache_file(organisms.HOST_TABLES[code]))
        assert not set(table.host_id) & set(other.host_id)


def test_the_slot_resolves_to_the_four_control_timepoints(table):
    """Filled, and filled by the CONTROL arm only."""
    declared = [s for s in slots.all_slots("Pf") if s.name == SLOT]
    assert len(declared) == 1
    assert tuple(declared[0].patterns) == COLUMNS
    assert set(COLUMNS) <= set(table.columns)
    # The series pairs every timepoint with a dsEcR knockdown. A knockdown is a perturbation and a
    # different question, so no column of this table may name it.
    assert not [c for c in table.columns if "ecr" in c.lower() or "knock" in c.lower()]


def test_every_value_is_a_share_of_the_whole_library_and_not_of_the_kept_genes(table):
    """Counts per million of ALL mosquito molecules, so the four timepoints compare.

    The protein-coding build does not carry every feature the submitters quantified, so these
    columns sum to 572,054 to 745,639 per million rather than to a million. That gap is the
    non-coding share of the library and it is left visible on purpose: renormalising to the kept
    genes would quietly turn "a quarter of this library is not in the gene build" into nothing.
    """
    totals = table[list(COLUMNS)].sum()
    assert (totals > 500_000).all() and (totals < 800_000).all()
    assert (table[list(COLUMNS)] >= 0).all().all()


def test_the_mapping_is_not_complete_and_the_shortfall_stays_missing(table):
    """40 genes of the Anopheles gene space have no value here, and missing is what they get.

    The deposit is keyed on the AGAP ids of the Cell Ranger reference the submitters built (13,845
    features) and the space on the Ensembl Metazoa AgamP4 protein-coding build (12,614 genes). The
    join is on the identifier itself, so there is no mapping table to be wrong -- but it is not
    onto, and a reader who assumed it was would read 40 absences as zeroes.
    """
    space = pd.read_parquet("results/host_gene_spaces_2026_11_09/ag_gene_nodes.parquet",
                            columns=["gene_id"])
    assert len(table) == 12_574 < len(space) == 12_614
    assert set(table.host_id) < set(space.gene_id)


def test_the_registry_entry_says_pseudobulk_and_cites_its_accession_and_paper():
    """GEO-derived columns ship with their accession and paper cited (Einar, 2026-11-13), and a
    pseudobulk says it is one: a reader must not take this for a bulk library."""
    entry = datasets.get("host_anopheles_midgut_transcriptome")
    assert entry.organism == organisms.ANOPHELES
    assert entry.accession == "GSE284537" and entry.pmid == "41125888"
    assert tuple(entry.columns) == COLUMNS
    assert "PSEUDOBULK" in entry.note and "not the same thing as a bulk library" in entry.note
    for column in COLUMNS:
        assert datasets.provenance(column, organisms.ANOPHELES).key == entry.key


def test_the_deposit_reader_returns_nothing_without_the_raw_table():
    """A machine without the dataset root gets an empty frame, never a half-built column."""
    assert deposits.anopheles_midgut_transcriptome("/nonexistent").empty
    assert organisms.ANOPHELES in host.HOST_TABLES
