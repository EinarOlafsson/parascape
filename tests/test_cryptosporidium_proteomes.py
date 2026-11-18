"""Instruction 41, 2026-11-18: two Cryptosporidium slots filled from PRIDE deposits, and what they are not.

`Cp_protein abundance · other life stages` <- PXD023935 (oocyst total proteome) and
`Cp_secretome / excreted` <- PXD028969 (proteins quantified in the medium of sporozoites on HCT-8 cells).
A column keyed to a public deposit ships as licence G only with its accession and its paper, and the
identifier family of each source must be declared rather than assumed: UniProt Iowa II accessions reach
`cgd` exactly, `CPATCC_` tags only through the 2,571-gene map.
"""
from __future__ import annotations

import importlib.util
import json
import os

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPACE = os.path.join(ROOT, "results", "cryptosporidium_space_2026_11_08")
COLUMNS = {"protein_oocyst_ppm": ("PXD023935", "36775128"),
           "secreted_sporozoite_detected": ("PXD028969", "41137792")}


@pytest.fixture(scope="module")
def nodes():
    return pd.read_parquet(os.path.join(SPACE, "cp_gene_nodes.parquet"))


@pytest.fixture(scope="module")
def meta():
    return json.load(open(os.path.join(SPACE, "column_metadata.json")))


def test_every_g_column_names_its_accession_and_its_paper(meta):
    for column, (accession, pmid) in COLUMNS.items():
        record = meta[column]
        assert record["licence"].startswith("G")
        assert record["accession"] == accession and record["pmid"] == pmid
        assert accession in record["citation"] and f"PMID {pmid}" in record["citation"]
    for column, record in meta.items():
        if record["licence"].startswith("G"):
            assert record.get("accession") and record.get("citation"), column


def test_the_table_still_has_one_row_per_gene_and_the_values_are_cgd_keyed(nodes):
    assert len(nodes) == 3923 and nodes.gene_id.is_unique
    assert nodes.gene_id.str.fullmatch(r"cgd\d+_\d+").all()
    assert nodes.protein_oocyst_ppm.notna().sum() == 2041
    assert nodes.secreted_sporozoite_detected.notna().sum() == 536


def test_the_oocyst_proteome_is_a_share_not_a_count_and_is_not_onto(nodes):
    values = nodes.protein_oocyst_ppm.dropna()
    assert (values > 0).all()
    # parts per million of the identified mass: groups that do not map keep their share of the
    # denominator, so the shipped values sum to LESS than a million and the shortfall is the unmapped mass.
    assert 800_000 < values.sum() < 1_000_000
    assert nodes.protein_oocyst_ppm.isna().sum() > 1500, "missing is unknown, never filled with a zero"


def test_the_secretome_is_a_detection_flag_with_no_invented_negatives(nodes):
    flag = nodes.secreted_sporozoite_detected.dropna()
    assert set(flag) == {1.0}
    assert nodes.secreted_sporozoite_detected.isna().sum() == len(nodes) - 536


def test_the_slots_are_answered_by_these_columns_only():
    import sys
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import generate_slot_table as G
    assert G.SPACE_PATTERNS["Cp"]["protein abundance · other life stages"] == ("protein_oocyst_ppm",)
    assert G.SPACE_PATTERNS["Cp"]["secretome / excreted"] == ("secreted_sporozoite_detected",)


def _script():
    spec = importlib.util.spec_from_file_location(
        "cp_deposits", os.path.join(ROOT, "scripts", "build_cryptosporidium_proteome_deposits.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_a_group_is_assigned_to_one_gene_or_to_none():
    """The mapping rule: every token that can be mapped must agree; nothing is mapped by proximity."""
    assign = _script()._assign
    uniprot = {"A3FPL8": "cgd7_10", "A3FPM1": "cgd7_20"}
    cpatcc = {"CPATCC_0005950": "cgd7_10", "CPATCC_0006170": "cgd7_99"}
    assert assign("CPATCC_0005950-RA-p1;A3FPL8", uniprot, cpatcc) == ("cgd7_10", "both agree")
    assert assign("A3FPL8;Q99999", uniprot, cpatcc) == ("cgd7_10", "UniProt")
    assert assign("CPATCC_0005950-RA-p1", uniprot, cpatcc) == ("cgd7_10", "CPATCC map")
    # the two families disagree about one group: dropped, not resolved in favour of either
    assert assign("CPATCC_0006170-RA-p1;A3FPM1", uniprot, cpatcc) == (None, "ambiguous")
    assert assign("A3FPL8;A3FPM1", uniprot, cpatcc) == (None, "ambiguous")
    assert assign("Q11111;CPATCC_9999999-RA-p1", uniprot, cpatcc) == (None, "unmapped")
