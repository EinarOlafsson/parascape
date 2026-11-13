"""Instruction 41: the correspondence between the two Cryptosporidium annotations.

`CPATCC_` (GCA_015245375.1, IOWA-ATCC) and `cgd` (GCA_000165345.1, Iowa II) are two independent
annotations of two assemblies of one organism, with no shared identifier and no NCBI Gene record to
join on. The map is built from identical protein sequences by `scripts/build_cpatcc_cgd_map.py`.

The tests that matter are not that the file exists. They are the two that could tell a real
correspondence from a plausible-looking one: that no gene is claimed twice, and that the map is
STRUCTURALLY coherent -- the ATCC numbering should walk the cgd chromosomes in long contiguous
runs, because that is what two annotations of the same genome look like. Noise that happened to
match a few thousand sequences would be scattered.
"""
import os
import re

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAP = os.path.join(ROOT, "results", "cryptosporidium_space_2026_11_08", "cpatcc_to_cgd.tsv")

pytestmark = pytest.mark.skipif(not os.path.exists(MAP), reason="the map has not been built")


@pytest.fixture(scope="module")
def pairs():
    m = pd.read_csv(MAP, sep="\t")
    m["atcc_n"] = m.cpatcc_locus_tag.str.extract(r"_(\d+)").astype(int)
    m["chrom"] = m.cgd_gene_id.str.extract(r"cgd(\d+)_").astype(int)
    return m.sort_values("atcc_n").reset_index(drop=True)


def test_no_gene_is_claimed_twice(pairs):
    """A one-to-one map or it is not a correspondence."""
    assert pairs.cpatcc_locus_tag.is_unique
    assert pairs.cgd_gene_id.is_unique


def test_every_identifier_has_the_shape_of_its_annotation(pairs):
    assert pairs.cpatcc_locus_tag.str.match(r"^CPATCC_\d+$").all()
    assert pairs.cgd_gene_id.str.match(r"^cgd\d+_\d+").all()


def test_the_map_walks_the_chromosomes_in_long_blocks_rather_than_scattering(pairs):
    """The structural check, and the one that separates a correspondence from a coincidence.

    Two annotations of one genome number their genes along the same chromosomes, so walking the
    ATCC numbering should cross chromosome boundaries only a handful of times. Scattered single
    genes are the signature of sequence matches that are real individually and wrong as a map --
    repeated domains, families, contamination.
    """
    runs = (pairs.chrom != pairs.chrom.shift()).cumsum()
    sizes = pairs.groupby(runs).size()
    assert len(sizes) <= 20, f"{len(sizes)} blocks is scattering, not a correspondence"
    assert (sizes == 1).sum() == 0, "isolated single-gene blocks mean the map is partly noise"
    assert sizes.median() > 100


def test_all_eight_chromosomes_are_represented(pairs):
    assert set(pairs.chrom) == set(range(1, 9)), sorted(set(pairs.chrom))


def test_the_coverage_is_recorded_rather_than_assumed_complete(pairs):
    """Two thirds, and the file must not pretend otherwise. 3,892 ATCC genes and 3,805 cgd genes
    exist; this maps 2,573 of them. A reader who assumes completeness would silently drop a third
    of the oocyst measurements."""
    assert 2400 <= len(pairs) <= 2700, len(pairs)
    assert len(pairs) < 3805, "a complete map would need alignment, not exact identity"
