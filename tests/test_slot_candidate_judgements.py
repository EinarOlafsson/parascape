"""Instruction 41: the judgements made on slot candidates, held against the deposits.

Acceptance for instruction 41 asks that every candidate "resolves, its title matches, and it names
the right organism -- assert it in a test that can be re-run offline against a cached response".
This is that test. The cached responses in ``tests/data/geo_candidate_cache/`` are the sample
records exactly as GEO served them on 2026-11-13, reduced to the declaration lines.

**Why a test and not only a note.** Instruction 41 records refusals so "the same searches are not
run again expecting a different answer". A note cannot enforce that; a search run next month returns
the same series with the same inviting title, and nothing stops it being accepted a second time. The
evidence that refused it is cached here so the refusal can be checked rather than remembered.

Every assertion below is a statement about what the SUBMITTER declared, not about what a title
promises. Three candidates in two days have promised a compartment the deposit does not contain.
"""
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "data", "geo_candidate_cache")


def records(accession: str) -> list:
    """(title, declarations) per sample, from the cached GEO reply."""
    out, title, said = [], None, []
    with open(os.path.join(CACHE, f"{accession}.txt")) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("!Sample_title"):
                if title is not None:
                    out.append((title, said))
                title, said = line.split("=", 1)[1].strip(), []
            elif title is not None and "=" in line:
                said.append(line.split("=", 1)[1].strip().lower())
    if title is not None:
        out.append((title, said))
    return out


def declared(accession: str, needle: str) -> int:
    """How many samples declare this, anywhere in their characteristics."""
    return sum(1 for _t, said in records(accession)
               if any(needle in value for value in said))


# ----------------------------------------------------------------- the two refusals


def test_the_anopheles_candidate_is_salivary_gland_and_never_midgut():
    """GSE284425 was accepted for `Pf_host transcriptome - Anopheles midgut` "pending a tissue
    check". The tissue check fails: every sample is salivary gland, which is the other end of the
    vector and a different parasite stage."""
    n = len(records("GSE284425"))
    assert n == 62, n
    assert declared("GSE284425", "salivary gland") == 62
    assert declared("GSE284425", "midgut") == 0, \
        "a midgut sample has appeared; the refusal needs revisiting"


def test_the_cryptosporidium_candidate_has_no_oocyst_sample():
    """GSE232438 was accepted for `Cp_transcription - oocyst`. The organism is right and it is a
    genuine lifecycle series, but no SAMPLE is an oocyst -- every mention of the word is in the
    growth protocol, where the oocyst is what was put on the cells."""
    rs = records("GSE232438")
    assert len(rs) == 22, len(rs)
    titles = {t for t, _ in rs}
    assert not any("oocyst" in t.lower() for t in titles), sorted(titles)
    # What it actually contains: a time course, an in vivo arm, male gametes, and a perturbation.
    assert any(t.startswith("crypto") for t in titles)
    assert any("MybM_OE" in t for t in titles)


# ----------------------------------------------------------------- the replacements found


def test_the_replacement_oocyst_candidate_measures_oocysts_rather_than_naming_them():
    """GSE271211 is the difference between the two: here the oocyst is the thing measured. Three
    `Inoculum` samples at time zero are undigested oocysts in triplicate."""
    rs = records("GSE271211")
    assert len(rs) == 15, len(rs)
    assert declared("GSE271211", "cryptosporidium parvum") == 15
    inoculum = [t for t, _ in rs if "inoculum" in t.lower()]
    assert len(inoculum) == 3, inoculum
    assert declared("GSE271211", "gastric effluents") == 6
    assert declared("GSE271211", "ileal effluents") == 6


def test_the_replacement_midgut_candidate_is_midgut_but_single_cell():
    """GSE284537 fits the slot on tissue, organisms and design. The obstacle is the readout, and it
    is recorded here so the decision is taken on what it is rather than rediscovered."""
    rs = records("GSE284537")
    assert len(rs) == 32, len(rs)
    assert declared("GSE284537", "midgut") == 32
    assert declared("GSE284537", "single cell") == 32, \
        "if this is no longer single-cell the slot decision changes"


# ----------------------------------------------------------------- the title trap itself


@pytest.mark.parametrize("accession,promised,actually", [
    ("GSE284425", "midgut", "salivary gland"),
    ("GSE271070", "midgut", "body without head"),
])
def test_a_title_promising_a_compartment_the_samples_do_not_contain(accession, promised, actually):
    """The recurring trap, pinned. GSE271070 is titled "...symbiotic bacteria in the mosquito
    midgut" and every sample is a whole body with the head removed. Judge from the samples."""
    assert declared(accession, actually) == len(records(accession))
    assert declared(accession, promised) == 0
