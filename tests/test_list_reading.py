"""Reading a LIST of genes (instruction 98, WP9), on a made-up table.

The rules these hold `parascape.list_reading` to, which are the rules of the module mould:
every share is reported against the share outside the list with the number both rest on,
absence is "unknown" and never a "no", percentages are whole numbers, and nothing handed in
is dropped without being named.
"""
import numpy as np
import pandas as pd
import pytest

from parascape import list_reading as LR


@pytest.fixture
def nodes():
    """600 genes: a compartment, a boolean, a product and a number.

    Built so that `dense granules` is rare overall (1 gene in 10) and so a list of them is
    unmistakably unusual, which is the thing the card exists to notice.
    """
    rng = np.random.default_rng(0)
    n = 600
    return pd.DataFrame({
        "gene_id": [f"TGME49_{i:06d}" for i in range(n)],
        "compartment": rng.choice(["nucleus - chromatin", "dense granules", "cytosol", ""],
                                  n, p=[0.4, 0.1, 0.3, 0.2]),
        "host_facing": rng.choice([True, False], n),
        "product": [f"hypothetical protein {i}" for i in range(n)],
        "fitness_score": rng.normal(size=n),
    })


@pytest.fixture
def granules(nodes):
    """Twenty dense-granule genes and five cytosolic ones: a list with something in common."""
    dg = nodes.index[nodes.compartment == "dense granules"][:20]
    other = nodes.index[nodes.compartment == "cytosol"][:5]
    return [str(g) for g in nodes.gene_id.iloc[list(dg) + list(other)]]


# ------------------------------------------------------------------ what is a list
@pytest.mark.parametrize("text,expected", [
    ("A, B, C", ["A", "B", "C"]),
    ("A\nB\nC", ["A", "B", "C"]),                 # a column pasted out of a spreadsheet
    ("A\tB; C  D", ["A", "B", "C", "D"]),
    ("A, a, A", ["A"]),                           # the same thing twice is one thing
    ("", []),
    ("   ", []),
])
def test_identifiers_are_split_however_they_were_typed(text, expected):
    assert LR.split_identifiers(text) == expected


def test_one_identifier_is_not_a_list():
    """What makes something a list is that it names several things. One is a search."""
    assert not LR.looks_like_a_list("TGME49_208830")
    assert LR.looks_like_a_list("TGME49_208830, TGME49_223920")


def test_what_names_no_gene_is_reported_and_not_dropped(nodes, granules):
    resolved = LR.resolve(nodes, granules + ["NOT_A_GENE", granules[0]])
    assert resolved.count == 25
    assert resolved.missing == ("NOT_A_GENE",)
    assert resolved.repeated == 1
    said = resolved.sentence()
    assert "NOT_A_GENE" in said and "not dropped" in said and "25 genes" in said


def test_identifiers_are_matched_ignoring_case_and_spelled_as_the_table_does(nodes):
    one = str(nodes.gene_id.iloc[3])
    assert LR.resolve(nodes, [one.lower()]).genes == (one,)


def test_a_list_of_nothing_reads_as_nothing(nodes):
    assert LR.resolve(nodes, []).count == 0
    assert LR.resolve(None, ["A"]).missing == ("A",)


# ------------------------------------------------------ what the list shares
def test_what_the_list_shares_is_reported_against_the_other_genes(nodes, granules):
    shared = LR.shared_labels(nodes, LR.resolve(nodes, granules).genes)
    assert shared, "a list of dense-granule genes shares being in the dense granules"
    top = shared[0]
    assert top.value == "dense granules" and top.column == "compartment"
    assert top.in_list == 20 and top.list_measured == 25
    assert top.list_percent == 80
    assert 0 < top.outside_percent < 20, "the comparison is with the rest of the genes"
    assert top.outside_measured > 100


def test_a_share_is_only_reported_with_what_it_is_a_share_of(nodes, granules):
    """Against chance, with n. Half a list in the nucleus means nothing until you know
    that half the genome is."""
    shared = LR.shared_labels(nodes, LR.resolve(nodes, granules).genes)
    said = LR.shared_sentence(shared[0])
    for part in ("20", "25", "80 %", "random"):
        assert part in said, part
    assert str(shared[0].outside_measured) in said.replace(",", "")
    assert f"{shared[0].outside_percent} %" in said


def test_a_label_no_more_common_in_the_list_than_outside_it_is_not_a_finding(nodes):
    """`host_facing` is a coin toss in this table, so no list is unusual in it."""
    genes = [str(g) for g in nodes.gene_id.iloc[:40]]
    shared = LR.shared_labels(nodes, genes)
    assert all(s.column != "host_facing" for s in shared), [s.value for s in shared]


def test_a_product_name_is_not_something_a_list_shares(nodes, granules):
    """Forty genes sharing forty different products is not a finding, and the column that
    would say so is kept out by how many distinct answers it holds."""
    assert "product" not in LR._label_columns(nodes)


def test_one_finding_said_by_several_measurements_is_one_row(nodes, granules):
    """Five measurements place a gene in the real table and they mostly agree; five rows
    saying "dense granules" is one finding said five times."""
    genes = LR.resolve(nodes, granules).genes
    both = nodes.assign(lopit_map=nodes["compartment"])
    shared = LR.shared_labels(both, genes)
    rows = [s for s in shared if s.value == "dense granules"]
    assert len(rows) == 1, [(s.column, s.value) for s in shared]
    assert rows[0].also, "the agreeing measurement must be named, not dropped"
    assert "lopit map" in LR.shared_sentence(rows[0]) or \
        "compartment" in LR.shared_sentence(rows[0])


def test_how_surprising_is_said_in_plain_words_not_as_a_number(nodes, granules):
    shared = LR.shared_labels(nodes, LR.resolve(nodes, granules).genes)
    words = LR.surprise_words(shared[0])
    assert "random lists" in words
    assert "p" not in words.split() and "hypergeom" not in words


def test_a_list_that_shares_nothing_says_so_rather_than_showing_an_empty_table(nodes):
    genes = [str(g) for g in nodes.gene_id.iloc[:40]]
    resolved = LR.resolve(nodes, genes)
    lead = LR.shared_lead(LR.shared_labels(nodes, genes), resolved)
    assert "not unusual" in lead and "not a failure" in lead


def test_the_bar_a_finding_has_to_clear_is_said_out_loud(nodes, granules):
    """A bar that hides things without saying so is a bar nobody can argue with."""
    resolved = LR.resolve(nodes, granules)
    lead = LR.shared_lead(LR.shared_labels(nodes, resolved.genes), resolved)
    assert f"{LR.MIN_TIMES:g} times" in lead and f"once in {LR.MIN_ONE_IN:,}" in lead


# ------------------------------------------------------ how much is measured
def test_coverage_counts_the_list_and_the_other_genes_both(nodes, granules):
    genes = LR.resolve(nodes, granules).genes
    covered = LR.coverage(nodes, genes)
    assert covered, "the table has columns of a known kind"
    place = next(c for c in covered if c.category == "place")
    assert place.listed == 25 and place.outside == 575
    assert place.in_list + place.unknown == place.listed
    assert 0 <= place.percent <= 100 and 0 <= place.outside_percent <= 100
    assert isinstance(place.percent, int), "whole numbers, as Einar asked"


def test_what_is_not_measured_is_unknown_and_never_a_no(nodes, granules):
    genes = LR.resolve(nodes, granules).genes
    thin = nodes.copy()
    thin.loc[thin.index[:300], "compartment"] = ""
    thin.loc[thin.gene_id.isin(genes), "compartment"] = ""
    covered = LR.coverage(thin.drop(columns=["host_facing"]), genes)
    place = next(c for c in covered if c.category == "place")
    assert place.in_list == 0 and place.mostly_unknown
    said = LR.unknown_sentence(place)
    assert "unknown" in said
    assert " no " not in said.replace("not ", "")


def test_the_emptiest_kinds_of_knowledge_come_first(nodes, granules):
    genes = LR.resolve(nodes, granules).genes
    thin = nodes.copy()
    thin.loc[thin.gene_id.isin(genes), "compartment"] = ""
    poor = LR.unknown_for_most(LR.coverage(thin.drop(columns=["host_facing"]), genes))
    assert poor and poor[0].category == "place"
    assert [c.percent for c in poor] == sorted(c.percent for c in poor)


def test_a_gap_in_the_data_is_told_apart_from_a_gap_in_the_list(nodes, granules):
    """A list nobody measured and a kind of knowledge nobody measured for anything are
    different problems with different answers."""
    genes = LR.resolve(nodes, granules).genes
    everywhere = LR.Coverage("place", "Where it is", 1, 2, 25, 40, 575)
    only_here = LR.Coverage("place", "Where it is", 1, 2, 25, 550, 575)
    assert "the gap is in the data" in LR.unknown_sentence(everywhere)
    assert "measured less than the genome" in LR.unknown_sentence(only_here)


# ------------------------------------------- where the list sits in the parasite
def test_where_the_list_sits_counts_every_compartment_it_is_found_in(nodes, granules):
    sp = LR.spread(nodes, LR.resolve(nodes, granules).genes, "compartment")
    assert sp.commonest == "dense granules"
    assert sp.count_of("dense granules") == 20 and sp.count_of("cytosol") == 5
    assert sp.measured == 25 and sp.listed == 25 and sp.unknown == 0
    assert sp.percent_of("dense granules") == 80
    assert sp.outside_percent_of("dense granules") < 20


def test_a_list_with_no_measured_compartment_says_unknown_and_draws_nothing(nodes):
    unplaced = [str(g) for g in nodes[nodes.compartment == ""].gene_id.iloc[:10]]
    sp = LR.spread(nodes, unplaced, "compartment")
    assert sp.measured == 0 and sp.unknown == 10 and not sp.commonest
    said = LR.spread_sentence(sp)
    assert "unknown" in said and "not the same as them being nowhere" in said


def test_the_drawing_says_how_much_of_the_list_it_is_speaking_for(nodes, granules):
    sp = LR.spread(nodes, LR.resolve(nodes, granules).genes, "compartment")
    trust = LR.spread_trust(sp)
    assert "25" in trust and "compartments" in trust and "dense granules" in trust


def test_the_spread_can_be_read_for_an_organelle_other_than_the_commonest(nodes, granules):
    sp = LR.spread(nodes, LR.resolve(nodes, granules).genes, "compartment")
    said = LR.spread_sentence(sp, "cytosol")
    assert "cytosol" in said and "5" in said


def test_a_column_that_is_not_in_the_table_is_not_an_error(nodes, granules):
    sp = LR.spread(nodes, granules, "no_such_column")
    assert sp.measured == 0 and sp.listed == 0
    assert LR.spread_trust(sp) == ""
