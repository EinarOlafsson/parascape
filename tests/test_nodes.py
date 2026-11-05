"""The node model behind the home screen: its shape, its rules and its honesty.

These are the checks that cost something if they fail. A node with no summary is a one-word card
nobody can act on; a claim-bearing node that does not carry its probability, uncertainty and source
tier makes trust a destination again, which is the thing instruction 74 removed a whole category to
avoid; and a node whose data is absent must say so rather than fail when it is clicked.

No Qt here: `parascape.nodes` imports none, and that is worth keeping true.
"""
import pytest

from parascape import nodes as ND
from parascape import question_registry as QR


@pytest.fixture(scope="module")
def catalogue():
    return ND.catalogue()


def test_level_one_is_the_three_shapes_of_the_matrix_and_nothing_else():
    """A row, a block of rows, a column. No `Sources` and no `Space`: provenance is an attribute of
    every value rather than a place, and a UMAP is a view of the rows rather than a kind of query."""
    assert [c.id for c in ND.CATEGORIES] == ["gene", "list", "label"]
    assert [c.title for c in ND.CATEGORIES] == ["Gene", "List", "Label"]
    assert all(c.purpose.strip() for c in ND.CATEGORIES)


def test_every_category_has_cards(catalogue):
    by_category = ND.by_category(catalogue)
    assert list(by_category) == list(ND.CATEGORY_IDS)
    for category_id, cards in by_category.items():
        assert cards, f"{category_id} has no cards"


def test_every_node_says_what_clicking_it_does_in_a_title_of_at_most_two_words(catalogue):
    """One-word titles were Einar's rule; the summary is where the explaining goes, and it is
    mandatory rather than decorative."""
    for node in catalogue:
        assert node.summary.strip(), f"{node.id} has no summary"
        assert len(node.title.split()) <= ND.MAX_TITLE_WORDS, f"{node.title!r} is too long a title"


def test_every_node_resolves_to_an_action(catalogue):
    for node in catalogue:
        action = ND.resolve(node)
        assert action.kind in ND.DESTINATION_KINDS
        assert action.node_id == node.id
        if action.kind == "unavailable":
            assert action.describe == node.unavailable_reason
        else:
            assert action.target, f"{node.id} resolves to no target"
            assert action.describe


def test_the_label_cards_are_in_the_order_einar_settled_on(catalogue):
    """Strategies is third, not an advanced afterthought: his first requirement was that nothing of
    the previous ParaScape be lost, which means all 39 strategies stay first-class."""
    assert [c.title for c in ND.by_category(catalogue)["label"]] == [
        "Catalogue", "Infer", "Strategies", "Scorecard", "Map", "Runs", "Fill", "Compare"]
    assert not [c for c in catalogue if c.advanced], "nothing is filed as advanced any more"


def test_an_advanced_node_still_renders_last_in_its_category():
    """The mechanism outlives the one node that used it, because the next one will want it."""
    ordinary = ND.Node("a", "A", "panel", "gene", "first", ND.Destination("dock", "evidence"))
    advanced = ND.Node("b", "B", "panel", "gene", "last", ND.Destination("dock", "evidence"),
                       advanced=True)
    assert [n.id for n in ND.by_category([advanced, ordinary])["gene"]] == ["a", "b"]


def test_nesting_is_three_deep_and_questions_hang_under_the_catalogue(catalogue):
    assert {n.level for n in catalogue} <= {2, 3}
    questions = ND.children("label_catalogue", catalogue)
    assert questions
    for node in questions:
        assert node.level == 3 and node.kind == "question" and node.question_id
        assert node.category == "label"
    assert all(n.level == 2 for n in catalogue if n.id not in {q.id for q in questions})


def test_a_question_node_exists_for_every_registered_question_except_the_retired_one(catalogue):
    ids = {n.question_id for n in ND.children("label_catalogue", catalogue)}
    expected = {q.id for q in QR.QUESTIONS.values() if q.status != "retired"}
    assert ids == expected
    assert "has_domain_ipr000719" not in ids       # a rejected question is a worked example


def test_a_question_carries_both_numbers_and_its_coverage_inline(catalogue):
    """Trust is integrated: the cross-validated estimate and the frozen test, side by side, on the
    node itself -- not a link to a scoreboard elsewhere."""
    node = ND.get("question_required_without_glucose", catalogue)
    board = node.scoreboard
    assert board is not None and board.scored
    line = board.inline()
    # One decimal, not whole percent: the difference between 89.1 % and 89.9 % decides whether a
    # question clears the bar, and rounding it away in the display hides that.
    assert f"{board.cross_validated_lower:.1%}" in line
    assert f"{board.frozen_test_lower:.1%}" in line
    # Both numbers side by side, in words that need no glossary: the rotated estimate and the one
    # measured on genes hidden from the model.
    assert "when the hidden genes rotate" in line and "on hidden genes" in line
    assert "not a negative" in board.coverage()
    assert str(ND.GENOME_GENES) in board.coverage().replace(",", "")


def test_a_question_that_cannot_be_scored_says_so_with_its_reason(catalogue):
    node = ND.get("question_favourable_immune_response_human", catalogue)
    assert not node.available
    assert "cannot be answered yet" in node.unavailable_reason
    assert ND.resolve(node).kind == "unavailable"


def test_a_question_below_the_frozen_test_threshold_has_no_frozen_number(catalogue):
    """318 labels cannot afford a partition held back from everything, and the node must not
    pretend otherwise."""
    node = ND.get("question_screen_egress_phenotype", catalogue)
    board = node.scoreboard
    assert board.frozen_test_lower is None
    assert not board.may_fill
    assert "not graded on hidden genes" in board.inline()
    assert "not enough known answers to grade" in board.inline()


def test_a_claim_bearing_node_declares_all_three_claim_fields(catalogue):
    showing = [n for n in catalogue if n.shows_claim]
    assert showing
    for node in showing:
        assert tuple(node.claim_fields) == ND.CLAIM_FIELDS
        assert "probability" in node.claim_line() and "source tier" in node.claim_line()


def test_a_claim_without_its_uncertainty_or_tier_is_refused():
    with pytest.raises(ValueError):
        ND.Node("x", "X", "gene", "gene", "shows a claim", ND.Destination("dock", "evidence"),
                shows_claim=True, claim_fields=("probability",))


@pytest.mark.parametrize("kwargs", [
    {"title": "A title that is far too long"},               # three words or more
    {"summary": "   "},                                      # no statement of what it does
    {"kind": "not a kind"},
    {"category": "trust"},                                   # the dropped category stays dropped
    {"available": False},                                    # unavailable with no reason
    {"level": 3},                                            # level 3 with no parent
])
def test_the_node_model_refuses_a_card_a_user_could_not_act_on(kwargs):
    base = dict(id="x", title="X", kind="panel", category="gene", summary="does a thing",
                destination=ND.Destination("dock", "evidence"))
    base.update(kwargs)
    with pytest.raises(ValueError):
        ND.Node(**base)


def test_an_unavailable_node_reports_a_reason_and_resolving_it_does_not_raise():
    missing = ND.Presence(information_map=False,
                          information_map_reason="needs the information map release")
    catalogue = ND.catalogue(presence=missing)
    facts = ND.get("gene_facts", catalogue)
    assert not facts.available
    assert "information map release" in facts.unavailable_reason
    assert ND.resolve(facts) == ND.Action("gene_facts", "unavailable",
                                          describe=facts.unavailable_reason)
    assert all(ND.resolve(n) for n in catalogue)


def test_the_list_nodes_wait_for_a_list_and_say_so(catalogue):
    for node_id in ("list_answer", "list_map", "list_export"):
        node = ND.get(node_id, catalogue)
        assert not node.available and "gene list" in node.unavailable_reason
    loaded = ND.catalogue(presence=ND.Presence(gene_list_loaded=True))
    assert ND.get("list_answer", loaded).available
    assert ND.get("list_drop", loaded).available      # the drop target never needs a list


def test_the_pulldown_node_names_the_background_it_is_waiting_for(catalogue):
    node = ND.get("list_pulldown", catalogue)
    assert not node.available
    assert "background table" in node.unavailable_reason
    assert "must not be used as a filter" in node.unavailable_reason


def test_nothing_may_be_filled_until_a_frozen_bound_clears_the_bar(catalogue):
    """`Fill` is unavailable precisely when no question has cleared the frozen-test bound, and the
    reason says it is that bound and not the cross-validated one."""
    node = ND.get("label_fill", catalogue)
    fillable = [n for n in ND.children("label_catalogue", catalogue)
                if n.scoreboard is not None and n.scoreboard.may_fill]
    assert node.available == bool(fillable)
    if not fillable:
        assert "frozen-test" in node.unavailable_reason


def test_the_catalogue_header_carries_effective_independent_axes():
    header = ND.catalogue_header()
    axes = ND.effective_independent_axes()
    assert "Effective independent axes" in header
    assert str(axes.effective) in header
    assert "upper bound" in header or axes.measured
    assert "not a negative" in header
    assert 0 < axes.answerable <= axes.questions


def test_search_filters_and_ranks_nodes(catalogue):
    assert [n.title for n in ND.rank("pulldown_background", catalogue)] == ["Pulldown"]
    assert ND.rank("glucose", catalogue)[0].question_id == "required_without_glucose"
    assert [n.title for n in ND.rank("partners", catalogue)][0] == "Partners"
    assert len(ND.rank("", catalogue)) == len(catalogue)
    assert ND.rank("no such thing anywhere", catalogue) == []


def test_nodes_extend_the_help_index_rather_than_duplicating_it(catalogue):
    from parascape import help_index as HI
    rows = ND.search_rows(catalogue)
    assert len(rows) == len(catalogue)
    assert {r.kind for r in rows} == {"node"}
    assert "node" in HI.KIND_ORDER and HI.KIND_LABEL["node"]
    assert {r.data["node"] for r in rows} == {n.id for n in catalogue}
    assert any(e.kind == "node" for e in HI.search(HI.static_entries(), "pulldown_background"))


def test_a_scoreboard_with_no_run_behind_it_says_it_is_unscored():
    board = ND.Scoreboard(question_id="x")
    assert not board.scored
    assert "not yet scored" in board.inline()
    assert "not a negative" in board.coverage()


def test_scoreboards_are_read_from_the_harness_results(tmp_path):
    import json
    run = tmp_path / "results" / "run_a"
    run.mkdir(parents=True)
    (run / "results.json").write_text(json.dumps({"q": {
        "question_id": "q", "n_labels": 1000, "bar": 0.9,
        "cross_validated": {"precision_lower90": 0.8, "coverage_at_bar": 0.1},
        "frozen_test": {"precision_lower90": 0.95, "reaches_bar": True, "coverage_at_bar": 0.2},
        "flags": []}}))
    boards = ND.scoreboards(tmp_path)
    assert boards["q"].may_fill and boards["q"].genes_at_bar == 200
    assert "200 genes answerable" in boards["q"].inline()
    assert ND.scoreboards(tmp_path / "nowhere") == {}


def test_the_node_model_needs_no_qt():
    import ast
    import os
    source = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "parascape", "nodes.py")
    tree = ast.parse(open(source).read())
    imported = [n.module or "" for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)]
    imported += [a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names]
    assert not [m for m in imported if "Qt" in m or "qt" in m.split(".")]


@pytest.mark.parametrize("node_id", ["gene_similar", "list_profile", "list_mechanism", "list_rank"])
def test_a_node_whose_view_does_not_exist_yet_says_so_rather_than_being_left_out(catalogue,
                                                                                 node_id):
    """Phase 1 is a navigation layer. A card for work the application cannot do yet is shown greyed
    with the reason, because a card left off the screen is a gap nobody notices."""
    node = ND.get(node_id, catalogue)
    assert not node.available
    assert "needs" in node.unavailable_reason
    assert ND.resolve(node).kind == "unavailable"


def test_the_inference_nodes_are_present_because_inference_is_what_this_application_does(catalogue):
    """The first draft of this catalogue browsed the map and omitted the strategies, the embedding
    and label inference, which is most of the program."""
    ids = {n.id for n in catalogue}
    assert {"gene_predicted", "gene_similar", "list_profile", "list_mechanism", "list_rank",
            "list_map", "label_infer", "label_strategies", "label_scorecard",
            "label_map"} <= ids


def test_the_strategy_card_promises_every_strategy_grouped_by_mechanism(catalogue):
    """Einar's first requirement for Label: lose nothing. The card says every way to infer, so it
    cannot go stale as the catalogue grows, and the strategies are reachable from it."""
    node = ND.get("label_strategies", catalogue)
    assert "every way to infer" in node.summary
    assert node.destination.target == "strategies"
    assert "how it was tested" in node.summary


def test_runs_is_a_node_even_though_no_panel_reads_the_ledger_yet(catalogue):
    node = ND.get("label_runs", catalogue)
    assert not node.available
    assert "harness_holdout_log" in node.unavailable_reason
    assert "re-runnable" in node.summary


def test_provenance_is_baked_into_the_gene_card_rather_than_given_a_category(catalogue):
    assert "sources" not in ND.CATEGORY_IDS
    card = ND.get("gene_card", catalogue)
    assert "datasets" in card.summary and "papers" in card.summary


def test_the_map_is_one_lens_under_all_three_categories(catalogue):
    """`Clusters` and `Project` undersold the embedding: it is one view entered with three contexts,
    and the embedding strategies stay under Label rather than becoming a browse category."""
    lens = [ND.get(node_id, catalogue) for node_id in ND.MAP_LENS]
    assert [n.title for n in lens] == ["Map", "Map", "Map"]
    assert [n.category for n in lens] == ["gene", "list", "label"]
    assert {n.destination.target for n in lens} == {"maps"}
    assert len({n.destination.section for n in lens}) == 3, "each lens needs its own context"
    assert not [n for n in catalogue if n.title in ("Project", "Clusters")]


def test_a_node_that_produces_an_inference_is_flagged_for_the_user_s_own_layer(catalogue):
    """The hook for Keep: an inference is written to the user's own layer, with its probability,
    uncertainty, strategy and date -- never to the shipped map, and never into Facts."""
    flagged = {n.id for n in catalogue if n.produces_inference}
    assert flagged == {"gene_predicted", "list_answer", "label_infer", "label_fill"}
    for node in catalogue:
        if node.produces_inference:
            assert node.shows_claim, f"{node.id} shows an inference without its uncertainty"
    assert not ND.get("gene_facts", catalogue).produces_inference


def test_facts_may_never_be_turned_into_an_inference_surface():
    base = dict(id="gene_facts", title="Facts", kind="gene", category="gene",
                summary="every trait recorded", destination=ND.Destination("dock", "information map"),
                shows_claim=True, claim_fields=ND.CLAIM_FIELDS, produces_inference=True)
    with pytest.raises(ValueError):
        ND.Node(**base)


def test_an_inference_without_its_uncertainty_is_refused():
    with pytest.raises(ValueError):
        ND.Node("x", "X", "gene", "gene", "an inference with no uncertainty",
                ND.Destination("dock", "discoveries"), produces_inference=True)


# --------------------------------------------------------------------------- the navigator
GENES = [("TGME49_208830", "GRA16 dense granule protein"),
         ("TGME49_233460", "MYR1 translocon component")]


def test_one_search_returns_genes_questions_and_index_rows_interleaved(catalogue):
    """The navigator Einar asked for: one bar for everything, every result a node."""
    from parascape import help_index as HI
    index = HI.static_entries()
    genes = ND.search("GRA16", catalogue, genes=GENES, index=index)
    assert genes and genes[0].kind == "gene"
    assert genes[0].destination.kind == "gene"
    assert genes[0].destination.target == "TGME49_208830"
    assert "dense granule" in genes[0].summary
    strategies = ND.search("hdbscan", catalogue, genes=GENES, index=index)
    assert strategies and all(n.destination.kind == "entry" for n in strategies[:2])
    assert len({n.id for n in strategies}) == len(strategies), "a result must not appear twice"
    mixed = ND.search("glucose", catalogue, genes=GENES, index=index)
    assert any(n.kind == "question" for n in mixed), "the question must still be reachable"
    assert any(n.from_search for n in mixed), "and the index rows beside it"


def test_a_product_name_finds_the_gene_not_only_its_accession(catalogue):
    found = ND.search("dense granule", catalogue, genes=GENES)
    assert [n.destination.target for n in found] == ["TGME49_208830"]


def test_a_searched_node_may_carry_the_name_of_what_it_found(catalogue):
    """The one-or-two-word rule is about what the CATALOGUE calls its nodes; a result is named by
    the thing it found."""
    node = ND.gene_node("TGME49_208830", "GRA16 dense granule protein")
    assert node.from_search and node.title == "TGME49_208830"
    for card in catalogue:
        assert not card.from_search, "a catalogue card is never a search result"
        assert len(card.title.split()) <= ND.MAX_TITLE_WORDS


def test_an_empty_or_unmatched_query_is_honest(catalogue):
    assert ND.search("", catalogue, genes=GENES) == list(catalogue)
    assert ND.search("zzzznothing", catalogue, genes=GENES) == []
    assert ND.genes_matching("x", None) == []
    assert ND.genes_matching("", GENES) == []


def test_a_broken_gene_lookup_does_not_break_the_search(catalogue):
    def explode(_query):
        raise RuntimeError("the table went away")
    assert ND.search("gra16", catalogue, genes=explode) == []


def test_the_shared_layout_scales_come_from_the_theme():
    """spaCR's spacing and radius scales, so the two programs read as one family."""
    from parascape import theme as TH
    assert TH.SPACING == {"xs": 4, "sm": 8, "md": 12, "lg": 16, "xl": 24, "xxl": 32}
    assert TH.RADIUS["tile"] == 16 and TH.RADIUS["card"] == 14
