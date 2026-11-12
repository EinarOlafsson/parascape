"""What comes next with the thing in hand (instruction 98, WP8).

spaCR's rule, which these hold ParaScape to: an offer either works, or is shown greyed
with the reason it cannot and what would fix it. A blocked offer is never dropped --
dropping it leaves somebody with no idea the module exists.
"""
from parascape import modules as MOD
from parascape import next_steps as NS
from parascape import ports as P


def test_every_module_that_reads_a_gene_is_offered():
    """Einar: "all the modules that take 1 gene and do something should come downstream
    of the gene card module"."""
    keys = [o.module for o in NS.offers(P.GENE, "TGME49_208830", without="gene")]
    assert set(keys) == set(P.consumers_of(P.GENE)) - {"gene"}


def test_ready_offers_come_before_blocked_ones():
    offered = NS.offers(P.GENE, "TGME49_208830", without="gene")
    assert [o.ok for o in offered] == sorted((o.ok for o in offered), reverse=True)


def test_a_blocked_offer_keeps_its_place_and_carries_the_reason_and_the_fix():
    offered = NS.offers(P.GENE, "TGME49_208830", without="gene")
    blocked = [o for o in offered if not o.ok]
    assert blocked, "Structure has no view built, so something must be blocked"
    for offer in blocked:
        assert offer.blocked and offer.fix
        assert offer.label.endswith(NS.NOT_READY)
        assert offer.blocked in offer.tip and offer.fix in offer.tip


def test_a_ready_offer_says_what_it_will_open_with():
    offer = next(o for o in NS.offers(P.GENE, "TGME49_208830", without="gene") if o.ok)
    assert "TGME49_208830" in offer.tip and offer.name in offer.tip
    assert not offer.label.endswith(NS.NOT_READY)


def test_a_window_that_cannot_open_a_module_greys_that_offer_with_its_own_reason():
    """A dock that failed to install greys its module with THAT reason, not a general one."""
    offered = NS.offers(P.GENE, "TGME49_208830", without="gene",
                        available={"partners": "the evidence view is not in this window"})
    partners = next(o for o in offered if o.module == "partners")
    assert not partners.ok and "not in this window" in partners.blocked


def test_offering_with_nothing_in_hand_still_lists_what_a_gene_would_let_you_do():
    """An empty gene card must still show what a gene is for, rather than an empty row."""
    offered = NS.offers(P.GENE, "", without="gene")
    assert len(offered) == len(MOD.downstream_of_subject(P.GENE, without="gene"))
    needs_one = next(o for o in offered if o.module == "partners")
    assert not needs_one.ok and "none is chosen" in needs_one.blocked


def test_the_row_draws_a_button_for_every_offer_and_disables_the_blocked_ones(qtbot):
    row = NS.NextSteps(P.GENE)
    qtbot.addWidget(row)
    offered = row.set_subject("TGME49_208830", without="gene")
    assert set(row.buttons) == {o.module for o in offered}
    for offer in offered:
        button = row.buttons[offer.module]
        assert button.isEnabled() is offer.ok, offer.module
        assert button.toolTip(), offer.module
    assert "TGME49_208830" in row.lead.text()


def test_clicking_a_ready_offer_hands_on_the_module_and_the_gene(qtbot):
    row = NS.NextSteps(P.GENE)
    qtbot.addWidget(row)
    row.set_subject("TGME49_208830", without="gene")
    seen = []
    row.chosen.connect(lambda key, kind, value: seen.append((key, kind, value)))
    row.buttons["partners"].click()
    assert seen == [("partners", "gene", "TGME49_208830")]


def test_the_row_says_how_many_can_run_with_nothing_chosen(qtbot):
    row = NS.NextSteps(P.GENE)
    qtbot.addWidget(row)
    row.set_subject("", without="gene")
    assert "can open without one" in row.lead.text()
    assert not row.buttons["partners"].isEnabled()


# ---------------------------------------- instruction 98, WP9: with a LIST in hand
LIST = "TGME49_208830,TGME49_223920,TGME49_214320"


def test_every_module_that_reads_a_list_is_offered():
    """Einar: "after clicking it all the moduals that come downstream of gene lists"."""
    keys = [o.module for o in NS.offers(P.GENE_LIST, LIST, without="list")]
    assert set(keys) == set(P.consumers_of(P.GENE_LIST)) - {"list"}


def test_a_list_is_said_as_a_count_rather_than_as_its_accessions():
    """A gene is its own name and reads well. A list of forty is "40 genes", and printing
    the forty accessions in a tooltip would be a wall of text where a count was wanted."""
    offer = next(o for o in NS.offers(P.GENE_LIST, LIST, without="list", shown="3 genes")
                 if o.ok)
    assert "3 genes" in offer.tip and "TGME49_208830," not in offer.tip
    assert offer.value == LIST, "the value itself is still what is handed over"
    assert offer.words == "3 genes"


def test_the_gene_card_is_ready_for_a_list_because_it_steps_through_one():
    """The port correction of WP9: the Gene card declared one gene as REQUIRED and a list
    as optional, so somebody holding a list was told it could not work -- in the one place
    a list is the thing in hand. Its own port says "a list to step through, one gene at a
    time", and it has a button for taking the gene list whole."""
    offered = {o.module: o for o in NS.offers(P.GENE_LIST, LIST, without="list")}
    for key in ("gene", "partners", "expression"):
        assert offered[key].ok, (key, offered[key].blocked)


def test_a_module_that_needs_something_else_is_greyed_with_its_own_reason():
    """Infer is asked a QUESTION, which a list is not, so it keeps its place and says so."""
    offered = {o.module: o for o in NS.offers(P.GENE_LIST, LIST, without="list")}
    assert not offered["infer"].ok
    assert "question" in offered["infer"].blocked and offered["infer"].fix


def test_with_nothing_in_hand_the_gene_card_still_says_a_gene_or_a_list_is_wanted():
    """The alternative is said out loud: "one gene is needed" hides half of its own fix."""
    readiness = P.check_ready("gene", P.Subject(P.GENE_LIST))
    assert not readiness.ok
    message = readiness.errors[0].message
    assert "one gene" in message and "a list of genes" in message
