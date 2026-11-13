"""The Navigator: one field over the animated field, answering with anything.

Each test below pins something that went wrong while building it, found by
rendering the page rather than reading the code.
"""
import pytest
from PySide6 import QtCore

from parascape import modules as MOD
from parascape import navigator as NV

GENES = [("TGME49_208830", "dense granule protein GRA16"),
         ("TGME49_223920", "rhoptry neck protein RON3"),
         ("TGME49_214320", "facilitative glucose transporter GT1")]



@pytest.fixture(autouse=True)
def _animation_on_by_default(monkeypatch):
    """These tests are about the animation, so they keep its shipped default even where CI turns it off."""
    monkeypatch.delenv("PARASCAPE_AMBIENT_DEFAULT", raising=False)

@pytest.fixture
def index():
    from parascape import help_index as HI
    return HI.static_entries()


# ------------------------------------------------------------------ the ranking


def test_a_module_named_outright_comes_first(index):
    assert NV.suggestions("map", index=index)[0] == NV.suggestions("map", index=index)[0]
    first = NV.suggestions("map", index=index)[0]
    assert (first.kind, first.key) == ("module", "map")


def test_a_gene_is_found_by_its_product(index):
    hits = NV.suggestions("GRA16", genes=GENES, index=index)
    assert ("gene", "TGME49_208830") in {(s.kind, s.key) for s in hits}


def test_a_gene_is_found_by_its_accession(index):
    hits = NV.suggestions("TGME49_223920", genes=GENES, index=index)
    assert any(s.kind == "gene" and s.key == "TGME49_223920" for s in hits)


def test_typing_a_gene_offers_the_gene_card_first():
    """Einar, 2026-11-12: "if i type a gene into the navigator i get a few moduals
    that pop up, one should be the gene card"."""
    hits = NV.suggestions("TGME49_223920", genes=GENES, index=None)
    assert hits[0].kind == "module" and hits[0].key == "gene"
    assert hits[0].subject == ("gene", "TGME49_223920")
    assert "TGME49_223920" in hits[0].detail


def test_the_modules_offered_for_a_gene_all_read_one_gene():
    from parascape import ports as P
    hits = NV.suggestions("TGME49_223920", genes=GENES, index=None)
    offered = [s.key for s in hits if s.kind == "module"]
    assert offered, "a typed gene offers no module at all"
    for key in offered:
        assert key in P.consumers_of(P.GENE), key
        assert dict(s.subject for s in hits if s.key == key)


def test_naming_a_module_still_wins_over_the_gene_modules():
    """Typing "map" means the Map module, even when a gene product mentions a map."""
    hits = NV.suggestions("map", genes=GENES, index=None)
    assert hits[0].kind == "module" and hits[0].key == "map"


def test_questions_are_still_searchable(index):
    """Dropping the old Gene/List/Label cards first dropped the question cards with
    them, and "glucose" stopped finding the glucose question."""
    kinds = [s.kind for s in NV.suggestions("glucose", index=index)]
    assert "question" in kinds


def test_methods_and_datasets_come_from_the_same_ranking(index):
    kinds = {s.key if s.kind == "group" else s.kind
             for s in NV.suggestions("conformal", index=index)}
    assert "strategy" in kinds


def test_map_gives_a_handful_of_cards_not_two_dozen(index):
    """Einar: "there are way too many options when I type map"."""
    cards = NV.suggestions("map", genes=GENES, index=index)
    assert len(cards) <= 8, [c.title for c in cards]
    assert not any(c.kind == "panel" and c.title.lower() == "map" for c in cards), \
        "a panel repeating the module's name"
    groups = [c for c in cards if c.kind == "group"]
    assert groups and all(len(g.items) > 1 for g in groups)
    # Nothing but modules that only mention maps, and the panel named like the module,
    # is dropped: every other result is a card or inside one.
    everything = [s for s in NV.all_matches("map", genes=GENES, index=index)
                  if s.kind != "module" and not (s.kind == "panel" and s.title.lower() == "map")]
    inside = sum(len(g.items) for g in groups) + sum(1 for c in cards
                                                     if c.kind not in ("group", "module"))
    assert inside == len(everything), "grouping lost results"


def test_a_short_answer_is_not_grouped(index):
    cards = NV.suggestions("TGME49_223920", genes=GENES, index=index)
    assert all(c.kind != "group" for c in cards)


def test_the_old_node_cards_do_not_come_back(index):
    """The modules replaced Gene/List/Label as the front door."""
    for s in NV.suggestions("map", index=index):
        assert not (s.kind == "gene" and s.title == "Map"), "the old Gene ▸ Map card resurfaced"


def test_an_empty_query_suggests_nothing():
    assert NV.suggestions("") == [] and NV.suggestions("   ") == []


def test_a_broken_index_does_not_break_the_field():
    def boom():
        raise RuntimeError("index unavailable")
    assert [s.kind for s in NV.suggestions("map", index=boom)][:1] == ["module"]


# ---------------------------------------------------------------- card text


def test_a_dataset_subtitle_drops_its_breadcrumb():
    """On a two-line card the slot-tree path used both lines and said nothing."""
    raw = "Slot tree ▸ Tg ▸ relational measurements ▸ host effect — host effect human K562; invasion"
    assert NV._subtitle(raw).startswith("host effect human K562")


def test_a_plain_subtitle_is_left_alone():
    assert NV._subtitle("facilitative glucose transporter GT1") == "facilitative glucose transporter GT1"


def test_clip_shortens_at_a_word_with_an_ellipsis():
    out = NV._clip("one two three four five six seven", 14)
    assert out.endswith("…") and len(out) <= 15 and " " not in out[-2:]


# --------------------------------------------------------------------- widget


@pytest.fixture
def page(qtbot, index):
    widget = NV.NavigatorPage(genes=GENES, index=index)
    qtbot.addWidget(widget)
    widget.resize(1200, 800)
    widget.show()
    return widget


def type_and_wait(page, qtbot, text):
    page.field.setText(text)
    page.refresh()


def test_typing_alone_shows_nothing_until_return(page, qtbot):
    """Einar: "when I type things in the navigator they instantly appear, they should wait
    until I hit enter"."""
    page.field.setText("rhoptry")
    page._edited("rhoptry")
    qtbot.wait(300)
    assert page._cards == [], "cards appeared before Return"
    assert "Return" in page.feedback.text()


def test_return_searches_and_a_second_return_opens(page, qtbot):
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    page.field.setText("map")
    page.submit()
    assert page._cards and seen == [], "the first Return searches, it does not open"
    page.submit()
    assert seen == [] and page._path == ["map"], "the second Return shows the module network"
    QtCore.QCoreApplication.processEvents()
    from PySide6.QtTest import QTest
    QTest.keyClick(page, QtCore.Qt.Key.Key_Return)
    assert seen == [("module", "map")], "Return in the network opens the highlighted module"


def test_one_click_shows_what_follows_and_a_second_opens(page, qtbot):
    """Einar: "the network ... of what is downstream of map. And a second click to open
    the module." spaCR: "Click once to open the module network. Click twice to open it"."""
    from parascape import modules as MOD
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    type_and_wait(page, qtbot, "map")
    page._cards[0].click()
    assert seen == [] and page._path == ["map"]
    assert page._columns[1] == list(MOD.downstream("map"))[:NV.FLOW_FIRST]
    assert not page.field.isVisible() and page.back.isVisible()
    page._flow_cards[(1, "network")].click()
    assert page._path == ["map", "network"] and seen == []
    assert page._columns[2] == list(MOD.downstream("network"))[:NV.FLOW_MORE]
    page._flow_cards[(1, "network")].click()
    assert seen == [("module", "network")]


def test_back_leaves_the_network_for_the_results(page, qtbot):
    type_and_wait(page, qtbot, "map")
    page._cards[0].click()
    page.back.click()
    assert page._path == [] and page.field.isVisible() and page._cards


def test_a_group_opens_to_its_members_and_back_returns(page, qtbot):
    type_and_wait(page, qtbot, "map")
    group = next(c for c in page._cards if c.suggestion.kind == "group")
    n = len(group.suggestion.items)
    group.click()
    assert len(page._cards) == n and page.back.isVisible()
    page.back.click()
    assert any(c.suggestion.kind == "group" for c in page._cards)


def test_reset_empties_everything(page, qtbot):
    """Einar: "navigator should be reset every time I go back to it"."""
    type_and_wait(page, qtbot, "map")
    page._cards[0].click()
    page.reset()
    assert page.field.text() == "" and page._cards == [] and page._path == []
    assert page.field.isVisible() and not page.back.isVisible()


def test_every_downstream_module_exists():
    """What follows a module is computed from the declared subjects now, so this is
    really a check that the two tables have not drifted: every key the graph answers
    with is a module the navigator can draw a card for."""
    from parascape import modules as MOD
    keys = {m.key for m in MOD.MODULES}
    for module in MOD.MODULES:
        following = MOD.downstream(module.key)
        assert set(following) <= keys, (module.key, set(following) - keys)
        assert module.key not in following, f"{module.key} follows itself"
        assert len(set(following)) == len(following), module.key


def test_editing_after_a_search_asks_for_return_again(page, qtbot):
    type_and_wait(page, qtbot, "rhoptry")
    page.field.setText("rhoptry neck")
    page._edited("rhoptry neck")
    assert page.feedback.text() == "Press Return to search."


def test_cards_are_spacrs_squares_with_an_icon_and_a_light_centred_name(page, qtbot):
    """Einar: the cards must look exactly like spaCR's. spaCR's are 224 px squares (smaller
    only when the window cannot fit them), an icon over an 18 px light name, nothing else."""
    from PySide6 import QtWidgets
    page.resize(1600, 900)
    type_and_wait(page, qtbot, "rhoptry")
    card = page._cards[0]
    final = card._final
    assert final.width() == final.height() == NV.CARD_SIDE
    assert card.objectName() == "EasySuggestionCard"
    labels = card.findChildren(QtWidgets.QLabel)
    assert any(l.pixmap() is not None and not l.pixmap().isNull() for l in labels), "no icon"
    assert "18px" in card.name_label.styleSheet() and "300" in card.name_label.styleSheet()
    assert len([l for l in labels if l.text()]) == 1, "spaCR cards carry the name and nothing else"


def test_the_cards_sit_above_the_field_and_the_details_below_it(page, qtbot):
    """spaCR's arrangement: the row of cards over the centred field, the description at the
    foot of the page."""
    page.resize(1600, 900)
    type_and_wait(page, qtbot, "rhoptry")
    field = page.field.geometry()
    assert page._cards[0]._final.bottom() < field.top()
    assert page.details.geometry().top() > field.bottom()


def test_a_gene_card_names_the_product_not_the_accession(page, qtbot):
    """Eight cards reading TGME49_... could not be told apart without hovering each.

    The gene cards come after the modules that read a gene now (instruction 98, WP8),
    so this looks for the gene card rather than taking the first one."""
    type_and_wait(page, qtbot, "TGME49_223920")
    card = next(c for c in page._cards if c.suggestion.kind == "gene")
    assert "RON3" in card.name_label.text()
    assert "TGME49_223920" in card.suggestion.detail


def test_hovering_a_card_describes_it_below(page, qtbot):
    type_and_wait(page, qtbot, "rhoptry")
    second = page._cards[1]
    page.explain(second)
    assert second.suggestion.title.split()[0] in page.details.toPlainText()


def test_every_card_icon_exists(index):
    """A missing PNG would draw an empty square, silently."""
    import os
    from parascape import modules as MOD
    for m in MOD.MODULES:
        assert os.path.isfile(NV.icon_path(NV.Suggestion("module", m.key, m.name))), m.key
    for kind in ("gene", "dataset", "question", "strategy", "guide", "panel"):
        assert os.path.isfile(NV.icon_path(NV.Suggestion(kind, "x", "x"))), kind


def test_return_opens_the_first_suggestion(page, qtbot):
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    type_and_wait(page, qtbot, "map")
    page.open_current()
    assert seen == [("module", "map")]


def test_return_before_the_debounce_still_opens_what_was_typed(page, qtbot):
    """Typing then pressing Return immediately must not open the PREVIOUS row."""
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    page.field.setText("map")                  # debounce running, row not built yet
    page.open_current()
    assert seen == [("module", "map")]


def test_an_unavailable_module_explains_rather_than_opens(page, qtbot):
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    type_and_wait(page, qtbot, "fill")
    assert page.current().key == "fill"
    page.open_current()
    assert seen == [], "Fill has nothing to write; opening it would show an empty view"
    assert "not available" in page.feedback.text()


def test_nothing_found_says_so(page, qtbot):
    type_and_wait(page, qtbot, "zzqxnotathing")
    assert not page._cards and "Nothing here matches" in page.feedback.text()


def test_the_ai_symbol_sits_inside_the_field_on_the_right(page, qtbot):
    """A child of a QLineEdit is drawn at (0, 0) unless something positions it; the
    first render put the assistant above the text at the top left."""
    type_and_wait(page, qtbot, "rhoptry")
    page._place_ai()
    geo = page.ai_button.geometry()
    assert geo.x() > page.field.width() // 2
    assert 0 <= geo.y() and geo.bottom() <= page.field.height()
    assert page.field.textMargins().right() >= geo.width(), "typed text stops before the symbol"


def test_the_ai_symbol_is_a_switch_that_is_off_until_clicked(page):
    """As in spaCR: a ring whose dot is drawn only while the assistant is on."""
    assert page.ai_button.isCheckable() and not page.ai_button.isChecked()
    assert "off" in page.ai_button.toolTip()


def test_a_search_with_no_answer_offers_the_assistant_instead_of_asking_it(page, qtbot):
    """With the assistant off nothing is sent anywhere, and the page says how to turn it on."""
    type_and_wait(page, qtbot, "zzqxnotathing")
    assert page._ai_worker is None
    assert "ring" in page.feedback.text()


def test_with_the_assistant_on_a_search_with_no_answer_is_answered_below(page, qtbot, monkeypatch):
    """The answer is streamed into the details panel, the ring shows it is busy, then settles."""
    from parascape import chat
    monkeypatch.setattr(chat, "stream", lambda *a, **k: iter(["The ", "answer ", "is 42."]))
    page._ai_provider = "parascape"
    page.ai_button.setChecked(True)
    page.field.setText("zzqxnotathing")
    page.refresh()
    assert page.ai_button.busy
    qtbot.waitUntil(lambda: not page.ai_button.busy, timeout=5000)
    assert "The answer is 42." in page.details.toPlainText()


def test_turning_the_assistant_off_cancels_an_answer_in_progress(page, qtbot, monkeypatch):
    import threading
    from parascape import chat
    release = threading.Event()

    def slow(*a, **k):
        release.wait(5)
        yield "late"
    monkeypatch.setattr(chat, "stream", slow)
    page._ai_provider = "parascape"
    page.ai_button.setChecked(True)
    page.ask_ai("anything")
    assert page.ai_button.busy
    page.ai_button.setChecked(False)
    page._toggle_ai(False)
    release.set()
    assert not page.ai_button.busy and page._ai_worker is None
    qtbot.waitUntil(lambda: not page._ai_retired, timeout=5000)


def test_the_first_click_with_no_provider_asks_which_one(page, monkeypatch):
    """Not guessing a provider: choosing is also when a person agrees to what is shared."""
    asked = []
    monkeypatch.setattr(page, "choose_provider", lambda *a, **k: asked.append(True) or "")
    page._ai_provider = ""
    page.ai_button.setChecked(True)
    page._toggle_ai(True)
    assert asked and not page.ai_button.isChecked()


def test_the_provider_sets_the_colour_of_the_dot(page):
    from parascape.ai_button import PROVIDER_COLORS
    for name, colour in PROVIDER_COLORS.items():
        page.ai_button.set_provider(name)
        assert page.ai_button.color().name().upper() == colour.upper()


def test_typing_pulses_the_dot_without_sending_anything(page, qtbot):
    page.ai_button.setChecked(True)
    page._ai_provider = "parascape"
    page.field.textEdited.emit("g")
    assert page.ai_button._pulse > 0 and page._ai_worker is None


def test_escape_clears_the_field(page, qtbot):
    page.field.setText("rhoptry")
    qtbot.keyClick(page.field, QtCore.Qt.Key.Key_Escape)
    assert page.field.text() == ""


# --------------------------------------------------------------------- landing


@pytest.fixture
def landing(qtbot, index, tmp_path, monkeypatch):
    from parascape.backdrops import _compat as AC
    path = str(tmp_path / "ambient.ini")
    monkeypatch.setattr(AC, "_settings",
                        lambda: QtCore.QSettings(path, QtCore.QSettings.Format.IniFormat))
    from parascape.landing import Landing
    widget = Landing(genes=GENES, index=index)
    qtbot.addWidget(widget)
    return widget


def test_the_landing_opens_on_the_navigator(landing):
    assert landing.page == "Navigator"


def test_the_switch_reaches_the_modules(landing):
    landing.show_page("Modules")
    assert landing.page == "Modules" and landing.stack.currentWidget() is landing.home


def test_one_backdrop_serves_both_pages(landing):
    """Installing one per page would restart the animation on every switch and run
    two timers for one picture."""
    assert landing.animated
    assert landing.home.animated is False, "the tiles page must not install its own"
    before = landing._ambient
    landing.show_page("Modules")
    landing.show_page("Navigator")
    assert landing._ambient is before


def test_a_module_from_either_page_arrives_the_same_way(landing):
    seen = []
    landing.module_chosen.connect(seen.append)
    landing.navigator.field.setText("map")
    landing.navigator.open_current()
    landing.home.open_module("gene")
    assert seen == ["map", "gene"]


def test_a_gene_from_the_navigator_opens_the_gene_card_holding_that_gene(landing):
    """Einar, 2026-11-12: "now when i put in a gene and pressed the modual that came up
    for the gene, i was brought to the old version of parascape." It opened the gene as a
    NODE, which raised the pre-navigator evidence dock. Return on a typed gene now opens
    the gene card module, and the gene travels with it."""
    seen = []
    landing.module_chosen.connect(lambda key: seen.append(key))
    landing.navigator.field.setText("TGME49_223920")
    landing.navigator.open_current()
    assert seen == ["gene"]
    assert landing.navigator.subject == ("gene", "TGME49_223920")


def test_a_gene_card_among_the_results_still_opens_as_an_item(landing):
    """The gene results themselves are still items; `open_item` is what routes them to
    the gene card now, so the kind on the signal has not changed."""
    seen = []
    landing.item_chosen.connect(lambda kind, key: seen.append((kind, key)))
    landing.navigator.field.setText("TGME49_223920")
    landing.navigator.refresh()
    card = next(c for c in landing.navigator._cards if c.suggestion.kind == "gene")
    landing.navigator._open(landing.navigator._cards.index(card))
    assert seen == [("gene", "TGME49_223920")]


def test_the_theme_can_be_changed_live(landing):
    from parascape.backdrops import engines as E
    landing.set_theme("blobs")
    assert landing.animated
    landing.set_theme(E.NO_ANIMATION)
    assert not landing.animated
    # And the page paints its own colour again. Left transparent with no backdrop
    # behind it, it would show through to whatever is underneath.
    assert landing.property("ambient") == "off"
    assert landing.home.property("ambient") == "off"
    assert not landing.stack.testAttribute(QtCore.Qt.WidgetAttribute.WA_TranslucentBackground)


# ------------------------------------ instruction 98, WP8: clicking the gene card
def test_clicking_the_gene_card_once_reveals_the_modules_downstream_of_one_gene(page, qtbot):
    """Einar, 2026-11-12: "when i click the gene card once i should see the moduals
    downstream of one gene (like the network in spacr). so all the modules that take 1
    gene and do something should come downstream of the gene card module"."""
    from parascape import modules as MOD
    from parascape import ports as P
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    type_and_wait(page, qtbot, "TGME49_223920")
    card = next(c for c in page._cards
                if c.suggestion.kind == "module" and c.suggestion.key == "gene")
    card.click()
    assert seen == [] and page._path == ["gene"], "the first click opens the network"
    assert set(page._columns[1]) == set(P.consumers_of(P.GENE)) - {"gene"}
    assert page._columns[1] == list(MOD.downstream_of_subject(P.GENE, without="gene"))


def test_every_module_revealed_carries_the_gene_and_says_why_it_cannot(page, qtbot):
    """Instruction 98: offered pre-filled, and an offer that cannot work shown greyed
    WITH THE REASON."""
    type_and_wait(page, qtbot, "TGME49_223920")
    next(c for c in page._cards
         if c.suggestion.kind == "module" and c.suggestion.key == "gene").click()
    assert page.subject == ("gene", "TGME49_223920")
    greyed = 0
    for (column, key), card in page._flow_cards.items():
        if column != 1:
            continue
        s = card.suggestion
        assert s.subject == ("gene", "TGME49_223920"), key
        if s.blocked:
            greyed += 1
            assert s.blocked.lower() in s.detail.lower(), key
            assert "Not ready" in s.detail, key
        else:
            assert "TGME49_223920" in s.detail, key
    assert greyed, "Structure has no view built, so one offer must be greyed with its reason"


def test_a_greyed_offer_refuses_to_open_and_says_what_is_missing(page, qtbot):
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    type_and_wait(page, qtbot, "TGME49_223920")
    next(c for c in page._cards
         if c.suggestion.kind == "module" and c.suggestion.key == "gene").click()
    blocked = next(k for (c, k) in page._flow_cards if c == 1
                   and page._flow_cards[(c, k)].suggestion.blocked)
    page._flow_cards[(1, blocked)].click()         # make it the step shown
    page._flow_cards[(1, blocked)].click()         # and try to open it
    assert seen == [], "a greyed offer must not open"
    assert "not available" in page.flow_hint.text() or "cannot" in page.flow_hint.text()


def test_the_subject_is_forgotten_when_the_navigator_is_reset(page, qtbot):
    type_and_wait(page, qtbot, "TGME49_223920")
    assert page.subject[0] == "gene"
    page.reset()
    assert page.subject == ("", "")


# ------------------------------- instruction 98, WP9: a typed LIST of genes
LIST_TEXT = "TGME49_208830, TGME49_223920, TGME49_214320"


def test_several_identifiers_that_name_genes_are_read_as_a_list():
    """What makes something a list is that it names SEVERAL genes. One is a search, and
    three words that name none stay a search too."""
    assert NV.genes_in_text(LIST_TEXT, GENES) == [g for g, _ in GENES]
    assert NV.genes_in_text("TGME49_208830", GENES) == []
    assert NV.genes_in_text("rhoptry neck protein", GENES) == []


def test_a_pasted_list_is_recognised_however_it_was_typed():
    for text in ("TGME49_208830 TGME49_223920", "TGME49_208830\nTGME49_223920",
                 "tgme49_208830;TGME49_223920"):
        assert len(NV.genes_in_text(text, GENES)) == 2, text
    # Spelled as the table spells them, never as they were typed.
    assert NV.genes_in_text("tgme49_208830, tgme49_223920", GENES) == \
        ["TGME49_208830", "TGME49_223920"]


def test_typing_a_list_offers_the_gene_list_card_first():
    """Einar, 2026-11-12: "if i type in a list of genes, the same logic but this time not
    gene card instead a gene list module"."""
    hits = NV.suggestions(LIST_TEXT, genes=GENES, index=None)
    assert hits[0].kind == "module" and hits[0].key == "list"
    assert hits[0].subject == ("gene_list", ",".join(g for g, _ in GENES))
    assert "3 genes" in hits[0].detail


def test_the_modules_offered_for_a_list_all_read_a_list():
    from parascape import ports as P
    hits = NV.suggestions(LIST_TEXT, genes=GENES, index=None)
    offered = [s.key for s in hits if s.kind == "module"]
    assert offered and offered[0] == "list"
    for key in offered:
        assert key in P.consumers_of(P.GENE_LIST), key
        card = next(s for s in hits if s.kind == "module" and s.key == key)
        assert card.subject[0] == P.GENE_LIST and card.subject[1], key


def test_the_genes_of_a_typed_list_are_one_group_card_not_three_hundred():
    hits = NV.suggestions(LIST_TEXT, genes=GENES, index=None)
    group = next(s for s in hits if s.kind == "group")
    assert group.key == "gene" and len(group.items) == 3
    assert all(m.kind == "gene" for m in group.items)


def test_a_list_is_said_as_a_count_never_as_its_accessions():
    """A gene is its own name; a list of forty is "40 genes"."""
    from parascape import ports as P
    assert NV.subject_words(P.GENE, "TGME49_208830") == "TGME49_208830"
    assert NV.subject_words(P.GENE_LIST, LIST_TEXT) == "3 genes"
    assert NV.subject_words(P.GENE_LIST, "") == ""


def test_clicking_the_gene_list_card_reveals_the_modules_downstream_of_a_list(page, qtbot):
    """Einar: "after clicking it all the moduals that come downstream of gene lists"."""
    from parascape import modules as MOD
    from parascape import ports as P
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    type_and_wait(page, qtbot, LIST_TEXT)
    card = next(c for c in page._cards
                if c.suggestion.kind == "module" and c.suggestion.key == "list")
    card.click()
    assert seen == [] and page._path == ["list"], "the first click opens the network"
    assert set(page._columns[1]) == set(P.consumers_of(P.GENE_LIST)) - {"list"}
    assert page._columns[1] == list(MOD.downstream_of_subject(P.GENE_LIST, without="list"))
    assert NV.FLOW_FIRST is None, "all of them, so a number here cannot go stale again"
    assert len(page._flow_cards) == 1 + len(page._columns[1]), "every one of them is drawn"


def test_every_module_revealed_after_a_list_carries_it_and_says_why_it_cannot(page, qtbot):
    """Offered pre-filled, and an offer that cannot work shown greyed WITH THE REASON."""
    type_and_wait(page, qtbot, LIST_TEXT)
    next(c for c in page._cards
         if c.suggestion.kind == "module" and c.suggestion.key == "list").click()
    assert page.subject == ("gene_list", ",".join(g for g, _ in GENES))
    greyed = 0
    for (column, key), card in page._flow_cards.items():
        if column != 1:
            continue
        s = card.suggestion
        assert s.subject[0] == "gene_list", key
        if s.blocked:
            greyed += 1
            assert s.blocked.lower() in s.detail.lower(), key
            assert "Not ready" in s.detail, key
        else:
            # The count, never the accessions.
            assert "3 genes" in s.detail, key
            assert "TGME49_208830," not in s.detail, key
    assert greyed, "Infer needs a question, so one offer must be greyed with its reason"


def test_the_first_column_follows_the_subject_in_hand_not_a_slice_of_everything(page, qtbot):
    """`modules.downstream("list")` is the readers of a list, of a trait and of one gene in
    one tuple; the first eight of that is not "what follows a list"."""
    from parascape import ports as P
    from parascape import modules as MOD
    type_and_wait(page, qtbot, LIST_TEXT)
    assert page.following("list", first=True) == \
        MOD.downstream_of_subject(P.GENE_LIST, without="list")
    assert page.following("list", first=False) == MOD.downstream("list")


def test_a_card_name_too_long_for_a_small_card_is_elided_not_clipped(page, qtbot):
    """Eight cards in the height of one is an 84 px square, and "Expression" was drawn as
    "xpression" -- clipped at both ends, so neither end said which module it was."""
    type_and_wait(page, qtbot, LIST_TEXT)
    card = next(c for c in page._cards
                if c.suggestion.kind == "module" and c.suggestion.key == "list")
    card.fit(84)
    assert card.name_label.text() != "Gene list" or len("Gene list") < 10
    long_one = next(c for c in page._cards if c.suggestion.kind == "module")
    long_one.fit(84)
    from PySide6 import QtGui
    font = QtGui.QFont(long_one.name_label.font())
    font.setPixelSize(max(11, min(18, 84 // 9)))
    assert QtGui.QFontMetrics(font).horizontalAdvance(long_one.name_label.text()) <= 84 - 20
    long_one.fit(224)
    assert long_one.name_label.text() == NV.card_name(long_one.suggestion)


# ----------------------------------------------------- what the search ranks on (WP4)
# Einar: "the navigator should have a search bar like spacr". Typing a gene and a gene
# list already worked; what WP4 added is module names, view names, what a module DOES,
# and the words a parasitologist would use rather than the internal name. The bands are
# declared in `modules.py` beside the reason for the order between them.


def test_a_name_typed_in_full_beats_everything_else():
    first = NV.suggestions("map", index=None)[0]
    assert (first.kind, first.key) == ("module", "map")
    assert MOD.match_score("map", "Map") == MOD.NAME_EXACT


def test_the_word_a_parasitologist_uses_finds_the_module():
    """"interactome" is the name of no tile and appears in no note, so before WP4 it
    found nothing at all. The words are declared on the catalogue row."""
    for typed, wanted in (("interactome", {"partners", "network"}),
                          ("LOPIT", {"cell"}),
                          ("CRISPR", {"screen"}),
                          ("rhoptry", {"cell"}),
                          ("AlphaFold", {"structure"}),
                          ("UMAP", {"map"})):
        keys = {m.key for m in MOD.search(typed)}
        assert wanted <= keys, (typed, keys)


def test_a_declared_word_never_outranks_a_name():
    """Somebody who typed a name meant that module. "network" is a declared word of
    Partners and the name of another module, and the other module wins."""
    assert MOD.search("network")[0].key == "network"


def test_what_a_module_does_is_matched_and_ranked_last():
    """The subjects a module reads come from the ports, so "trait" finds every module
    that reads one without the word being written on any row."""
    from parascape import ports as P
    keys = {m.key for m in MOD.search("trait")}
    assert set(P.consumers_of(P.TRAIT)) <= keys
    scorecard = MOD.get("scorecard")
    assert MOD.match_score("trait", scorecard.name, scorecard.words,
                           MOD.what_it_does(scorecard)) == MOD.WHAT_IT_DOES
    assert MOD.match_score("trait", scorecard.name) is None, \
        "the word is not on the row; it comes from the declared subjects"


def test_a_view_behind_a_module_is_findable_by_its_own_name():
    """WP10 gave every panel behind the navigator a name; WP4 makes the field answer to
    it. Before this a view could be reached only through the help index of a live
    window, and then only as a row called "panel"."""
    cards = NV.suggestions("coexpression", index=None)
    assert cards[0].kind == "view" and cards[0].key == "Coexpression"
    assert cards[0].owner == "expression"
    assert "Expression" in cards[0].detail


def test_every_declared_view_can_be_found_by_typing_its_name():
    for module, view in MOD.views():
        found = [s for s in NV.all_matches(view.name, index=None)
                 if s.kind == "view" and s.key == view.name]
        assert found, view.name
        assert found[0].owner == module.key


def test_a_module_wins_its_own_name_back_from_its_views():
    """"expression" names the module; "coexpression" names a tab of the panel it shows."""
    assert NV.suggestions("expression", index=None)[0].key == "expression"


def test_a_view_that_matched_is_never_silently_dropped(index):
    """Only four cards sit at the top. A view that does not win one of them falls into
    the groups below rather than being lost, which is what "grouping lost results"
    checks for the other kinds."""
    cards = NV.suggestions("expression", genes=GENES, index=index)
    reachable = set()
    for card in cards:
        reachable |= ({(card.kind, card.key)} if card.kind != "group"
                      else {(m.kind, m.key) for m in card.items})
    for s in NV.all_matches("expression", genes=GENES, index=index):
        if s.kind == "view":
            assert ("view", s.key) in reachable, s.key


def test_a_module_card_says_what_it_takes_and_what_follows_it():
    card = next(s for s in NV.all_matches("partners", index=None) if s.kind == "module")
    assert "Takes one gene" in card.detail
    assert "Leads to" in card.detail


def test_the_words_of_a_query_narrow_rather_than_widen():
    both = {m.key for m in MOD.search("gene list")}
    assert "list" in both
    assert both <= {m.key for m in MOD.search("gene")}


def test_a_search_with_matches_shows_cards_and_never_asks_the_assistant(page, qtbot,
                                                                       monkeypatch):
    """Deliberately different from spaCR and recorded as such: only a search with NO
    match goes to the assistant."""
    asked = []
    monkeypatch.setattr(page, "ask_ai", lambda text: asked.append(text))
    page.ai_button.setChecked(True)
    type_and_wait(page, qtbot, "interactome")
    assert page._cards and not asked
    type_and_wait(page, qtbot, "zzzzzznothinghere")
    assert not page._cards and asked == ["zzzzzznothinghere"]


def test_opening_a_view_from_the_navigator_names_its_module(page, qtbot):
    chosen = []
    page.chosen.connect(lambda kind, key: chosen.append((kind, key)))
    type_and_wait(page, qtbot, "coexpression")
    page._open(0)
    assert chosen == [("view", "Coexpression")]
    assert MOD.find_view("Coexpression")[0].key == "expression"
