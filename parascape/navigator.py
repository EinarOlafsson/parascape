"""The Navigator: one field, centred over the animated field, that takes anything.

PORTED AND MODIFIED FROM spaCR
==============================
The form is spaCR's Navigator (``spacr/qt/widgets/home.py``: ``_EasyHomeField``,
``_EasyPathEdit``, the suggestion cards and the details browser) -- "a single
centered editor over the window's unchanged ambient field".

    Copyright (c) 2025-2026, Einar Birnir Olafsson and the spaCR contributors
    BSD 3-Clause with additional noncommercial terms
    SPDX-License-Identifier: LicenseRef-spaCR-BSD-3-Clause-Additional-Terms
    Full text: LICENSES/spaCR-BSD-3-Clause-Additional-Terms.txt

**Modifications made in porting.** spaCR's Navigator inspects a *path* -- an image
folder, a plate -- and proposes a pipeline for it. ParaScape's input is not a
path: it is a gene, a question, a module, a dataset or a method. So the field is
kept and what answers it is replaced. Suggestions come from
:func:`parascape.modules.search` and :func:`parascape.nodes.search`, which ranks
genes, questions, strategies, datasets, panels and guide pages in one call, so the
best answer wins whatever kind of thing it is. spaCR's path inspection, protocol
editor and workflow builder are not ported; its AI symbol is, drawn and
behaving as in spaCR (`parascape.ai_button`): a ring-and-dot switch in the field, an answer
written below it, a right-click to choose the provider.

Kept from spaCR: the 64 px pill field with its 28 px radius and 16 px type; the
row of suggestion cards beneath it with faded edges; a details pane for the card
under the cursor; keyboard first (Down into the cards, arrows across them, Return
opens, Escape clears); and the rule that it works with no AI provider at all --
every suggestion is local.
"""
from __future__ import annotations

from dataclasses import dataclass
from html import escape

from PySide6 import QtCore, QtGui, QtWidgets

from . import iconset
from . import modules as MOD
from . import theme as TH
from .ai_button import AiButton, ProviderDialog

#: Field geometry, spaCR's: a tall pill, so the one thing on the page reads as the
#: thing to use.
FIELD_H = 64
FIELD_RADIUS = 28
FIELD_MIN_W = 720             # spaCR: min(max(720, placeholder width), width - 48)
CARD_SIDE = 224               # spaCR: square cards, side min(224, room)
CARD_GAP = 20
ICON_PX = 96
GROW_MS = 220                 # the cards grow into place, as spaCR's do
GROW_STAGGER_MS = 35
DETAILS_H = 150
MAX_CARDS = 24
FLOW_SIDE = 170               # a next-step node in the module network; the root keeps CARD_SIDE
#: How many next steps the first column after the clicked module shows. Six, not four,
#: because Einar asked for ALL the modules that take one gene to come downstream of the
#: gene card, and there are six of them: "so all the modules that take 1 gene and do
#: something should come downstream of the gene card module."
FLOW_FIRST = 6
#: How many each column after that shows: fewer, because by then it is a path being
#: traced rather than the question that was asked.
FLOW_MORE = 4
FLOW_GAP_X = 110              # room between columns, for the edges
FLOW_GAP_Y = 22
AI_ROOM = 60                  # right margin of the field: the AI symbol and a gap


@dataclass(frozen=True)
class Suggestion:
    """One card under the field.

    :param kind: ``module``, ``gene``, ``question``, ``strategy``, ``dataset``,
        ``panel`` or ``guide`` -- what the card opens.
    :param key: what opening it is told: a module key, a gene id, or an index title.
    :param title: the card's heading.
    :param subtitle: one line under it.
    :param detail: HTML for the details pane.
    """

    kind: str
    key: str
    title: str
    subtitle: str = ""
    detail: str = ""
    items: tuple = ()             # a ``group`` card's members
    #: Why this offer cannot be taken with what is in hand, or "" when it can.
    #: Instruction 98: "an offer that cannot work shown greyed WITH THE REASON".
    #: Separate from `modules.unavailable_reason`, which is about the window rather
    #: than about what the user is holding; a card is greyed for either.
    blocked: str = ""
    #: What is being carried into the module: ``(subject kind, value)``, so an offer
    #: opens pre-filled rather than empty.
    subject: tuple = ()


#: The word shown on a card for each kind, so a gene and a dataset with similar
#: names are never confused.
KIND_WORD = {"module": "Module", "gene": "Gene", "question": "Question",
             "strategy": "Method", "dataset": "Dataset", "panel": "Panel",
             "guide": "Guide", "node": "Card"}


#: How a group card names its members, singular and plural.
GROUP_WORDS = {"gene": ("gene", "genes"), "question": ("question", "questions"),
               "strategy": ("method", "methods"), "dataset": ("dataset", "datasets"),
               "panel": ("setting or panel", "settings and panels"),
               "guide": ("guide page", "guide pages"), "node": ("card", "cards")}
#: The order groups appear in, after the modules.
GROUP_ORDER = ("gene", "question", "strategy", "dataset", "panel", "guide", "node")
#: Modules shown before the groups: the ones the text names, at most this many.
MAX_MODULES = 4
#: Below this many results in all, they are shown one card each rather than grouped.
GROUP_WHEN_MORE_THAN = 6


def all_matches(text: str, genes=None, index=None, limit: int = MAX_CARDS) -> list:
    """Everything that matches `text`, one card each, best first (before grouping).

    Modules first when one is named outright -- typing "map" means the Map module --
    then the shared ranking over genes and the help index. The node catalogue is
    passed empty on purpose: those are the old Gene/List/Label cards, and the
    modules have replaced them as the front door.

    :param genes: `(gene_id, product)` pairs, or a callable taking the query.
    :param index: `help_index` rows, or a callable returning them.
    """
    q = (text or "").strip()
    if not q:
        return []
    out = []
    for m in MOD.search(q)[:6]:
        reason = MOD.unavailable_reason(m)
        detail = (f"<h3>{escape(m.name)}</h3><p>{escape(m.note)}.</p>"
                  f"<p><i>{escape(m.band)}</i></p>")
        if reason:
            detail += f"<p><b>Not available:</b> {escape(reason)}.</p>"
        out.append(Suggestion("module", m.key, m.name, m.note, detail))
    try:
        from . import nodes as N
        # The node catalogue minus its old Gene/List/Label cards: the modules replaced
        # those as the front door, but the QUESTION cards are still the only place the
        # registered questions are searchable. Passing an empty catalogue to drop the
        # cards also dropped every question -- "glucose" stopped finding the glucose
        # question -- so keep exactly those.
        questions = [n for n in N.catalogue() if getattr(n, "kind", "") == "question"]
        found = N.search(q, nodes=questions, genes=genes, index=index, limit=limit)
    except Exception:                       # a bad index must not break the field
        found = []
    for node in found:
        kind = getattr(node, "kind", "node") or "node"
        title = str(getattr(node, "title", ""))
        summary = str(getattr(node, "summary", "") or "")
        key = title if kind != "gene" else title.split()[0]
        detail = f"<h3>{escape(title)}</h3>"
        if summary:
            detail += f"<p>{escape(summary)}</p>"
        detail += f"<p><i>{escape(KIND_WORD.get(kind, kind.title()))}</i></p>"
        out.append(Suggestion(kind, key, title, summary, detail))
    seen, unique = set(), []
    for s in out:
        if (s.kind, s.key) not in seen:
            seen.add((s.kind, s.key))
            unique.append(s)
    return unique[:limit]


def gene_in(found) -> str:
    """The gene a list of results is about, or ``""`` when it is about no gene.

    :param found: suggestions from :func:`all_matches`.
    :returns: the identifier of the best gene among them.
    """
    for s in found:
        if s.kind == "gene":
            return s.key
    return ""


def module_offers_for_genes(found) -> list:
    """The modules to offer when what was typed found genes, the gene card first.

    Each carries the gene as its subject, so opening one fills it in rather than
    leaving it to be typed a second time. Ready ones first and a blocked one last is
    :func:`parascape.modules.downstream_of_subject`, which is spaCR's own order in its
    chaining strip.

    :param found: suggestions from :func:`all_matches`.
    :returns: module suggestions, at most :data:`MAX_MODULES`, or [] when no gene matched.
    """
    gene = gene_in(found)
    if not gene:
        return []
    from . import ports as PORTS
    keys = ["gene"] + [k for k in MOD.downstream_of_subject(PORTS.GENE, without="gene")]
    out = []
    for key in keys[:MAX_MODULES]:
        m = MOD.get(key)
        reason = MOD.unavailable_reason(m)
        detail = (f"<h3>{escape(m.name)}</h3><p>{escape(m.note)}.</p>"
                  f"<p>Opens with <b>{escape(gene)}</b> already in it.</p>"
                  f"<p><i>{escape(m.band)}</i></p>")
        if reason:
            detail += f"<p><b>Not available:</b> {escape(reason)}.</p>"
        out.append(Suggestion("module", key, m.name, m.note, detail,
                              blocked=reason, subject=(PORTS.GENE, gene)))
    return out


def suggestions(text: str, genes=None, index=None, limit: int = MAX_CARDS) -> list:
    """The cards the field shows for `text`: few, and each a different thing.

    Einar: "there are way too many options when I type map". "map" gave 24 cards: the
    Map module, three modules that merely mention maps, a panel also called Map, eight
    colour-map settings, four methods, a question and seven genes. Now: the modules the
    text names (at most four; ones that only mention it in passing only when nothing is
    named), a panel that repeats a module's name dropped, and everything else one
    **group card per kind** -- "7 genes", "4 methods" -- that opens to its members. A
    short answer (six results or fewer) is shown as it is.
    """
    found = all_matches(text, genes=genes, index=index, limit=limit)
    if not found:
        return []
    modules = [s for s in found if s.kind == "module"]
    named = [s for s in modules if (text or "").strip().lower() in s.title.lower()]
    if not named:
        # Einar, 2026-11-12: "if i type a gene into the navigator i get a few moduals that
        # pop up, one should be the gene card". Typing a gene used to offer gene CARDS --
        # nodes, which opened the pre-navigator evidence dock -- and no module at all. The
        # modules that read one gene are offered first now, the gene card first of them.
        #
        # Only when the text does not NAME a module: typing "map" means the Map module, and
        # burying it under gene modules because some gene product mentions a map would be a
        # worse answer than the one this replaces.
        modules = module_offers_for_genes(found) + modules
    modules = (named or modules)[:MAX_MODULES]
    module_names = {s.title.lower() for s in modules}
    rest = [s for s in found if s.kind != "module"
            and not (s.kind == "panel" and s.title.lower() in module_names)]
    if len(modules) + len(rest) <= GROUP_WHEN_MORE_THAN:
        return modules + rest
    out = list(modules)
    by_kind: dict = {}
    for s in rest:
        by_kind.setdefault(s.kind, []).append(s)
    for kind in sorted(by_kind, key=lambda k: GROUP_ORDER.index(k) if k in GROUP_ORDER else 99):
        members = by_kind[kind]
        if len(members) == 1:
            out.append(members[0])
            continue
        one, many = GROUP_WORDS.get(kind, (kind, kind + "s"))
        names = "".join(f"<li>{escape(card_name(m))}</li>" for m in members[:12])
        more = f"<p>…and {len(members) - 12} more.</p>" if len(members) > 12 else ""
        detail = (f"<h3>{len(members)} {escape(many)}</h3><ul>{names}</ul>{more}"
                  "<p><i>Click to see them one card each.</i></p>")
        out.append(Suggestion("group", kind, f"{len(members)} {many}",
                              ", ".join(card_name(m) for m in members[:3]), detail,
                              items=tuple(members)))
    return out


def _subtitle(text: str) -> str:
    """The informative part of an index summary.

    Dataset summaries open with their place in the slot tree -- "Slot tree ▸ Tg ▸
    relational measurements ▸ parasite-host ▸ host effect — host effect human K562;
    Toxoplasma invasion assay…". On a two-line card the breadcrumb used both lines
    and the card said nothing. What follows the dash is what the dataset is.
    """
    text = " ".join(str(text or "").split())
    if " — " in text and ("▸" in text.split(" — ", 1)[0]):
        text = text.split(" — ", 1)[1]
    return text


def _clip(text: str, n: int) -> str:
    """Shorten to about `n` characters at a word boundary, with an ellipsis."""
    text = " ".join(str(text or "").split())
    if len(text) <= n:
        return text
    cut = text[:n].rsplit(" ", 1)[0]
    return (cut or text[:n]).rstrip(" ·,;:—-") + "…"


def icon_path(suggestion: Suggestion) -> str:
    """The spaCR artwork for a card, as a file path, from the one icon provider.

    Kept as a function of the card so a caller never needs to know whether the card is a module
    or a kind of result; the tables live in :mod:`parascape.iconset` and nowhere else.

    :param suggestion: the card.
    """
    if suggestion.kind == "module":
        return iconset.module_icon_file(suggestion.key)
    if suggestion.kind == "group":
        return iconset.kind_icon_file(suggestion.key)
    return iconset.kind_icon_file(suggestion.kind)


def card_name(suggestion: Suggestion) -> str:
    """What the card says, in spaCR's one large light line.

    spaCR's cards show a module name and nothing else. ParaScape's results are not all
    modules: a search for "rhoptry" returns eight genes, and eight cards reading
    TGME49_... would be indistinguishable without hovering each. So a gene card shows
    its product, the thing that tells them apart, and the accession goes to the
    details panel with everything else.
    """
    if suggestion.kind == "gene" and suggestion.subtitle:
        return _clip(suggestion.subtitle, 60)
    return _clip(suggestion.title, 60)


#: spaCR's maturity colours, used for the cards' hover rim and fill.
STAGE_HOVER = {"stable": "#3B82F6", "beta": "#FF00FF", "alpha": "#00CEC8"}


def _rgba(hex_colour: str, alpha: float) -> str:
    c = QtGui.QColor(hex_colour)
    return f"rgba({c.red()}, {c.green()}, {c.blue()}, {alpha})"


def navigator_qss(palette: dict) -> str:
    """spaCR's Navigator stylesheet, value for value (`_EasyHomePage.refresh`)."""
    P = palette
    surface = "rgb(25, 25, 25)"               # spaCR's own field and card fill, not a token
    details = "rgba(25, 25, 25, 0.92)"
    stage_rules = "\n".join(
        f"""QPushButton#EasySuggestionCard[stage="{stage}"]:hover {{
    background-color: {_rgba(hue, 0.22)}; border: 1px solid {hue}; }}
QPushButton#EasySuggestionCard[stage="{stage}"]:pressed {{
    background-color: {_rgba(hue, 0.40)}; border: 1px solid {hue}; }}"""
        for stage, hue in STAGE_HOVER.items())
    return f"""
QWidget#NavigatorPage {{ background: transparent; }}
QWidget#NavigatorPage QWidget {{ font-family: "Open Sans"; font-weight: 300; }}
QLineEdit#NavigatorField {{
    background: {surface}; color: {P["fg"]}; border: none; border-radius: {FIELD_RADIUS}px;
    font-size: 16px; padding: 0px 24px; placeholder-text-color: {P["fg_muted"]};
    selection-background-color: {P["accent"]};
}}
QPushButton#EasySuggestionCard {{
    background: {surface}; color: {P["fg"]}; border: 1px solid {P["border"]};
    border-radius: 16px; padding: 12px; min-width: 0px; min-height: 0px; text-align: left;
}}
QPushButton#EasySuggestionCard QLabel {{ color: {P["fg"]}; background: transparent;
    border: none; font-weight: 300; }}
QPushButton#EasySuggestionCard[current="true"] {{ border: 1px solid {STAGE_HOVER["stable"]}; }}
{stage_rules}
QTextBrowser#NavigatorDetails {{
    color: {P["fg"]}; background: {details}; border: 1px solid {P["border"]};
    border-radius: 10px; padding: 8px; font-size: 16px;
}}
QTextBrowser#NavigatorDetails QScrollBar:vertical {{ background: transparent; width: 8px; }}
QTextBrowser#NavigatorDetails QScrollBar::handle:vertical {{ background: {P["border"]};
    border-radius: 4px; min-height: 24px; }}
QTextBrowser#NavigatorDetails QScrollBar::add-line:vertical,
QTextBrowser#NavigatorDetails QScrollBar::sub-line:vertical {{ height: 0px; }}
QPushButton#NavigatorBack {{ color: {P["fg"]}; background: {surface};
    border: 1px solid {P["border"]}; border-radius: 14px; padding: 5px 16px; font-size: 14px; }}
QPushButton#NavigatorBack:hover {{ border-color: {STAGE_HOVER["stable"]}; }}
QLabel#NavigatorFeedback {{ color: {P["fg_muted"]}; background: transparent; font-size: 14px; }}
QScrollBar#NavigatorRowScroll:horizontal {{ background: transparent; height: 10px; }}
QScrollBar#NavigatorRowScroll::handle:horizontal {{ background: {P["border"]};
    border-radius: 5px; min-width: 40px; }}
QScrollBar#NavigatorRowScroll::add-line:horizontal,
QScrollBar#NavigatorRowScroll::sub-line:horizontal {{ width: 0px; }}
"""


class SuggestionCard(QtWidgets.QPushButton):
    """One result, as spaCR draws a module: its icon, and under it its name, large and light.

    No description on the card, as in spaCR: hovering or focusing it puts the description
    in the panel at the foot of the page (`NavigatorPage.explain`).
    """

    hovered = QtCore.Signal(object)

    def __init__(self, suggestion: Suggestion, parent=None):
        super().__init__(parent)
        self.suggestion = suggestion
        self.setObjectName("EasySuggestionCard")
        unavailable = bool(suggestion.blocked) or (
            suggestion.kind == "module"
            and bool(MOD.unavailable_reason(MOD.get(suggestion.key))))
        self.setProperty("stage", "" if unavailable else "stable")
        self.setProperty("current", "false")
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.setAccessibleName(card_name(suggestion))
        self.setAccessibleDescription(suggestion.subtitle or suggestion.title)
        self._final = QtCore.QRect()
        self._animation = None

        col = QtWidgets.QVBoxLayout(self)
        col.setContentsMargins(12, 12, 12, 12)
        col.setSpacing(10)
        col.addStretch(1)
        icon = QtWidgets.QLabel()
        icon.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        icon.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        ratio = self.devicePixelRatioF() or 1.0
        icon.setPixmap(iconset.pixmap(icon_path(suggestion), ICON_PX, "dark", ratio))
        icon.setFixedHeight(ICON_PX)
        if unavailable:
            fade = QtWidgets.QGraphicsOpacityEffect(icon)
            fade.setOpacity(0.35)
            icon.setGraphicsEffect(fade)
        col.addWidget(icon)
        name = QtWidgets.QLabel(card_name(suggestion))
        name.setObjectName("CardName")
        name.setWordWrap(True)
        name.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter | QtCore.Qt.AlignmentFlag.AlignTop)
        name.setStyleSheet("font-size: 18px; font-weight: 300;")
        name.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        col.addWidget(name)
        col.addStretch(1)
        self.name_label = name
        self.icon_label = icon

    def fit(self, side: int) -> None:
        """Shrink the symbol so the name has room on a card smaller than a spaCR square.

        The module network squeezes a column of six into the height of one, and at that
        size a 96 px symbol and an 18 px name do not both fit: the name was drawn across
        the middle of the symbol, so "Information" read as "formatio". The symbol gives way,
        because the name is the part that tells two cards apart.

        :param side: the edge of the card, in device-independent pixels.
        """
        px = max(24, min(ICON_PX, side - 72))
        if self.icon_label.height() == px:
            return
        ratio = self.devicePixelRatioF() or 1.0
        self.icon_label.setPixmap(iconset.pixmap(icon_path(self.suggestion), px, "dark", ratio))
        self.icon_label.setFixedHeight(px)
        self.name_label.setStyleSheet(
            f"font-size: {max(11, min(18, side // 9))}px; font-weight: 300;")

    def set_current(self, on: bool) -> None:
        """Mark this card as the one Return would open. A property, so the stylesheet
        draws the rim; repolished, because Qt does not restyle on its own."""
        self.setProperty("current", "true" if on else "false")
        self.style().unpolish(self)
        self.style().polish(self)

    def enterEvent(self, event):                          # noqa: N802
        """Describe this card in the panel below, as spaCR's `explain_card` does."""
        self.hovered.emit(self)
        super().enterEvent(event)

    def focusInEvent(self, event):                        # noqa: N802
        """Keyboard focus describes the card too, so arrowing across reads each one."""
        self.hovered.emit(self)
        super().focusInEvent(event)


class NavigatorPage(QtWidgets.QWidget):
    """spaCR's Navigator: one centred field, the cards above it, the details at the foot.

    Placed by hand in `resizeEvent`, as spaCR's `_EasyHomeField._arrange` does, not by a
    layout -- which is also what lets the cards grow into place.

    Searching happens on Return, not on every keystroke: typing changes nothing below the
    field until it is submitted, and a second Return opens the highlighted card.

    `chosen(kind, key)` fires when a card is opened. The assistant is the ring-and-dot switch at
    the end of the field, as in spaCR: when it is on and a search finds nothing, the question is
    answered in the panel below the field.
    """

    chosen = QtCore.Signal(str, str)

    PROMPT = "Type a gene, a question, a module or a dataset, then press Return."

    def __init__(self, genes=None, index=None, parent=None):
        super().__init__(parent)
        self.setObjectName("NavigatorPage")
        self._genes, self._index = genes, index
        self._cards: list = []
        self._current = -1
        self._searched = None
        self.setStyleSheet(navigator_qss(TH.palette_for("dark")))

        self.field = QtWidgets.QLineEdit(self)
        self.field.setObjectName("NavigatorField")
        self.field.setFixedHeight(FIELD_H)
        self.field.setPlaceholderText("A gene, a question, a module or a dataset  ·  Return to search  ·  ◉ asks the assistant")
        self.field.setAccessibleName("Search ParaScape")
        self.field.setToolTip("Type a gene, a product, a question, a module or a dataset and "
                              "press Return. Everything here works with no assistant at all.")
        # No clear button: Einar does not want the small cross anywhere. Escape empties the field
        # (the key handling below), so the field needs no helper here.
        self.field.setClearButtonEnabled(False)
        self.field.setTextMargins(0, 0, AI_ROOM, 0)
        self.ai_button = AiButton(self.field)
        self._ai_state_init()

        self.feedback = QtWidgets.QLabel(self.PROMPT, self)
        self.feedback.setObjectName("NavigatorFeedback")
        self.feedback.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        self.feedback.setWordWrap(True)

        self.details = QtWidgets.QTextBrowser(self)
        self.details.setObjectName("NavigatorDetails")
        self.details.setOpenLinks(False)
        self.details.hide()

        self._scroll = QtWidgets.QScrollBar(QtCore.Qt.Orientation.Horizontal, self)
        self._scroll.setObjectName("NavigatorRowScroll")
        self._scroll.valueChanged.connect(lambda _v: self._arrange(animate=False))
        self._scroll.hide()

        # The module network (instruction 91): one click on a module shows it with the
        # modules it leads to; a click on one of those shows what follows it; a second
        # click on a highlighted module opens it. `_path` is the chain clicked so far.
        self._path: list = []
        self._flow_cards: dict = {}
        self._columns: list = []
        self._top: list = []                    # the search results, to come back to
        self._expanded = None                   # the group card being shown open, if any
        self.back = QtWidgets.QPushButton("← Back", self)
        self.back.setObjectName("NavigatorBack")
        self.back.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.back.setToolTip("Back to the search results (Escape).")
        self.back.clicked.connect(self.go_back)
        self.back.hide()
        self.flow_hint = QtWidgets.QLabel("", self)
        self.flow_hint.setObjectName("NavigatorFeedback")
        self.flow_hint.hide()

        self.field.textEdited.connect(self._edited)
        self.field.textEdited.connect(self.ai_button.typed)
        self.field.returnPressed.connect(self.submit)
        self.field.installEventFilter(self)

    # ------------------------------------------------------------------ layout
    def _field_rect(self) -> QtCore.QRect:
        w, h = self.width(), self.height()
        width = min(max(FIELD_MIN_W, 0), max(180, w - 48))
        y = max(16, (h - FIELD_H) // 2)
        return QtCore.QRect((w - width) // 2, y, width, FIELD_H)

    def resizeEvent(self, event):                          # noqa: N802
        """Place the field, the cards above it and the details below, as spaCR does."""
        super().resizeEvent(event)
        self._arrange(animate=False)

    def _arrange(self, animate: bool) -> None:
        if self._path:
            self._arrange_flow(animate)
            return
        w, h = self.width(), self.height()
        field = self._field_rect()
        self.field.setGeometry(field)
        self._place_ai()
        self.feedback.setGeometry(field.left(), field.bottom() + 10, field.width(), 44)

        details_w = field.width()
        self.details.setGeometry((w - details_w) // 2, h - DETAILS_H - 16, details_w, DETAILS_H)

        room = max(64, field.top() - 40)
        side = min(CARD_SIDE, room, max(96, w - 64))
        top = max(16, field.top() - side - 24)
        n = len(self._cards)
        total = n * side + max(0, n - 1) * CARD_GAP
        left = 16
        span = max(1, w - 32)
        overflow = max(0, total - span)
        self._scroll.blockSignals(True)
        self._scroll.setRange(0, overflow)
        self._scroll.setPageStep(span)
        self._scroll.setSingleStep(side + CARD_GAP)
        self._scroll.blockSignals(False)
        self._scroll.setGeometry(left + 8, top + side + 5, max(1, span - 16), 10)
        self._scroll.setVisible(bool(overflow))
        start = left + max(0, (span - total) // 2) - self._scroll.value()
        for i, card in enumerate(self._cards):
            card._final = QtCore.QRect(start + i * (side + CARD_GAP), top, side, side)
            if animate and _motion_allowed():
                self._grow(card, i)
            elif card._animation is None:
                card.setGeometry(card._final)
            visible = card._final.right() >= left and card._final.left() <= left + span
            card.setVisible(visible)
        self._fade_edges(left, span)
        self.details.raise_()

    def _grow(self, card, i: int) -> None:
        """spaCR's arrival: each card grows from its centre, one after another."""
        final = card._final
        card.setGeometry(QtCore.QRect(final.center(), QtCore.QSize(1, 1)))
        animation = QtCore.QPropertyAnimation(card, b"geometry", card)
        animation.setDuration(GROW_MS)
        animation.setStartValue(QtCore.QRect(final.center(), QtCore.QSize(1, 1)))
        animation.setEndValue(final)
        animation.setEasingCurve(QtCore.QEasingCurve.Type.OutCubic)

        def done(c=card):
            c._animation = None
            c.setGeometry(c._final)
        animation.finished.connect(done)
        card._animation = animation
        QtCore.QTimer.singleShot(i * GROW_STAGGER_MS, card, animation.start)

    def _fade_edges(self, left: int, span: int) -> None:
        """Fade cards at an edge the row continues past, as spaCR's opacity mask does."""
        for card in self._cards:
            r = card._final
            fade = 1.0
            if self._scroll.value() and r.left() < left + 44:
                fade = max(0.25, (r.right() - left) / max(1, r.width()))
            if self._scroll.value() < self._scroll.maximum() and r.right() > left + span - 44:
                fade = min(fade, max(0.25, (left + span - r.left()) / max(1, r.width())))
            effect = card.graphicsEffect()
            if fade >= 0.999:
                if effect is not None:
                    card.setGraphicsEffect(None)
            else:
                if not isinstance(effect, QtWidgets.QGraphicsOpacityEffect):
                    effect = QtWidgets.QGraphicsOpacityEffect(card)
                    card.setGraphicsEffect(effect)
                effect.setOpacity(fade)

    def _place_ai(self) -> None:
        """Put the AI symbol inside the right end of the field, centred, as spaCR does."""
        f = self.field
        self.ai_button.move(f.width() - self.ai_button.width() - 12,
                            (f.height() - self.ai_button.height()) // 2)
        self.ai_button.raise_()

    def wheelEvent(self, event):                           # noqa: N802
        """Scroll the card row sideways with the mouse wheel when it overflows."""
        if self._scroll.isVisible():
            delta = event.angleDelta().y() or event.angleDelta().x()
            self._scroll.setValue(self._scroll.value() - delta)
            event.accept()
            return
        super().wheelEvent(event)

    # ------------------------------------------------------------------ search
    def _edited(self, _text: str) -> None:
        """Typing changes nothing below the field until Return -- Einar's request."""
        text = self.field.text().strip()
        if not text:
            self._cancel_ai()
            self._clear()
            self.feedback.setText(self.PROMPT)
        elif text != self._searched:
            self.feedback.setText("Press Return to search.")

    def submit(self) -> None:
        """Return: search what is typed, or, when that is already on screen, open the card."""
        text = self.field.text().strip()
        if not text:
            return
        if text == self._searched and self._cards:
            self._clicked(max(self._current, 0))
        else:
            self.refresh()

    def open_current(self) -> None:
        """Open the highlighted card, searching first if the field has not been searched.

        Without that, a fast typist pressing Return would open the card for what they had
        typed a moment ago rather than what is in the field.
        """
        text = self.field.text().strip()
        if text and text != self._searched:
            self.refresh()
        if self._cards:
            self._open(max(self._current, 0))

    def _clear(self) -> None:
        for card in self._cards:
            card.setParent(None)
            card.deleteLater()
        self._cards = []
        self._current = -1
        self._searched = None
        self.details.hide()
        self._scroll.hide()

    def refresh(self) -> None:
        """Search for what is in the field and grow the cards into place."""
        text = self.field.text().strip()
        self._cancel_ai()
        self._leave_flow()
        self._clear()
        if not text:
            self.feedback.setText(self.PROMPT)
            return
        self._searched = text
        self.set_subject()
        found = suggestions(text, genes=self._genes, index=self._index)
        # What is in hand, for everything opened from here on: the gene the search is
        # about. A module card carries it on itself, but a card the user reaches another
        # way -- Return, the arrow keys, a group opened and closed -- must find it too.
        for s in found:
            if s.subject:
                self.set_subject(*s.subject)
                break
        if not found:
            if self.ai_button.isChecked():
                self.feedback.setText(f"Nothing here matches “{text}”. Asking the assistant…")
                self.ask_ai(text)
            else:
                self.feedback.setText(
                    f"Nothing here matches “{text}”. Try a gene id or a product name, or "
                    "click the ring at the end of the field to turn the assistant on and ask.")
            return
        self._top = found
        self._expanded = None
        self._show(found)
        counts: dict = {}
        for s in found:
            kind = s.key if s.kind == "group" else s.kind
            counts[kind] = counts.get(kind, 0) + (len(s.items) if s.kind == "group" else 1)
        words = []
        for kind, n in counts.items():
            one, many = GROUP_WORDS.get(kind, (KIND_WORD.get(kind, kind).lower(),
                                               KIND_WORD.get(kind, kind).lower() + "s"))
            words.append(f"{n} {one if n == 1 else many}")
        self.feedback.setText(" · ".join(words) + " — click a module to see what follows it")

    def _show(self, found: list) -> None:
        """Put one card per suggestion in the row and grow them into place."""
        for card in self._cards:
            card.setParent(None)
            card.deleteLater()
        self._cards = []
        self._current = -1
        for i, s in enumerate(found):
            card = SuggestionCard(s, self)
            card.clicked.connect(lambda _=False, k=i: self._clicked(k))
            card.hovered.connect(self.explain)
            card.installEventFilter(self)
            card.show()
            self._cards.append(card)
        self._scroll.setValue(0)
        self._arrange(animate=True)
        self._set_current(0)

    def _clicked(self, i: int) -> None:
        """A click on a result: a group opens, a module shows its network, the rest open."""
        s = self._cards[i].suggestion
        if s.kind == "group":
            self.expand(s)
        elif s.kind == "module":
            self.explore(s.key, s.subject)
        else:
            self._open(i)

    def expand(self, group: Suggestion) -> None:
        """Show a group card's members one card each; Back returns to the results."""
        self._expanded = group
        self._show(list(group.items))
        self.feedback.setText(f"{group.title} for “{self._searched}” — Escape goes back")
        self.back.show()
        self._place_back()

    def go_back(self) -> None:
        """One step back: out of the module network, then out of an opened group."""
        if self._path:
            self._leave_flow()
            self._arrange(animate=False)
            return
        if self._expanded is not None:
            self._expanded = None
            self.back.hide()
            self._show(self._top)
            self.feedback.setText("Click a module to see what follows it.")

    def reset(self) -> None:
        """Empty the Navigator, as it is every time Home is shown (Einar's request)."""
        self._leave_flow()
        self._cancel_ai()
        self._clear()
        self._top, self._expanded = [], None
        self.set_subject()
        self.back.hide()
        self.field.clear()
        self.feedback.setText(self.PROMPT)
        self._arrange(animate=False)

    # ------------------------------------------------------------ the assistant
    def _ai_state_init(self) -> None:
        """Set up the state of the assistant: no provider yet, nothing running.

        The provider is remembered between runs, as spaCR remembers it, so the first click of a
        later session turns the assistant on without asking again.
        """
        self._ai_provider = self._saved_provider()
        self._ai_worker = None
        self._ai_retired: list = []
        self._ai_text = ""
        self._ai_timer = QtCore.QTimer(self)
        self._ai_timer.setSingleShot(True)
        self._ai_timer.setInterval(50)
        self._ai_timer.timeout.connect(self._render_ai)
        self.ai_button.set_provider(self._ai_provider or "parascape")
        self.ai_button.clicked.connect(self._toggle_ai)
        self.ai_button.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.ai_button.customContextMenuRequested.connect(self.choose_provider)
        self._update_ai_tooltip()

    @staticmethod
    def _saved_provider() -> str:
        """The provider chosen in an earlier session, or an empty string when none was."""
        from . import chat
        value = str(QtCore.QSettings("parascape", "parascape").value("navigator/ai_provider", "")
                    or "")
        return value if value in {p.name for p in chat.PROVIDERS} else ""

    def _update_ai_tooltip(self) -> None:
        """Say what the AI symbol is doing now and how to change it, as spaCR does."""
        from . import chat
        label = next((p.label for p in chat.PROVIDERS if p.name == self._ai_provider), "")
        if self.ai_button.isChecked():
            self.ai_button.setToolTip(f"AI: {label}. Click to turn it off; right-click to change "
                                      "the provider. A search that finds nothing is answered "
                                      "below the field.")
        else:
            self.ai_button.setToolTip("The assistant is off. Click to turn it on; right-click to "
                                      "choose which assistant answers.")

    def _toggle_ai(self, on: bool) -> None:
        """Turn the assistant on or off from a click on the symbol.

        Turning it off cancels an answer in progress. Turning it on with no provider chosen yet
        asks which one, rather than guessing, and leaves it off if the question is dismissed.

        :param on: the state the symbol now shows.
        """
        if not on:
            was_running = self._ai_worker is not None
            self._cancel_ai()
            if was_running:
                self.feedback.setText("Cancelled.")
        elif not self._ai_provider:
            self.ai_button.setChecked(False)
            self.choose_provider(turn_on=True)
            return
        self.ai_button.set_provider(self._ai_provider or "parascape")
        self._update_ai_tooltip()

    def choose_provider(self, *_args, turn_on: bool = True) -> str:
        """Ask which assistant should answer, and turn the symbol on when one is chosen.

        :param turn_on: switch the assistant on after a choice, as spaCR does.
        :returns: the provider chosen, or an empty string when the question was dismissed.
        """
        from . import chat
        dialog = ProviderDialog(chat.PROVIDERS, self._ai_provider, self)
        try:
            from . import glass
            glass.dress(dialog)
        except Exception:
            pass
        chosen = ""
        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            chosen = dialog.selected()
            self._cancel_ai()
            self._ai_provider = chosen
            QtCore.QSettings("parascape", "parascape").setValue("navigator/ai_provider", chosen)
            self.ai_button.set_provider(chosen)
            if turn_on:
                self.ai_button.setChecked(True)
        dialog.deleteLater()
        self._update_ai_tooltip()
        return chosen

    @staticmethod
    def _ai_system() -> str:
        """The standing instructions for the assistant in the Navigator.

        The caveats of the chat panel, plus the list of modules so that an answer can point to the
        module that does what was asked, and the rule that the Navigator only advises.
        """
        from . import chat
        modules = "; ".join(f"{m.name} ({m.note})" for m in MOD.MODULES)
        return (f"{chat.GROUNDING}\n\nThe Navigator is an advice interface. Answer in text, and "
                "name the ParaScape module that would do what is asked. Do not claim to run "
                f"anything.\n\nThe modules are: {modules}.")

    def ask_ai(self, text: str) -> None:
        """Send one question to the chosen assistant and write the answer below the field.

        Nothing is sent unless the symbol is on and a provider is chosen. A provider that is not
        installed says how to install it, in the panel, instead of failing silently.

        :param text: the question.
        """
        from . import chat
        provider = next((p for p in chat.PROVIDERS if p.name == self._ai_provider), None)
        if provider is None:
            self.feedback.setText("No assistant chosen. Right-click the ring to choose one.")
            return
        self._cancel_ai()
        self.details.show()
        self.details.raise_()
        if not provider.available():
            self.details.setMarkdown(
                f"**{provider.label}** is not installed.\n\nInstall it with\n\n"
                f"    {provider.install_hint}\n\nthen sign in with\n\n    {provider.login_command}"
                "\n\nor right-click the ring and choose the free ParaScape AI, which needs nothing.")
            return
        self._ai_text = ""
        self.details.setMarkdown("Thinking…")
        worker = chat._Worker(provider, text, self._ai_system(), self)
        worker.chunk.connect(self._ai_chunk)
        worker.done.connect(self._ai_done)
        self._ai_worker = worker
        self.ai_button.set_busy(True)
        worker.start()

    def _ai_chunk(self, text: str) -> None:
        """Take a piece of the answer and redraw the panel, at most every fifty milliseconds."""
        if self.sender() is not self._ai_worker:
            return
        self._ai_text += text
        if not self._ai_timer.isActive():
            self._ai_timer.start()

    def _render_ai(self) -> None:
        """Draw the answer so far in the details panel, as formatted text."""
        self.details.setMarkdown(self._ai_text or "Thinking…")
        self.details.show()
        self.details.raise_()

    def _ai_done(self) -> None:
        """The answer is complete: stop the waiting motion and draw the whole of it."""
        if self.sender() is not self._ai_worker:
            return
        self._ai_timer.stop()
        self._retire(self._ai_worker)
        self._ai_worker = None
        self.ai_button.set_busy(False)
        self._ai_text = self._ai_text.strip() or "(no output; try again or choose another assistant)"
        self._render_ai()

    def _retire(self, worker) -> None:
        """Keep a stopped worker alive until its thread has really finished.

        Destroying a `QThread` that is still running aborts the interpreter, so a cancelled worker
        is kept here and released by its own `finished` signal.

        :param worker: the worker to let go of.
        """
        if worker is None:
            return
        for signal in (worker.chunk, worker.done):
            try:
                signal.disconnect()
            except (RuntimeError, TypeError):
                pass
        self._ai_retired.append(worker)

        def release(w=worker):
            """Drop the finished worker from the list and let Qt delete it."""
            if w in self._ai_retired:
                self._ai_retired.remove(w)
            w.deleteLater()
        if worker.isFinished():
            release()
        else:
            worker.finished.connect(release)

    def _cancel_ai(self) -> None:
        """Stop any answer in progress, quietly, without waiting for the thread."""
        self._ai_timer.stop()
        worker, self._ai_worker = self._ai_worker, None
        if worker is not None:
            worker.stop()
            self._retire(worker)
        self.ai_button.set_busy(False)

    def hideEvent(self, event):                            # noqa: N802
        """Stop an answer that nobody is looking at when another page comes on screen."""
        self._cancel_ai()
        super().hideEvent(event)

    def closeEvent(self, event):                           # noqa: N802
        """Never leave a running thread behind when the page goes away."""
        self._cancel_ai()
        for worker in list(self._ai_retired):
            worker.wait(2000)
        super().closeEvent(event)

    # ------------------------------------------------------------ module network
    @property
    def subject(self) -> tuple:
        """What is in hand, as ``(subject kind, value)``; ``("", "")`` when nothing is.

        Set by a search that found a gene, and carried into whatever module is opened, so
        a module reached through the Navigator opens already reading it. Einar's journey --
        type a gene, click the card it offers -- is this value travelling.
        """
        return getattr(self, "_subject", ("", ""))

    def set_subject(self, kind: str = "", value: str = "") -> None:
        """Say what is in hand. Empty arguments mean nothing is.

        :param kind: a subject from :mod:`parascape.ports`, such as ``"gene"``.
        :param value: what is in it.
        """
        self._subject = (str(kind or ""), str(value or ""))

    def explore(self, key: str, subject: tuple = ()) -> None:
        """First click on a module: show it with the modules it leads to.

        spaCR: "Click once to open the module network. Click twice to open the module."

        :param key: the module clicked.
        :param subject: ``(kind, value)`` to carry into it, from the card that was clicked.
        """
        if subject:
            self.set_subject(*subject)
        self._path = [key]
        self._build_flow()

    def _flow_clicked(self, column: int, key: str) -> None:
        """A node clicked: the highlighted one opens; another becomes the next step shown."""
        if column < len(self._path) and self._path[column] == key:
            self._open_module(key)
            return
        self._path = self._path[:column] + [key]
        self._build_flow()

    def _open_module(self, key: str) -> None:
        """Open a module from the network, unless something stops it -- and then say what.

        Two different things can stop it and they are reported apart: the window has no
        such view (the catalogue knows), or the module cannot work on what is in hand (the
        declared subjects know, through `ports.check_ready`).
        """
        reason = MOD.unavailable_reason(MOD.get(key))
        if reason:
            self.flow_hint.setText(f"{MOD.get(key).name} is not available yet: {reason}.")
            return
        blocked, fix = self._readiness(key)
        if blocked:
            self.flow_hint.setText(f"{MOD.get(key).name} cannot work on this: {blocked}. {fix}")
            return
        self.chosen.emit("module", key)

    def _readiness(self, key: str) -> tuple:
        """``(why it cannot work on what is in hand, what to do)``; ``("", "")`` when it can.

        :param key: the module key being offered.
        """
        kind, value = self.subject
        if not kind:
            return "", ""
        from . import next_steps as NS
        for offer in NS.offers(kind, value, without=""):
            if offer.module == key:
                return offer.blocked, offer.fix
        return "", ""

    def _build_flow(self) -> None:
        """Make the cards for the chain clicked so far and what follows its last step.

        The first column after the clicked module shows everything that follows it, up to
        :data:`FLOW_FIRST`; later columns show fewer. Each card is checked against what is
        in hand and a card that cannot work is drawn greyed, carrying the reason and the fix
        in the panel at the foot -- instruction 98, and spaCR's own rule for its chaining
        strip: an offer either works, or is shown greyed with the reason it cannot.
        """
        old = dict(self._flow_cards)
        self._columns = [[self._path[0]]]
        for step in self._path:
            limit = FLOW_FIRST if len(self._columns) == 1 else FLOW_MORE
            following = list(MOD.downstream(step))[:limit]
            if not following:
                break
            self._columns.append(following)
            if len(self._columns) > len(self._path):
                break
        self._flow_cards = {}
        fresh = []
        kind, value = self.subject
        for c, column in enumerate(self._columns):
            for key in column:
                ident = (c, key)
                card = old.pop(ident, None)
                if card is None:
                    m = MOD.get(key)
                    detail = (f"<h3>{escape(m.name)}</h3><p>{escape(m.note)}.</p>"
                              f"<p><i>{escape(m.band)}</i></p>")
                    blocked, fix = self._readiness(key)
                    blocked = blocked or MOD.unavailable_reason(m)
                    if blocked:
                        detail += (f"<p><b>Not ready:</b> {escape(blocked)}.</p>"
                                   + (f"<p>{escape(fix)}</p>" if fix else ""))
                    elif value:
                        detail += (f"<p>Opens with <b>{escape(value)}</b> already in it.</p>")
                    nxt = MOD.downstream(key)
                    if nxt:
                        detail += ("<p>Leads to: " + ", ".join(escape(MOD.get(k).name)
                                                               for k in nxt) + ".</p>")
                    card = SuggestionCard(
                        Suggestion("module", key, m.name, m.note, detail, blocked=blocked,
                                   subject=(kind, value) if kind else ()), self)
                    card.clicked.connect(lambda _=False, col=c, k=key: self._flow_clicked(col, k))
                    card.hovered.connect(self.explain)
                    fresh.append(card)
                card.set_current(c < len(self._path) and self._path[c] == key)
                card.show()
                self._flow_cards[ident] = card
        for card in old.values():
            card.setParent(None)
            card.deleteLater()
        for card in self._cards:
            card.hide()
        for widget in (self.field, self.feedback, self._scroll):
            widget.hide()
        self.back.show()
        self.flow_hint.show()
        last = MOD.get(self._path[-1])
        self.flow_hint.setText(
            f"{last.name}: click a highlighted module again to open it, or a next step to "
            "see what follows it.")
        self.explain(self._flow_cards[(len(self._path) - 1, self._path[-1])])
        for card in fresh:
            card._animation = None
        self._arrange_flow(animate=True, fresh=fresh)
        self.update()

    def _leave_flow(self) -> None:
        if not self._path and not self._flow_cards:
            return
        self._path = []
        self._columns = []
        for card in self._flow_cards.values():
            card.setParent(None)
            card.deleteLater()
        self._flow_cards = {}
        self.flow_hint.hide()
        self.back.setVisible(self._expanded is not None)
        self.field.show()
        self.feedback.show()
        for card in self._cards:
            card.show()
        self.update()

    def _arrange_flow(self, animate: bool, fresh=()) -> None:
        """Columns left to right: the module, what it leads to, what the chosen one leads to."""
        w, h = self.width(), self.height()
        self._place_back()
        self.flow_hint.setGeometry(140, 22, max(100, w - 160), 28)
        top = 70
        bottom = h - DETAILS_H - 32
        room = max(120, bottom - top)
        xs, sides = [], []
        x = 32
        for c, column in enumerate(self._columns):
            n = len(column)
            side = CARD_SIDE if c == 0 else FLOW_SIDE
            side = int(min(side, (room - (n - 1) * FLOW_GAP_Y) / max(1, n)))
            xs.append(x)
            sides.append(side)
            x += side + FLOW_GAP_X
        overflow = max(0, x - FLOW_GAP_X - (w - 32))
        for c, column in enumerate(self._columns):
            n, side = len(column), sides[c]
            height = n * side + (n - 1) * FLOW_GAP_Y
            y0 = top + max(0, (room - height) // 2)
            for i, key in enumerate(column):
                card = self._flow_cards.get((c, key))
                if card is None:
                    continue
                card._final = QtCore.QRect(xs[c] - overflow, y0 + i * (side + FLOW_GAP_Y),
                                           side, side)
                card.fit(side)
                card.setVisible(card._final.right() > 0)
                if animate and card in fresh and _motion_allowed():
                    self._grow(card, i)
                elif card._animation is None:
                    card.setGeometry(card._final)
        self.details.setGeometry((w - self._field_rect().width()) // 2, h - DETAILS_H - 16,
                                 self._field_rect().width(), DETAILS_H)
        self.details.raise_()
        self.update()

    def _place_back(self) -> None:
        hint = self.back.sizeHint()
        self.back.setGeometry(24, 20, max(96, hint.width()), hint.height())
        self.back.raise_()

    def paintEvent(self, event):                           # noqa: N802
        """The edges of the module network, from each step to what it leads to."""
        super().paintEvent(event)
        if not self._path or len(self._columns) < 2:
            return
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        hue = QtGui.QColor(STAGE_HOVER["stable"])
        for c in range(1, len(self._columns)):
            parent = self._flow_cards.get((c - 1, self._path[c - 1]))
            if parent is None:
                continue
            start = QtCore.QPointF(parent._final.right(), parent._final.center().y())
            for key in self._columns[c]:
                child = self._flow_cards.get((c, key))
                if child is None:
                    continue
                end = QtCore.QPointF(child._final.left(), child._final.center().y())
                on_path = c < len(self._path) and self._path[c] == key
                colour = QtGui.QColor(hue)
                colour.setAlphaF(0.95 if on_path else 0.45)
                painter.setPen(QtGui.QPen(colour, 2.4 if on_path else 1.4))
                bend = (end.x() - start.x()) * 0.5
                path = QtGui.QPainterPath(start)
                path.cubicTo(QtCore.QPointF(start.x() + bend, start.y()),
                             QtCore.QPointF(end.x() - bend, end.y()), end)
                painter.drawPath(path)
        painter.end()

    def explain(self, card) -> None:
        """Show a card's description in the panel at the foot of the page."""
        self.details.setHtml(card.suggestion.detail)
        self.details.show()
        self.details.raise_()

    def _set_current(self, i: int) -> None:
        if not self._cards:
            return
        i = max(0, min(len(self._cards) - 1, i))
        if 0 <= self._current < len(self._cards):
            self._cards[self._current].set_current(False)
        self._current = i
        card = self._cards[i]
        card.set_current(True)
        self._reveal(card)
        self.explain(card)

    def _reveal(self, card) -> None:
        """Scroll the row so a card arrowed onto is inside it."""
        span = max(1, self.width() - 32)
        r = card._final
        if r.left() < 16:
            self._scroll.setValue(self._scroll.value() - (16 - r.left()))
        elif r.right() > 16 + span:
            self._scroll.setValue(self._scroll.value() + (r.right() - 16 - span))

    def current(self):
        """The suggestion Return would open, or None."""
        return self._cards[self._current].suggestion if 0 <= self._current < len(self._cards) else None

    def _open(self, i: int) -> None:
        """Open the card at `i`, refusing with a reason when it cannot be opened."""
        s = self._cards[i].suggestion
        if s.kind == "group":
            self.expand(s)
            return
        if s.subject:
            self.set_subject(*s.subject)
        if s.kind == "module" and MOD.unavailable_reason(MOD.get(s.key)):
            self._set_current(i)
            self.feedback.setText(f"{s.title} is not available yet: "
                                  f"{MOD.unavailable_reason(MOD.get(s.key))}.")
            return
        self.chosen.emit(s.kind, s.key)

    # ----------------------------------------------------------------- keyboard
    def keyPressEvent(self, event):                        # noqa: N802
        """In the module network, Escape goes back and Return opens the highlighted module."""
        K = QtCore.Qt.Key
        if self._path and event.key() == K.Key_Escape:
            self.go_back()
            return
        if self._path and event.key() in (K.Key_Return, K.Key_Enter):
            self._open_module(self._path[-1])
            return
        super().keyPressEvent(event)

    def eventFilter(self, obj, event):                    # noqa: N802
        """Up or Down from the field enters the cards; arrows move; Escape clears or returns."""
        if obj is self.field and event.type() in (QtCore.QEvent.Type.Resize,
                                                  QtCore.QEvent.Type.Show):
            self._place_ai()
            return False
        if event.type() != QtCore.QEvent.Type.KeyPress:
            return False
        key = event.key()
        K = QtCore.Qt.Key
        if obj is self.field:
            if key in (K.Key_Up, K.Key_Down) and self._cards:
                self._cards[max(self._current, 0)].setFocus()
                return True
            if key == K.Key_Escape:
                if self._expanded is not None:
                    self.go_back()
                    return True
                self.field.clear()
                self._edited("")
                return True
            return False
        if obj in self._cards:
            if key in (K.Key_Right, K.Key_Tab):
                self._set_current(self._current + 1)
                self._cards[self._current].setFocus()
                return True
            if key in (K.Key_Left, K.Key_Backtab):
                self._set_current(self._current - 1)
                self._cards[self._current].setFocus()
                return True
            if key in (K.Key_Up, K.Key_Down, K.Key_Escape):
                self.field.setFocus()
                return True
            if key in (K.Key_Return, K.Key_Enter):
                self._open(self._current)
                return True
        return False


def _motion_allowed() -> bool:
    """spaCR grows its cards only when animation is on and Reduce motion is off."""
    try:
        from .backdrops import _compat as AC
        return AC.get_ambient_enabled() and not AC.get_reduce_motion()
    except Exception:
        return False
