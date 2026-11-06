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
editor and workflow builder are not ported; its AI button is, routed to
ParaScape's assistant.

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

from . import modules as MOD
from . import theme as TH

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
FLOW_GAP_X = 110              # room between columns, for the edges
FLOW_GAP_Y = 22


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


#: spaCR's card artwork, by ParaScape module and by the kind of a result. The images
#: are spaCR's own (white line art on transparency, under spaCR's licence), so a
#: ParaScape card carries the same marks as a spaCR one rather than look-alikes.
_MODULE_ICON = {
    "gene": "profiler", "structure": "layer_viewer", "partners": "graph_builder",
    "expression": "pca", "map": "umap", "network": "graph_builder", "cell": "image_scatter",
    "gallery": "trellis", "questions": "experiment_design", "infer": "ml_analyze",
    "strategies": "model_zoo", "scorecard": "model_compare", "runs": "run_history",
    "fill": "regression", "drop": "download", "screen": "hit_list", "list": "tabulate",
    "export": "methods_export", "datasets": "data_manager", "information": "feature_dict",
    "jobs": "queue", "console": "report",
}
_KIND_ICON = {"gene": "profiler", "dataset": "db_browser", "question": "experiment_design",
              "strategy": "model_zoo", "guide": "report", "panel": "settings", "node": "settings"}


def icon_path(suggestion: Suggestion) -> str:
    """The spaCR artwork for a card, as a file path."""
    from pathlib import Path
    if suggestion.kind == "module":
        name = _MODULE_ICON.get(suggestion.key)
    elif suggestion.kind == "group":
        name = _KIND_ICON.get(suggestion.key)
    else:
        name = _KIND_ICON.get(suggestion.kind)
    name = name or "settings"
    return str(Path(__file__).with_name("data") / "icons" / "spacr" / f"{name}.png")


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
QPushButton#NavigatorAsk {{ color: {P["fg"]}; background: transparent; border: none;
    padding: 0px; font-size: 16px; }}
QPushButton#NavigatorAsk:hover {{ color: {STAGE_HOVER["stable"]}; }}
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
        unavailable = (suggestion.kind == "module"
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
        pixmap = QtGui.QPixmap(icon_path(suggestion))
        if not pixmap.isNull():
            ratio = self.devicePixelRatioF() or 1.0
            pixmap = pixmap.scaled(int(ICON_PX * ratio), int(ICON_PX * ratio),
                                   QtCore.Qt.AspectRatioMode.KeepAspectRatio,
                                   QtCore.Qt.TransformationMode.SmoothTransformation)
            pixmap.setDevicePixelRatio(ratio)
            icon.setPixmap(pixmap)
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

    `chosen(kind, key)` fires when a card is opened; `ask(text)` when the assistant is.
    """

    chosen = QtCore.Signal(str, str)
    ask = QtCore.Signal(str)

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
        self.field.setPlaceholderText("A gene, a question, a module or a dataset  ·  Return to search")
        self.field.setAccessibleName("Search ParaScape")
        self.field.setClearButtonEnabled(True)
        self.field.setTextMargins(0, 0, 64, 0)
        self.ask_button = QtWidgets.QPushButton("Ask", self.field)
        self.ask_button.setObjectName("NavigatorAsk")
        self.ask_button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.ask_button.setToolTip("Ask the assistant instead of picking a card. Everything else "
                                   "on this page works with no assistant at all.")
        self.ask_button.clicked.connect(lambda: self.ask.emit(self.field.text().strip()))

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
        self._place_ask()
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

    def _place_ask(self) -> None:
        """Inset Ask inside the field's right margin, left of the clear button."""
        f = self.field
        hint = self.ask_button.sizeHint()
        clear = 34 if f.isClearButtonEnabled() and f.text() else 18
        self.ask_button.setGeometry(f.width() - hint.width() - clear,
                                    (f.height() - hint.height()) // 2,
                                    hint.width(), hint.height())

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
        self._place_ask()
        text = self.field.text().strip()
        if not text:
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
        self._leave_flow()
        self._clear()
        if not text:
            self.feedback.setText(self.PROMPT)
            return
        self._searched = text
        found = suggestions(text, genes=self._genes, index=self._index)
        if not found:
            self.feedback.setText(f"Nothing here matches “{text}”. Ask the assistant, "
                                  "or try a gene id or a product name.")
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
            self.explore(s.key)
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
        self._clear()
        self._top, self._expanded = [], None
        self.back.hide()
        self.field.clear()
        self.feedback.setText(self.PROMPT)
        self._arrange(animate=False)

    # ------------------------------------------------------------ module network
    def explore(self, key: str) -> None:
        """First click on a module: show it with the modules it leads to.

        spaCR: "Click once to open the module network. Click twice to open the module."
        """
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
        reason = MOD.unavailable_reason(MOD.get(key))
        if reason:
            self.flow_hint.setText(f"{MOD.get(key).name} is not available yet: {reason}.")
            return
        self.chosen.emit("module", key)

    def _build_flow(self) -> None:
        """Make the cards for the chain clicked so far and what follows its last step."""
        old = dict(self._flow_cards)
        self._columns = [[self._path[0]]]
        for step in self._path:
            following = list(MOD.downstream(step))[:4]
            if not following:
                break
            self._columns.append(following)
            if len(self._columns) > len(self._path):
                break
        self._flow_cards = {}
        fresh = []
        for c, column in enumerate(self._columns):
            for key in column:
                ident = (c, key)
                card = old.pop(ident, None)
                if card is None:
                    m = MOD.get(key)
                    detail = (f"<h3>{escape(m.name)}</h3><p>{escape(m.note)}.</p>"
                              f"<p><i>{escape(m.band)}</i></p>")
                    nxt = MOD.downstream(key)
                    if nxt:
                        detail += ("<p>Leads to: " + ", ".join(escape(MOD.get(k).name)
                                                               for k in nxt) + ".</p>")
                    card = SuggestionCard(Suggestion("module", key, m.name, m.note, detail), self)
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
        s = self._cards[i].suggestion
        if s.kind == "group":
            self.expand(s)
            return
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
            self._place_ask()
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
