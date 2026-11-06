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


#: The word shown on a card for each kind, so a gene and a dataset with similar
#: names are never confused.
KIND_WORD = {"module": "Module", "gene": "Gene", "question": "Question",
             "strategy": "Method", "dataset": "Dataset", "panel": "Panel",
             "guide": "Guide", "node": "Card"}


def suggestions(text: str, genes=None, index=None, limit: int = MAX_CARDS) -> list:
    """Everything the field should offer for `text`, best first.

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
    name = (_MODULE_ICON.get(suggestion.key) if suggestion.kind == "module"
            else _KIND_ICON.get(suggestion.kind)) or "settings"
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
            self._open(max(self._current, 0))
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
        for i, s in enumerate(found):
            card = SuggestionCard(s, self)
            card.clicked.connect(lambda _=False, k=i: self._open(k))
            card.hovered.connect(self.explain)
            card.installEventFilter(self)
            card.show()
            self._cards.append(card)
        kinds = {}
        for s in found:
            word = KIND_WORD.get(s.kind, s.kind).lower()
            kinds[word] = kinds.get(word, 0) + 1
        self.feedback.setText(" · ".join(f"{n} {k}{'s' if n != 1 else ''}"
                                         for k, n in kinds.items())
                              + " — Return opens the highlighted card")
        self._scroll.setValue(0)
        self._arrange(animate=True)
        self._set_current(0)

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
        if s.kind == "module" and MOD.unavailable_reason(MOD.get(s.key)):
            self._set_current(i)
            self.feedback.setText(f"{s.title} is not available yet: "
                                  f"{MOD.unavailable_reason(MOD.get(s.key))}.")
            return
        self.chosen.emit(s.kind, s.key)

    # ----------------------------------------------------------------- keyboard
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
