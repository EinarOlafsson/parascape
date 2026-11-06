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
FIELD_MAX_W = 760
CARD_W = 220
CARD_H = 104
#: How long typing must pause before suggestions refresh. spaCR waits 600 ms because
#: it inspects a path on disk; ParaScape's search is in memory, so it can answer
#: sooner without searching on every keystroke of a fast typist.
DEBOUNCE_MS = 180
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


def navigator_qss(palette: dict) -> str:
    """spaCR's Navigator rules, against ParaScape's palette."""
    P = palette
    return f"""
QWidget#NavigatorPage {{ background: transparent; }}
QLineEdit#NavigatorField {{
    background: {P["surface"]};
    color: {P["fg"]};
    border: 1px solid {P["border_soft"]};
    border-radius: {FIELD_RADIUS}px;
    font-family: {TH.FONT_FAMILY};
    font-size: 16px;
    padding: 0px 24px;
    placeholder-text-color: {P["fg_muted"]};
    selection-background-color: {P["accent"]};
}}
QLineEdit#NavigatorField:focus {{ border: 1px solid {P["accent"]}; }}
QPushButton#NavigatorAsk {{
    color: {P["fg_muted"]}; background: transparent; border: none;
    font-size: 15px; padding: 0px 6px;
}}
QPushButton#NavigatorAsk:hover {{ color: {P["accent"]}; }}
QPushButton#SuggestionCard {{
    background: {P["surface"]};
    color: {P["fg"]};
    border: 1px solid {P["border_soft"]};
    border-radius: 16px;
    padding: 10px 12px;
    text-align: left;
}}
QPushButton#SuggestionCard:hover, QPushButton#SuggestionCard:focus {{
    border: 1px solid {P["accent"]};
    background: {P["surface_hi"]};
}}
QPushButton#SuggestionCard[current="true"] {{ border: 1px solid {P["accent"]}; }}
QLabel#CardKind {{ color: {P["accent"]}; background: transparent; font-size: 11px;
                   font-weight: 600; letter-spacing: 0.4px; }}
QLabel#CardTitle {{ color: {P["fg"]}; background: transparent; font-size: 14px; font-weight: 500; }}
QLabel#CardSub {{ color: {P["fg_muted"]}; background: transparent; font-size: 11px; font-weight: 300; }}
QTextBrowser#NavigatorDetails {{
    color: {P["fg"]}; background: {P["surface"]};
    border: 1px solid {P["border_soft"]}; border-radius: 10px;
    padding: 8px 12px; font-size: 14px;
}}
QTextBrowser#NavigatorDetails QScrollBar:vertical {{ background: transparent; width: 8px; }}
QTextBrowser#NavigatorDetails QScrollBar::handle:vertical {{ background: {P["border"]};
    border-radius: 4px; min-height: 24px; }}
QTextBrowser#NavigatorDetails QScrollBar::add-line:vertical,
QTextBrowser#NavigatorDetails QScrollBar::sub-line:vertical {{ height: 0px; }}
QLabel#NavigatorFeedback {{ color: {P["fg_muted"]}; background: transparent; font-size: 13px; }}
QLabel#NavigatorHint {{ color: {P["fg_muted"]}; background: transparent; font-size: 12px; }}
QScrollArea#NavigatorRow, QWidget#NavigatorRowBody {{ background: transparent; border: none; }}
QScrollArea#NavigatorRow QScrollBar:horizontal {{ background: transparent; height: 8px; }}
QScrollArea#NavigatorRow QScrollBar::handle:horizontal {{ background: {P["border"]};
    border-radius: 4px; min-width: 32px; }}
QScrollArea#NavigatorRow QScrollBar::add-line:horizontal,
QScrollArea#NavigatorRow QScrollBar::sub-line:horizontal {{ width: 0px; }}
"""


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


class SuggestionCard(QtWidgets.QPushButton):
    """One answer: what kind of thing it is, its name, one line."""

    def __init__(self, suggestion: Suggestion, parent=None):
        super().__init__(parent)
        self.suggestion = suggestion
        self.setObjectName("SuggestionCard")
        self.setProperty("current", "false")
        self.setFixedSize(CARD_W, CARD_H)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.setAccessibleName(f"{KIND_WORD.get(suggestion.kind, suggestion.kind)}: {suggestion.title}")
        self.setToolTip(suggestion.subtitle or suggestion.title)
        col = QtWidgets.QVBoxLayout(self)
        col.setContentsMargins(12, 8, 12, 8)
        col.setSpacing(2)
        kind = QtWidgets.QLabel(KIND_WORD.get(suggestion.kind, suggestion.kind.title()).upper())
        kind.setObjectName("CardKind")
        title = QtWidgets.QLabel(suggestion.title)
        title.setObjectName("CardTitle")
        # About two lines at 11 px across a 196 px card. A QLabel centres vertically, so
        # text taller than its box was cropped at the TOP, into the title; aligning to
        # the top makes any overflow fall off the bottom instead.
        sub = QtWidgets.QLabel(_clip(_subtitle(suggestion.subtitle), 58))
        sub.setObjectName("CardSub")
        sub.setWordWrap(True)
        sub.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignTop)
        # Two lines at most: a dataset summary can run to hundreds of characters, and
        # the first render showed it pouring out of the bottom of the card.
        sub.setMaximumHeight(2 * QtGui.QFontMetrics(sub.font()).lineSpacing() + 2)
        for label in (kind, title, sub):
            label.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents)
            col.addWidget(label)
        col.addStretch(1)
        # Elide a long title rather than letting it widen the card past its row.
        metrics = QtGui.QFontMetrics(title.font())
        title.setText(metrics.elidedText(suggestion.title, QtCore.Qt.TextElideMode.ElideRight,
                                         CARD_W - 28))

    def set_current(self, on: bool) -> None:
        """Mark this card as the one Return would open. A property, so the stylesheet
        draws the accent rim; repolished, because Qt does not restyle on its own."""
        self.setProperty("current", "true" if on else "false")
        self.style().unpolish(self)
        self.style().polish(self)


class NavigatorPage(QtWidgets.QWidget):
    """The centred field, its suggestion row and the details pane.

    `chosen(kind, key)` fires when a suggestion is opened. `ask(text)` fires when the
    person asks the assistant rather than picking a suggestion.

    :param genes: `(gene_id, product)` pairs, or a callable, for gene suggestions.
    :param index: `help_index` rows, or a callable, for everything else.
    """

    chosen = QtCore.Signal(str, str)
    ask = QtCore.Signal(str)

    def __init__(self, genes=None, index=None, parent=None):
        super().__init__(parent)
        self.setObjectName("NavigatorPage")
        self._genes, self._index = genes, index
        self._cards: list = []
        self._current = -1
        self.setStyleSheet(navigator_qss(TH.palette_for("dark")))

        col = QtWidgets.QVBoxLayout(self)
        col.setContentsMargins(TH.SPACING["xl"], 0, TH.SPACING["xl"], TH.SPACING["xl"])
        col.setSpacing(TH.SPACING["md"])
        col.addStretch(3)

        title = QtWidgets.QLabel("ParaScape")
        title.setObjectName("HeroTitle")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        col.addWidget(title)

        row = QtWidgets.QHBoxLayout()
        row.addStretch(1)
        self.field = QtWidgets.QLineEdit()
        self.field.setObjectName("NavigatorField")
        self.field.setFixedHeight(FIELD_H)
        self.field.setMaximumWidth(FIELD_MAX_W)
        self.field.setMinimumWidth(420)
        self.field.setPlaceholderText("A gene, a question, a module or a dataset")
        self.field.setAccessibleName("Search ParaScape")
        self.field.setClearButtonEnabled(True)
        self.ask_button = QtWidgets.QPushButton("Ask", self.field)
        self.ask_button.setObjectName("NavigatorAsk")
        self.ask_button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.ask_button.setToolTip("Ask the assistant instead of picking a suggestion. "
                                   "Everything else on this page works with no assistant at all.")
        self.ask_button.clicked.connect(lambda: self.ask.emit(self.field.text().strip()))
        self.field.setTextMargins(0, 0, 64, 0)
        self.field.installEventFilter(self)
        row.addWidget(self.field, 4)
        row.addStretch(1)
        col.addLayout(row)

        self.feedback = QtWidgets.QLabel(
            "Try GRA16, glucose, rhoptry, conformal, or the name of a module.")
        self.feedback.setObjectName("NavigatorHint")
        self.feedback.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        col.addWidget(self.feedback)

        self.row = QtWidgets.QScrollArea()
        self.row.setObjectName("NavigatorRow")
        self.row.setWidgetResizable(True)
        self.row.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.row.setFixedHeight(CARD_H + 18)
        self.row.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        body = QtWidgets.QWidget()
        body.setObjectName("NavigatorRowBody")
        self._row_layout = QtWidgets.QHBoxLayout(body)
        self._row_layout.setContentsMargins(0, 0, 0, 0)
        self._row_layout.setSpacing(TH.SPACING["md"])
        # A stretch on BOTH sides centres the cards under the field while they fit; once
        # they overflow, the row's width grows and both stretches collapse, so it simply
        # scrolls. With one stretch, a single result sat at the far left edge of the
        # window, disconnected from the centred field it answered.
        self._row_layout.addStretch(1)
        self._row_layout.addStretch(1)
        self.row.setWidget(body)
        self.row.viewport().setAutoFillBackground(False)
        self.row.hide()
        col.addWidget(self.row)

        self.details = QtWidgets.QTextBrowser()
        self.details.setObjectName("NavigatorDetails")
        self.details.setOpenLinks(False)
        self.details.setMaximumHeight(160)
        self.details.setMaximumWidth(FIELD_MAX_W)
        self.details.hide()
        drow = QtWidgets.QHBoxLayout()
        drow.addStretch(1)
        drow.addWidget(self.details, 4)
        drow.addStretch(1)
        col.addLayout(drow)
        col.addStretch(4)

        self._debounce = QtCore.QTimer(self)
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(DEBOUNCE_MS)
        self._debounce.timeout.connect(self.refresh)
        self.field.textChanged.connect(lambda _t: (self._debounce.start(), self._place_ask()))
        self.field.returnPressed.connect(self.open_current)
        self.field.installEventFilter(self)

    # ------------------------------------------------------------------ search
    def refresh(self) -> None:
        """Rebuild the suggestion row for what is in the field now."""
        text = self.field.text().strip()
        for card in self._cards:
            card.setParent(None)
            card.deleteLater()
        self._cards = []
        self._current = -1
        if not text:
            self.row.hide()
            self.details.hide()
            self.feedback.setText("Try GRA16, glucose, rhoptry, conformal, or the name of a module.")
            return
        found = suggestions(text, genes=self._genes, index=self._index)
        if not found:
            self.row.hide()
            self.details.hide()
            self.feedback.setText(f"Nothing here matches “{text}”. "
                                  "Ask the assistant, or try a gene id or a product name.")
            return
        for i, s in enumerate(found):
            card = SuggestionCard(s)
            card.clicked.connect(lambda _=False, k=i: self._open(k))
            card.installEventFilter(self)
            self._row_layout.insertWidget(self._row_layout.count() - 1, card)
            self._cards.append(card)
        kinds = {}
        for s in found:
            kinds[KIND_WORD.get(s.kind, s.kind)] = kinds.get(KIND_WORD.get(s.kind, s.kind), 0) + 1
        self.feedback.setText(" · ".join(f"{n} {k.lower()}{'s' if n != 1 else ''}"
                                         for k, n in kinds.items()) + " — Return opens the first")
        self.row.show()
        self._set_current(0)

    def _set_current(self, i: int) -> None:
        if not self._cards:
            return
        i = max(0, min(len(self._cards) - 1, i))
        if 0 <= self._current < len(self._cards):
            self._cards[self._current].set_current(False)
        self._current = i
        card = self._cards[i]
        card.set_current(True)
        self.row.ensureWidgetVisible(card, 24, 0)
        self.details.setHtml(card.suggestion.detail)
        self.details.show()

    def current(self):
        """The suggestion Return would open, or None."""
        return self._cards[self._current].suggestion if 0 <= self._current < len(self._cards) else None

    def open_current(self) -> None:
        """Open the highlighted suggestion -- what Return does.

        If typing has not paused long enough for the row to catch up, the row is
        rebuilt first: otherwise a fast typist pressing Return opens the suggestion
        for what they had typed a moment ago, not for what is in the field.
        """
        if self._debounce.isActive():         # Return before the row caught up with typing
            self._debounce.stop()
            self.refresh()
        if self._cards:
            self._open(max(self._current, 0))

    def _open(self, i: int) -> None:
        s = self._cards[i].suggestion
        if s.kind == "module" and MOD.unavailable_reason(MOD.get(s.key)):
            self._set_current(i)
            self.feedback.setText(f"{s.title} is not available yet: "
                                  f"{MOD.unavailable_reason(MOD.get(s.key))}.")
            return
        self.chosen.emit(s.kind, s.key)

    # ----------------------------------------------------------------- keyboard
    def _place_ask(self) -> None:
        """Inset Ask inside the field's right margin, left of the clear button.

        A child of a QLineEdit is drawn at its top-left corner unless something
        positions it -- the first render put Ask above the text. spaCR re-arranges
        its inset buttons on every resize and style change for the same reason.
        """
        f = self.field
        hint = self.ask_button.sizeHint()
        clear = 34 if f.isClearButtonEnabled() and f.text() else 14
        self.ask_button.setGeometry(f.width() - hint.width() - clear,
                                    (f.height() - hint.height()) // 2,
                                    hint.width(), hint.height())

    def eventFilter(self, obj, event):                    # noqa: N802
        """Down from the field enters the cards; arrows move; Escape clears or returns."""
        if obj is self.field and event.type() in (QtCore.QEvent.Type.Resize,
                                                  QtCore.QEvent.Type.Show):
            self._place_ask()
            return False
        if event.type() != QtCore.QEvent.Type.KeyPress:
            return False
        key = event.key()
        K = QtCore.Qt.Key
        if obj is self.field:
            if key == K.Key_Down and self._cards:
                self._cards[max(self._current, 0)].setFocus()
                return True
            if key == K.Key_Escape:
                self.field.clear()
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
            if key in (K.Key_Up, K.Key_Escape):
                self.field.setFocus()
                return True
            if key in (K.Key_Return, K.Key_Enter):
                self._open(self._current)
                return True
        return False
