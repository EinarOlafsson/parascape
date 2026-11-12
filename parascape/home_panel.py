"""The home screen: the three things ParaScape can be asked about, as cards.

Instruction 74 phase 1, and the answer to the oldest complaint in instruction 70 -- "as it is now, I
get lost right away." The application already had every panel on this screen; what it did not have
was a place that says what can be asked at all. So this dock is the entry point, and it is a
navigation layer: **every card opens a panel that already exists**. Nothing was removed to make room
for it, which is why no workflow breaks (phase 5 retires a dock only once each has a node here).

The screen follows the data model rather than a guess at the user's mood. The map is a matrix, so a
query is one row (`Gene`), a block of rows (`List`) or one column (`Label`), and there is no fourth
section: provenance is an attribute of every value rather than a place, and a UMAP is a view of the
rows rather than a kind of query. `parascape.nodes` holds that model, the availability probes and
the resolver; this module draws it and nothing else, which is what lets the whole screen be built
and clicked in a test with no window at all.

Three things here are deliberate and not decoration:

* **Trust is integrated, not visited.** There is no "how reliable is this" section, because a
  separate room for reliability makes the default experience untrustworthy. A card whose
  destination shows a claim says on its face that the claim arrives with its probability, its
  uncertainty and its source tier; a question card carries its own scoreboard inline -- the
  cross-validated estimate beside the frozen-test result -- its coverage, and whether it may fill
  anything at all. The label catalogue's header carries effective independent axes, where it
  qualifies every row beneath it.
* **An unavailable card is shown, not hidden.** Greyed, with the reason on it ("needs a background frequency
  background", "needs the mechanism view"), and clicking it says the reason rather than failing. A
  card left off the screen is a gap nobody notices; a greyed one with a reason is a decision.
* **One search bar is the navigator.** The field at the top returns genes, questions, strategies,
  datasets, panels and guide pages together, ranked and interleaved by `help_index.search` -- the
  same code behind the field beside the Help menu, extended with node kinds rather than replaced.
  Every result is a card, and clicking a gene lands on *that gene's* card with the gene selected.
  Type, Down, Return, as the Help field already behaves.
* **The layout is spaCR's, not a new one.** `theme.RADIUS` and `theme.SPACING` are spaCR's scales
  already; the cards follow its `QPushButton#AppTile` rule -- a launcher sitting on the home pane
  earns an edge only when it is interactive, because a resting rim on every tile made the pane read
  as a ruled table -- and the type sizes follow its `TYPOGRAPHY` (a bold header for the title, body
  for the line under it, small for the honesty lines).
* **A dropped file is a gene list.** The window accepts a drop anywhere and routes it to `List ▸
  Drop`, whose view shows the file's name, sheet, sniffed separator and identifier column -- all
  changeable -- the resolution report with the unresolved identifiers **listed rather than
  discarded**, and the extra per-gene columns with the treatment their *numbers* suggest. The
  reading is `parascape.importer` through `parascape.gene_drop`; reading is slow enough for a
  workbook that it happens on a `QThread`.
"""
from __future__ import annotations

import os

from PySide6 import QtCore, QtGui, QtWidgets

from . import gene_drop as GD
from . import no_cross
from . import nodes as ND
from . import theme as TH

#: Every control's tooltip, in one place: the panel builds controls in several methods, and a
#: release once failed its checks for a single missing one.
TIPS = {
    "search": "One field for everything: a gene by accession or product, a question, a strategy, a "
              "dataset, a panel, a page of the guide. Results are ranked and interleaved by the "
              "same search that sits beside the Help menu, and every one of them is a card -- "
              "clicking a gene lands on that gene's card with the gene selected. Down moves into "
              "the results, Return opens one; empty the box to see the whole home screen again.",
    "found": "What the search found beyond the home cards: genes, strategies, datasets, panels and "
             "guide pages, best match first.",
    "card": "Open this. The line under the title says exactly what that does.",
    "unavailable_card": "Not available yet. The card says why; clicking it repeats the reason and "
                        "changes nothing.",
    "refresh": "Read the harness results and the released information map again, and redraw the "
               "cards. Use it after a scoring run has finished.",
    "back": "Back to the four query units.",
    "catalogue": "Every registered question as a row: how many labels it has, the cross-validated "
                 "estimate and the frozen-test result side by side, how many genes it can answer "
                 "at the 90 % bar, and its state. Click a row to open that question.",
    "catalogue_header": "What qualifies the whole list: how many questions are registered, how "
                        "many have labels, and how many INDEPENDENT axes those questions really "
                        "hold. A list of questions is not a list of facts.",
    "question": "This question's scoreboard, both numbers side by side, its label coverage, and "
                "the evidence families it may and may not see.",
    "choose_file": "Choose a file of genes: csv, tsv, txt, tab, xlsx, xls or parquet. You can also "
                   "drop one anywhere on the window.",
    "sheet": "Which sheet of the workbook holds the table. Published supplements routinely put it "
             "behind a legend sheet, so this is a choice rather than a guess.",
    "column": "Which column holds the identifiers. The suggestion is the column whose values look "
              "most like gene accessions; change it if the suggestion is wrong.",
    "report": "Every identifier in the file and what became of it: resolved, ambiguous, or "
              "unresolved with the reason. Nothing is silently dropped -- an unresolved identifier "
              "is listed here rather than being missing from a shorter list.",
    "extras": "The other per-gene columns found in the file, each with its range and the "
              "quantification its NUMBERS suggest, for you to confirm. The range decides, never "
              "the filename: a column reaching 16,520 has not been logged and one stopping at 9.7 "
              "has, whatever the file is called.",
    "list_action": "Do this with the list that is loaded.",
    "export": "Write the list with every identifier's fate and the facts looked up for it, so the "
              "unresolved ones travel with the file rather than vanishing from it.",
}

#: Card geometry. Wide enough for a one-line summary to wrap in two, on spaCR's spacing grid.
CARD_WIDTH = 248
CARD_HEIGHT = 128
#: The radius of a home-screen launcher, from the shared scale: spaCR's `AppTile`.
CARD_RADIUS = TH.RADIUS["tile"]
#: spaCR's spacing scale, so the sections sit on the same grid as a spaCR pane.
GAP = TH.SPACING
#: Cards per row before wrapping.
CARDS_PER_ROW = 3
#: Rows the drop report shows before scrolling; the rest are there, in the table.
REPORT_ROWS = 400

#: The three pages: the cards, the label catalogue, the dropped list.
PAGE_HOME, PAGE_CATALOGUE, PAGE_DROP = 0, 1, 2


def _palette() -> dict:
    """The colours of the theme in force, for the painted cards."""
    app = QtWidgets.QApplication.instance()
    name = str(app.property("parascapeTheme") or "dark") if app else "dark"
    try:
        return TH.palette_for(name)
    except Exception:                                  # an unknown theme must not stop the screen
        return TH.DARK


class _Catalogue(QtCore.QThread):
    """Rebuild the node catalogue off the main thread: it reads the results and the map release."""

    ready = QtCore.Signal(object)
    failed = QtCore.Signal(str)

    def __init__(self, gene_list_loaded: bool, parent=None):
        super().__init__(parent)
        self.gene_list_loaded = gene_list_loaded

    def run(self):
        """Read the scoreboards and the presence probes, and hand back the catalogue."""
        try:
            self.ready.emit(ND.catalogue(gene_list_loaded=self.gene_list_loaded))
        except Exception as error:                     # a bad results file must not kill the home dock
            self.failed.emit(f"{type(error).__name__}: {error}")


class _DropReader(QtCore.QThread):
    """Read a dropped file off the main thread: a workbook takes seconds, and the window must live."""

    ready = QtCore.Signal(object)
    failed = QtCore.Signal(str)

    def __init__(self, path: str, sheet=None, known_ids=None, gene_column: str = "", parent=None):
        super().__init__(parent)
        self.path, self.sheet = path, sheet
        self.known_ids, self.gene_column = known_ids, gene_column

    def run(self):
        """Read the file and hand back a `gene_drop.Drop`, or the reason it could not be read."""
        try:
            self.ready.emit(GD.read_drop(self.path, sheet=self.sheet, known_ids=self.known_ids,
                                         gene_column=self.gene_column))
        except GD.UnsupportedFormat as error:
            self.failed.emit(str(error))
        except (OSError, ValueError, KeyError, TypeError, ImportError) as error:
            self.failed.emit(f"{os.path.basename(self.path)} could not be read: {error}")


class NodeCard(QtWidgets.QAbstractButton):
    """One node as a card: its one-word title, what clicking it does, and its honesty lines.

    Painted rather than styled, for `guided_panel.Choice`'s reason: the card reads the same on all
    four themes, and the three extra lines -- the claim line, the scoreboard, the reason a card is
    unavailable -- have to be legible without a stylesheet rule per theme. An unavailable card is
    drawn dimmed with a dashed border, so it cannot be mistaken for one that works.
    """

    def __init__(self, node: ND.Node, parent=None):
        super().__init__(parent)
        self.node = node
        self.setObjectName(f"NodeCard_{node.id}")
        self.setText(node.title)
        self.setAccessibleName(f"{node.title}: {node.summary}")
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.setMinimumSize(CARD_WIDTH, CARD_HEIGHT)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                           QtWidgets.QSizePolicy.Policy.Fixed)
        self.setToolTip(TH.tip(self.tip_text()))

    def tip_text(self) -> str:
        """What the tooltip says: the summary, then whatever else the card is obliged to state."""
        parts = [f"{self.node.title} -- {self.node.summary}."]
        if not self.node.available:
            parts.append(f"Not available: {self.node.unavailable_reason}.")
            parts.append(TIPS["unavailable_card"])
            return " ".join(parts)
        if self.node.shows_claim:
            parts.append("Every claim it shows arrives with its probability, its uncertainty and "
                         "the tier of the source behind it, in the same view.")
        if self.node.scoreboard is not None:
            parts.append(self.node.scoreboard.inline() + ".")
            parts.append(self.node.scoreboard.coverage() + ".")
        parts.append(TIPS["card"])
        return " ".join(parts)

    def lines(self) -> list:
        """The lines under the title, in the order they are drawn."""
        out = [self.node.summary]
        if not self.node.available:
            out.append("unavailable: " + self.node.unavailable_reason)
            return out
        if self.node.scoreboard is not None:
            out.append(self.node.scoreboard.inline())
        if self.node.shows_claim:
            out.append(self.node.claim_line())
        return out

    def sizeHint(self):
        """A fixed card; the grid gives it the width it has."""
        return QtCore.QSize(CARD_WIDTH, CARD_HEIGHT)

    def paintEvent(self, event):
        p = _palette()
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        rect = QtCore.QRectF(self.rect()).adjusted(1.5, 1.5, -1.5, -1.5)
        live = self.node.available
        # spaCR's rule for a launcher on the home pane: no resting rim, because a rim on every tile
        # made the pane read as a ruled table; only an interactive state earns an edge. An
        # unavailable tile is the exception -- it carries a dashed one, so it cannot be mistaken
        # for a tile that works.
        active = (self.hasFocus() or self.underMouse()) and live
        if active:
            painter.setPen(QtGui.QPen(QtGui.QColor(p["accent"]), 1.6))
        elif not live:
            pen = QtGui.QPen(QtGui.QColor(p["border"]), 1.2)
            pen.setStyle(QtCore.Qt.PenStyle.DashLine)
            painter.setPen(pen)
        else:
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QColor(p["accent_soft"] if self.isDown() and live
                                      else p["surface_hi"] if active
                                      else p["surface_alt"] if live
                                      else p["surface"]))
        painter.drawRoundedRect(rect, CARD_RADIUS, CARD_RADIUS)
        title = QtGui.QFont(self.font())
        title.setBold(True)
        if title.pointSizeF() > 0:
            title.setPointSizeF(title.pointSizeF() * 1.18)
        painter.setFont(title)
        painter.setPen(QtGui.QColor(p["fg"] if live else p["fg_muted"]))
        box = rect.adjusted(GAP["md"], GAP["sm"], -GAP["md"], -GAP["sm"])
        height = QtGui.QFontMetricsF(title).height()
        painter.drawText(QtCore.QRectF(box.left(), box.top(), box.width(), height),
                         int(QtCore.Qt.AlignmentFlag.AlignLeft
                             | QtCore.Qt.AlignmentFlag.AlignVCenter), self.node.title)
        small = QtGui.QFont(self.font())
        if small.pointSizeF() > 0:
            small.setPointSizeF(max(small.pointSizeF() * 0.86, 6.0))
        painter.setFont(small)
        top = box.top() + height + 2
        for i, line in enumerate(self.lines()):
            painter.setPen(QtGui.QColor(p["fg_muted"] if i == 0 and live else p["fg_dim"]))
            area = QtCore.QRectF(box.left(), top, box.width(), box.bottom() - top)
            if area.height() <= 2:
                break
            flags = int(QtCore.Qt.TextFlag.TextWordWrap
                        | QtCore.Qt.AlignmentFlag.AlignTop | QtCore.Qt.AlignmentFlag.AlignLeft)
            used = painter.boundingRect(area, flags, line)
            painter.drawText(area, flags, line)
            top = used.bottom() + 3
        painter.end()

    def keyPressEvent(self, event):
        K = QtCore.Qt.Key
        if event.key() in (K.Key_Return, K.Key_Enter):
            self.click()
            event.accept()
            return
        super().keyPressEvent(event)


class QuestionCatalogue(QtWidgets.QWidget):
    """Every registered question as a row, with both numbers and its state, under one header.

    The header is the point of the page as much as the rows are: it carries effective independent
    axes, so a list of questions cannot read as a list of independent facts. A question with too
    few labels is a row in a state, not a page somewhere else.
    """

    #: Emitted with the level-3 question node of the row that was chosen.
    chosen = QtCore.Signal(object)

    COLUMNS = ("label", "question", "answer type", "known answers", "right when rotated",
               "right on hidden genes",
               "genes at the bar", "state")

    def __init__(self, question_nodes, parent=None):
        super().__init__(parent)
        self.setObjectName("QuestionCatalogue")
        layout = QtWidgets.QVBoxLayout(self)
        self.header = QtWidgets.QLabel(ND.catalogue_header())
        self.header.setObjectName("QuestionCatalogueHeader")
        self.header.setWordWrap(True)
        self.header.setToolTip(TH.tip(TIPS["catalogue_header"]))
        layout.addWidget(self.header)
        self.table = QtWidgets.QTableWidget(0, len(self.COLUMNS))
        self.table.setObjectName("QuestionCatalogueTable")
        self.table.setHorizontalHeaderLabels(list(self.COLUMNS))
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setToolTip(TH.tip(TIPS["catalogue"]))
        self.table.cellClicked.connect(self._chosen)
        layout.addWidget(self.table, 1)
        self.detail = QtWidgets.QTextBrowser()
        self.detail.setObjectName("QuestionCatalogueDetail")
        self.detail.setToolTip(TH.tip(TIPS["question"]))
        self.detail.setMinimumHeight(140)
        layout.addWidget(self.detail)
        self.nodes = []
        self.show_questions(question_nodes)

    def show_questions(self, question_nodes) -> None:
        """Fill the table from the level-3 question nodes."""
        from . import question_registry as QR
        self.nodes = list(question_nodes)
        self.header.setText(ND.catalogue_header())
        self.table.setRowCount(len(self.nodes))
        for row, node in enumerate(self.nodes):
            board = node.scoreboard
            state = ("cannot be answered yet" if not node.available else
                     "answerable" if board is not None and board.scored else "not yet scored")
            if node.available and board is not None and board.scored \
                    and board.frozen_test_lower is None:
                state = "ranks genes only: not enough known answers to grade"
            question = QR.QUESTIONS.get(node.question_id)
            values = (node.title,
                      question.text if question is not None else node.summary,
                      question.answer_type if question is not None else "",
                      f"{board.n_labels:,}" if board and board.n_labels else "not counted",
                      "" if not board or board.cross_validated_lower is None
                      else f"≥ {board.cross_validated_lower:.3f}",
                      "" if not board or board.frozen_test_lower is None
                      else f"≥ {board.frozen_test_lower:.3f}",
                      f"{board.genes_at_bar:,}" if board and board.may_fill else "0",
                      state)
            for column, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(str(value))
                item.setToolTip(f"{node.summary}\n{'' if not board else board.coverage()}")
                self.table.setItem(row, column, item)
        self.table.resizeColumnsToContents()
        if self.nodes:
            self.show_detail(0)

    def show_detail(self, row: int) -> None:
        """Write one question's scoreboard, coverage and evidence rules under the table."""
        from html import escape
        if not 0 <= row < len(self.nodes):
            return
        node = self.nodes[row]
        board = node.scoreboard
        from . import question_registry as QR
        try:
            question = QR.get(node.question_id)
        except KeyError:
            question = None
        body = [f"<h3>{escape(node.title)}</h3>",
                f"<p>{escape(question.text if question else node.summary)}</p>"]
        if board is not None:
            body.append(f"<p><b>Scoreboard</b><br>{escape(board.inline())}</p>")
            body.append(f"<p><b>Coverage</b><br>{escape(board.coverage())}</p>")
        if not node.available:
            body.append(f"<p><b>State</b><br>{escape(node.unavailable_reason)}</p>")
        if question is not None:
            body.append(f"<p><b>Labels</b><br>{escape(question.label_source)} "
                        f"({escape(question.label_tier)})</p>")
            body.append("<p><b>Forbidden evidence</b><br>"
                        + escape(", ".join(question.evidence_forbidden) or "none") + "</p>")
        self.detail.setHtml("".join(body))

    def _chosen(self, row: int, _column: int) -> None:
        self.show_detail(row)
        if 0 <= row < len(self.nodes):
            self.chosen.emit(self.nodes[row])


class DropView(QtWidgets.QWidget):
    """A dropped file of genes: what was read, what resolved, and what else was in it.

    The file's own decisions -- sheet, identifier column -- are controls rather than a message,
    because a dropped file whose decisions were made silently is a list nobody can defend. The
    resolution report keeps the unresolved identifiers on screen: a list that quietly shrinks is how
    a published screen lost every one of its rows once already.
    """

    #: Emitted with a node id when one of the list actions is asked for.
    action = QtCore.Signal(str)
    #: Emitted with a path to read, when the user chooses or changes something.
    reread = QtCore.Signal(str, object, str)

    def __init__(self, catalogue=None, parent=None):
        super().__init__(parent)
        self.setObjectName("DropView")
        self.drop = None
        self.nodes = tuple(catalogue) if catalogue is not None else ND.catalogue()
        layout = QtWidgets.QVBoxLayout(self)
        self.headline = QtWidgets.QLabel(
            "Drop a file of genes anywhere on this window, or choose one. "
            f"Accepted: {', '.join(GD.READABLE)}.")
        self.headline.setWordWrap(True)
        layout.addWidget(self.headline)
        row = QtWidgets.QHBoxLayout()
        self.choose = QtWidgets.QPushButton("Choose a file…")
        self.choose.setToolTip(TH.tip(TIPS["choose_file"]))
        self.choose.clicked.connect(self._choose)
        self.sheet = QtWidgets.QComboBox()
        self.sheet.setToolTip(TH.tip(TIPS["sheet"]))
        self.sheet.setEnabled(False)
        self.sheet.currentTextChanged.connect(self._sheet_changed)
        self.column = QtWidgets.QComboBox()
        self.column.setToolTip(TH.tip(TIPS["column"]))
        self.column.setEnabled(False)
        self.column.currentTextChanged.connect(self._column_changed)
        for widget, label in ((self.choose, ""), (self.sheet, "sheet"), (self.column, "identifier")):
            if label:
                tag = QtWidgets.QLabel(label)
                tag.setToolTip(widget.toolTip())
                row.addWidget(tag)
            row.addWidget(widget)
        row.addStretch(1)
        layout.addLayout(row)
        self.summary = QtWidgets.QLabel("No list loaded.")
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.report = QtWidgets.QTableWidget(0, 4)
        self.report.setObjectName("DropReport")
        self.report.setHorizontalHeaderLabels(["identifier in the file", "gene", "state", "why"])
        self.report.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.report.setToolTip(TH.tip(TIPS["report"]))
        layout.addWidget(self.report, 1)
        self.extras = QtWidgets.QListWidget()
        self.extras.setObjectName("DropExtras")
        self.extras.setToolTip(TH.tip(TIPS["extras"]))
        self.extras.setMaximumHeight(120)
        layout.addWidget(self.extras)
        # One button per `List` node except Drop itself, built from the catalogue rather than
        # listed here, so a node added to the category gets its button without an edit.
        self.buttons = {}
        actions = QtWidgets.QHBoxLayout()
        for node in [n for n in self.nodes
                     if n.category == "list" and n.level == 2 and n.id != "list_drop"]:
            button = QtWidgets.QPushButton(node.title)
            button.setObjectName(f"DropAction_{node.id}")
            button.setToolTip(TH.tip(
                f"Not available: {node.unavailable_reason}" if not node.available
                else TIPS["export"] if node.id == "list_export"
                else f"{TIPS['list_action']} {node.title}: {node.summary}"))
            button.setEnabled(False)
            button.clicked.connect(lambda _c=False, n=node.id: self.action.emit(n))
            actions.addWidget(button)
            self.buttons[node.id] = button
        actions.addStretch(1)
        layout.addLayout(actions)

    def _choose(self) -> None:
        pattern = " ".join(f"*{ext}" for ext in GD.READABLE)
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open a file of genes", "",
                                                        f"Gene lists ({pattern})")
        if path:
            self.reread.emit(path, None, "")

    def _sheet_changed(self, sheet: str) -> None:
        if self.drop is not None and sheet and sheet != self.drop.sheet:
            self.reread.emit(self.drop.path, sheet, "")

    def _column_changed(self, column: str) -> None:
        if self.drop is not None and column and column != self.drop.gene_column:
            self.reread.emit(self.drop.path, self.drop.sheet or None, column)

    def show_drop(self, drop) -> None:
        """Show a file that has been read: its facts, its report and its extra columns."""
        self.drop = drop
        self.headline.setText(drop.headline())
        for box, values, current in ((self.sheet, drop.sheets, drop.sheet),
                                     (self.column, [str(c) for c in drop.table.columns],
                                      drop.gene_column)):
            box.blockSignals(True)
            box.clear()
            box.addItems([str(v) for v in values])
            if current:
                box.setCurrentText(str(current))
            box.setEnabled(bool(values))
            box.blockSignals(False)
        self.summary.setText(drop.report.summary())
        rows = drop.report.rows[:REPORT_ROWS]
        self.report.setRowCount(len(rows))
        for i, resolution in enumerate(rows):
            for j, value in enumerate((resolution.raw, resolution.gene_id, resolution.state,
                                       resolution.reason
                                       or (", ".join(resolution.candidates) if resolution.candidates
                                           else ""))):
                item = QtWidgets.QTableWidgetItem(str(value))
                item.setToolTip(TH.tip(TIPS["report"]))
                self.report.setItem(i, j, item)
        self.report.resizeColumnsToContents()
        self.extras.clear()
        for extra in drop.extras:
            item = QtWidgets.QListWidgetItem(extra.line())
            item.setToolTip(TH.tip(TIPS["extras"]))
            self.extras.addItem(item)
        if not drop.extras:
            item = QtWidgets.QListWidgetItem("no per-gene data columns in this file: a list of "
                                             "identifiers only")
            item.setToolTip(TH.tip(TIPS["extras"]))
            self.extras.addItem(item)
        for node_id, button in self.buttons.items():
            node = ND.get(node_id, self.nodes)
            # A list is loaded, so what still greys a button is the node itself: no background frequency
            # background, no set-enrichment view. The reason stays on the button.
            waiting_only_for_a_list = node.unavailable_reason == ND.NEEDS_LIST
            live = bool(drop.genes) and (node.available or waiting_only_for_a_list)
            button.setEnabled(live)
            button.setToolTip(TH.tip(
                f"{TIPS['list_action']} {node.title}: {node.summary}" if live
                else f"Not available: {node.unavailable_reason}"))

    def show_problem(self, message: str) -> None:
        """Say why a file could not be read, with what is accepted, and keep the view usable."""
        self.summary.setText(message)
        self.headline.setText(f"{message} Accepted: {', '.join(GD.READABLE)}.")


class HomePanel(QtWidgets.QWidget):
    """The three query units as sections of cards, with search, the catalogue and the drop view.

    Three pages in a stack, which is the three-deep nesting and no deeper: the cards, then a level-2
    destination drawn here (the label catalogue, the dropped list), then the panel a node opens.
    """

    #: Emitted with a `nodes.Node` when a card is chosen. `install` turns it into an action.
    node_chosen = QtCore.Signal(object)
    #: Emitted with a line for the status bar.
    status = QtCore.Signal(str)

    def __init__(self, catalogue=None, known_ids=None, gene_lookup=None, index=None, parent=None,
                 cards: bool = True):
        super().__init__(parent)
        # `cards=False` is how the window installs it (instruction 85): the landing is the only
        # home, so the node-card board and its search are not drawn, and the panel opens on the
        # gene list. The catalogue and the drop page are unchanged; Back from the catalogue
        # returns to the gene list instead of to a board that is not there.
        self.with_cards = cards
        self._home_page = PAGE_HOME if cards else PAGE_DROP
        self._backs: list = []
        self.setObjectName("HomePanel")
        self.setAcceptDrops(True)
        self.known_ids = known_ids
        # The navigator's two extra sources: the window's gene table and the help index. Passed in
        # rather than read here, so the panel is still buildable without a window.
        self.gene_lookup = gene_lookup
        self.index = index
        self.found: list = []
        self.nodes = tuple(catalogue) if catalogue is not None else ND.catalogue()
        self.cards: list = []
        self._reader = None
        self._loader = None
        self.drop = None

        outer = QtWidgets.QVBoxLayout(self)
        top = QtWidgets.QHBoxLayout()
        self.search = QtWidgets.QLineEdit()
        self.search.setObjectName("HomeSearch")
        no_cross.clear_on_escape(self.search)     # no small cross; Escape empties the field
        self.search.setPlaceholderText("What do you want to know?")
        self.search.setToolTip(TH.tip(TIPS["search"]))
        self.search.textChanged.connect(self.filter_cards)
        self.search.installEventFilter(self)
        self.refresh_button = QtWidgets.QPushButton("Refresh")
        self.refresh_button.setToolTip(TH.tip(TIPS["refresh"]))
        self.refresh_button.clicked.connect(self.reload)
        top.addWidget(self.search, 1)
        top.addWidget(self.refresh_button)
        outer.addLayout(top)

        self.pages = QtWidgets.QStackedWidget()
        outer.addWidget(self.pages, 1)

        self.board = QtWidgets.QWidget()
        self.board_layout = QtWidgets.QVBoxLayout(self.board)
        self.board_layout.setContentsMargins(0, 0, 0, 0)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setWidget(self.board)
        self.pages.insertWidget(PAGE_HOME, scroll)

        self.catalogue_view = QuestionCatalogue(ND.children("label_catalogue", self.nodes))
        self.catalogue_view.chosen.connect(self._question_chosen)
        self.pages.insertWidget(PAGE_CATALOGUE, self._with_back(self.catalogue_view))

        self.drop_view = DropView(self.nodes)
        self.drop_view.action.connect(self._list_action)
        self.drop_view.reread.connect(self.load_file)
        self.pages.insertWidget(PAGE_DROP, self._with_back(self.drop_view))

        self.note = QtWidgets.QLabel("")
        self.note.setObjectName("HomeNote")
        self.note.setWordWrap(True)
        outer.addWidget(self.note)
        self.build_cards()
        if not cards:
            self.search.hide()
            self.refresh_button.hide()
            self._backs[PAGE_DROP - 1].hide()          # the gene list is this panel's first page
            self.show_page(PAGE_DROP)

    # ------------------------------------------------------------------ building
    def _with_back(self, widget: QtWidgets.QWidget) -> QtWidgets.QWidget:
        """Wrap a page with a Back button, so no destination here is a one-way door."""
        holder = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(holder)
        layout.setContentsMargins(0, 0, 0, 0)
        back = QtWidgets.QPushButton("← Back")
        back.setObjectName("HomeBack")
        back.setToolTip(TH.tip(TIPS["back"]))
        back.clicked.connect(lambda: self.show_page(self._home_page))
        self._backs.append(back)
        bar = QtWidgets.QHBoxLayout()
        bar.addWidget(back)
        bar.addStretch(1)
        layout.addLayout(bar)
        layout.addWidget(widget, 1)
        return holder

    def build_cards(self) -> None:
        """Draw a section per category and a card per level-2 node, advanced ones last."""
        while self.board_layout.count():
            item = self.board_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()
            elif item.layout() is not None:
                self._drop_layout(item.layout())
        self.cards = []
        if self.found:
            self.board_layout.addLayout(self._section(
                "Found", f"what the search found beyond the home cards ({len(self.found)})",
                self.found, tip=TIPS["found"]))
        for category_id, cards in ND.by_category(self.nodes).items():
            category = ND.category(category_id)
            self.board_layout.addLayout(self._section(
                category.title, category.purpose, cards,
                tip=f"{category.title}: {category.purpose}. Every card below opens a panel that "
                    "already exists."))
        self.board_layout.addStretch(1)
        self.catalogue_view.show_questions(ND.children("label_catalogue", self.nodes))

    def _section(self, title: str, purpose: str, nodes, tip: str = "") -> QtWidgets.QVBoxLayout:
        """One titled section and its row of cards, laid out on spaCR's spacing grid."""
        column = QtWidgets.QVBoxLayout()
        column.setSpacing(GAP["sm"])
        column.setContentsMargins(0, GAP["sm"], 0, GAP["sm"])
        label = QtWidgets.QLabel(f"{title} — {purpose}")
        label.setObjectName("HomeSection")
        label.setWordWrap(True)
        font = QtGui.QFont(label.font())
        font.setBold(True)
        label.setFont(font)
        label.setToolTip(TH.tip(tip or f"{title}: {purpose}"))
        column.addWidget(label)
        grid = QtWidgets.QGridLayout()
        grid.setSpacing(GAP["sm"])
        for position, node in enumerate(nodes):
            card = NodeCard(node, self.board)
            card.clicked.connect(lambda _c=False, n=node: self.choose(n))
            grid.addWidget(card, position // CARDS_PER_ROW, position % CARDS_PER_ROW)
            self.cards.append(card)
        column.addLayout(grid)
        return column

    def _drop_layout(self, layout: QtWidgets.QLayout) -> None:
        """Delete a nested layout and everything in it, so a rebuild leaves nothing behind."""
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()
            elif item.layout() is not None:
                self._drop_layout(item.layout())

    def reload(self) -> None:
        """Read the results and the map release again on a thread, then redraw the cards."""
        if self._loader is not None:
            return
        self.note.setText("reading the scoreboards and the released map again…")
        loader = _Catalogue(self.drop is not None and bool(self.drop.genes), self)
        self._loader = loader
        loader.ready.connect(self._reloaded)
        loader.failed.connect(lambda message: self.note.setText(message))
        loader.finished.connect(self._loader_done)
        loader.start()

    def _reloaded(self, catalogue) -> None:
        self.nodes = tuple(catalogue)
        self.drop_view.nodes = self.nodes
        self.build_cards()
        self.filter_cards(self.search.text())
        self.note.setText("")

    def _loader_done(self) -> None:
        loader, self._loader = self._loader, None
        if loader is not None:
            loader.deleteLater()

    # ------------------------------------------------------------------ search
    def filter_cards(self, query: str = "") -> list:
        """The navigator: rank everything for `query`, draw what it found, and return the nodes.

        One field answers with genes, questions, strategies, datasets, panels and guide pages at
        once, interleaved by `help_index`'s ranking -- the same code behind the field beside the
        Help menu -- rather than with three lists stapled together. The home cards it matches keep
        their place in their own category, because the sections are the thing a user is learning;
        everything else it found is drawn as a `Found` section above them.
        """
        text = str(query or "").strip()
        ranked = ND.search(text, self.nodes, genes=self.gene_lookup, index=self.index)
        found = [n for n in ranked if n.from_search]
        if [n.id for n in found] != [n.id for n in self.found]:
            self.found = found
            self.build_cards()
        keep = {node.id for node in ranked}
        order = {node.id: i for i, node in enumerate(ranked)}
        for card in self.cards:
            card.setVisible(not text or card.node.id in keep)
        # `isHidden`, not `isVisibleTo`: the board lives in a page of the stack, so every card is
        # invisible-to-the-panel whenever the catalogue or the drop page is up, and a filter that
        # read it that way reported no results at all while the drop page was showing.
        shown = [card.node for card in self.cards if not card.isHidden()]
        shown.sort(key=lambda n: order.get(n.id, len(order)))
        if text and not shown:
            questions = [n for n in ranked if n.level == 3]
            self.note.setText(
                f"Nothing called “{text}”."
                + (f" {len(questions)} question(s) match it: open Label ▸ Catalogue." if questions
                   else " Try a gene accession, a product name, a question, a strategy or a "
                        "dataset."))
        elif text:
            self.note.setText(f"{len(shown)} result(s) for “{text}”"
                              + (f", {len(found)} of them found by search." if found else "."))
        else:
            self.note.setText("")
        return shown

    def eventFilter(self, obj, event):
        """Down in the search box moves to the first result, as the Help field behaves."""
        if obj is self.search and event.type() == QtCore.QEvent.Type.KeyPress \
                and event.key() in (QtCore.Qt.Key.Key_Down, QtCore.Qt.Key.Key_Tab):
            cards = self.shown_cards_widgets()
            if cards:
                cards[0].setFocus(QtCore.Qt.FocusReason.TabFocusReason)
                return True
        return super().eventFilter(obj, event)

    def shown_cards_widgets(self) -> list:
        """The cards the filter is showing, as widgets, in drawing order."""
        return [card for card in self.cards if not card.isHidden()]

    def shown_cards(self) -> list:
        """The cards the filter is showing, in drawing order."""
        return self.shown_cards_widgets()

    # ------------------------------------------------------------------ choosing
    def show_page(self, page: int) -> None:
        """Show one of the three pages."""
        self.pages.setCurrentIndex(page)

    def choose(self, node: ND.Node) -> str:
        """Act on a node: draw it here when the home screen owns it, else ask for its panel.

        Returns the line that went to the status bar. Never raises: an unavailable node says why,
        which is the whole reason a node carries a reason.
        """
        action = ND.resolve(node)
        if action.kind == "unavailable":
            message = f"{node.title}: {node.unavailable_reason}"
            self.note.setText(message)
            self.status.emit(message)
            return message
        if action.kind == "catalogue":
            self.show_page(PAGE_CATALOGUE)
            if action.target == "fillable":
                self.catalogue_view.show_questions(
                    [n for n in ND.children("label_catalogue", self.nodes)
                     if n.scoreboard is not None and n.scoreboard.may_fill])
            else:
                self.catalogue_view.show_questions(ND.children("label_catalogue", self.nodes))
            message = f"{node.title}: {ND.catalogue_header()}"
            self.status.emit(message)
            return message
        if action.kind == "question":
            self.show_page(PAGE_CATALOGUE)
            for row, candidate in enumerate(self.catalogue_view.nodes):
                if candidate.question_id == action.target:
                    self.catalogue_view.table.selectRow(row)
                    self.catalogue_view.show_detail(row)
                    break
            message = (f"{node.title}: {node.scoreboard.inline()}" if node.scoreboard is not None
                       else f"{node.title}")
            self.status.emit(message)
            return message
        if action.kind in ("gene", "entry"):
            self.node_chosen.emit(node)
            self.status.emit(action.describe)
            return action.describe
        if action.kind == "drop":
            self.show_page(PAGE_DROP)
            message = f"{node.title}: {node.summary}"
            self.note.setText(message)
            self.status.emit(message)
            return message
        self.node_chosen.emit(node)
        self.status.emit(action.describe)
        return action.describe

    def _question_chosen(self, node: ND.Node) -> None:
        self.status.emit(f"{node.title}: {node.scoreboard.inline()}"
                         if node.scoreboard is not None else node.title)

    def _list_action(self, node_id: str) -> None:
        self.choose(ND.get(node_id, self.nodes))

    # ------------------------------------------------------------------ a dropped file
    def load_file(self, path: str, sheet=None, gene_column: str = "") -> bool:
        """Read a file of genes on a thread and show it on the drop page. False if it was refused."""
        self.show_page(PAGE_DROP)
        if not GD.accepts(path):
            problem = str(GD.UnsupportedFormat(path))
            self.drop_view.show_problem(problem)
            self.note.setText(problem)
            self.status.emit(problem)
            return False
        if self._reader is not None:
            self.status.emit("still reading the previous file")
            return False
        self.note.setText(f"reading {os.path.basename(path)}…")
        reader = _DropReader(path, sheet, self.known_ids, gene_column, self)
        self._reader = reader
        reader.ready.connect(self._file_read)
        reader.failed.connect(self._file_failed)
        reader.finished.connect(self._reader_done)
        reader.start()
        return True

    def wait_for_file(self, ms: int = 30000) -> bool:
        """Block until a file being read is finished. For a caller that cannot wait on a signal."""
        reader = self._reader
        if reader is None:
            return True
        ok = reader.wait(ms)
        QtWidgets.QApplication.processEvents()
        return ok

    def _file_read(self, drop) -> None:
        self.drop = drop
        self.drop_view.show_drop(drop)
        message = f"{drop.name}: {drop.report.summary()}"
        self.note.setText(message)
        self.status.emit(message)
        self.reload()

    def _file_failed(self, message: str) -> None:
        self.drop_view.show_problem(message)
        self.note.setText(message)
        self.status.emit(message)

    def _reader_done(self) -> None:
        reader, self._reader = self._reader, None
        if reader is not None:
            reader.deleteLater()

    def genes(self) -> list:
        """The genes of the loaded list, or an empty list when none has been dropped."""
        return list(self.drop.genes) if self.drop is not None else []

    # ------------------------------------------------------------------ drag and drop
    @staticmethod
    def dropped_paths(event) -> list:
        """The local file paths a drag event carries."""
        data = event.mimeData()
        if data is None or not data.hasUrls():
            return []
        return [url.toLocalFile() for url in data.urls() if url.isLocalFile()]

    def dragEnterEvent(self, event):
        if self.dropped_paths(event):
            event.acceptProposedAction()

    def dragMoveEvent(self, event):
        if self.dropped_paths(event):
            event.acceptProposedAction()

    def dropEvent(self, event):
        paths = self.dropped_paths(event)
        if not paths:
            return
        event.acceptProposedAction()
        self.load_file(paths[0])

    def closeEvent(self, event):
        # A QThread destroyed while it is still running aborts the process, and the dock can be
        # closed in the middle of reading a workbook. `information_map_panel` learnt this one.
        for thread in (self._reader, self._loader):
            if thread is not None and thread.isRunning():
                thread.wait(5000)
        super().closeEvent(event)

    def keyPressEvent(self, event):
        K = QtCore.Qt.Key
        if event.key() in (K.Key_Left, K.Key_Up):
            self.focusPreviousChild()
            event.accept()
            return
        if event.key() in (K.Key_Right, K.Key_Down):
            self.focusNextChild()
            event.accept()
            return
        super().keyPressEvent(event)


class _WindowDrops(QtCore.QObject):
    """Catches a file dropped anywhere on the window and routes it to the gene-list dock.

    An event filter rather than a subclass, because the window is `app.Window` and this is a panel
    being installed into it: the filter is removed with the panel and the window needs no knowledge
    of it.
    """

    def __init__(self, window, panel: HomePanel):
        super().__init__(window)
        self.window, self.panel = window, panel
        window.setAcceptDrops(True)
        window.installEventFilter(self)

    def eventFilter(self, obj, event):
        """Accept a dragged file and hand a dropped one to the home panel."""
        T = QtCore.QEvent.Type
        if obj is self.window and event.type() in (T.DragEnter, T.DragMove, T.Drop):
            paths = HomePanel.dropped_paths(event)
            if not paths:
                return False
            event.acceptProposedAction()
            if event.type() == T.Drop:
                # A file dropped on the landing must not be read behind it.
                from .landing_install import leave_landing
                leave_landing(self.window)
                dock = getattr(self.window, "home_dock", None)
                if dock is not None:
                    dock.show()
                    dock.raise_()
                self.panel.choose(ND.get("list_drop", self.panel.nodes))
                self.panel.load_file(paths[0])
            return True
        return False


def gene_lookup(window):
    """A callable the navigator can ask for the genes matching a query.

    Reads the window's own node table -- identifier, product and whatever alias column it carries --
    so a search for "GRA16" or for "dense granule" reaches the gene. Returns None when the window
    has no table, and the search then simply offers no genes.
    """
    table = getattr(window, "nodes", None)
    if table is None or "gene_id" not in getattr(table, "columns", ()):
        return None
    columns = [c for c in ("gene_id", "product", "symbol", "alias", "gene_name")
               if c in table.columns]

    def matching(query: str) -> list:
        """The genes whose identifier, product or alias contains `query`, best (shortest) first."""
        text = str(query or "").strip().lower()
        if len(text) < 2:
            return []
        frame = table[columns]
        hit = None
        for column in columns:
            found = frame[column].astype(str).str.lower().str.contains(text, regex=False,
                                                                       na=False)
            hit = found if hit is None else (hit | found)
        rows = frame[hit].head(ND.GENE_LIMIT) if hit is not None else frame.head(0)
        product = "product" if "product" in columns else columns[-1]
        return [(str(row.gene_id), str(getattr(row, product, "")))
                for row in rows.itertuples(index=False)]

    return matching


def _index(window) -> list:
    """The help index for this window, built once and kept, for the navigator to rank.

    The same index the search field beside Help uses: built from the live window when there is one,
    so a panel or a command that is not in this window is not offered.
    """
    from . import help_index as HI
    from . import help_search as HS
    field = HS.field_of(window)
    if field is not None:
        return field.index()
    cached = getattr(window, "_home_index", None)
    if cached is None:
        cached = HI.build_index(window)
        window._home_index = cached
    return cached


def open_node(window, node: ND.Node) -> str:
    """Open a node's destination in the main window, and say in the status bar what happened.

    The only place an `Action` is acted on. Every branch reports rather than raising: a destination
    that has gone (a dock that failed to install) must leave the home screen usable.
    """
    action = ND.resolve(node)
    message = action.describe or node.title
    try:
        if action.kind == "unavailable":
            message = f"{node.title}: {node.unavailable_reason}"
        elif action.kind == "dock":
            dock = next((d for d in window.findChildren(QtWidgets.QDockWidget)
                         if d.windowTitle() == action.target), None)
            if dock is None:
                message = f"the {action.target} panel is not available in this window"
            else:
                dock.show()
                dock.raise_()
                if action.tab:
                    tabs = dock.findChild(QtWidgets.QTabWidget)
                    found = -1
                    if tabs is not None:
                        for index in range(tabs.count()):
                            if tabs.tabText(index) == action.tab:
                                tabs.setCurrentIndex(index)
                                found = index
                                break
                    message = (f"{action.target} ▸ {action.tab}" if found >= 0
                               else f"{action.target}: {action.tab} is not a tab here any more")
                else:
                    message = f"{action.target}" + (f" ▸ {action.section}" if action.section else "")
        elif action.kind == "window":
            method = getattr(window, action.target, None)
            if not callable(method):
                message = f"{node.title} cannot be opened in this window ({action.target})"
            elif action.target == "run_my_list":
                panel = getattr(window, "home", None)
                genes = panel.genes() if panel is not None else []
                message = (f"answering for {len(genes)} genes" if genes
                           else "no gene list loaded: drop or paste one first")
                if genes:
                    method(genes)
            else:
                method()
        elif action.kind == "dialog":
            from .data_download_dialog import DataDownloadDialog
            dialog = DataDownloadDialog(window)
            window._home_dialog = dialog              # kept alive; a local would be collected
            dialog.show()
        elif action.kind == "workflow":
            window._open_workflows(int(action.target))
        elif action.kind == "gene":
            # Select the gene first, then show its card: the card reads the selection, so the other
            # order shows the previous gene's card for a moment.
            chooser = getattr(window, "_workflow_gene", None)
            if callable(chooser):
                chooser(action.target)
            dock = getattr(window, "right_dock", None)
            if dock is not None:
                dock.show()
                dock.raise_()
            message = f"{action.target}: its card"
        elif action.kind == "entry":
            from . import help_search as HS
            message = HS.open_entry(window, node.entry) or f"opened {node.title}"
    except Exception as error:                        # a broken destination must not take the window
        message = f"could not open {node.title}: {type(error).__name__}: {error}"
    try:
        window.statusBar().showMessage(message, 6000)
    except Exception:
        pass
    return message


#: The dock's title. It was "home" until instruction 85: the landing is the only home now, and
#: what this dock still holds is the gene list a dropped file becomes, and the question catalogue.
DOCK_TITLE = "gene list"


def install(window) -> HomePanel | None:
    """Add the gene-list dock and make the whole window a target for a dropped file of genes.

    Until instruction 85 this was the home dock: a board of node cards, raised at startup. The
    landing replaced it as the home, and Einar found the old board still waiting behind it -- open
    Map from the Navigator and the right-hand dock group showed "the old home screen with tiles".
    So the board is not drawn (`cards=False`) and the dock is not raised: it is tabbed in with the
    others and comes up when a file is dropped, when the Drop module is opened, or when a search
    result it owns (a question, the catalogue) is chosen. Its attributes keep their names,
    `window.home` and `window.home_dock`, because the drop routing and the help search use them.
    """
    nodes_table = getattr(window, "nodes", None)
    known = (nodes_table["gene_id"].astype(str).tolist()
             if nodes_table is not None and "gene_id" in getattr(nodes_table, "columns", ())
             else None)
    panel = HomePanel(known_ids=known, gene_lookup=gene_lookup(window),
                      index=lambda: _index(window), cards=False)
    panel.status.connect(lambda message: window.statusBar().showMessage(message, 6000))
    panel.node_chosen.connect(lambda node: open_node(window, node))
    dock = QtWidgets.QDockWidget(DOCK_TITLE, window)
    dock.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
    scroll = QtWidgets.QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
    scroll.setMinimumWidth(120)
    scroll.setWidget(panel)
    dock.setWidget(scroll)
    window.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, dock)
    anchor = (getattr(window, "information_map_dock", None)
              or getattr(window, "discoveries_dock", None)
              or getattr(window, "strategies_dock", None)
              or getattr(window, "analysis_dock", None))
    if anchor is not None:
        window.tabifyDockWidget(anchor, dock)
        anchor.raise_()                               # tabifying puts the new tab on top
    window.home, window.home_dock = panel, dock
    window._home_drops = _WindowDrops(window, panel)
    return panel
