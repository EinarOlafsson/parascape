"""The GENE LIST CARD: what a list of genes shares that other genes do not.

Instruction 98, WP9. Einar asked for it in the same breath as the gene card: *"if i type in a
list of genes, the same logic but this time not gene card instead a gene list module and after
clicking it all the moduals that come downstream of gene lists."*

**The one sentence the mould asks for before anything is laid out.** This module lets you take
a list of genes and see what the genes of it share that other genes do not -- and how much of
the list is measured at all -- so that you can tell what the list points to, and which of its
gaps are worth inferring.

WHY THE CARD SHOWS WHAT IT SHOWS
================================
A list card that is a table of identifiers fails the first goal of the mould: it makes the
information no more accessible than the file the list came in, and synthesises nothing. The
list of identifiers is therefore the FOLDED section at the bottom, and what opens first is the
only thing about a list that cannot be read off the genes one at a time:

1. **What the genes of the list share that other genes do not.** Every label any gene can
   carry, counted inside the list and outside it, kept where the list holds more of it than
   the rest of the genome, with how surprising that is -- how rarely a list of the same size
   picked at random would hold as many. A share with nothing to compare it against cannot be
   read: half the list in the nucleus means nothing until you know what half the genome does.
2. **Where the genes of the list are in the parasite**, on the same drawing the gene card uses
   (:class:`parascape.celldiagram.CellDiagram`), coloured for one compartment at a time, with
   how much of the list is in it against how much of the rest of the genome is. The drawing is
   the thing Einar named as something the restructure must keep, and a list has a place in the
   cell exactly as one gene does.
3. **How much of the list is measured**, per kind of knowledge, against the rest of the genes.
   A list that nobody has measured and a kind of knowledge nobody has measured for anything
   look identical until both numbers are on screen; they have different answers.
4. **What is unknown for most of the list**, which is where inference is the only way forward
   and is therefore what this card hands to Infer -- the second goal of the mould.

Every sentence here is computed by :mod:`parascape.list_reading`, which has no Qt in it, so
each one is tested on a made-up table rather than by looking at a window.

WHAT FOLLOWS A LIST
===================
The row at the foot is :mod:`parascape.next_steps` -- the same strip, from the same function,
that the gene card uses -- asked about :data:`parascape.ports.GENE_LIST` instead of one gene.
The eight modules that read a list are offered pre-filled with it, each checked by
``ports.check_ready``, and one that cannot work keeps its place greyed with the reason.
"""
from __future__ import annotations

import os
from html import escape

import numpy as np
import pandas as pd
from PySide6 import QtCore, QtWidgets

from . import gene_reading as GR
from . import list_reading as LR
from . import ports as PORTS
from . import theme as TH
from . import organisms
from .recipes_panel import widget_setting
from .shell import ModuleShell, solid_dropdown


def _muted(text: str) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


def _browser() -> QtWidgets.QTextBrowser:
    view = QtWidgets.QTextBrowser()
    view.setOpenExternalLinks(True)
    return view


def _table(headers: list) -> QtWidgets.QTableWidget:
    """A read-only table with these headings, sized to its contents."""
    table = QtWidgets.QTableWidget(0, len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.verticalHeader().hide()
    table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
    table.horizontalHeader().setStretchLastSection(True)
    return table


class ListModule(QtCore.QObject):
    """The gene list card for one window."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.resolved = LR.Resolved()
        self.shared: tuple = ()
        self.coverage: tuple = ()
        self.source = ""
        self._label_columns = None                     # computed once; see `label_columns`
        self.shell = shell = ModuleShell.for_module("list")
        shell.module = self
        shell.setAcceptDrops(True)
        shell.installEventFilter(self)
        self._build_settings()
        self._build_runtime()
        shell.home_requested.connect(lambda: _go_home(window))

    # ------------------------------------------------------------------ the list itself
    @property
    def genes(self) -> list:
        """The genes of the list that are in the gene table of this window."""
        return list(self.resolved.genes)

    def label_columns(self) -> list:
        """The columns worth asking what the list shares, computed once for this window.

        Cached because it does not depend on the list -- it is a property of the table --
        and because deciding it means counting the distinct values of every column, which
        is seconds of work on a table of eight thousand genes and several hundred columns.
        Recomputing it on every keystroke would make the card feel broken.
        """
        if self._label_columns is None:
            try:
                self._label_columns = LR._label_columns(self._nodes(), self._organism())
            except Exception:                          # a card must open on any table
                self._label_columns = []
        return self._label_columns

    # ------------------------------------------------------------------ left
    def _build_settings(self) -> None:
        shell = self.shell
        mine = QtWidgets.QWidget(shell)
        v = QtWidgets.QVBoxLayout(mine)
        v.setContentsMargins(6, 2, 6, 6)
        self.box = QtWidgets.QPlainTextEdit()
        self.box.setPlaceholderText("Gene ids, separated by commas, spaces or new lines")
        self.box.setToolTip(TH.tip(
            "Paste a column out of a spreadsheet or type the ids with commas between them -- "
            "either works. An identifier that names no gene in this table is reported below "
            "rather than dropped, so a list of 300 that reaches 250 genes says which 50 did "
            "not land. Ctrl+Return reads the list."))
        self.box.setMaximumHeight(120)
        # Ctrl+Return reads the list. A bare Return cannot: this is a box with several lines
        # in it, and taking Return would stop anybody pasting a column.
        self.box.installEventFilter(self)
        v.addWidget(self.box)
        row = QtWidgets.QHBoxLayout()
        read = QtWidgets.QPushButton("Read this list")
        read.setToolTip(TH.tip("Match what is typed above against the genes of this table and "
                               "read the list. Ctrl+Return does the same."))
        read.clicked.connect(lambda _checked=False: self._typed())
        choose = QtWidgets.QPushButton("Choose a list…")
        choose.setToolTip(TH.tip("A table of gene ids: CSV, TSV or Excel. You can also drop "
                                 "the file on this module."))
        choose.clicked.connect(self._choose_list)
        row.addWidget(read)
        row.addWidget(choose)
        v.addLayout(row)
        row2 = QtWidgets.QHBoxLayout()
        use = QtWidgets.QPushButton("Use the gene list")
        use.setToolTip(TH.tip("The list last dropped on ParaScape (the gene list panel)."))
        use.clicked.connect(self._use_gene_list)
        gated = QtWidgets.QPushButton("Use the genes gated on the map")
        gated.setToolTip(TH.tip("The genes inside the last gate drawn on the map. Nothing "
                                "gated yet means nothing to take."))
        gated.clicked.connect(self._use_gated)
        clear = QtWidgets.QPushButton("Clear")
        clear.setToolTip(TH.tip("Forget the list and empty the card."))
        # Through a lambda: `clicked` carries a checked flag, and handing it straight to a
        # slot that takes no arguments is a call that only fails when somebody presses it.
        clear.clicked.connect(lambda _checked=False: self.clear())
        row2.addWidget(use)
        row2.addWidget(gated)
        row2.addWidget(clear)
        v.addLayout(row2)
        self.list_note = _muted("No genes yet.")
        self.list_note.setToolTip(TH.tip(
            "How many genes the list reached, how many identifiers reached none, and how many "
            "were repeated. What reached nothing is listed, never silently dropped."))
        v.addWidget(self.list_note)
        shell.add_setting("Your genes", mine, subject=PORTS.GENE_LIST)

        lens = QtWidgets.QWidget(shell)
        ll = QtWidgets.QVBoxLayout(lens)
        ll.setContentsMargins(6, 2, 6, 6)
        self.place_box = QtWidgets.QComboBox()
        self.place_box.setToolTip(TH.tip(
            "Which measurement is read as 'where they are' when the parasite is drawn. "
            "Several measurements place a gene in the cell -- the curated compartment and the "
            "individual estimates behind it -- and they do not always agree. This is the same "
            "choice as 'color by' in the find panel whenever the two can be the same."))
        self.place_box.addItems(self._placement_columns())
        solid_dropdown(self.place_box)
        self.place_box.setCurrentText(GR.preferred_placement_column(
            self._placement_columns(), self._window_colour()))
        self.place_box.currentTextChanged.connect(self._place_chosen)
        ll.addWidget(self.place_box)
        self.place_note = _muted("")
        ll.addWidget(self.place_note)
        shell.add_setting("Where they are, read by", lens)
        self._follow_window_colour()

        kinds = QtWidgets.QWidget(shell)
        kl = QtWidgets.QVBoxLayout(kinds)
        kl.setContentsMargins(6, 2, 6, 6)
        kl.addWidget(_muted("What to read about them:"))
        self.include: dict = {}
        for category, title in GR.ORDER:
            if category in LR.SKIP_CATEGORIES:
                continue
            box = QtWidgets.QCheckBox(title)
            box.setChecked(True)
            box.setToolTip(TH.tip(f"Count {title.lower()} when asking what the list shares, "
                                  "and report how much of the list it is measured for."))
            box.toggled.connect(lambda *_: self.render())
            kl.addWidget(box)
            self.include[category] = box
        shell.add_setting("What to read", kinds)

        advanced = QtWidgets.QWidget(shell)
        al = QtWidgets.QVBoxLayout(advanced)
        al.setContentsMargins(6, 2, 6, 6)
        al.addWidget(_muted("What counts as something the list shares:"))
        self.min_in_list = QtWidgets.QSpinBox()
        self.min_in_list.setRange(2, 500)
        self.min_in_list.setValue(LR.MIN_IN_LIST)
        self.min_in_list.setPrefix("at least ")
        self.min_in_list.setSuffix(" genes of the list carry it")
        self.min_in_list.setToolTip(TH.tip(
            "A label one gene carries is a fact about that gene, not about the list, so two "
            "is the lowest this goes. Raise it for a long list, where two genes sharing "
            "something is not news."))
        self.min_in_list.valueChanged.connect(lambda *_: self.render())
        al.addWidget(self.min_in_list)
        self.how_many = QtWidgets.QSpinBox()
        self.how_many.setRange(3, 60)
        self.how_many.setValue(12)
        self.how_many.setPrefix("show the top ")
        self.how_many.setToolTip(TH.tip("How many shared labels to show, most surprising "
                                        "first."))
        self.how_many.valueChanged.connect(lambda *_: self.render())
        al.addWidget(self.how_many)
        shell.add_setting("Advanced", advanced, expanded=False)

        keep = QtWidgets.QWidget(shell)
        kp = QtWidgets.QVBoxLayout(keep)
        kp.setContentsMargins(6, 2, 6, 6)
        kp.addWidget(_muted("Keep this reading as files you can hand to somebody else, "
                            "or read one back."))
        save = QtWidgets.QPushButton("Save everything as a bundle…")
        save.setToolTip(TH.tip("One folder: the genes of the list, what they share with the "
                               "numbers behind it, how much of the list is measured, the "
                               "methods text with references, and a manifest whose hashes "
                               "show if anything was edited later."))
        save.clicked.connect(lambda _checked=False: self.save_bundle())
        reopen = QtWidgets.QPushButton("Open a saved bundle…")
        reopen.setToolTip(TH.tip("Read a bundle saved here back: its files are checked against "
                                 "their hashes and its list is restored."))
        reopen.clicked.connect(lambda _checked=False: self.open_bundle())
        kp.addWidget(save)
        kp.addWidget(reopen)
        shell.add_setting("Save", keep, expanded=False)
        shell.add_recipes(self.recipe_settings(), after_apply=self.render)

    def recipe_settings(self) -> list:
        """The settings of the Gene list card that a recipe keeps (instruction 98, WP11).

        THE LIST ITSELF IS NOT ONE: a list of genes is the subject of this module, so a recipe
        is how the card READS a list and is then applied to whichever list is in it. The two
        bars that decide what counts as something the list shares are the settings worth
        keeping under a name, because they are the choice somebody argues about and will want
        to use again on the next list.
        """
        out = [widget_setting("placement_column", self.place_box, "where they are, read by"),
               widget_setting("min_in_list", self.min_in_list,
                              "how many of the list must carry it"),
               widget_setting("how_many", self.how_many, "how many shared labels to show")]
        for category, box in self.include.items():
            out.append(widget_setting(f"read_{category}", box, f"what to read: {category}"))
        return out

    # ------------------------------------------------------------------ right
    def _build_runtime(self) -> None:
        shell, window = self.shell, self.window
        self._build_shared()
        self._build_parasite()

        covered = QtWidgets.QWidget(shell)
        cv = QtWidgets.QVBoxLayout(covered)
        cv.setContentsMargins(8, 4, 8, 8)
        self.coverage_lead = _muted("")
        self.coverage_lead.setToolTip(TH.tip(
            "How many kinds of knowledge this list is measured for at all, and how many are "
            "unknown for most of it."))
        cv.addWidget(self.coverage_lead)
        # A TABLE rather than sentences, which is what looking at it decided. Fifteen kinds
        # of knowledge as fifteen paragraphs showed two of them in the height the pane has,
        # with the other thirteen behind a scrollbar; as rows, the same height shows six and
        # they can be read down the column instead of one at a time.
        self.coverage_view = _table(["kind of knowledge", "of your list", "of other genes",
                                     "unknown for"])
        self.coverage_view.setToolTip(
            "How much of the list is measured for each kind of knowledge, and how much of the "
            "rest of the genes is. A list measured less than the genome, and a kind of "
            "knowledge nobody has measured for anything, are different problems with "
            "different answers; the two columns tell them apart.")
        cv.addWidget(self.coverage_view, 1)
        shell.add_section("How much of the list is measured", covered,
                          extent=380, stretch=2)

        self.unknown_view = _browser()
        shell.link_style(self.unknown_view)
        self.unknown_view.setToolTip(
            "What is unknown for most of the list. These are the places where inference is "
            "the only way forward, which is what the Predict button below is for. Unknown "
            "never means no.")
        shell.add_section("What is unknown for most of the list, and what can be inferred",
                          self.unknown_view, extent=240, stretch=1)

        self.members = _table(["gene", "product", "where it is"])
        self.members.setToolTip("The genes of the list themselves. Click one to read it on "
                                "its own in the gene card.")
        self.members.cellDoubleClicked.connect(lambda r, _c: self._open_member(r))
        shell.add_section("The genes of the list", self.members, extent=220, stretch=1,
                          folded=True)

        from .next_steps import NextSteps
        self.next_steps = NextSteps(PORTS.GENE_LIST, shell,
                                    lead="Next, with this list:")
        self.next_steps.chosen.connect(lambda key, _kind, value: self._open(key, value))
        shell.add_section("What you can do with this list next", self.next_steps,
                          extent=110, stretch=0)
        self._offer_next_steps()

        shell.add_action("Put the list on the map", self._on_map, primary=True,
                         tip="Open the Map module with these genes marked, and read what "
                             "surrounds them.")
        shell.add_action("Predict what is missing", lambda: self._open("infer"),
                         tip="Answer a question about these genes from the evidence, with how "
                             "often such answers are right (Infer). Infer is asked a question "
                             "rather than a list, which is why it is here as well as among "
                             "the next steps.")
        shell.add_action("Read one gene", self._read_one,
                         tip="Open the gene card on the gene selected below, or on the first "
                             "of the list.")
        shell.add_action("Export the list", lambda: self._open("export"),
                         tip="Write these genes out, with the columns you choose.")
        shell.borrow("Console and AI chat", [("Console", getattr(window, "console_dock", None)),
                                             ("AI chat", getattr(window, "chat_dock", None))],
                     folded=True)
        self.render()

    def _build_shared(self) -> None:
        """What the genes of the list share, which is the thing this card is about.

        First in the runtime column, large, with the focus: the mould asks the right-hand
        column to open with the real current data rather than with a method or a button, and
        for a list the real data is the comparison with every other gene. Nothing else here
        is a fact about the LIST -- the rest can be read off the genes one at a time.
        """
        shell = self.shell
        holder = QtWidgets.QWidget(shell)
        v = QtWidgets.QVBoxLayout(holder)
        v.setContentsMargins(8, 4, 8, 8)
        self.shared_lead = _muted("")
        self.shared_lead.setToolTip(TH.tip(
            "What the list was compared with: every other gene in the table that is measured "
            "for the same thing. A share of the list on its own cannot be read -- half of it "
            "in the nucleus is only a finding if half the genome is not."))
        v.addWidget(self.shared_lead)
        self.shared_view = _table(["what they share", "measurement", "of your list",
                                   "of other genes", "how surprising"])
        self.shared_view.setToolTip(
            "Every label the genes of the list carry more often than the other genes do, most "
            "surprising first. 'How surprising' is how rarely a list of this size picked at "
            "random would hold as many. Click a row that names a compartment to colour it in "
            "the parasite below.")
        self.shared_view.cellClicked.connect(lambda r, _c: self._shared_clicked(r))
        v.addWidget(self.shared_view, 1)
        shell.add_section("What the genes of the list share", holder, extent=620, stretch=5,
                          focus=True)

    def _build_parasite(self) -> None:
        """The parasite drawn, with the organelle most of the list is found in filled in.

        The same drawing as the gene card, for the same reason: Einar named the schematic as
        something the restructure must keep, and a list sits somewhere in the cell exactly as
        one gene does. A SECOND INSTANCE of the widget rather than the one from the `find`
        dock moved here -- the rule of the mould, "move, do not copy", is about CONTROLS
        whose methods must keep reading them, and the `find` dock is lent whole to the Cell
        module, which could not then borrow a widget this card had taken.
        """
        shell = self.shell
        from .celldiagram import CellDiagram, available as diagram_available
        holder = QtWidgets.QWidget(shell)
        # The caption BESIDE the drawing, never under it. The drawing is portrait, so its
        # width follows its height; under a caption in a column a thousand pixels wide it is
        # a thumbnail with nine hundred pixels of nothing on either side, and in a pane too
        # short for both Qt draws the caption ON TOP of the parasite. Both of those were
        # found by looking at WP8 and neither was caught by a test.
        layout = QtWidgets.QHBoxLayout(holder)
        layout.setContentsMargins(8, 4, 8, 8)
        layout.setSpacing(16)
        self.parasite = CellDiagram(holder) if diagram_available() else None
        if self.parasite is not None:
            self.parasite.compartment_clicked.connect(self._organelle_clicked)
            self.parasite.setMinimumHeight(0)
            layout.addWidget(self.parasite, 2)
        words = QtWidgets.QWidget(holder)
        wl = QtWidgets.QVBoxLayout(words)
        wl.setContentsMargins(0, 0, 0, 0)
        wl.addStretch(1)
        self.parasite_note = _muted("Type or choose a list on the left to see where in the "
                                    "cell the genes of it are.")
        self.parasite_note.setToolTip(TH.tip(
            "How much of the list is in the organelle coloured in, against how much of the "
            "rest of the genes is. Click another organelle to read that one instead."))
        wl.addWidget(self.parasite_note)
        self.parasite_trust = _muted("")
        self.parasite_trust.setToolTip(TH.tip(
            "How much of the list has a measured compartment at all, how many compartments "
            "the placed ones are spread over, and that the unmeasured ones are unknown "
            "rather than nowhere."))
        wl.addWidget(self.parasite_trust)
        wl.addStretch(1)
        layout.addWidget(words, 3)
        shell.add_section("Where the genes of the list are in the parasite", holder,
                          extent=400, stretch=4)
        self.compartment = ""                          # the organelle currently coloured

    # ------------------------------------------------------------------ the window
    def _nodes(self):
        return getattr(self.window, "nodes", None)

    def _organism(self):
        return getattr(getattr(self.window, "panel", None), "organism", None)

    def _reasons(self) -> dict:
        """Which modules this window cannot actually open, by the docks it really built.

        Asked of the window rather than of the catalogue so a greyed offer says what is true
        HERE: a dock that failed to install greys its module with that reason.
        """
        try:
            from .landing_install import _reasons
            return _reasons(self.window)
        except Exception:
            return {}

    def _placement_columns(self) -> list:
        """The measurements of the gene table that place a gene in an organelle the drawing has."""
        from .celldiagram import COMPARTMENT_SL
        try:
            return GR.placement_columns(self._nodes(), COMPARTMENT_SL)
        except Exception:                              # a card must open on any table
            return []

    def _window_colour(self) -> str:
        """What the window is colouring by right now, or "" when it has no such control."""
        box = getattr(self.window, "category_box", None)
        return box.currentText() if box is not None else ""

    def _follow_window_colour(self) -> None:
        """Keep the lens in step with the colour-by control of the window, both ways.

        One source of truth, the same rule the gene card is written to: the card and the map
        must never show two different compartments for one gene. The window wins when its
        choice places genes; when it is something else -- a run, a binned number -- the card
        keeps the placement it had, because there is nothing sensible to follow it to.
        """
        box = getattr(self.window, "category_box", None)
        if box is not None:
            box.currentTextChanged.connect(self._colour_changed_elsewhere)

    def _colour_changed_elsewhere(self, text: str) -> None:
        """The window was recoloured: follow it when the new colouring places genes."""
        if text and text in self._placement_columns() and self.place_box.currentText() != text:
            self.place_box.setCurrentText(text)

    def _place_chosen(self, text: str) -> None:
        """The lens changed: redraw the parasite, and recolour the window to match."""
        box = getattr(self.window, "category_box", None)
        if box is not None and text and box.currentText() != text:
            try:
                if text in self.window.color_sources():
                    box.setCurrentText(text)
            except Exception:                          # recolouring is not worth a crash
                pass
        self.compartment = ""
        self._render_parasite()

    # ------------------------------------------------------------------ choosing a list
    def set_genes(self, ids, source: str = "") -> LR.Resolved:
        """Take a list of identifiers, resolve it against the table, and read it.

        :param ids: the identifiers, in any order; strings.
        :param source: where the list came from, for the line under the box.
        :returns: the :class:`parascape.list_reading.Resolved` the card is now showing.
        """
        self.resolved = LR.resolve(self._nodes(), list(ids or []))
        self.source = source
        self.compartment = ""
        self.render()
        return self.resolved

    def _typed(self) -> None:
        """Read whatever is in the box as a list of identifiers."""
        self.set_genes(LR.split_identifiers(self.box.toPlainText()), "what you typed")

    def load_list(self, path: str) -> None:
        """Read a file of genes and show the list it holds.

        :param path: a CSV, TSV or Excel file with a column of identifiers.
        """
        from . import gene_drop as GD
        nodes = self._nodes()
        known = nodes["gene_id"].astype(str).tolist() if nodes is not None else None
        try:
            drop = GD.read_drop(path, known_ids=known)
        except Exception as error:
            self.shell.say(f"{os.path.basename(path)} could not be read: {error}")
            return
        self.box.setPlainText("\n".join(drop.genes))
        self.set_genes(drop.genes, os.path.basename(path))

    def _choose_list(self) -> None:
        from . import gene_drop as GD
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self.shell, "Choose a list of genes", "",
            "Tables (" + " ".join(f"*{e}" for e in GD.READABLE) + ")")
        if path:
            self.load_list(path)

    def _use_gene_list(self) -> None:
        panel = getattr(self.window, "home", None)
        genes = panel.genes() if panel is not None else []
        if not genes:
            self.shell.say("No gene list loaded: drop a file on this module or choose one.")
            return
        self.box.setPlainText("\n".join(str(g) for g in genes))
        self.set_genes(genes, "the gene list")

    def _use_gated(self) -> None:
        """Take the genes inside the last gate drawn on the map.

        The gate is how the application already says "this set of genes is the one in hand",
        and it is what `Export gated selection` writes, so the list card reads the same set
        rather than asking somebody to export it and load the file back.
        """
        nodes, gated = self._nodes(), getattr(self.window, "gated", None)
        if nodes is None or gated is None or not len(gated):
            self.shell.say("Nothing is gated on the map: drag a gate over some genes first.")
            return
        genes = [str(g) for g in nodes["gene_id"].astype(str).to_numpy()[np.asarray(gated, int)]]
        self.box.setPlainText("\n".join(genes))
        self.set_genes(genes, "the gate on the map")

    def clear(self) -> None:
        """Forget the list and empty the card.

        The mould asks every module for a way back to nothing. Without it the only way to
        stop looking at a list was to paste another one over it.
        """
        self.box.clear()
        self.set_genes([], "")
        self.shell.say("")

    def eventFilter(self, obj, event):                 # noqa: N802 - Qt name
        """Accept a file dragged onto the module and read its genes as the list."""
        T = QtCore.QEvent.Type
        if obj is self.shell and event.type() in (T.DragEnter, T.DragMove, T.Drop):
            from .home_panel import HomePanel
            paths = HomePanel.dropped_paths(event)
            if not paths:
                return False
            event.acceptProposedAction()
            if event.type() == T.Drop:
                self.load_list(paths[0])
            return True
        # `getattr`, because the filter is installed on the shell in `__init__` and events
        # arrive while the settings are still being built: `self.box` does not exist yet, and
        # naming it directly raised an AttributeError inside a Qt override, which Qt reports
        # as three hundred nested "Error calling Python override" lines.
        if (obj is getattr(self, "box", None) and event.type() == T.KeyPress
                and event.key() in (QtCore.Qt.Key.Key_Return, QtCore.Qt.Key.Key_Enter)
                and event.modifiers() & QtCore.Qt.KeyboardModifier.ControlModifier):
            self._typed()
            return True
        return False

    # ------------------------------------------------------------------ reading it
    def wanted(self, category: str) -> bool:
        """Whether this kind of knowledge is ticked in 'What to read'."""
        box = self.include.get(category)
        return box is None or box.isChecked()

    def wanted_columns(self) -> list:
        """The label columns to ask what the list shares, after the ticks on the left.

        The ticks are kinds of knowledge and the comparison is per column, so the choice has
        to be turned into columns somewhere; here, once, rather than in each sentence.
        """
        return [c for c in self.label_columns() if self.wanted(GR.category_of(c, self._organism()))]

    def spread(self):
        """Where the genes of the list sit, by the measurement the lens names.

        :returns: a :class:`parascape.list_reading.Spread`; its ``measured`` is 0 when
            nothing of the list has a compartment.
        """
        return LR.spread(self._nodes(), self.genes, self.place_box.currentText())

    def render(self) -> None:
        """Read the list again and rewrite every section of the card.

        One entry point, called by every control, because a card where each control updated
        its own corner is a card that can show two readings of two different lists at once.
        """
        self.list_note.setText(self.resolved.sentence()
                               + (f" From {escape(self.source)}." if self.source else ""))
        nodes = self._nodes()
        self.coverage = LR.coverage(nodes, self.genes, self._organism()) if self.genes else ()
        self.shared = (LR.shared_labels(nodes, self.genes, columns=self.wanted_columns(),
                                        min_in_list=self.min_in_list.value(),
                                        limit=self.how_many.value(),
                                        organism=self._organism())
                       if self.genes else ())
        self._render_shared()
        self._render_parasite()
        self._render_coverage()
        self._render_unknown()
        self._render_members()
        self._offer_next_steps()

    def _offer_next_steps(self) -> None:
        """Offer every module that reads a list, pre-filled with this one."""
        self.next_steps.set_subject(",".join(self.genes), available=self._reasons(),
                                    without="list", shown=self.in_hand())

    def in_hand(self) -> str:
        """How to SAY what is in hand: "40 genes", not forty accessions.

        A list is the first subject in ParaScape whose value does not read as a name. Every
        offer and every tooltip that would have printed the value prints this instead.
        """
        n = self.resolved.count
        return f"{n:,} gene{'' if n == 1 else 's'}" if n else ""

    def _render_shared(self) -> None:
        """Write what the genes of the list share, most surprising first."""
        self.shared_lead.setText(LR.shared_lead(self.shared, self.resolved))
        table = self.shared_view
        table.setRowCount(len(self.shared))
        for i, s in enumerate(self.shared):
            measurement = str(s.column).replace("_", " ")
            if s.also:
                measurement += f" (and {len(s.also):,} more that agree)"
            cells = (GR.value_text(s.column, s.value),
                     measurement,
                     f"{s.in_list:,} of {s.list_measured:,} ({s.list_percent} %)",
                     f"{s.outside_percent} % of {s.outside_measured:,}",
                     LR.surprise_words(s))
            for j, text in enumerate(cells):
                item = QtWidgets.QTableWidgetItem(text)
                item.setToolTip(LR.shared_sentence(s))
                table.setItem(i, j, item)
        table.resizeColumnsToContents()
        table.horizontalHeader().setStretchLastSection(True)

    def _shared_clicked(self, row: int) -> None:
        """A shared label clicked: colour that compartment in the parasite, when it is one."""
        if not 0 <= row < len(self.shared):
            return
        s = self.shared[row]
        if s.column == self.place_box.currentText() or s.column in self._placement_columns():
            if s.column != self.place_box.currentText():
                self.place_box.setCurrentText(s.column)
            self.compartment = s.value
            self._render_parasite()
            self.shell.say(f"{s.value} coloured in the parasite")
        else:
            self.shell.say(LR.surprise_words(s).capitalize() + f": {s.value}")

    def _adopt_lens(self, compartment: str) -> bool:
        """Colour the window by the lens when it has no colour for the organelle being shown.

        ONE PALETTE, which is the rule `celldiagram` is built on and the gene card follows:
        the fills come from `window.color_of`, the same dict the points on the map are drawn
        from, so the drawing and the map cannot disagree about what a colour means. The cost
        of that rule is this method -- when the window is coloured by something else, there
        is no colour for "dense granules" and the organelle comes out unfilled.

        Rather than mix in a colour of its own, the card moves the window onto the
        measurement it is reading, and says so.

        :param compartment: the organelle about to be coloured.
        :returns: whether the window was recoloured.
        """
        colours = getattr(self.window, "color_of", {}) or {}
        box = getattr(self.window, "category_box", None)
        column = self.place_box.currentText()
        if not compartment or compartment in colours or box is None or not column:
            return False
        try:
            if column not in self.window.color_sources():
                return False
            box.setCurrentText(column)
        except Exception:                              # colour is not worth a crash
            return False
        self.shell.say(f"The map and the find panel are now coloured by "
                       f"{column.replace('_', ' ')}, so the drawing and the map mean the "
                       f"same thing by a colour.")
        return True

    def _render_parasite(self) -> None:
        """Fill the organelle the list mostly sits in, and say what that rests on."""
        sp = self.spread()
        name = self.compartment or sp.commonest
        self.place_note.setText(f"Reading {str(sp.column).replace('_', ' ')}."
                                if sp.column else "")
        self._adopt_lens(name)
        if self.parasite is not None:
            # The window's own palette, never a second one: two palettes are two claims
            # about what a colour means, and one of them drifts.
            self.parasite.set_palette(getattr(self.window, "color_of", {}) or {}, name)
        self.parasite_note.setText(LR.spread_sentence(sp, name))
        self.parasite_trust.setText(LR.spread_trust(sp))

    def _organelle_clicked(self, name: str) -> None:
        """An organelle clicked: read the list for that compartment, and tick it in the window."""
        self.compartment = name
        self._render_parasite()
        picker = getattr(self.window, "select_compartment", None)
        if callable(picker):
            picker(name)
        self.shell.say(f"{name}: {self.spread().count_of(name):,} genes of the list")

    def _render_coverage(self) -> None:
        """Write how much of the list is measured, kind of knowledge by kind of knowledge.

        Least measured first, because the emptiest kinds are the ones that decide what to do
        next -- they are what the Predict button is for -- and the ones measured for all of
        it need no decision at all.
        """
        shown = sorted((c for c in self.coverage if self.wanted(c.category)),
                       key=lambda c: (c.percent, -c.unknown))
        table = self.coverage_view
        table.setRowCount(len(shown))
        for i, c in enumerate(shown):
            cells = (c.title,
                     f"{c.in_list:,} of {c.listed:,} ({c.percent} %)",
                     f"{c.outside_percent} % of {c.outside:,}",
                     f"{c.unknown:,}" if c.unknown else "none of it")
            for j, text in enumerate(cells):
                item = QtWidgets.QTableWidgetItem(text)
                item.setToolTip(LR.coverage_sentence(c))
                if j == 3 and c.unknown:
                    item.setForeground(QtCore.Qt.GlobalColor.gray)
                table.setItem(i, j, item)
        table.resizeColumnsToContents()
        table.horizontalHeader().setStretchLastSection(True)
        if not shown:
            self.coverage_lead.setText(
                "Choose a list on the left. This says how much of it is measured for each "
                "kind of knowledge, against how much of the rest of the genes is.")
            return
        poor = [c for c in shown if c.mostly_unknown]
        self.coverage_lead.setText(
            f"Of the {len(shown):,} kinds of knowledge this table holds, "
            f"{len(shown) - len(poor):,} are measured for more than half of your "
            f"{self.resolved.count:,} genes and {len(poor):,} are unknown for most of them.")

    def _render_unknown(self) -> None:
        """Write what is unknown for most of the list, which is what inference is for."""
        poor = [c for c in LR.unknown_for_most(self.coverage) if self.wanted(c.category)]
        if not self.genes:
            self.unknown_view.setHtml(
                "<p>Choose a list on the left. This says what is unknown for most of it, "
                "which is where inference is the only way forward.</p>")
            return
        if not poor:
            self.unknown_view.setHtml(
                "<p>Every kind of knowledge chosen is measured for more than half of this "
                "list. There is no gap here big enough to need inference.</p>")
            return
        body = "".join(f"<p>{LR.unknown_sentence(c)}</p>" for c in poor)
        self.unknown_view.setHtml(
            body + "<p style='color:#999'>Predict what is missing, below, asks a question "
                   "about these genes and answers it with how often such an answer is "
                   "right.</p>")

    def _render_members(self) -> None:
        """Write the genes of the list themselves, with what each is and where it is."""
        nodes, table = self._nodes(), self.members
        genes = self.genes
        table.setRowCount(len(genes))
        column = self.place_box.currentText()
        for i, gene in enumerate(genes):
            product, where = "", ""
            if nodes is not None:
                hit = nodes[nodes["gene_id"].astype(str) == gene]
                if len(hit):
                    product = str(hit.iloc[0].get("product", "") or "")
                    where = GR.placed_at(nodes, gene, column) if column else ""
            for j, text in enumerate((gene, product[:70], where or "unknown")):
                item = QtWidgets.QTableWidgetItem(text)
                item.setToolTip("Double-click to read this gene on its own in the gene card.")
                if j == 2 and not where:
                    item.setForeground(QtCore.Qt.GlobalColor.gray)
                table.setItem(i, j, item)
        table.resizeColumnsToContents()
        table.horizontalHeader().setStretchLastSection(True)

    def _open_member(self, row: int) -> str:
        """Open one gene of the list in the gene card; returns what was opened."""
        if not 0 <= row < len(self.genes):
            return ""
        from .landing_install import open_module
        return open_module(self.window, "gene", subject=PORTS.GENE, value=self.genes[row])

    def _read_one(self) -> str:
        """Open the gene card on the selected gene of the list, or on the first of it."""
        row = self.members.currentRow()
        return self._open_member(row if row >= 0 else 0)

    # ------------------------------------------------------------------ actions
    def _on_map(self) -> None:
        """Open the Map module with the genes of the list marked on it."""
        nodes = self._nodes()
        map_shell = getattr(self.window, "map_module", None)
        if map_shell is None or nodes is None or not self.genes:
            return
        ids = nodes["gene_id"].astype(str).to_numpy()
        index = [int(i) for i in np.flatnonzero(np.isin(ids, np.array(self.genes)))]
        from .map_module import open_in
        open_in(self.window)
        if index:
            map_shell.module.select(index, "")

    def _open(self, key: str, value: str = "") -> str:
        """Open another module with this list already in it; returns what was opened.

        The pre-filling half of an offer: a next step that opened an empty module would
        leave the list to be pasted again, which is the whole thing the chaining strip
        exists to stop. The list goes through `landing_install.open_module`, which is the
        one route every module is opened by.
        """
        from .landing_install import open_module
        return open_module(self.window, key, subject=PORTS.GENE_LIST,
                           value=value or ",".join(self.genes))

    def save_bundle(self, folder: str = "") -> str:
        """Write the current list and its readings as an export bundle; returns the folder."""
        from . import export_bundle
        if not self.genes:
            self.shell.say("choose a list first; the bundle holds its reading")
            return ""
        folder = folder or export_bundle.ask_folder(self.shell)
        if not folder:
            return ""
        members = pd.DataFrame({"gene_id": self.genes})
        shared = pd.DataFrame([{
            "measurement": s.column, "label": s.value, "of the list": s.in_list,
            "the list measured": s.list_measured, "of other genes": s.outside_with,
            "other genes measured": s.outside_measured,
            "how surprising": LR.surprise_words(s)} for s in self.shared])
        covered = pd.DataFrame([{
            "kind of knowledge": c.title, "measured for": c.in_list, "of the list": c.listed,
            "share of the list": c.percent, "share of other genes": c.outside_percent}
            for c in self.coverage])
        tables = {"the genes of the list": members}
        kinds = {"the genes of the list": "measured"}
        if len(shared):
            tables["what the list shares"], kinds["what the list shares"] = shared, "computed"
        if len(covered):
            tables["how much is measured"], kinds["how much is measured"] = covered, "computed"
        if self.resolved.missing:
            tables["identifiers that named no gene"] = pd.DataFrame(
                {"identifier": list(self.resolved.missing)})
            kinds["identifiers that named no gene"] = "measured"
        export_bundle.write(
            folder,
            organism=getattr(self.window.panel, "organism", organisms.TOXOPLASMA),
            # The identifiers themselves, not a count: a bundle that said "40 genes" could
            # be read but not reopened, and reopening it is what `open_bundle` is for.
            query={"module": "list", "genes": list(self.genes), "source": self.source},
            tables=tables, kinds=kinds,
            columns=[s.column for s in self.shared if s.column in (self._nodes().columns
                                                                   if self._nodes() is not None
                                                                   else ())],
            nodes=self._nodes())
        self.shell.say(f"bundle written: {folder}")
        return folder

    def open_bundle(self, folder: str = "") -> dict:
        """Read a bundle saved by this module back and restore the list it held."""
        from . import export_bundle
        return export_bundle.reopen(self, folder)


def build(window) -> ModuleShell | None:
    """The gene list module for `window`, or None when the window has no gene table."""
    if getattr(window, "nodes", None) is None:
        return None
    return ListModule(window).shell


def open_in(window, genes: str = "") -> str:
    """Show the gene list module over the window, reading `genes` when given.

    :param window: the main window.
    :param genes: the list, as identifiers separated by commas, spaces or new lines -- the
        same text the Navigator resolved, so nothing has to be typed a second time. A single
        identifier is a list of one and is read as such rather than refused.
    """
    from .landing_install import enter_landing
    landing = getattr(window, "landing", None)
    shell = getattr(window, "list_module", None)
    if landing is None or shell is None:
        return "the gene list module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    ids = LR.split_identifiers(genes)
    if ids:
        shell.module.box.setPlainText("\n".join(ids))
        shell.module.set_genes(ids, "the navigator")
    return "Gene list"


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)
