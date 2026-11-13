"""The GENE CARD: everything recorded about one gene, synthesised, with its sources.

Built with the parascape-module mould (instruction 96); Einar chose "synthesis first"
(2026-11-06). The module lets you take a gene -- or step through a list -- and read what is
known about it in sentences (where it is, what it binds, how much the parasite needs it, when
it is made...) each with its sources, then the evidence itself grouped the same way, then what
is NOT known and what can be inferred, with how often such a call is right.

The facts come from the published information map (`information_map.PublishedMap`, loaded by
the information-map panel) and the window's own gene table; `gene_reading` turns them into the
sentences, so every sentence is testable without a window.

WHAT WORK PACKAGE 8 OF INSTRUCTION 98 CHANGED, AND WHY
=====================================================
Einar, 2026-11-12: "if i type a gene into the navigator i get a few moduals that pop up, one
should be the gene card (overhauled to mach the spacr module structure but with the same
information as before including the parasite schematic), when i click the gene card once i
should see the moduals downstream of one gene (like the network in spacr)."

Three things followed, and nothing was taken away to make room for them:

1. **The parasite schematic is the thing itself.** The mould asks the right-hand column to
   open with the real current data, large, never a button away. For a gene card that is the
   parasite drawn with the organelle the gene is found in filled in -- the schematic Einar
   named by name. It was in the `find` dock and nowhere else, so somebody reading a gene
   could not see where it sat without leaving the card. The drawing is
   :class:`parascape.celldiagram.CellDiagram`, fed from the window's own ``color_of``, so
   the gene card and the map cannot disagree about what a colour means. Clicking an
   organelle ticks that compartment in the window, as it does beside the map.

   A SECOND INSTANCE of the widget, not the window's one moved here. The mould says move
   rather than copy, and that rule is about CONTROLS whose methods must keep reading them.
   The drawing is a stateless SVG view -- its state is the palette and the selection, both
   passed in on every redraw -- and the `find` dock is lent whole to the Cell module, which
   could not then borrow a widget this card had taken.

2. **The lens**, which the mould asks every module for: which measurement is read as "where
   it is". There are six placement columns in the gene table; the card follows the colour-by
   control of the window when that is one of them, and says which it is reading.

3. **What follows a gene** is offered at the foot of the card, as spaCR offers successors
   above its Run button (:mod:`parascape.next_steps`): every module that reads one gene,
   pre-filled with this one, each checked with ``ports.check_ready``, and one that cannot
   work shown greyed with the reason rather than dropped.
"""
from __future__ import annotations

import os
from html import escape

import numpy as np
import pandas as pd
from PySide6 import QtCore, QtWidgets

from . import gene_reading as GR
from . import ports as PORTS
from . import theme as TH
from .recipes_panel import widget_setting
from .shell import ModuleShell, solid_dropdown
from . import organisms


def _muted(text: str) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


def _browser() -> QtWidgets.QTextBrowser:
    view = QtWidgets.QTextBrowser()
    view.setOpenExternalLinks(True)
    return view


class GeneModule(QtCore.QObject):
    """The Gene module for one window."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.gene = ""
        self.genes: list = []                          # a list being stepped through
        self.reading = None
        self.shell = shell = ModuleShell.for_module("gene")
        shell.module = self
        shell.setAcceptDrops(True)
        shell.installEventFilter(self)
        self._build_settings()
        self._build_runtime()
        shell.home_requested.connect(lambda: _go_home(window))
        self._waiting = QtCore.QTimer(self)
        self._waiting.setInterval(300)
        self._waiting.timeout.connect(self._map_ready)

    # ------------------------------------------------------------------ left
    def _build_settings(self) -> None:
        shell = self.shell
        mine = QtWidgets.QWidget(shell)
        v = QtWidgets.QVBoxLayout(mine)
        v.setContentsMargins(6, 2, 6, 6)
        self.field = QtWidgets.QLineEdit()
        self.field.setPlaceholderText("A gene id or a word of its product, then Return")
        self.field.setToolTip(TH.tip("TGME49_208830, or GRA16, or 'rhoptry neck'. Several "
                                     "matches are listed below to pick from."))
        self.field.returnPressed.connect(self._typed)
        v.addWidget(self.field)
        row = QtWidgets.QHBoxLayout()
        choose = QtWidgets.QPushButton("Choose a list…")
        choose.setToolTip(TH.tip("A table of gene ids; step through it with the list below. "
                                 "You can also drop the file on this module."))
        choose.clicked.connect(self._choose_list)
        use = QtWidgets.QPushButton("Use the gene list")
        use.setToolTip(TH.tip("The list last dropped on ParaScape (the gene list panel)."))
        use.clicked.connect(self._use_gene_list)
        clear = QtWidgets.QPushButton("Clear")
        clear.setToolTip(TH.tip("Forget the gene and the list, and empty the card."))
        # Through a lambda: `clicked` carries a checked flag, and handing it straight to a
        # slot that takes no arguments is a call that only fails when somebody presses it.
        clear.clicked.connect(lambda _checked=False: self.clear())
        row.addWidget(choose)
        row.addWidget(use)
        row.addWidget(clear)
        v.addLayout(row)
        self.matches = QtWidgets.QListWidget()
        self.matches.setToolTip("Click a gene to read it. Up and Down step through the list.")
        self.matches.setMaximumHeight(170)
        self.matches.currentItemChanged.connect(
            lambda item, _old: self.show(item.data(QtCore.Qt.ItemDataRole.UserRole)) if item else None)
        v.addWidget(self.matches)
        self.match_note = _muted("")
        v.addWidget(self.match_note)
        shell.add_setting("Your gene", mine, subject=PORTS.GENE)

        lens = QtWidgets.QWidget(shell)
        ll = QtWidgets.QVBoxLayout(lens)
        ll.setContentsMargins(6, 2, 6, 6)
        self.place_box = QtWidgets.QComboBox()
        self.place_box.setToolTip(TH.tip(
            "Which measurement is read as 'where it is' when the parasite is drawn. Several "
            "measurements place a gene in the cell -- the curated compartment, and the "
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
        shell.add_setting("Where it is, read by", lens)
        self._follow_window_colour()

        kinds = QtWidgets.QWidget(shell)
        kl = QtWidgets.QVBoxLayout(kinds)
        kl.setContentsMargins(6, 2, 6, 6)
        kl.addWidget(_muted("What to read about it:"))
        self.include: dict = {}
        for category, title in GR.ORDER:
            box = QtWidgets.QCheckBox(title)
            box.setChecked(category != "sequence")
            box.setToolTip(TH.tip(f"Include {title.lower()} in the summary and the evidence."))
            box.toggled.connect(lambda *_: self.render())
            kl.addWidget(box)
            self.include[category] = box
        shell.add_setting("What to read", kinds)

        advanced = QtWidgets.QWidget(shell)
        al = QtWidgets.QVBoxLayout(advanced)
        al.setContentsMargins(6, 2, 6, 6)
        self.show_predicted = QtWidgets.QCheckBox("Show predicted values (never mixed with facts)")
        self.show_predicted.setChecked(True)
        self.show_predicted.setToolTip(TH.tip("Predictions appear only under 'What can be "
                                              "inferred', with their certainty."))
        self.show_predicted.toggled.connect(lambda *_: self.render())
        al.addWidget(self.show_predicted)
        shell.add_setting("Advanced", advanced, expanded=False)

        keep = QtWidgets.QWidget(shell)
        kp = QtWidgets.QVBoxLayout(keep)
        kp.setContentsMargins(6, 2, 6, 6)
        kp.addWidget(_muted("Keep this reading as files you can hand to somebody else, "
                            "or read one back."))
        save = QtWidgets.QPushButton("Save everything as a bundle…")
        save.setToolTip(TH.tip("One folder: what is measured about the gene, what is unusual, "
                               "your kept claims (marked inferred), the methods text with "
                               "references, and a manifest whose hashes show if anything was "
                               "edited later."))
        save.clicked.connect(lambda: self.save_bundle())
        reopen = QtWidgets.QPushButton("Open a saved bundle…")
        reopen.setToolTip(TH.tip("Read a bundle saved here back: its files are checked against "
                                 "their hashes and its view is restored."))
        reopen.clicked.connect(lambda: self.open_bundle())
        kp.addWidget(save)
        kp.addWidget(reopen)
        shell.add_setting("Save", keep, expanded=False)
        shell.add_recipes(self.recipe_settings(), after_apply=self.render)

    def recipe_settings(self) -> list:
        """The settings of the Gene card that a recipe keeps (instruction 98, WP11).

        THE GENE ITSELF IS NOT ONE. Einar settled that a gene is the SUBJECT of this module
        and not a setting of it, so the field, the list and the matches are left out: a recipe
        says how the card READS a gene, and is then applied to whichever gene is in it. What is
        left is a genuine set of choices -- which measurement is read as "where it is", which
        kinds of knowledge are summarised, and whether predictions are shown beside the facts --
        and before this there was no way to keep a working set of them under a name while
        every other module had one.

        Each kind of knowledge is its own setting, named by its category, so a category a data
        release drops is reported by name when an old recipe is applied rather than silently
        ignored.
        """
        out = [widget_setting("placement_column", self.place_box, "where it is, read by"),
               widget_setting("show_predicted", self.show_predicted,
                              "show predicted values")]
        for category, box in self.include.items():
            out.append(widget_setting(f"read_{category}", box, f"what to read: {category}"))
        return out

    # ------------------------------------------------------------------ right
    def _build_runtime(self) -> None:
        shell, window = self.shell, self.window
        self._build_parasite()
        self.known = _browser()
        shell.link_style(self.known)
        self.known.setToolTip("What is known about the gene, one sentence per kind of "
                              "knowledge, each with its sources.")
        shell.add_section("What is known", self.known, extent=300, stretch=3)

        self.unusual = _browser()
        shell.link_style(self.unusual)
        self.unusual.setToolTip("Of everything measured on this gene, what is furthest from other "
                                "genes: each as a percentile among every gene measured, with its "
                                "source. Study effort and identifiers are left out.")
        shell.add_section("What is unusual about it", self.unusual, extent=170, stretch=1)

        holder = QtWidgets.QWidget(shell)
        hl = QtWidgets.QVBoxLayout(holder)
        hl.setContentsMargins(0, 0, 0, 0)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.evidence_page = QtWidgets.QWidget()
        self.evidence_layout = QtWidgets.QVBoxLayout(self.evidence_page)
        self.evidence_layout.setContentsMargins(8, 4, 8, 8)
        self.evidence_layout.addStretch(1)
        scroll.setWidget(self.evidence_page)
        hl.addWidget(scroll)
        shell.add_section("The evidence", holder, extent=260, stretch=2)

        self.links = _browser()
        shell.link_style(self.links)
        self.links.setOpenLinks(False)
        self.links.anchorClicked.connect(lambda url: self.open_counterpart(url.toString()))
        self.links.setToolTip("Verified links from this gene to genes of other organisms: measured "
                              "interactions and orthology, each kept as its own kind. Click the other "
                              "end to open it in its own organism.")
        shell.add_section("Links to other organisms", self.links, extent=150, stretch=1, folded=True)

        self.inferred = _browser()
        shell.link_style(self.inferred)
        self.inferred.setOpenLinks(False)
        self.inferred.anchorClicked.connect(self._inferred_link)
        self.inferred.setToolTip("What is not known about the gene, and what the data suggest, "
                                 "each with how often such a call is right.")
        shell.add_section("What is not known, and what can be inferred", self.inferred,
                          extent=170, stretch=1)

        self.methods = _browser()
        shell.link_style(self.methods)
        self.methods.setOpenLinks(False)
        self.methods.anchorClicked.connect(lambda url: self.open_strategy(url.toString().split(":", 1)[-1]))
        self.methods.setToolTip("Every method in the catalogue: what it said about this gene when "
                                "the recorded answer for it was hidden, and how its calibration "
                                "graded it.")
        shell.add_section("What each method said about it, held out", self.methods, extent=200,
                          stretch=1, folded=True)

        from .next_steps import NextSteps
        self.next_steps = NextSteps(PORTS.GENE, shell)
        self.next_steps.chosen.connect(lambda key, _kind, gene: self._open(key, gene))
        shell.add_section("What you can do with this gene next", self.next_steps,
                          extent=110, stretch=0)
        self.next_steps.set_subject("", available=self._reasons(), without="gene")

        shell.add_action("Put it on the map", self._on_map, primary=True,
                         tip="Open the Map module with this gene marked and its neighbours read.")
        shell.add_action("Predict what is missing", lambda: self._open("infer"),
                         tip="Answer a question about this gene from the evidence, with how often "
                             "such answers are right (Infer). Infer is asked a question rather "
                             "than a gene, which is why it is here and not among the next steps.")
        shell.add_action("Full evidence panel", self._full_panel,
                         tip="The detailed evidence panel of the main window for this gene.")
        shell.add_action("Previous", lambda: self._step(-1), tip="The previous gene of the list.")
        shell.add_action("Next", lambda: self._step(1), tip="The next gene of the list.")
        shell.borrow("Console and AI chat", [("Console", getattr(window, "console_dock", None)),
                                             ("AI chat", getattr(window, "chat_dock", None))],
                     folded=True)
        self.known.setHtml("<p>Type a gene on the left to read what is known about it.</p>")

    # ------------------------------------------------------- where it is in the parasite
    def _build_parasite(self) -> None:
        """The parasite drawn, with the organelle this gene is found in filled in.

        The thing the card is about, so it is first in the runtime column, large and showing
        the real reading rather than a sample. Einar named the schematic as something the
        restructure must keep; this is where it is kept.
        """
        shell = self.shell
        from .celldiagram import CellDiagram, available as diagram_available
        holder = QtWidgets.QWidget(shell)
        # The caption BESIDE the drawing, not under it. The drawing is portrait, so its
        # width follows its height; under a caption in a column a thousand pixels wide it
        # was a thumbnail with nine hundred pixels of empty space on either side of it.
        layout = QtWidgets.QHBoxLayout(holder)
        layout.setContentsMargins(8, 4, 8, 8)
        layout.setSpacing(16)
        self.parasite = CellDiagram(holder) if diagram_available() else None
        if self.parasite is not None:
            self.parasite.compartment_clicked.connect(self._organelle_clicked)
            # NO minimum height. Asking for 260 px in a section the splitter opened at 265
            # left the layout unable to fit the caption, and Qt laid the caption ON TOP of
            # the drawing -- the parasite was drawn through the sentence describing it.
            self.parasite.setMinimumHeight(0)
            layout.addWidget(self.parasite, 2)
        words = QtWidgets.QWidget(holder)
        wl = QtWidgets.QVBoxLayout(words)
        wl.setContentsMargins(0, 0, 0, 0)
        wl.addStretch(1)
        self.parasite_note = _muted("Type a gene on the left to see where in the cell it is.")
        self.parasite_note.setToolTip(TH.tip(
            "Where the gene is found, which measurement says so, and -- when nothing is "
            "measured -- that it is unknown, which is not the same as the gene being nowhere."))
        wl.addWidget(self.parasite_note)
        self.parasite_trust = _muted("")
        self.parasite_trust.setToolTip(TH.tip(
            "How many genes this measurement covers, and how often a guess with no evidence "
            "would be right -- the number any suggestion has to beat."))
        wl.addWidget(self.parasite_trust)
        wl.addStretch(1)
        layout.addWidget(words, 3)
        # The largest share of the runtime column, which is what "the thing itself, large"
        # means when ten sections compete for the height. The splitter opens each pane at
        # its EXTENT, scaled so they fill the column -- stretch only decides who grows when
        # the window is resized -- so the number that makes the parasite large is this one.
        # At 420 against the 1,650 the other nine ask for, it opened 186 px tall.
        shell.add_section("Where it is in the parasite", holder, extent=900, stretch=6,
                          focus=True)

    def _placement_columns(self) -> list:
        """The measurements of the gene table that place a gene in an organelle the drawing has."""
        from .celldiagram import COMPARTMENT_SL
        try:
            return GR.placement_columns(self._nodes(), COMPARTMENT_SL)
        except Exception:                                  # a card must open on any table
            return []

    def _window_colour(self) -> str:
        """What the window is colouring by right now, or "" when it has no such control."""
        box = getattr(self.window, "category_box", None)
        return box.currentText() if box is not None else ""

    def _follow_window_colour(self) -> None:
        """Keep the lens in step with the colour-by control of the window, both ways.

        One source of truth: the card and the map must never show two different compartments
        for one gene. The window wins when its choice is a placement; when it is something
        else entirely -- a run, a binned number -- the card keeps the placement it had, because
        there is nothing sensible to follow it to.
        """
        box = getattr(self.window, "category_box", None)
        if box is None:
            return
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
            except Exception:                              # recolouring is not worth a crash
                pass
        self._render_parasite()

    def placement(self):
        """Where the current gene is placed, with the chance baseline; None with no gene.

        :returns: a :class:`parascape.gene_reading.Placement`, or None.
        """
        column = self.place_box.currentText()
        if not self.gene or not column:
            return None
        nodes = self._nodes()
        source = "compartment_source" if nodes is not None and \
            "compartment_source" in nodes.columns else ""
        return GR.placement(nodes, self.gene, column, source)

    def _adopt_lens(self, place) -> bool:
        """Colour the window by the lens when it has no colour for what the card is showing.

        ONE PALETTE, which is the rule `celldiagram` is built on: the fills come from
        `window.color_of`, the same dict the points on the map are drawn from, so the
        drawing and the map cannot disagree about what a colour means. The cost of that
        rule is this method: when the window is coloured by something else -- it opens on
        `compartment_best`, whose groupings the artwork has no shape for -- there is no
        colour for "dense granules" and the organelle came out unfilled.

        Rather than mix in a colour of its own, the card moves the window onto the
        measurement it is reading, and says so. Nothing is invented and the two views go on
        meaning the same thing.

        :param place: the reading from :meth:`placement`.
        :returns: whether the window was recoloured.
        """
        colours = getattr(self.window, "color_of", {}) or {}
        box = getattr(self.window, "category_box", None)
        if not place.known or place.compartment in colours or box is None:
            return False
        try:
            if place.column not in self.window.color_sources():
                return False
            box.setCurrentText(place.column)
        except Exception:                                  # colour is not worth a crash
            return False
        self.shell.say(f"The map and the find panel are now coloured by "
                       f"{place.column.replace('_', ' ')}, so the drawing and the map "
                       f"mean the same thing by a colour.")
        return True

    def _render_parasite(self) -> None:
        """Fill the organelle the gene is found in, and say what that rests on."""
        place = self.placement()
        if place is None:
            self.parasite_note.setText("Type a gene on the left to see where in the cell it is.")
            self.parasite_trust.setText("")
            if self.parasite is not None:
                self.parasite.set_palette(getattr(self.window, "color_of", {}) or {}, "")
            return
        self.place_note.setText(f"Reading {place.column.replace('_', ' ')}.")
        self._adopt_lens(place)
        if self.parasite is not None:
            # The window's own palette, never a second one: two palettes are two claims about
            # what a colour means and one of them drifts.
            self.parasite.set_palette(getattr(self.window, "color_of", {}) or {},
                                      place.compartment)
            shared = self.parasite.showing()
        else:
            shared = ""
        note = GR.placement_sentence(place)
        if shared:
            note += f" <span style='color:#999'>({escape(shared)})</span>"
        self.parasite_note.setText(note)
        self.parasite_trust.setText(GR.placement_trust(place))

    def _organelle_clicked(self, name: str) -> None:
        """An organelle clicked: tick that compartment in the window, as the find panel does."""
        picker = getattr(self.window, "select_compartment", None)
        if callable(picker):
            picker(name)
        self.shell.say(f"{name} ticked in the find panel")

    # ------------------------------------------------------------------ choosing a gene
    def _nodes(self):
        return getattr(self.window, "nodes", None)

    def _reasons(self) -> dict:
        """Which modules this window cannot actually open, by the docks it really built.

        Asked of the window rather than of the catalogue so a greyed offer says what is
        true HERE -- a dock that failed to install greys its module with that reason, not
        with a general one.
        """
        try:
            from .landing_install import _reasons
            return _reasons(self.window)
        except Exception:
            return {}

    def clear(self) -> None:
        """Forget the gene and the list, and empty the card.

        The mould asks every module for a way back to nothing. Without it the only way to
        stop looking at a gene was to type another one, which left somebody who wanted a
        clean start reading a card they had finished with.
        """
        self.gene = ""
        self.genes = []
        self.reading = None
        self.field.clear()
        self.matches.blockSignals(True)
        self.matches.clear()
        self.matches.blockSignals(False)
        self.match_note.setText("")
        self.known.setHtml("<p>Type a gene on the left to read what is known about it.</p>")
        for browser in (self.unusual, self.links, self.inferred, self.methods):
            browser.setHtml("")
        self._render_parasite()
        self.next_steps.set_subject("", available=self._reasons(), without="gene")
        self.shell.say("")

    def _product(self, gene: str) -> str:
        nodes = self._nodes()
        if nodes is None or "product" not in nodes:
            return ""
        hit = nodes[nodes["gene_id"].astype(str) == gene]
        return str(hit["product"].iloc[0]) if len(hit) else ""

    def matching(self, text: str) -> list:
        """Gene ids matching typed text: an exact id first, then ids and products containing it."""
        nodes = self._nodes()
        text = (text or "").strip()
        if nodes is None or not text:
            return []
        ids = nodes["gene_id"].astype(str)
        low = text.lower()
        exact = ids[ids.str.lower() == low].tolist()
        contains = ids[ids.str.lower().str.contains(low, regex=False)].tolist()
        if "product" in nodes:
            prod = nodes["product"].astype(str).str.lower()
            contains += ids[prod.str.contains(low, regex=False)].tolist()
        return list(dict.fromkeys(exact + contains))[:200]

    def _fill_matches(self, genes: list, note: str) -> None:
        self.matches.blockSignals(True)
        self.matches.clear()
        for gene in genes:
            item = QtWidgets.QListWidgetItem(f"{gene} — {self._product(gene)[:60]}")
            item.setData(QtCore.Qt.ItemDataRole.UserRole, gene)
            self.matches.addItem(item)
        self.matches.blockSignals(False)
        self.match_note.setText(note)
        if genes:
            self.matches.setCurrentRow(0)

    def _typed(self) -> None:
        genes = self.matching(self.field.text())
        self.genes = genes
        self._fill_matches(genes, f"{len(genes)} matching gene{'s' if len(genes) != 1 else ''}"
                           if genes else "No gene matches that.")

    def load_list(self, path: str) -> None:
        """Read a file of genes and step through them."""
        from . import gene_drop as GD
        nodes = self._nodes()
        known = nodes["gene_id"].astype(str).tolist() if nodes is not None else None
        try:
            drop = GD.read_drop(path, known_ids=known)
        except Exception as error:
            self.shell.say(f"{os.path.basename(path)} could not be read: {error}")
            return
        genes = [g for g in drop.genes if known is None or g in set(known)]
        self.genes = genes
        self._fill_matches(genes, f"{len(genes)} genes from {os.path.basename(path)}")

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
        self.genes = list(genes)
        self._fill_matches(self.genes, f"{len(genes)} genes from the gene list")

    def show_list(self, genes) -> list:
        """Step through a list of genes, starting at the first; returns the ones kept.

        The other half of the port the card has always declared -- "a list to step through,
        one gene at a time". It is how a list carried in from the gene list card arrives
        (instruction 98, WP9): the identifiers that name a gene of this table fill the list
        on the left, the first is read, and Previous and Next walk the rest. An identifier
        that names nothing is left out of the walk and said in the line under the list,
        never silently counted as a gene.

        :param genes: the identifiers, in the order they were given.
        :returns: the ones that are genes of this table.
        """
        from . import list_reading as LR
        resolved = LR.resolve(self._nodes(), list(genes or []))
        self.genes = list(resolved.genes)
        self._fill_matches(self.genes, resolved.sentence())
        return self.genes

    def _step(self, by: int) -> None:
        row = self.matches.currentRow()
        if self.matches.count():
            self.matches.setCurrentRow(max(0, min(self.matches.count() - 1, row + by)))

    def eventFilter(self, obj, event):                    # noqa: N802 - Qt name
        """Accept a file dragged onto the module and step through its genes."""
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
        return False

    # ------------------------------------------------------------------ reading it
    def _published(self):
        panel = getattr(self.window, "information_map", None)
        if panel is None:
            return None
        if panel.map is None:
            panel.load()
        return panel.map

    def _map_ready(self) -> None:
        if self._published() is not None:
            self._waiting.stop()
            if self.gene:
                self.show(self.gene)

    def show(self, gene: str) -> None:
        """Read one gene and show it."""
        if not gene:
            return
        self.gene = gene
        published = self._published()
        if published is None:
            self._waiting.start()
            self.known.setHtml(f"<p>Reading the information map for {escape(gene)}…</p>")
        nodes = self._nodes()
        row = None
        if nodes is not None:
            hit = nodes[nodes["gene_id"].astype(str) == gene]
            row = hit.iloc[0] if len(hit) else None
        self.reading = GR.read(published, gene, self._product(gene), node_row=row)
        self.render()

    def _wanted(self, category: str) -> bool:
        box = self.include.get(category)
        return box is None or box.isChecked()

    def render(self) -> None:
        """Write the synthesis, the evidence sections and the inferences for the current gene."""
        r = self.reading
        if r is None:
            return
        sentences = [(c, s) for c, s in GR.synthesis(r) if self._wanted(c)]
        head = (f"<h2>{escape(r.gene)}</h2><p style='color:#999'>{escape(r.product)}</p>"
                if r.product else f"<h2>{escape(r.gene)}</h2>")
        body = "".join(f"<p>{s}</p>" for _c, s in sentences) or \
            "<p>Nothing recorded about this gene in the kinds of knowledge chosen.</p>"
        self.known.setHtml(head + body)
        self._render_parasite()
        self.next_steps.set_subject(r.gene, available=self._reasons(), without="gene")
        self._render_unusual(r)
        self._render_methods(r)
        self._render_evidence(r)
        self._render_links(r)
        self._render_inferences(r)

    def open_bundle(self, folder: str = "") -> dict:
        """Read a bundle saved by this module back and restore its view (64.39)."""
        from . import export_bundle
        return export_bundle.reopen(self, folder)

    def save_bundle(self, folder: str = "") -> str:
        """Write the current gene as an export bundle; returns the folder."""
        from . import export_bundle
        from . import kept_claims
        r = self.reading
        if r is None:
            self.shell.say("choose a gene first; the bundle holds its reading")
            return ""
        folder = folder or export_bundle.ask_folder(self.shell)
        if not folder:
            return ""
        frames = [f.assign(category=c) for c, f in (r.facts or {}).items() if f is not None and len(f)]
        facts = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        tables, kinds = {"what is measured": facts}, {"what is measured": "measured"}
        unusual = GR.unusual(self._nodes(), r.gene)
        if len(unusual):
            tables["what is unusual"], kinds["what is unusual"] = unusual, "computed"
        if r.predicted is not None and len(r.predicted):
            tables["predicted cells"], kinds["predicted cells"] = r.predicted, "inferred"
        kept = kept_claims.claims(r.gene)
        if len(kept):
            tables["your kept claims"], kinds["your kept claims"] = kept, "inferred"
        columns = [str(c) for c in facts.get("column", pd.Series(dtype=str)).dropna().unique()
                   if c in self._nodes().columns] if len(facts) else []
        export_bundle.write(folder, organism=getattr(self.window.panel, "organism", organisms.TOXOPLASMA),
                            query={"module": "gene", "gene": r.gene, "product": r.product},
                            tables=tables, kinds=kinds, columns=columns, nodes=self._nodes())
        self.shell.say(f"bundle written: {folder}")
        return folder

    def _render_methods(self, r) -> None:
        """Every catalogue method and its held-out calls on this gene (instruction 64.30)."""
        try:
            profile = GR.strategy_profile(r.gene)
        except Exception as error:
            self.methods.setHtml(f"<p>The hold-out record could not be read: {escape(str(error))}</p>")
            return
        spoke = profile[profile["status"] == "held-out calls"]
        head = (f"<p><b>{len(spoke)}</b> of the {len(profile)} methods were tested on this gene with its "
                "answer hidden; the others are listed with why.</p>")
        try:
            stale = GR.stale_records(r.gene, self._nodes())
            if stale:
                head += ("<p><b>Stale:</b> the hold-out record was built on an older answer for this gene -- "
                         + escape("; ".join(stale)) + ". Its right/wrong verdicts on those labels are out of "
                         "date until the record is rebuilt.</p>")
        except Exception:
            pass
        try:
            head += "".join(f"<p>{line}</p>" for line in GR.consensus_lines(r.gene))
        except Exception:
            pass
        try:
            frozen = GR.prospective_lines(r.gene)
            if frozen:
                head += "".join(f"<p>{line}</p>" for line in frozen)
        except Exception:
            pass
        try:
            meta = GR.meta_lines(r.gene)
            if meta:
                head += "<p><i>Combined (a model that learned whom to trust):</i></p>" + \
                    "".join(f"<p>{line}</p>" for line in meta)
        except Exception:
            pass
        rows = "".join(
            f"<tr><td><a href='strategy:{escape(p.key)}'>{escape(p.strategy)}</a></td><td>{escape(p.grade)}</td>"
            f"<td>{escape(p.calls or p.status)}</td></tr>" for p in profile.itertuples(index=False))
        self.methods.setHtml(head + "<table cellspacing=4><tr><th align=left>method</th>"
                             "<th align=left>graded</th><th align=left>held out, it said</th></tr>"
                             + rows + "</table>")

    def open_strategy(self, key: str) -> bool:
        """One click from a method to its card and scorecard: the Strategies module, that strategy shown."""
        panel = getattr(self.window, "strategy_panel", None)
        if panel is None or not key:
            return False
        try:
            panel.select(key)
        except KeyError:
            return False
        from .dock_module import open_in
        open_in(self.window, "strategies")
        return True

    def _render_unusual(self, r) -> None:
        """The gene's most extreme measurements, among the kinds of knowledge chosen."""
        frame = GR.unusual(self._nodes(), r.gene, top=30)
        frame = frame[[self._wanted(c) for c in frame["category"]]].head(10) if len(frame) else frame
        lines = GR.unusual_sentences(frame)
        self.unusual.setHtml("".join(f"<p>{s}</p>" for s in lines) or
                             "<p>Nothing measured on this gene stands out from other genes, or too "
                             "little is measured on it to say.</p>")

    def _render_evidence(self, r) -> None:
        while self.evidence_layout.count() > 1:
            item = self.evidence_layout.takeAt(0)
            if item.widget() is not None:
                item.widget().deleteLater()
        for category, title in GR.ORDER:
            rows = r.facts.get(category)
            if rows is None or not len(rows) or not self._wanted(category):
                continue
            rows = GR.informative(rows) if category not in ("place",) else rows
            if not len(rows):
                continue
            table = QtWidgets.QTableWidget(len(rows), 4)
            table.setHorizontalHeaderLabels(["measurement", "value", "source", "how sure"])
            table.verticalHeader().hide()
            table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
            for i, row in enumerate(rows.itertuples(index=False)):
                source = GR._text(getattr(row, "reference", "")) or \
                    "; ".join(r.sources_for(str(row.column), 2)) or GR._text(getattr(row, "dataset", "")) \
                    or GR._text(row.source_tier)
                p = getattr(row, "probability", float("nan"))
                dataset = GR._text(getattr(row, "dataset", ""))
                what = (f"{dataset} · {str(row.column).replace('_', ' ')}" if dataset
                        else GR._text(row.description) or str(row.column))
                cells = (what,
                         GR.value_text(str(row.column), row.value), source,
                         "" if p != p else f"{p:.0%}")
                for j, text in enumerate(cells):
                    table.setItem(i, j, QtWidgets.QTableWidgetItem(text))
            table.resizeColumnsToContents()
            table.horizontalHeader().setStretchLastSection(True)
            table.setMinimumHeight(min(260, 30 + 24 * len(rows)))
            section = self.shell.section(f"{title} ({len(rows)})", table,
                                         expanded=category in ("place", "interaction", "fitness"))
            self.evidence_layout.insertWidget(self.evidence_layout.count() - 1, section)

    def _render_links(self, r) -> None:
        """The gene's cross-organism links (instruction 64.24), the other end resolved in its own space."""
        from . import bridges
        try:
            links = bridges.links_of(r.gene)
        except Exception as error:
            self.links.setHtml(f"<p>The links could not be read: {escape(str(error))}</p>")
            return
        self._counterparts = {}
        lines = []
        for i, link in enumerate(links.head(40).to_dict("records")):
            c = bridges.counterpart(link, from_gene=r.gene)
            self._counterparts[str(i)] = c
            other = escape(f"{c['organism']} {c['id']}")
            target = f"<a href='open:{i}'>{other}</a>" if c["space"] else other
            lines.append(f"<p>{escape(link['relation'])}: {target} -- {escape(c['status'])}"
                         + (f" <span style='color:#999'>({escape(str(link['reference']))[:80]})</span>"
                            if str(link.get("reference", "")) not in ("", "nan") else "") + "</p>")
        more = f"<p style='color:#999'>and {len(links) - 40} more</p>" if len(links) > 40 else ""
        self.links.setHtml("".join(lines) + more or
                           "<p>No verified link from this gene to another organism.</p>")

    def open_counterpart(self, anchor: str) -> str:
        """Open the other end of a link in its own organism; returns what was opened (or why not).

        A parasite gene opens in that organism map and its Gene module; a host gene is shown from its own
        host gene space, with that space's measurements and gene count.
        """
        c = getattr(self, "_counterparts", {}).get(anchor.split(":", 1)[-1])
        if not c or not c["space"]:
            return ""
        from . import organisms
        space = organisms.SPACES.get(c["space"])
        if space is not None and space.distribution == "wheel":
            window = self.window
            if getattr(window, "species", None) != space.species and hasattr(window, "open_species"):
                window = window.open_species(space.species)
            open_in(window, c["genes"][0])
            return f"{c['space']}:{c['genes'][0]}"
        rows = c["rows"]
        body = "".join(
            "<tr>" + "".join(f"<td>{escape(str(v))}</td>" for v in row) + "</tr>"
            for row in rows.head(5).itertuples(index=False))
        self.links.setHtml(
            f"<p><b>{escape(c['organism'])} {escape(c['id'])}</b> -- {escape(c['status'])}. "
            "<a href='back:'>back to the links</a></p>"
            "<table cellspacing=4><tr>" + "".join(f"<th align=left>{escape(str(k))}</th>" for k in rows.columns)
            + "</tr>" + body + "</table>")
        self.links.anchorClicked.disconnect()
        self.links.anchorClicked.connect(self._back_to_links)
        return f"{c['space']}:{','.join(c['genes'])}"

    def _back_to_links(self, _url) -> None:
        self.links.anchorClicked.disconnect()
        self.links.anchorClicked.connect(lambda url: self.open_counterpart(url.toString()))
        if self.reading is not None:
            self._render_links(self.reading)

    def _render_inferences(self, r) -> None:
        lines = [s for s in GR.inferences(r)
                 if self.show_predicted.isChecked() or not s.startswith("Predicted")]
        suggestion = self._map_suggestion(r.gene)
        if suggestion:
            lines.append(suggestion + " <a href='keep:map'>keep this</a>")
        lines.extend(self._kept_lines(r.gene))
        self.inferred.setHtml("".join(f"<p>{s}</p>" for s in lines)
                              or "<p>Every key trait of this gene is measured.</p>")

    def _kept_lines(self, gene: str) -> list:
        """The user's own kept claims for this gene -- shown as kept, never as facts -- and any
        that a release now measures, said loudly (instruction 78)."""
        from . import kept_claims as KC
        kept = KC.claims(gene)
        if not len(kept):
            return []
        superseded = {}
        published = self._published()
        if published is not None:
            checked = KC.check_against(published)
            superseded = {row.claim_id: row for row in checked.itertuples(index=False)}
        out = []
        for c in kept.itertuples(index=False):
            p = "" if c.probability != c.probability else f", {float(c.probability):.0%}"
            line = (f"<b>Kept by you</b>: {escape(str(c.trait))} = <b>{escape(str(c.value))}</b> "
                    f"({escape(str(c.method))}{p}; release {escape(str(c.release))}, "
                    f"{escape(str(c.kept_at)[:10])}) <a href='forget:{c.claim_id}'>forget</a>")
            hit = superseded.get(c.claim_id)
            if hit is not None:
                line += (f"<br><b style='color:#d9a441'>Now measured</b>: {escape(hit.measured_value)}"
                         f" -- the kept claim {'agrees' if hit.agrees else 'DISAGREES'}; the "
                         f"measurement replaces it.")
            out.append(line)
        return out

    def _inferred_link(self, url) -> None:
        from PySide6 import QtGui
        text = url.toString()
        from . import kept_claims as KC
        if text.startswith("forget:"):
            KC.revert(text.split(":", 1)[1])
            self.render()
        elif text == "keep:map":
            self.keep_map_suggestion()
        else:
            QtGui.QDesktopServices.openUrl(url)

    def keep_map_suggestion(self):
        """Keep what the gene's map neighbours suggest, with its hit rate, in your own layer."""
        from . import kept_claims as KC
        from . import map_reading as MR
        map_module = getattr(getattr(self.window, "map_module", None), "module", None)
        nodes = self._nodes()
        if map_module is None or nodes is None or not self.gene:
            return None
        index = np.flatnonzero(nodes["gene_id"].astype(str).to_numpy() == self.gene)
        if not len(index):
            return None
        hood = MR.neighbours(self.window.xyz, map_module.labels(), [int(index[0])],
                             getattr(self.window, "placed", None))
        if not hood.top:
            return None
        label, n, of = hood.top
        agree = map_module.agreement()
        call = agree.for_label(label) if agree is not None else None
        claim = KC.keep(KC.KeptClaim(
            gene=self.gene, trait=map_module.label_name(), value=label,
            probability=call[0] if call else float("nan"),
            uncertainty=(f"right {call[0]:.0%} of {call[1]} such calls on this map; chance "
                         f"{agree.chance:.0%}") if call else "no record on this map",
            method="map neighbours", method_version=f"k={MR.K}",
            question=f"{map_module.label_name()} from its {of} nearest labelled genes on the map",
            accepted_score=call[0] if call else float("nan")))
        self.shell.say(f"Kept: {self.gene} {claim.trait} = {claim.value}")
        self.render()
        return claim

    def _map_suggestion(self, gene: str) -> str:
        """What the gene's neighbours on the current map suggest, with how often that is right."""
        map_module = getattr(getattr(self.window, "map_module", None), "module", None)
        nodes = self._nodes()
        if map_module is None or nodes is None or getattr(self.window, "xyz", None) is None:
            return ""
        index = np.flatnonzero(nodes["gene_id"].astype(str).to_numpy() == gene)
        if not len(index):
            return ""
        from . import map_reading as MR
        labels = map_module.labels()
        hood = MR.neighbours(self.window.xyz, labels, [int(index[0])],
                             getattr(self.window, "placed", None))
        if not hood.top:
            return ""
        label, n, of = hood.top
        name = map_module.label_name()
        text = (f"On the current map, {n} of its {of} nearest genes of known {escape(name)} are "
                f"<b>{escape(label)}</b>")
        agree = map_module.agreement()
        call = agree.for_label(label) if agree is not None else None
        if call is not None:
            text += (f"; when neighbours say {escape(label)} they are right {call[0]:.0%} of the "
                     f"time on this map ({agree.chance:.0%} by chance)")
        return text + "."

    # ------------------------------------------------------------------ actions
    def _on_map(self) -> None:
        nodes = self._nodes()
        map_shell = getattr(self.window, "map_module", None)
        if map_shell is None or nodes is None or not self.gene:
            return
        index = np.flatnonzero(nodes["gene_id"].astype(str).to_numpy() == self.gene)
        from .map_module import open_in
        open_in(self.window)
        if len(index):
            map_shell.module.select([int(index[0])], "")

    def _open(self, key: str, gene: str = "") -> str:
        """Open another module with this gene already in it; returns what was opened.

        The second half of the offer: a next step that opened an empty module would leave
        the gene to be typed again, which is the thing the chaining strip exists to stop.
        The gene is handed over through `landing_install.open_module`, which is the one
        route every module is opened by.
        """
        from .landing_install import open_module
        gene = gene or self.gene
        return open_module(self.window, key, subject=PORTS.GENE, value=gene)

    def _full_panel(self) -> None:
        from .landing_install import leave_landing
        leave_landing(self.window)
        if self.gene and hasattr(self.window, "_workflow_gene"):
            self.window._workflow_gene(self.gene)
        dock = getattr(self.window, "right_dock", None)
        if dock is not None:
            dock.show()
            dock.raise_()


def build(window) -> ModuleShell | None:
    """The Gene module for `window`."""
    if getattr(window, "nodes", None) is None:
        return None
    return GeneModule(window).shell


def open_in(window, gene: str = "") -> str:
    """Show the Gene module over the window, reading `gene` when given.

    :param window: the main window.
    :param gene: one gene, or SEVERAL separated by commas -- a list, which the card steps
        through one gene at a time. Its own port says so ("a list to step through, one gene
        at a time") and there is a button for taking the gene list whole, but until
        instruction 98 WP9 a list handed in here went into the search field as one string
        and matched nothing: somebody who came from the gene list card and clicked Gene got
        an empty card. See `list_reading.split_identifiers`.
    """
    from .landing_install import enter_landing
    from . import list_reading as LR
    landing = getattr(window, "landing", None)
    shell = getattr(window, "gene_module", None)
    if landing is None or shell is None:
        return "the Gene module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    ids = LR.split_identifiers(gene)
    if len(ids) > 1:
        shell.module.show_list(ids)
    elif gene:
        shell.module.field.setText(gene)
        shell.module._typed()
    return "Gene"


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)
