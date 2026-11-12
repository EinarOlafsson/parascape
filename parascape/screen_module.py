"""The Screen module: your measured table joined to everything already recorded.

WHERE IT CAME FROM (instruction 98, WP10)
=========================================
This was the "Compare a screen" tab of `workflows.py`, the pre-navigator window Einar
landed in. Nothing else in ParaScape calls `prioritization.compare_screen`, so retiring
that window without moving this would have deleted the only way to bring a screen of your
own in. The controls are MOVED here: the same column choosers, the same source hash, the
same comparison and the same export record.

What was wrong before this, and is now right: the Screen tile opened the `compare` dock,
which holds the PUBLISHED expression studies. A tile called "join your screen to
everything already recorded" that opened the published RNA-seq of other people is the
"everything behind the navigator is unchanged" complaint in one click. The Screen module
reads a screen now.

WHAT IT IS FOR
--------------
It lets you take your own table of genes and values and see each row joined to what is
already recorded about that gene -- including the rows that did not join -- so that you
can decide which hits are worth following and carry them on to the rest of ParaScape.

WHAT IT REFUSES TO DO
---------------------
An identifier that did not resolve stays in the table and in the export, as does a gene
with no measurement. Duplicate gene rows are refused rather than averaged, because an
average of two measurements of two different things is neither of them.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from PySide6 import QtCore, QtWidgets

from . import ports as PORTS
from . import theme as TH
from .shell import ModuleShell, solid_dropdown

#: Rows shown in the table. The export carries every one.
SHOWN = 1000
#: The columns of a comparison worth reading, in the order they are read.
COMPARISON_COLUMNS = ("input_gene_id", "gene_id", "mapping_status", "screen_effect",
                      "measurement_status")


def _muted(text: str) -> QtWidgets.QLabel:
    """A dimmed line of explanation, styled by the shell."""
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


def _table():
    """A read-only table: reading a screen must not edit it."""
    table = QtWidgets.QTableWidget()
    table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
    table.setAlternatingRowColors(True)
    table.horizontalHeader().setStretchLastSection(True)
    table.verticalHeader().setVisible(False)
    return table


def fill(table, frame, limit: int = SHOWN) -> None:
    """Show at most `limit` rows of `frame`, leaving the frame itself alone."""
    frame = frame.head(limit)
    table.setSortingEnabled(False)
    table.clear()
    table.setColumnCount(len(frame.columns))
    table.setRowCount(len(frame))
    table.setHorizontalHeaderLabels([str(c) for c in frame.columns])
    for i, row in enumerate(frame.itertuples(index=False, name=None)):
        for j, value in enumerate(row):
            text = "" if value is None or (pd.api.types.is_scalar(value)
                                           and pd.isna(value)) else str(value)
            item = QtWidgets.QTableWidgetItem(text)
            item.setToolTip(text)
            table.setItem(i, j, item)
    table.resizeColumnsToContents()
    table.setSortingEnabled(True)


class ScreenModule(QtCore.QObject):
    """The Screen module for one window: its shell and what each control does."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.nodes = window.nodes
        self.screen: pd.DataFrame | None = None
        self.source: dict = {}
        self.comparison: pd.DataFrame | None = None
        self.shell = shell = ModuleShell(
            "screen", "Screen", "join your screen to everything already recorded",
            instruction="Open your measured table, say which columns it uses, then compare.")
        shell.module = self
        self._build_settings()
        self._build_runtime()
        shell.home_requested.connect(lambda: _go_home(window))

    # ------------------------------------------------------------------ the left column
    def _build_settings(self) -> None:
        """The settings column: the table, which columns it uses, and the export."""
        shell = self.shell

        mine = QtWidgets.QWidget(shell)
        ml = QtWidgets.QVBoxLayout(mine)
        ml.setContentsMargins(6, 2, 6, 6)
        row = QtWidgets.QHBoxLayout()
        open_button = QtWidgets.QPushButton("Open a table…")
        open_button.setToolTip(TH.tip(
            "A gene-level CSV or TSV from spaCR or from anywhere else: one row per gene, "
            "with a column of identifiers and a column of measured values."))
        open_button.clicked.connect(self._choose)
        clear = QtWidgets.QPushButton("Clear")
        clear.setToolTip(TH.tip("Forget the table and the comparison made from it."))
        clear.clicked.connect(lambda _checked=False: self.clear())
        row.addWidget(open_button)
        row.addWidget(clear)
        ml.addLayout(row)
        self.file_note = _muted("No table open. Nothing is read from disk until you open one.")
        ml.addWidget(self.file_note)
        shell.add_setting("Your screen", mine)

        which = QtWidgets.QWidget(shell)
        wl = QtWidgets.QVBoxLayout(which)
        wl.setContentsMargins(6, 2, 6, 6)
        self.gene_column = QtWidgets.QComboBox()
        self.gene_column.setToolTip(TH.tip(
            "The column holding one gene identifier per row. Rows whose identifier does "
            "not resolve stay in the table, marked, rather than being dropped."))
        solid_dropdown(self.gene_column)
        self.value_column = QtWidgets.QComboBox()
        self.value_column.setToolTip(TH.tip(
            "The measured column: an effect size, a count, a score. A missing value stays "
            "missing and never becomes a zero."))
        solid_dropdown(self.value_column)
        form = QtWidgets.QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.addRow("Gene identifier", self.gene_column)
        form.addRow("Measured value", self.value_column)
        holder = QtWidgets.QWidget(which)
        holder.setLayout(form)
        wl.addWidget(holder)
        wl.addWidget(_muted("Both are guessed from the file and both can be changed. The "
                            "guess is a starting point, never a decision."))
        shell.add_setting("Which columns it uses", which)

        advanced = QtWidgets.QWidget(shell)
        al = QtWidgets.QVBoxLayout(advanced)
        al.setContentsMargins(6, 2, 6, 6)
        self.hash_note = _muted("Open a table and its name and content hash appear here.")
        self.hash_note.setToolTip(TH.tip(
            "The SHA-256 of the file as it was read. It travels with every export, so a "
            "comparison can be traced to the exact file it was made from."))
        al.addWidget(self.hash_note)
        shell.add_setting("Advanced", advanced, expanded=False)

        save = QtWidgets.QWidget(shell)
        sl = QtWidgets.QVBoxLayout(save)
        sl.setContentsMargins(6, 2, 6, 6)
        self.export_button = QtWidgets.QPushButton("Save the comparison…")
        self.export_button.setEnabled(False)
        self.export_button.setToolTip(TH.tip(
            "Every row, including the ones that did not resolve and the ones with no "
            "measurement, beside a record of the file and the columns chosen."))
        self.export_button.clicked.connect(self._export)
        sl.addWidget(self.export_button)
        shell.add_setting("Save", save, expanded=False)

    # ----------------------------------------------------------------- the right column
    def _build_runtime(self) -> None:
        """The runtime column: the joined table, how much of it joined, and what is next."""
        shell = self.shell

        joined = QtWidgets.QWidget(shell)
        jl = QtWidgets.QVBoxLayout(joined)
        jl.setContentsMargins(8, 4, 8, 8)
        self.caption = _muted("Open a table on the left. It is shown here as it is read, and "
                              "then again joined to what is recorded once you compare.")
        jl.addWidget(self.caption)
        self.table = _table()
        self.table.setToolTip(TH.tip("Double-click a row to select that gene everywhere else "
                                     "in ParaScape."))
        self.table.cellDoubleClicked.connect(self._row_picked)
        jl.addWidget(self.table, 1)
        shell.add_section("Your screen, joined", joined, extent=400, stretch=3, focus=True)

        self.reading = QtWidgets.QTextBrowser(shell)
        shell.link_style(self.reading)
        self.reading.setHtml("<p>Compare and this says how many rows resolved to a gene here, "
                             "how many carry a usable number, and what became of the rest.</p>")
        shell.add_section("How much of it joined", self.reading, extent=180, stretch=1)

        from .next_steps import NextSteps
        self.next_steps = NextSteps(PORTS.GENE_LIST, shell,
                                    lead="Next, with the genes of your screen:")
        self.next_steps.chosen.connect(lambda key, _kind, _value: self._open(key))
        shell.add_section("What you can do with these genes next", self.next_steps,
                          extent=110, stretch=0)
        self.next_steps.set_subject("", available=self._reasons(), without="screen")

        self.compare_button = shell.add_action(
            "Compare with what is recorded", self.compare, primary=True,
            tip="Join every row to the gene table: which identifiers resolved, which rows "
                "carry a usable number, and what the gene table already holds for each.")
        shell.add_action("Put these genes on the map", self._on_map,
                         tip="Open the Map module with the genes of your screen marked and "
                             "their neighbours read.")
        shell.add_action("Test these genes", self._test,
                         tip="Hand the genes to the strategies and ask what the rest of the "
                             "data would have said about them.")
        shell.add_action("Read a gene in Gene", self._to_gene,
                         tip="Open the selected gene in the Gene module: what is known about "
                             "it, with sources.")
        shell.borrow("Console and AI chat",
                     [("Console", getattr(self.window, "console_dock", None)),
                      ("AI chat", getattr(self.window, "chat_dock", None))],
                     folded=True)

    # ------------------------------------------------------------------- reading a file
    def _choose(self) -> None:
        """Ask for a table and read it, saying in the module why one cannot be read."""
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self.shell, "Your measured gene table", "", "Tables (*.csv *.tsv *.txt)")
        if not path:
            return
        try:
            self.load(path)
        except Exception as error:
            self.shell.say(f"cannot read that table: {error}")

    def load(self, path) -> pd.DataFrame:
        """Read a measured table and keep the hash of the file it came from.

        Public so a test can hand it a file without a dialog, and so that the hash is
        taken at the moment the bytes are read rather than at export, when the file on
        disk may no longer be the one that was compared.
        """
        path = Path(path)
        separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        self.screen = pd.read_csv(path, sep=separator)
        self.source = {"filename": path.name,
                       "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        self.comparison = None
        self.export_button.setEnabled(False)
        for box in (self.gene_column, self.value_column):
            box.clear()
            box.addItems([str(c) for c in self.screen.columns])
        if "gene_id" in self.screen.columns:
            self.gene_column.setCurrentText("gene_id")
        numeric = list(self.screen.select_dtypes(include=np.number))
        if numeric:
            self.value_column.setCurrentText(str(numeric[0]))
        fill(self.table, self.screen)
        self.caption.setText(f"{len(self.screen):,} rows read from {path.name}, as they are "
                             f"in the file. Check the two columns on the left, then compare.")
        self.file_note.setText(f"{path.name} · {len(self.screen):,} rows · "
                               f"{len(self.screen.columns)} columns")
        self.hash_note.setText(f"{path.name}\nSHA-256 {self.source['sha256'][:16]}…")
        self.shell.say(f"{os.path.basename(str(path))} read")
        return self.screen

    def clear(self) -> None:
        """Forget the table and everything read from it."""
        self.screen = None
        self.comparison = None
        self.source = {}
        self.gene_column.clear()
        self.value_column.clear()
        self.table.setRowCount(0)
        self.table.setColumnCount(0)
        self.export_button.setEnabled(False)
        self.caption.setText("Open a table on the left. It is shown here as it is read, and "
                             "then again joined to what is recorded once you compare.")
        self.file_note.setText("No table open. Nothing is read from disk until you open one.")
        self.hash_note.setText("Open a table and its name and content hash appear here.")
        self.reading.setHtml("<p>Compare and this says how many rows resolved to a gene here, "
                             "how many carry a measurement, and what became of the rest.</p>")
        self.next_steps.set_subject("", available=self._reasons(), without="screen")

    # ------------------------------------------------------------------- comparing it
    def compare(self) -> pd.DataFrame | None:
        """Join every row to the gene table, keeping the rows that did not join.

        Returns the comparison, or None with the reason said in the module. Duplicate
        gene rows raise inside `prioritization.compare_screen` and the message is shown
        as it is: an aggregation chosen here would be a decision made on somebody behalf.
        """
        from .prioritization import compare_screen
        if self.screen is None:
            self.shell.say("open a measured table first")
            return None
        gene_column = self.gene_column.currentText()
        value_column = self.value_column.currentText()
        try:
            self.comparison = compare_screen(self.nodes, self.screen, value_column, gene_column)
        except Exception as error:
            self.comparison = None
            self.export_button.setEnabled(False)
            self.shell.say(f"the comparison stopped: {error}")
            self.reading.setHtml(f"<p>The comparison stopped: {error}</p>")
            return None
        self.source.update(gene_column=gene_column, value_column=value_column)
        columns = [c for c in COMPARISON_COLUMNS if c in self.comparison.columns]
        if "product" in self.comparison.columns:
            columns.append("product")
        fill(self.table, self.comparison[columns])
        self.export_button.setEnabled(True)
        self._read_comparison()
        self.next_steps.set_subject(f"{len(self.genes()):,} genes",
                                    available=self._reasons(), without="screen")
        return self.comparison

    def genes(self) -> list:
        """The genes of the screen that resolved to a gene in this table."""
        if self.comparison is None:
            return []
        got = self.comparison[self.comparison.mapping_status.eq("mapped")]
        return [str(g) for g in got.gene_id.dropna().astype(str).unique()]

    def _read_comparison(self) -> None:
        """Say in sentences how much of the screen joined, and what became of the rest."""
        frame = self.comparison
        if frame is None:
            return
        total = len(frame)
        mapped = int(frame.mapping_status.eq("mapped").sum())
        unmapped = total - mapped
        # `measurement_status` is about YOUR column, not about the gene table: it says
        # whether the value in the row is a usable number. Reading it as "the gene is
        # already measured here" would be the wrong sentence over the right number.
        usable = int(frame.measurement_status.eq("measured").sum()) \
            if "measurement_status" in frame.columns else 0
        share = (100.0 * mapped / total) if total else 0.0
        self.caption.setText(
            f"{total:,} rows of your screen, joined to the gene table. Showing at most "
            f"{SHOWN:,}; the saved comparison has every row, the unresolved ones included.")
        self.reading.setHtml(
            f"<p><b>{mapped:,} of {total:,} rows</b> resolved to a gene here "
            f"({share:.1f} %).</p>"
            f"<ul><li>{unmapped:,} did not resolve. They stay in the table and in the "
            f"export, marked, because a row dropped is a row nobody knows to look for.</li>"
            f"<li>{usable:,} rows carry a usable number in your measured column; "
            f"{total - usable:,} were blank or were not a number, and stay missing.</li>"
            f"<li>{len(self.genes()):,} distinct genes, which is the list the next steps "
            f"below are offered on.</li></ul>"
            f"<p>A value in the gene table is often an aggregate: it does not recover the "
            f"replicate, the time point or the uncertainty when those were not kept.</p>")

    # ------------------------------------------------------------------------- going on
    def _row_picked(self, row: int, _column: int) -> None:
        """Select the gene of a clicked row everywhere else in ParaScape."""
        column = 1 if self.comparison is not None else 0
        item = self.table.item(row, column) or self.table.item(row, 0)
        if item is not None and item.text() and hasattr(self.window, "_workflow_gene"):
            self.window._workflow_gene(item.text())
            self.shell.say(f"{item.text()} selected")

    def _to_gene(self) -> None:
        """Open the selected gene in the Gene module."""
        item = self.table.currentItem()
        text = item.text() if item is not None else ""
        from .gene_module import open_in
        open_in(self.window, text)

    def _on_map(self) -> None:
        """Open the Map module with the genes of the screen marked."""
        genes = self.genes()
        if not genes:
            self.shell.say("compare first: the map is given the genes that resolved")
            return
        from .map_module import open_in
        open_in(self.window)
        shell = getattr(self.window, "map_module", None)
        if shell is not None and getattr(shell, "module", None) is not None:
            index = self.nodes.index[self.nodes.gene_id.astype(str).isin(set(genes))]
            shell.module.select([int(i) for i in index], f"{len(genes):,} genes of your screen")

    def _test(self) -> None:
        """Ask the strategies what the rest of the data would have said about these genes."""
        genes = self.genes()
        if not genes:
            self.shell.say("compare first: the strategies are given the genes that resolved")
            return
        run = getattr(self.window, "run_my_list", None)
        if callable(run):
            run(genes)
            self.shell.say(f"asking the strategies about {len(genes):,} genes")
        else:
            self.shell.say("the strategies are not available in this window")

    def _open(self, key: str) -> str:
        """Open another module, which is how the screen becomes the next step."""
        from .landing_install import open_module
        return open_module(self.window, key)

    def _reasons(self) -> dict:
        """Which modules this window cannot open, so a blocked offer says why."""
        from .landing_install import _reasons
        try:
            return _reasons(self.window)
        except Exception:
            return {}

    def _export(self) -> None:
        """Write the comparison and the record of where it came from, beside each other."""
        if self.comparison is None:
            self.shell.say("nothing to save: compare first")
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self.shell, "Save the comparison", "screen_comparison.csv", "CSV (*.csv)")
        if not path:
            return
        try:
            self.comparison.to_csv(path, index=False)
            Path(path).with_suffix(".json").write_text(
                json.dumps(self.source, indent=2) + "\n")
            self.shell.say(f"saved the comparison and its source record to {path}")
        except (OSError, ValueError) as error:
            self.shell.say(f"could not save: {error}")


def build(window) -> ModuleShell | None:
    """The Screen module for `window`, or None when there is no gene table to join to."""
    if getattr(window, "nodes", None) is None:
        return None
    return ScreenModule(window).shell


def open_in(window) -> str:
    """Show the Screen module over the window. Returns what the status bar was told."""
    from .landing_install import enter_landing
    landing = getattr(window, "landing", None)
    shell = getattr(window, "screen_module", None)
    if landing is None or shell is None:
        return "the Screen module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    return "Screen"


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)
