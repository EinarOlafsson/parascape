"""The Reach module: which cluster a treatment reaches, and which it could not have been seen in.

WHAT IT IS FOR
==============
A single-cell experiment gives clusters. A bulk experiment -- a drug, a knockout, a pH switch, a
time point -- gives one number per gene. Nothing joins them except the gene space, and the question
they answer together is the one the field needs answering: *if a compound only moves the genes of
one bradyzoite subtype, it is only going to clear that subtype.* That analysis was done by hand for
the bradyzoite subtypes of the mouse brain against the stored drug and stress contrasts. This module
is that analysis, for any cluster block and any contrast, so it can be done again by somebody who
did not write the script.

WHAT IT REFUSES TO DO
---------------------
Three refusals, each of which is the whole point and each of which lives in
:mod:`parascape.cluster_axis` so it cannot be clicked away:

* the baseline is a matched random gene set, never the number of genes in the cluster;
* a cluster the contrast barely measured is reported UNKNOWN and never as no effect;
* every tested cluster carries the smallest effect its own comparison could have told from its
  own baseline, so "nothing there" and "could not have seen it" read differently.

A fourth refusal is reported rather than enforced: where the matched sets cannot be centred on the
cluster in expression -- a cluster of the brightest genes in the genome has nothing brighter to be
compared with -- the module says so and declines to call it an effect.

WHERE THE TWO SIDES COME FROM
-----------------------------
The clusters come from a declared dataset whose columns are one value per group
(:func:`parascape.cluster_axis.cluster_blocks`), or from a table of your own with one column per
cluster. The contrasts come from the expression traits release: a shrunken log2 fold change per
gene computed inside one study against that study's own control, which is the only form in which
one study's treatment can be read beside another's.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from PySide6 import QtCore, QtWidgets

from . import cluster_axis as CA
from . import ports as PORTS
from . import theme as TH
from .recipes_panel import widget_setting
from .shell import ModuleShell, solid_dropdown

#: The columns of a result worth reading, in the order they are read.
SHOWN_COLUMNS = ("n_markers", "n_measured", "observed", "null_mean", "detection_limit",
                 "over_detection_limit", "q", "verdict")

#: What a row of the result says, in one word each, so the table does not need a legend.
REACHES = "reaches it"
NO_EFFECT = "no effect above what could be seen"
UNKNOWN = "n.m."
UNMATCHED = "cannot be told from expression"


def _muted(text: str) -> QtWidgets.QLabel:
    """A dimmed line of explanation, styled by the shell."""
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


def _table() -> QtWidgets.QTableWidget:
    """A read-only table: a result is read, not edited."""
    table = QtWidgets.QTableWidget()
    table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
    table.setAlternatingRowColors(True)
    table.horizontalHeader().setStretchLastSection(True)
    table.verticalHeader().setVisible(False)
    return table


def verdicts(result: CA.AxisResult) -> pd.Series:
    """One word per cluster for what the comparison concluded.

    Four outcomes and not two, because the three the engine separates are exactly the three a
    reader collapses: reached, measured and unmoved, not measured, and measured but not
    separable from expression.
    """
    table = result.table
    out = []
    for name, row in table.iterrows():
        if not bool(row["testable"]):
            out.append(UNKNOWN)
        elif not bool(row.get("expression_matched", True)):
            out.append(UNMATCHED)
        elif name in result.reached:
            out.append(REACHES)
        else:
            out.append(NO_EFFECT)
    return pd.Series(out, index=table.index, name="verdict")


def as_shown(result: CA.AxisResult) -> pd.DataFrame:
    """The result as the table shows it: the cluster, the numbers that matter, the verdict."""
    table = result.table.copy()
    table["verdict"] = verdicts(result)
    columns = [c for c in SHOWN_COLUMNS if c in table.columns]
    out = table[columns].reset_index().rename(columns={"index": "cluster"})
    return out


def fill(table: QtWidgets.QTableWidget, frame: pd.DataFrame) -> None:
    """Show `frame`, leaving the frame itself alone. A missing number stays blank, never a zero."""
    table.setSortingEnabled(False)
    table.clear()
    table.setColumnCount(len(frame.columns))
    table.setRowCount(len(frame))
    table.setHorizontalHeaderLabels([str(c) for c in frame.columns])
    for i, row in enumerate(frame.itertuples(index=False, name=None)):
        for j, value in enumerate(row):
            if value is None or (pd.api.types.is_scalar(value) and pd.isna(value)):
                text = ""
            elif isinstance(value, float):
                text = f"{value:.3g}"
            else:
                text = str(value)
            item = QtWidgets.QTableWidgetItem(text)
            item.setToolTip(text)
            table.setItem(i, j, item)
    table.resizeColumnsToContents()
    table.setSortingEnabled(True)


def reading(result: CA.AxisResult) -> str:
    """What the result does and does not say, in prose, for the lower section."""
    table = result.table
    reached, unknown = result.reached, result.unknown
    lines = [f"<p><b>{result.axis_name}</b>, against {len(table)} clusters, "
             f"{result.n_genes_in_axis:,} genes measured, {result.n_null:,} matched random sets "
             f"behind each cluster.</p>"]
    if reached:
        rows = []
        for name in reached:
            row = table.loc[name]
            rows.append(f"<li><b>{name}</b> &mdash; {row['over_detection_limit']:.1f} times the "
                        f"smallest effect this comparison could have seen, q = {row['q']:.3g}, "
                        f"on {int(row['n_measured'])} of {int(row['n_markers'])} markers</li>")
        lines.append("<p>Reached:</p><ul>" + "".join(rows) + "</ul>")
    else:
        lines.append("<p>No cluster is moved more than genes like its markers are moved anyway. "
                     "That is a result about this contrast, not about the clusters.</p>")
    untouched = [c for c in table.index if c not in reached and c not in unknown]
    if untouched:
        limits = ", ".join(f"{c} (could have seen {table.loc[c, 'detection_limit']:.2f})"
                           for c in untouched)
        lines.append(f"<p>Measured and not moved, with the smallest effect each comparison could "
                     f"have told from its own baseline: {limits}.</p>")
    if unknown:
        lines.append(f"<p>Not measured well enough to test, and therefore <b>{UNKNOWN}</b> rather "
                     f"than no effect: {', '.join(unknown)}.</p>")
    for note in result.notes:
        lines.append(f"<p>{note[0].upper() + note[1:]}.</p>")
    lines.append("<p>The baseline is a random gene set matched to each cluster in size and in "
                 "expression, not the number of genes in the cluster. Matching on expression is "
                 "what stops a cluster of highly expressed markers being scored against genes "
                 "nobody could detect.</p>")
    return "".join(lines)


class ReachModule(QtCore.QObject):
    """The Reach module for one window: its shell and what each control does."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.nodes = window.nodes
        self.blocks = CA.cluster_blocks(self.nodes)
        self.catalogue = CA.axis_catalogue()
        self.mine: pd.DataFrame | None = None
        self.source: dict = {}
        self.result: CA.AxisResult | None = None
        self.shell = shell = ModuleShell.for_module("reach")
        shell.module = self
        self._build_settings()
        self._build_runtime()
        self._fill_axes()
        shell.home_requested.connect(lambda: _go_home(window))

    # ------------------------------------------------------------------ the left column
    def _build_settings(self) -> None:
        """The settings column: where the clusters come from, which contrast, what it refuses."""
        shell = self.shell

        treatment = QtWidgets.QWidget(shell)
        tl = QtWidgets.QVBoxLayout(treatment)
        tl.setContentsMargins(6, 2, 6, 6)
        self.study = QtWidgets.QComboBox()
        self.study.setToolTip(TH.tip(
            "Narrow the contrasts to one study. Every fold change is computed inside one study "
            "against that study's own control, never pooled across studies."))
        solid_dropdown(self.study)
        self.study.currentIndexChanged.connect(self._fill_contrasts)
        self.contrast = QtWidgets.QComboBox()
        self.contrast.setToolTip(TH.tip(
            "The treatment, as one log2 fold change per gene. A drug, a knockout, a pH switch "
            "and a time point are the same shape of measurement here."))
        solid_dropdown(self.contrast)
        self.contrast.currentIndexChanged.connect(self._describe_axis)
        form = QtWidgets.QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.addRow("Study", self.study)
        form.addRow("Contrast", self.contrast)
        holder = QtWidgets.QWidget(treatment)
        holder.setLayout(form)
        tl.addWidget(holder)
        self.axis_note = _muted("Choose a contrast to see how many genes it measures.")
        tl.addWidget(self.axis_note)
        shell.add_setting("The treatment", treatment, subject=PORTS.TRAIT)

        clusters = QtWidgets.QWidget(shell)
        cl = QtWidgets.QVBoxLayout(clusters)
        cl.setContentsMargins(6, 2, 6, 6)
        self.block = QtWidgets.QComboBox()
        self.block.setToolTip(TH.tip(
            "A dataset whose columns are one value per group: the five bradyzoite subtypes of a "
            "brain cyst, for instance. Or a table of your own with one column per cluster."))
        solid_dropdown(self.block)
        for key, found in self.blocks.items():
            self.block.addItem(f"{found['name']} ({len(found['columns'])})", key)
        self.block.addItem("A table of my own…", "")
        self.block.currentIndexChanged.connect(self._describe_block)
        cl.addWidget(self.block)
        open_button = QtWidgets.QPushButton("Open a table…")
        open_button.setToolTip(TH.tip(
            "A CSV or TSV with a column of gene identifiers and one numeric column per cluster: "
            "the average expression per cluster that a single-cell paper publishes."))
        open_button.clicked.connect(self._choose)
        cl.addWidget(open_button)
        self.block_note = _muted("")
        cl.addWidget(self.block_note)
        shell.add_setting("Where the clusters come from", clusters)

        marking = QtWidgets.QWidget(shell)
        ml = QtWidgets.QVBoxLayout(marking)
        ml.setContentsMargins(6, 2, 6, 6)
        self.top = QtWidgets.QSpinBox()
        self.top.setRange(10, 2000)
        self.top.setValue(CA.TOP_MARKERS)
        self.top.setToolTip(TH.tip(
            "How many markers a cluster keeps, most specific first. A long list of weak markers "
            "is dominated by genes that are merely slightly higher in one cluster."))
        self.margin = QtWidgets.QDoubleSpinBox()
        self.margin.setRange(1.0, 10.0)
        self.margin.setSingleStep(0.05)
        self.margin.setValue(CA.MIN_MARGIN)
        self.margin.setToolTip(TH.tip(
            "How much higher than the runner-up cluster a gene must be to mark its own. At 1.0 "
            "every gene marks whichever cluster happens to be highest, which marks nothing."))
        form = QtWidgets.QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.addRow("Markers per cluster", self.top)
        form.addRow("Over the runner-up", self.margin)
        holder = QtWidgets.QWidget(marking)
        holder.setLayout(form)
        ml.addWidget(holder)
        ml.addWidget(_muted("A gene marks the cluster it is highest in, and only when it clears "
                            "the runner-up by this much. A gene expressed evenly marks nothing."))
        shell.add_setting("How a marker is chosen", marking)

        refuses = QtWidgets.QWidget(shell)
        rl = QtWidgets.QVBoxLayout(refuses)
        rl.setContentsMargins(6, 2, 6, 6)
        self.min_measured = QtWidgets.QSpinBox()
        self.min_measured.setRange(3, 500)
        self.min_measured.setValue(CA.MIN_MEASURED)
        self.min_measured.setToolTip(TH.tip(
            "Below this many markers measured by the contrast, the cluster is reported unknown "
            "and no test is run. An empty cell that means nobody looked must not read as an "
            "empty cell that means nothing there."))
        self.n_null = QtWidgets.QSpinBox()
        self.n_null.setRange(100, 20000)
        self.n_null.setSingleStep(100)
        self.n_null.setValue(CA.N_NULL)
        self.n_null.setToolTip(TH.tip(
            "How many matched random gene sets stand behind each cluster. This is the baseline, "
            "and it is never the number of genes in the cluster."))
        self.signed = QtWidgets.QCheckBox("Ask about direction instead of size")
        self.signed.setToolTip(TH.tip(
            "Off, a treatment that pushes half a programme up and half down has still reached "
            "it. On, the question becomes whether the cluster moves up or down on balance."))
        form = QtWidgets.QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.addRow("Markers measured, at least", self.min_measured)
        form.addRow("Matched random sets", self.n_null)
        holder = QtWidgets.QWidget(refuses)
        holder.setLayout(form)
        rl.addWidget(holder)
        rl.addWidget(self.signed)
        shell.add_setting("What it refuses to call an effect", refuses, expanded=False)

        save = QtWidgets.QWidget(shell)
        sl = QtWidgets.QVBoxLayout(save)
        sl.setContentsMargins(6, 2, 6, 6)
        self.export_button = QtWidgets.QPushButton("Save the result…")
        self.export_button.setEnabled(False)
        self.export_button.setToolTip(TH.tip(
            "Every cluster, including the ones reported unknown, beside a record of the contrast, "
            "the cluster block and the settings the result was made with."))
        self.export_button.clicked.connect(self._export)
        sl.addWidget(self.export_button)
        shell.add_setting("Save", save, expanded=False)
        shell.add_recipes(self.recipe_settings())
        self._describe_block()

    def recipe_settings(self) -> list:
        """The settings a recipe keeps. The contrast is the SUBJECT and is not one of them."""
        return [
            widget_setting("block", self.block, "where the clusters come from"),
            widget_setting("top", self.top, "markers per cluster"),
            widget_setting("margin", self.margin, "margin over the runner-up"),
            widget_setting("min_measured", self.min_measured, "markers measured, at least"),
            widget_setting("n_null", self.n_null, "matched random sets"),
            widget_setting("signed", self.signed, "ask about direction instead of size"),
        ]

    # ----------------------------------------------------------------- the right column
    def _build_runtime(self) -> None:
        """The runtime column: the verdict per cluster, what it does not say, and what is next."""
        shell = self.shell

        held = QtWidgets.QWidget(shell)
        hl = QtWidgets.QVBoxLayout(held)
        hl.setContentsMargins(8, 4, 8, 8)
        self.caption = _muted("Choose a contrast and a set of clusters, then test. One row per "
                              "cluster, including the ones nothing can be said about.")
        hl.addWidget(self.caption)
        self.table = _table()
        self.table.setToolTip(TH.tip("Double-click a cluster to carry its markers on as a gene "
                                     "list."))
        self.table.cellDoubleClicked.connect(self._row_picked)
        hl.addWidget(self.table, 1)
        shell.add_section("Which cluster it reaches", held, extent=360, stretch=3, focus=True)

        self.prose = QtWidgets.QTextBrowser(shell)
        shell.link_style(self.prose)
        self.prose.setHtml("<p>This says what the result does and does not support: what was "
                           "reached, what was measured and not moved, with the smallest effect "
                           "each comparison could have seen, and what was not measured at all.</p>")
        shell.add_section("What this does and does not say", self.prose, extent=200, stretch=2)

        from .next_steps import NextSteps
        self.next_steps = NextSteps(PORTS.GENE_LIST, shell,
                                    lead="Next, with the markers of the chosen cluster:")
        self.next_steps.chosen.connect(lambda key, _kind, _value: self._open(key))
        shell.add_section("What you can do with these genes next", self.next_steps,
                          extent=110, stretch=0)
        self.next_steps.set_subject("", available=self._reasons(), without="cluster")

        self.test_button = shell.add_action(
            "Test every cluster", self.run, primary=True,
            tip="Mark each cluster, then ask whether the contrast moves its markers more than "
                "matched random genes. Each cluster is corrected across the clusters tested.")
        shell.add_action("Put the markers on the map", self._on_map,
                         tip="Open the Map module with the markers of the chosen cluster marked.")
        shell.add_action("Read a gene in Gene", self._to_gene,
                         tip="Open the selected gene in the Gene module: what is known about it, "
                             "with sources.")
        shell.borrow("Console and AI chat",
                     [("Console", getattr(self.window, "console_dock", None)),
                      ("AI chat", getattr(self.window, "chat_dock", None))],
                     folded=True)

    # ------------------------------------------------------------------ what is on offer
    def _fill_axes(self) -> None:
        """Fill the study and contrast choosers from the stored contrasts."""
        self.study.blockSignals(True)
        self.study.clear()
        if self.catalogue.empty:
            self.study.addItem("no stored contrasts in this release", "")
            self.study.setEnabled(False)
            self.contrast.setEnabled(False)
            self.axis_note.setText("No expression traits release is present, so there is no "
                                   "contrast to test against. Build one, or bring your own axis.")
            self.study.blockSignals(False)
            return
        self.study.addItem(f"every study ({len(self.catalogue)} contrasts)", "")
        for name, rows in sorted(self.catalogue.groupby("study")):
            self.study.addItem(f"{name} ({len(rows)})", str(name))
        self.study.blockSignals(False)
        self._fill_contrasts()

    def _fill_contrasts(self) -> None:
        """The contrasts of the chosen study, longest-measured first."""
        if self.catalogue.empty:
            return
        wanted = self.study.currentData()
        rows = self.catalogue
        if wanted:
            rows = rows[rows["study"].astype(str) == wanted]
        rows = rows.sort_values("n_genes", ascending=False)
        self.contrast.blockSignals(True)
        self.contrast.clear()
        for column, row in rows.iterrows():
            label = str(row.get("contrast_label") or row.get("contrast") or column)
            self.contrast.addItem(f"{label}  [{row['study']}]", str(column))
        self.contrast.blockSignals(False)
        self._describe_axis()

    def _describe_axis(self) -> None:
        """How many genes the chosen contrast measures, which is its own limit on the answer."""
        column = self.contrast.currentData()
        if not column or self.catalogue.empty or column not in self.catalogue.index:
            self.axis_note.setText("Choose a contrast.")
            return
        row = self.catalogue.loc[column]
        self.axis_note.setText(
            f"{int(row['n_genes']):,} genes measured, from {row['release']}. A cluster whose "
            f"markers are not among them cannot be tested, and is reported {UNKNOWN}.")

    def _describe_block(self) -> None:
        """Which columns the chosen cluster block uses, named rather than implied."""
        key = self.block.currentData()
        if key:
            found = self.blocks[key]
            names = ", ".join(found["columns"].values())
            self.block_note.setText(f"{len(found['columns'])} clusters, from {names}. "
                                    f"{found['provides']}.")
        elif self.mine is not None:
            self.block_note.setText(
                f"{self.source.get('name', 'your table')}: {len(self.mine):,} genes, "
                f"{len(self.source.get('clusters', ()))} clusters "
                f"({', '.join(self.source.get('clusters', ()))}).")
        else:
            self.block_note.setText("Open a table with one column of gene identifiers and one "
                                    "numeric column per cluster.")

    # ------------------------------------------------------------------- reading a table
    def _choose(self) -> None:
        """Ask for a cluster table and read it, saying in the module why one cannot be read."""
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self.shell, "One column per cluster", "", "Tables (*.csv *.tsv *.txt)")
        if not path:
            return
        try:
            self.load(path)
        except (OSError, ValueError, KeyError) as error:
            self.shell.say(f"could not read that table: {error}")

    def load(self, path: str) -> pd.DataFrame:
        """Read a table of one average expression per cluster. Raises rather than guessing."""
        sep = "\t" if str(path).lower().endswith((".tsv", ".txt")) else ","
        frame = pd.read_csv(path, sep=sep)
        genes = next((c for c in frame.columns if str(c).lower() in
                      ("gene_id", "gene", "geneid", "gene id", "id")), None)
        if genes is None:
            raise ValueError("no column of gene identifiers; name one gene_id")
        numeric = [c for c in frame.columns
                   if c != genes and pd.api.types.is_numeric_dtype(frame[c])]
        if len(numeric) < 2:
            raise ValueError(f"{len(numeric)} numeric columns; two or more clusters are needed "
                             f"because a margin needs a runner-up")
        self.mine = frame.rename(columns={genes: "gene_id"})
        self.source = {"name": Path(path).name, "path": str(path), "clusters": tuple(map(str, numeric))}
        self.block.setCurrentIndex(self.block.count() - 1)
        self._describe_block()
        self.shell.say(f"read {len(frame):,} genes and {len(numeric)} clusters")
        return self.mine

    # ------------------------------------------------------------------------- the answer
    def frame_and_columns(self):
        """The gene table and the cluster columns to mark from, whichever source is chosen."""
        key = self.block.currentData()
        if key:
            return self.nodes, dict(self.blocks[key]["columns"])
        if self.mine is None:
            raise ValueError("open a table of your own, or choose a dataset block")
        return self.mine, {c: c for c in self.source["clusters"]}

    def markers(self) -> dict:
        """The marker set of each cluster, by the rule the settings name."""
        frame, columns = self.frame_and_columns()
        return CA.markers_from_blocks(frame, columns, top=self.top.value(),
                                      min_margin=self.margin.value())

    def run(self) -> CA.AxisResult | None:
        """Mark the clusters, test each against the chosen contrast, and show the result."""
        column = self.contrast.currentData()
        if not column:
            self.shell.say("choose a contrast: there is nothing to test the clusters against")
            return None
        try:
            frame, columns = self.frame_and_columns()
            marks = CA.markers_from_blocks(frame, columns, top=self.top.value(),
                                           min_margin=self.margin.value())
            axis = CA.axis_values(str(column))
        except (ValueError, KeyError, OSError) as error:
            self.shell.say(f"cannot test: {error}")
            return None
        empty = [name for name, genes in marks.items() if not genes]
        if len(empty) == len(marks):
            self.shell.say("no cluster has a gene of its own at this margin; lower it")
            return None
        # The expression used for matching is the cluster block itself, averaged over the
        # clusters: it is what the single-cell experiment measured, and matching on anything
        # else would match on a quantity the clusters were not defined by.
        expression = frame[list(columns.values())].apply(
            pd.to_numeric, errors="coerce").mean(axis=1)
        expression.index = frame["gene_id"].astype(str)
        label = self.contrast.currentText().split("  [")[0]
        self.result = CA.test_clusters_against_axis(
            {k: v for k, v in marks.items() if v}, axis, expression.dropna(),
            axis_name=label, n_null=self.n_null.value(),
            min_measured=self.min_measured.value(), absolute=not self.signed.isChecked())
        self.marks = marks
        fill(self.table, as_shown(self.result))
        self.prose.setHtml(reading(self.result))
        reached = self.result.reached
        self.caption.setText(
            f"{len(self.result.table)} clusters tested against {label}; "
            f"{len(reached)} reached" + (f": {', '.join(reached)}" if reached else ""))
        self.export_button.setEnabled(True)
        if empty:
            self.shell.say(f"no gene marks {', '.join(empty)} at this margin, so "
                           f"{'it is' if len(empty) == 1 else 'they are'} not tested")
        else:
            self.shell.say(self.result.summary())
        self._offer_the_chosen_cluster()
        return self.result

    def chosen_cluster(self) -> str:
        """The cluster whose row is selected, or the first one reached, or the first row."""
        if self.result is None:
            return ""
        row = self.table.currentRow()
        if row >= 0 and self.table.item(row, 0) is not None:
            return self.table.item(row, 0).text()
        reached = self.result.reached
        return reached[0] if reached else str(self.result.table.index[0])

    def genes(self) -> list:
        """The markers of the chosen cluster, which is what the next steps act on."""
        return list(getattr(self, "marks", {}).get(self.chosen_cluster(), ()))

    def _offer_the_chosen_cluster(self) -> None:
        # `set_subject` writes the lead itself, from the subject and how many offers can run;
        # assigning to `next_steps.lead` would overwrite the LABEL with a string, which is how
        # this first failed.
        self.next_steps.set_subject(self.chosen_cluster(), available=self._reasons(),
                                    without="cluster")

    # --------------------------------------------------------------------- what comes next
    def _row_picked(self, row: int, _column: int) -> None:
        self._offer_the_chosen_cluster()
        item = self.table.item(row, 0)
        if item is not None:
            self.shell.say(f"{item.text()}: {len(self.genes())} markers carried on")

    def _to_gene(self) -> None:
        """Open the first marker of the chosen cluster in the Gene module."""
        genes = self.genes()
        if not genes:
            self.shell.say("test first: the gene card is given a marker of a tested cluster")
            return
        from .gene_module import open_in
        open_in(self.window, genes[0])

    def _on_map(self) -> None:
        """Open the Map module with the markers of the chosen cluster marked."""
        genes = self.genes()
        if not genes:
            self.shell.say("test first: the map is given the markers of a tested cluster")
            return
        from .map_module import open_in
        open_in(self.window)
        shell = getattr(self.window, "map_module", None)
        if shell is not None and getattr(shell, "module", None) is not None:
            index = self.nodes.index[self.nodes.gene_id.astype(str).isin(set(genes))]
            shell.module.select([int(i) for i in index],
                               f"{len(genes):,} markers of {self.chosen_cluster()}")

    def _open(self, key: str) -> str:
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
        """Write the result and the record of what it was made from, beside each other."""
        if self.result is None:
            self.shell.say("nothing to save: test first")
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self.shell, "Save the result", "cluster_reach.csv", "CSV (*.csv)")
        if not path:
            return
        record = {"contrast": self.contrast.currentData(),
                  "contrast_label": self.contrast.currentText(),
                  "clusters_from": self.block.currentData() or self.source,
                  "markers_per_cluster": self.top.value(),
                  "margin_over_runner_up": self.margin.value(),
                  "markers_measured_at_least": self.min_measured.value(),
                  "matched_random_sets": self.n_null.value(),
                  "asked_about_direction": self.signed.isChecked(),
                  "genes_in_the_contrast": self.result.n_genes_in_axis,
                  "notes": list(self.result.notes)}
        try:
            self.result.table.assign(verdict=verdicts(self.result)).to_csv(path)
            Path(path).with_suffix(".json").write_text(json.dumps(record, indent=2) + "\n")
            self.shell.say(f"saved the result and what it was made from to {path}")
        except (OSError, ValueError) as error:
            self.shell.say(f"could not save: {error}")


def build(window) -> ModuleShell | None:
    """The Reach module for `window`, or None when there is no gene table to mark clusters in."""
    if getattr(window, "nodes", None) is None:
        return None
    return ReachModule(window).shell


def open_in(window) -> str:
    """Show the Reach module over the window. Returns what the status bar was told."""
    from .landing_install import enter_landing
    landing = getattr(window, "landing", None)
    shell = getattr(window, "reach_module", None)
    if landing is None or shell is None:
        return "the Reach module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    return "Reach"


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)
