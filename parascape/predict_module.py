"""The Predict module: predict any measured column, and say how often that is right.

WHERE IT CAME FROM, AND WHY IT IS A MODULE NOW (instruction 98, WP10)
=====================================================================
This was the "Predict a trait" tab of `workflows.py` -- the window Einar landed in when
he typed a gene and clicked the module offered: *"now when i put in a gene and pressed
the modual that came up for the gene, i was brought to the old version of parascape."*
WP10 retires that window. It could not simply be deleted: nothing else in ParaScape runs
`prediction.run`, so the eight methods, the holdout by family, the held-out metrics and
the full run export existed in one place and one place only. The controls are MOVED here
rather than rewritten, so `prediction.TaskSpec` is built from the same widgets and the
same background runner as before.

WHAT IT IS FOR
--------------
It lets you take one measured column and see which genes that carry no measurement the
data would call, each with the support behind it -- so that you can decide whether a call
of that kind is worth believing before anything is written into the map.

HOW IT DIFFERS FROM INFER
-------------------------
Infer answers a REGISTERED question: the holdout, the control and the leakage closure
arrive already chosen and already checked. Here you choose the column, the method and
what to hold out, which is more freedom and less protection. The module says so on screen
rather than leaving the two looking interchangeable.

WHAT IT REFUSES TO DO
---------------------
A hypothesis is never drawn as an annotation. Every call carries whether it is supported,
why it abstained when it did, and whether the probabilities behind it are calibrated; the
sentence above the table says in words that these are model calls and not measurements.
"""
from __future__ import annotations

import pandas as pd
from PySide6 import QtCore, QtWidgets

from . import ports as PORTS
from . import theme as TH
from .shell import ModuleShell, solid_dropdown

#: Metrics named in the trust sentence, in the order they are read out.
HEADLINE = ("accuracy", "macro_f1", "macro_average_precision", "rmse", "r2", "coverage")
#: How many candidate rows the table shows. The export carries every one.
SHOWN = 1000


def _muted(text: str) -> QtWidgets.QLabel:
    """A dimmed line of explanation, styled by the shell rather than by a colour here."""
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


def _combo(items, tip):
    """A choice whose stored value is independent of the words on screen.

    Kept from `workflows.py`: the method names on screen are plain words and the values
    behind them are what `prediction.TaskSpec` takes, and letting the two be the same
    string is how a rewording of a label silently changes a run.
    """
    box = QtWidgets.QComboBox()
    for label, value in items:
        box.addItem(label, value)
    box.setToolTip(TH.tip(tip))
    solid_dropdown(box)
    return box


def _table():
    """A read-only table, so reading a call cannot edit it into a measurement."""
    table = QtWidgets.QTableWidget()
    table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
    table.setAlternatingRowColors(True)
    table.horizontalHeader().setStretchLastSection(True)
    table.verticalHeader().setVisible(False)
    return table


def fill(table, frame, limit: int = SHOWN) -> None:
    """Show at most `limit` rows of `frame` in `table`, leaving the frame itself alone.

    Bounded on purpose: the whole frame goes to the export, and a table with forty
    thousand rows in it is slower to draw than the model was to fit.
    """
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


class PredictModule(QtCore.QObject):
    """The Predict module for one window: its shell, its settings and what each does."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.nodes = window.nodes
        self.runner = getattr(window, "jobs", None)
        self.result = None
        self.job = None
        self.mine: list = []
        self.shell = shell = ModuleShell(
            "predict", "Predict",
            "predict any measured column, and measure how often that is right",
            instruction="Choose a column to predict, a method, and what to hold out.")
        shell.module = self
        if self.runner is not None:
            self.runner.finished.connect(self._finished)
            self.runner.progress.connect(self._progress)
        self._build_settings()
        self._build_runtime()
        shell.home_requested.connect(lambda: _go_home(window))

    # ------------------------------------------------------------------ the left column
    def _build_settings(self) -> None:
        """The settings column: what to predict, whose genes, how, and what to hold out."""
        shell = self.shell
        columns = [c for c in self.nodes.columns if c != "gene_id"]

        what = QtWidgets.QWidget(shell)
        wl = QtWidgets.QVBoxLayout(what)
        wl.setContentsMargins(6, 2, 6, 6)
        self.target = _combo(
            [(c, c) for c in columns],
            "The measured column to predict. Genes with no measurement for it are never "
            "treated as a measured no; they are the ones a call is made about.")
        self.target.setEditable(True)
        self.target.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)
        if "compartment" in columns:
            self.target.setCurrentText("compartment")
        wl.addWidget(self.target)
        self.kind = _combo(
            [("A category", "classification"), ("A number", "regression")],
            "A category is predicted as a label with a probability. A number keeps the "
            "measurement as a number and reports how far out it was, with an interval.")
        wl.addWidget(self.kind)
        wl.addWidget(_muted("Genes with no measurement for this column are the ones a call "
                            "is made about."))
        shell.add_setting("What to predict", what)

        mine = QtWidgets.QWidget(shell)
        ml = QtWidgets.QVBoxLayout(mine)
        ml.setContentsMargins(6, 2, 6, 6)
        self.field = QtWidgets.QLineEdit()
        self.field.setPlaceholderText("Gene ids, separated by commas or spaces, then Return")
        self.field.setToolTip(TH.tip("Your own genes. The calls made about them are read out "
                                     "under the table, so a screen hit list can be checked "
                                     "against the answer without searching the table."))
        self.field.returnPressed.connect(self._genes_typed)
        ml.addWidget(self.field)
        row = QtWidgets.QHBoxLayout()
        choose = QtWidgets.QPushButton("Choose a list…")
        choose.setToolTip(TH.tip("A CSV, TSV, XLS or XLSX with gene ids in a column."))
        choose.clicked.connect(self._choose_list)
        use = QtWidgets.QPushButton("Use the gene list")
        use.setToolTip(TH.tip("The list last dropped on ParaScape (the gene list panel)."))
        use.clicked.connect(self._use_gene_list)
        clear = QtWidgets.QPushButton("Clear")
        clear.setToolTip(TH.tip("Forget your genes; the run itself is unchanged."))
        clear.clicked.connect(lambda _checked=False: self.set_genes([]))
        for button in (choose, use, clear):
            row.addWidget(button)
        ml.addLayout(row)
        self.mine_note = _muted("No genes of your own. The whole table is predicted either way.")
        ml.addWidget(self.mine_note)
        shell.add_setting("Your genes", mine)

        how = QtWidgets.QWidget(shell)
        hl = QtWidgets.QVBoxLayout(how)
        hl.setContentsMargins(6, 2, 6, 6)
        self.method = _combo(
            [("Straight line through the measurements", "linear"),
             ("Boosted trees", "boosted"),
             ("Nearest genes by their measurements", "neighbors"),
             ("Nearest genes after PCA", "pca"),
             ("Nearest genes on a UMAP", "umap"),
             ("Masked factors", "multiview"),
             ("The link network, as it is", "network"),
             ("Guess the commonest answer", "prior")],
            "Every method is scored on the same held-out groups, so they can be compared. "
            "Guess the commonest answer is the chance baseline: a method that does not beat "
            "it has found nothing. The map on screen is never used to predict.")
        self.method.currentIndexChanged.connect(self._method_changed)
        hl.addWidget(self.method)
        shell.add_setting("How to predict", how)

        hold = QtWidgets.QWidget(shell)
        hol = QtWidgets.QVBoxLayout(hold)
        hol.setContentsMargins(6, 2, 6, 6)
        groups = [(c, c) for c in columns
                  if not pd.api.types.is_numeric_dtype(self.nodes[c])]
        self.group = _combo(
            [("Separate genes, chosen at random", None)] + groups,
            "What is kept together, so that a gene is never scored on a group the method "
            "has already seen. Orthogroup keeps related proteins together, which is the "
            "honest question: splitting a family across the holdout answers an easier one.")
        if "orthogroup" in self.nodes.columns:
            self.group.setCurrentText("orthogroup")
        hol.addWidget(self.group)
        self.folds = QtWidgets.QSpinBox()
        self.folds.setRange(2, 10)
        self.folds.setValue(3)
        self.folds.setToolTip(TH.tip("How many times the data is split. Certainty is fitted "
                                     "inside each training split, on groups of its own."))
        form = QtWidgets.QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.addRow("Times to split", self.folds)
        holder = QtWidgets.QWidget(hold)
        holder.setLayout(form)
        hol.addWidget(holder)
        shell.add_setting("What to hold out", hold)

        advanced = QtWidgets.QWidget(shell)
        al = QtWidgets.QVBoxLayout(advanced)
        al.setContentsMargins(6, 2, 6, 6)
        self.threshold = QtWidgets.QDoubleSpinBox()
        self.threshold.setRange(0.0, 1.0)
        self.threshold.setSingleStep(0.05)
        self.threshold.setToolTip(TH.tip(
            "Below this probability the method says nothing rather than guessing. Zero "
            "applies only the rule about how much of a gene was measured. Check whether "
            "the probabilities are calibrated before reading one as a chance of being right."))
        aform = QtWidgets.QFormLayout()
        aform.setContentsMargins(0, 0, 0, 0)
        aform.addRow("Say nothing below", self.threshold)
        aholder = QtWidgets.QWidget(advanced)
        aholder.setLayout(aform)
        al.addWidget(aholder)
        self.sequence = QtWidgets.QCheckBox("Use the frozen protein sequence features")
        self.sequence.setToolTip(TH.tip(
            "The bundled ESM-2 table, joined by gene id. It folds nothing and downloads "
            "nothing. Sequence evidence is not automatically an improvement; it is scored "
            "on the same held-out genes as everything else."))
        self.sequence.setEnabled(self._sequence_possible())
        al.addWidget(self.sequence)
        shell.add_setting("Advanced", advanced, expanded=False)

        save = QtWidgets.QWidget(shell)
        sl = QtWidgets.QVBoxLayout(save)
        sl.setContentsMargins(6, 2, 6, 6)
        self.export_button = QtWidgets.QPushButton("Save the whole run…")
        self.export_button.setEnabled(False)
        self.export_button.setToolTip(TH.tip(
            "Every row, the per-class results, which genes were in which split, and the "
            "provenance -- so the run can be read again without this window."))
        self.export_button.clicked.connect(self._export)
        sl.addWidget(self.export_button)
        shell.add_setting("Save", save, expanded=False)

    # ----------------------------------------------------------------- the right column
    def _build_runtime(self) -> None:
        """The runtime column: the calls, how far to trust them, and what to do next."""
        shell = self.shell

        calls = QtWidgets.QWidget(shell)
        cl = QtWidgets.QVBoxLayout(calls)
        cl.setContentsMargins(8, 4, 8, 8)
        self.caption = _muted(
            "Nothing run yet. Choose a column on the left and press Predict. What appears "
            "here is what the method CALLS for genes that carry no measurement -- model "
            "calls, never measurements.")
        cl.addWidget(self.caption)
        self.predictions = _table()
        self.predictions.setToolTip(TH.tip(
            "One row per gene with no measurement for the chosen column. Double-click a "
            "row to select that gene everywhere else in ParaScape."))
        self.predictions.cellDoubleClicked.connect(self._row_picked)
        cl.addWidget(self.predictions, 1)
        shell.add_section("What the data would call", calls, extent=380, stretch=3, focus=True)

        self.trust = QtWidgets.QTextBrowser(shell)
        shell.link_style(self.trust)
        self.trust.setHtml("<p>Run something and this says how often a call of that kind was "
                           "right on the groups the method was not trained on, against the "
                           "chance baseline.</p>")
        shell.add_section("How far to trust it", self.trust, extent=200, stretch=1)

        self.yours = QtWidgets.QTextBrowser(shell)
        shell.link_style(self.yours)
        self.yours.setHtml("<p>Choose genes of your own on the left to read the calls made "
                           "about them here.</p>")
        shell.add_section("Your genes in the answer", self.yours, extent=150, stretch=1,
                          folded=True)

        from .next_steps import NextSteps
        self.next_steps = NextSteps(PORTS.CLAIM, shell,
                                    lead="Next, with what these calls say:")
        self.next_steps.chosen.connect(lambda key, _kind, _value: self._open(key))
        shell.add_section("What you can do with this next", self.next_steps,
                          extent=110, stretch=0)
        self.next_steps.set_subject("", available=self._reasons(), without="predict")

        self.run_button = shell.add_action(
            "Predict", self.run, primary=True,
            tip="Score the chosen method fold by fold on groups it was not trained on, then "
                "call the genes that carry no measurement.")
        self.cancel_button = shell.add_action(
            "Stop", self.cancel,
            tip="Stop between model fits, so nothing is killed half way through one.")
        self.cancel_button.setEnabled(False)
        shell.add_action("Read a gene in Gene", self._to_gene,
                         tip="Open the selected gene in the Gene module: what is known about "
                             "it, with sources.")
        shell.add_action("Ask a registered question instead", lambda: self._open("infer"),
                         tip="Infer answers a question whose holdout, control and leakage "
                             "closure are already chosen and checked. Less freedom, more "
                             "protection.")
        shell.borrow("Console and AI chat",
                     [("Console", getattr(self.window, "console_dock", None)),
                      ("AI chat", getattr(self.window, "chat_dock", None))],
                     folded=True)

    # ------------------------------------------------------------------------- settings
    def _sequence_possible(self) -> bool:
        """Whether the bundled sequence table can be joined to this gene table at all.

        Only the Toxoplasma accessions are in it, so the box is disabled rather than
        offered and then failing inside the job.
        """
        try:
            return bool(self.nodes.gene_id.astype(str).str.startswith("TGME49_").any())
        except Exception:
            return False

    def _method_changed(self) -> None:
        """Turn the sequence features off for the method that uses the link network.

        The network method works on the fixed graph rather than on per-gene features, so
        sequence features would be loaded and then ignored.
        """
        self.sequence.setEnabled(self.method.currentData() != "network"
                                 and self._sequence_possible())

    def spec(self):
        """The `prediction.TaskSpec` the settings describe, or None with a reason said.

        A method of its own so a test can read the settings without starting a job, which
        is how the eight methods are checked against the widgets that choose them.
        """
        from .prediction import TaskSpec
        target = self.target.currentText().strip()
        if target not in self.nodes.columns:
            self.shell.say(f"{target or 'that'} is not a column of this table.")
            return None
        try:
            return TaskSpec(target, kind=self.kind.currentData(),
                            method=self.method.currentData(),
                            group_column=self.group.currentData(),
                            folds=self.folds.value(),
                            min_probability=self.threshold.value())
        except ValueError as error:
            self.shell.say(str(error))
            return None

    def set_genes(self, genes) -> None:
        """Hold a list of your own genes, and read the answer against it."""
        self.mine = [str(g) for g in genes if str(g).strip()]
        if self.mine:
            self.mine_note.setText(f"{len(self.mine):,} genes of your own.")
        else:
            self.mine_note.setText("No genes of your own. The whole table is predicted "
                                   "either way.")
        self._show_mine()

    def _genes_typed(self) -> None:
        """Take the gene ids typed into the field as your list."""
        text = self.field.text().replace(",", " ")
        self.set_genes(text.split())

    def _choose_list(self) -> None:
        """Read a file of gene ids as your list, saying so when it cannot be read."""
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self.shell, "A table with gene ids", "",
            "Tables (*.csv *.tsv *.txt *.xls *.xlsx)")
        if not path:
            return
        try:
            from . import gene_drop as GD
            known = set(self.nodes.gene_id.astype(str))
            drop = GD.read_drop(path, known_ids=known)
            self.set_genes([g for g in drop.genes if g in known])
        except Exception as error:
            self.shell.say(f"could not read {path}: {error}")

    def _use_gene_list(self) -> None:
        """Take the list already dropped on ParaScape as your list."""
        panel = getattr(self.window, "home", None)
        genes = panel.genes() if panel is not None else []
        if genes:
            self.set_genes(genes)
        else:
            self.shell.say("no gene list is loaded: drop or paste one first")

    # ----------------------------------------------------------------------- running it
    def run(self) -> None:
        """Score the chosen method on held-out groups, then call the unmeasured genes.

        The fitting happens in the window shared job runner rather than here, so it shows
        up in the Jobs panel with everything else and can be stopped from there as well as
        from this module.
        """
        if self.job is not None and getattr(self.job, "active", False):
            self.shell.say("a run is already going; stop it first")
            return
        spec = self.spec()
        if spec is None:
            return
        if self.runner is None:
            self.shell.say("no job runner in this window, so nothing can be fitted")
            return
        nodes = self.nodes.copy()
        add_sequence = self.sequence.isChecked() and self.sequence.isEnabled()
        self.result = None
        self.export_button.setEnabled(False)
        self.predictions.setRowCount(0)
        self.caption.setText(f"Scoring {spec.method} on held-out groups of "
                             f"{spec.target}…")
        self.run_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.job = self.runner.submit(_work(nodes, spec, add_sequence, self.runner),
                                      f"Predict {spec.target} ({spec.method})")

    def cancel(self) -> None:
        """Ask the run to stop at the next boundary between model fits."""
        if self.job is not None:
            self.job.cancel()
        self.shell.say("Stop asked for. The fit on screen finishes before it stops, so "
                       "nothing is left half done.")

    def _progress(self, job_id, _percent, note) -> None:
        """Say what the run is doing, and only for the run of this module."""
        if self.job is not None and job_id == self.job.id:
            self.shell.say(str(note))

    def _finished(self, job_id, ok) -> None:
        """Put the controls back, and show a run only when it actually finished."""
        if self.job is None or job_id != self.job.id:
            return
        self.run_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        if not ok:
            reason = self.job.error or self.job.note or self.job.state
            self.caption.setText(f"No result: {reason}")
            self.shell.say("The job record keeps the error or the cancellation.")
            return
        self.result = self.job.result
        self.export_button.setEnabled(True)
        self.show_result(self.result)

    def show_result(self, result) -> None:
        """Draw a finished run: the calls, the trust sentence and your genes in it.

        Public so a test can hand it a result without fitting anything, which is how the
        wording about calibration and support is checked.
        """
        candidates = result.predictions[result.predictions.role.eq("unlabelled_candidate")]
        supported = int(candidates.supported.sum()) if len(candidates) else 0
        fixed = result.spec.method == "network"
        self.caption.setText(
            f"{len(candidates):,} genes carry no measurement for {result.spec.target}. "
            f"{supported:,} of them have enough evidence behind the call to be worth "
            f"reading. These are calls the method makes, not measurements. Showing at "
            f"most {SHOWN:,}; the saved run has every row.")
        columns = [c for c in ("gene_id", "prediction", "score", "feature_coverage",
                               "supported", "abstention_reason", "calibration_status",
                               "interval_lower", "interval_upper") if c in candidates]
        fill(self.predictions, candidates[columns])
        self.trust.setHtml(self._trust_html(result, fixed))
        self._show_mine()

    def _trust_html(self, result, fixed: bool) -> str:
        """How far to trust this run, in sentences, against the chance baseline."""
        from html import escape
        metrics = result.metrics
        lines = []
        for key in HEADLINE:
            if key in metrics:
                lines.append(f"<li>{escape(key.replace('_', ' '))}: "
                             f"<b>{metrics[key]:.3f}</b></li>")
            if key == "accuracy" and "accuracy" in metrics:
                # The chance baseline named in words, because an accuracy on its own says
                # nothing: the commonest answer already scores the share it covers.
                lines.append("<li>Compare it with the <i>Guess the commonest answer</i> "
                             "method on the same splits -- that is what chance is here.</li>")
        how = ("The link network is used as it is, so a gene is placed using links to "
               "genes in the training groups: the numbers describe this graph and do not "
               "transfer to a gene outside it."
               if fixed else
               "Nothing about the held-out groups reached the fitting, so the numbers are "
               "what to expect on genes like them.")
        group = result.spec.group_column or "nothing: the genes were split at random"
        return (f"<p>{escape(result.spec.target)} by <b>{escape(result.spec.method)}</b>, "
                f"held out by {escape(str(group))}, {result.spec.folds} splits.</p>"
                f"<ul>{''.join(lines)}</ul><p>{how}</p>")

    def _show_mine(self) -> None:
        """Read the calls made about your own genes, when there are both."""
        from html import escape
        if not self.mine:
            self.yours.setHtml("<p>Choose genes of your own on the left to read the calls "
                               "made about them here.</p>")
            return
        if self.result is None:
            self.yours.setHtml(f"<p>{len(self.mine):,} genes of your own. Predict something "
                               f"and the calls made about them appear here.</p>")
            return
        frame = self.result.predictions
        got = frame[frame.gene_id.astype(str).isin(set(self.mine))]
        if not len(got):
            self.yours.setHtml(f"<p>None of your {len(self.mine):,} genes is in this run.</p>")
            return
        rows = []
        for row in got.head(50).itertuples(index=False):
            call = getattr(row, "prediction", "")
            role = getattr(row, "role", "")
            rows.append(f"<li><b>{escape(str(row.gene_id))}</b>: {escape(str(call))} "
                        f"<i>({escape(str(role))})</i></li>")
        more = len(got) - 50
        tail = f"<p>…and {more:,} more.</p>" if more > 0 else ""
        self.yours.setHtml(
            f"<p>{len(got):,} of your {len(self.mine):,} genes are in this run.</p>"
            f"<ul>{''.join(rows)}</ul>{tail}")

    # ------------------------------------------------------------------------- going on
    def _row_picked(self, row: int, _column: int) -> None:
        """Select the gene of a clicked row everywhere else, without calling it measured."""
        item = self.predictions.item(row, 0)
        if item is not None and hasattr(self.window, "_workflow_gene"):
            self.window._workflow_gene(item.text())
            self.shell.say(f"{item.text()} selected")

    def _to_gene(self) -> None:
        """Open the selected gene in the Gene module."""
        item = self.predictions.currentItem()
        row = item.row() if item is not None else -1
        cell = self.predictions.item(row, 0) if row >= 0 else None
        from .gene_module import open_in
        open_in(self.window, cell.text() if cell is not None else "")

    def _open(self, key: str) -> str:
        """Open another module, which is how a call becomes the next step."""
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
        """Write the whole evaluated run to a folder, saying so when the write fails."""
        if self.result is None:
            self.shell.say("nothing to save: predict something first")
            return
        directory = QtWidgets.QFileDialog.getExistingDirectory(self.shell, "Save the whole run")
        if not directory:
            return
        try:
            self.result.save(directory)
            self.shell.say(f"saved the whole run to {directory}")
        except (OSError, ValueError) as error:
            self.shell.say(f"could not save: {error}")


def _work(nodes, spec, add_sequence: bool, runner):
    """The fitting, as a function for the job runner to call off the GUI thread.

    Built outside the module so the job holds no reference to a widget: a job that
    outlived its module and then touched one is the crash this shape avoids.
    """
    def work(job):
        """Fit with the numerical libraries bounded, reporting between model fits."""
        from .jobs import Stopped
        from .prediction import run

        def log(message):
            """Publish progress, and stop between fits when stopping was asked for."""
            if job.cancelled:
                raise Stopped("stopped between model fits")
            job.note = message
            runner._sig.progress.emit(job.id, -1, message)

        frame = nodes
        if add_sequence and spec.method != "network":
            log("Loading the frozen sequence features")
            from . import paths
            sequence = pd.read_parquet(paths.cache_file("esm_features.parquet"))
            frame = frame.merge(sequence, on="gene_id", how="left", validate="one_to_one")
        log("Preparing the splits")
        from threadpoolctl import threadpool_limits
        with threadpool_limits(limits=2):
            if spec.method == "network":
                from .network_prediction import bundled_network, run_network
                graph, sources = bundled_network(frame)
                return run_network(frame, spec, graph, list(sources), sources, log=log)
            return run(frame, spec, log=log)
    return work


def build(window) -> ModuleShell | None:
    """The Predict module for `window`, or None when there is no gene table to predict."""
    if getattr(window, "nodes", None) is None:
        return None
    return PredictModule(window).shell


def open_in(window) -> str:
    """Show the Predict module over the window. Returns what the status bar was told."""
    from .landing_install import enter_landing
    landing = getattr(window, "landing", None)
    shell = getattr(window, "predict_module", None)
    if landing is None or shell is None:
        return "the Predict module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    return "Predict"


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)
