"""The Infer module: ask a question of the data, and get genes named with how far to trust them.

Built with the parascape-module mould (instruction 96). It lets you take a question -- "which
unlabelled genes are host-facing?" -- and see the genes the data names for it, each with the
purity of the group it came from and whether an independent control agrees, so that you can
decide what to believe and carry the candidates on to the Gene module.

The controls are the analysis panel's own "8 · Questions" tab, MOVED here (not copied): the
question list, its description, the answer table and the report, so `check_question`,
`run_question` and `save_question_pdf` run exactly as before. New here: your own genes, read
against the answer -- which of them were named, and for what.
"""
from __future__ import annotations

import os
from html import escape

import numpy as np
import pandas as pd
from PySide6 import QtCore, QtWidgets

from . import theme as TH
from .shell import ModuleShell, solid_dropdown

QUESTIONS_TAB = "8 · Questions"


def _muted(text: str) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


class InferModule(QtCore.QObject):
    """The Infer module for one window."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.panel = panel = window.panel
        self.mine: list = []                           # gene ids the user brought
        self.shell = shell = ModuleShell(
            "infer", "Infer", "ask the data a question; get genes named, with how far to trust them",
            instruction="Choose a question, check what it may use, then answer it.")
        shell.module = self

        # ---------------------------------------------------------------- left
        q = QtWidgets.QWidget(shell)
        ql = QtWidgets.QVBoxLayout(q)
        ql.setContentsMargins(6, 2, 6, 6)
        solid_dropdown(panel.question_choice)
        # The questions are whole sentences; the box shows the start and the full question is on
        # the right. Sized by its contents it pushed the column past the window and hid buttons.
        panel.question_choice.setMinimumWidth(0)
        panel.question_choice.setSizeAdjustPolicy(
            QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        panel.question_choice.setMinimumContentsLength(28)
        panel.question_choice.view().setTextElideMode(QtCore.Qt.TextElideMode.ElideNone)
        panel.question_choice.view().setMinimumWidth(640)
        ql.addWidget(panel.question_choice)
        ql.addWidget(_muted("Each question comes with what it may use, the label it must recover "
                            "and an independent control. One is shipped because it is refused."))
        shell.add_setting("Your question", q)

        uses = QtWidgets.QWidget(shell)
        self.uses_layout = QtWidgets.QVBoxLayout(uses)
        self.uses_layout.setContentsMargins(6, 2, 6, 6)
        self.uses_note = _muted("")
        self.uses_layout.addWidget(self.uses_note)
        self.uses_boxes: dict = {}
        shell.add_setting("What the answer may use", uses)
        panel.question_choice.currentIndexChanged.connect(lambda *_: self.refresh_uses())

        mine = QtWidgets.QWidget(shell)
        ml = QtWidgets.QVBoxLayout(mine)
        ml.setContentsMargins(6, 2, 6, 6)
        self.field = QtWidgets.QLineEdit()
        self.field.setPlaceholderText("Gene ids, separated by commas (optional), then Return")
        self.field.setToolTip(TH.tip("Genes you care about: after the answer, you are told which "
                                     "of them were named and for what."))
        self.field.returnPressed.connect(self._typed)
        ml.addWidget(self.field)
        row = QtWidgets.QHBoxLayout()
        choose = QtWidgets.QPushButton("Choose a list…")
        choose.setToolTip("A table of gene ids to read against the answer.")
        choose.clicked.connect(self._choose_list)
        use = QtWidgets.QPushButton("Use the gene list")
        use.setToolTip("The list last dropped on ParaScape (the gene list panel).")
        use.clicked.connect(self._use_gene_list)
        row.addWidget(choose)
        row.addWidget(use)
        ml.addLayout(row)
        self.mine_note = _muted("No genes of yours: the whole answer is shown.")
        ml.addWidget(self.mine_note)
        shell.add_setting("Your genes", mine)

        save = QtWidgets.QWidget(shell)
        sl = QtWidgets.QVBoxLayout(save)
        sl.setContentsMargins(6, 2, 6, 6)
        sl.addWidget(panel.question_pdf)
        export = QtWidgets.QPushButton("Save the named genes as a table…")
        export.setToolTip("The answer table as CSV: every named gene with its purity and the "
                          "control's verdict.")
        export.clicked.connect(self._export)
        sl.addWidget(export)
        methods = QtWidgets.QPushButton("Copy the methods text")
        methods.setToolTip("A methods paragraph, the numbered references and BibTeX for this "
                           "question, from the sources of the columns it may use. A source whose "
                           "publication is not recorded is flagged, never cited.")
        methods.clicked.connect(panel.copy_methods_text)
        sl.addWidget(methods)
        reopen = QtWidgets.QPushButton("Open a saved bundle…")
        reopen.setToolTip("Read a bundle saved here back: its files are checked against their hashes, "
                          "its question is selected again and its named genes shown, with the left-out "
                          "evidence restored.")
        reopen.clicked.connect(lambda: self.open_bundle())
        sl.addWidget(reopen)
        bundle = QtWidgets.QPushButton("Save everything as a bundle…")
        bundle.setToolTip("One folder for this answer: the named genes (marked inferred), the "
                          "recovery scores, the map as a vector figure, the methods text with its "
                          "references, and a manifest whose hashes show if anything was edited later.")
        bundle.clicked.connect(self.save_bundle)
        sl.addWidget(bundle)
        kept = QtWidgets.QPushButton("Export your kept claims…")
        kept.setToolTip("Every claim you kept, with its provenance, as CSV.")
        kept.clicked.connect(self._export_kept)
        sl.addWidget(kept)
        shell.add_setting("Save the answer", save, expanded=False)

        # ---------------------------------------------------------------- right
        shell.link_style(panel.question_detail)
        shell.add_section("The question", panel.question_detail, extent=200, stretch=1)
        answer = QtWidgets.QWidget(shell)
        al = QtWidgets.QVBoxLayout(answer)
        al.setContentsMargins(8, 4, 8, 8)
        self.answer_note = _muted("Answer the question to see the genes it names.")
        al.addWidget(self.answer_note)
        panel.question_table.setMinimumHeight(140)
        panel.question_table.cellDoubleClicked.connect(lambda r, _c: self._read_row(r))
        al.addWidget(panel.question_table, 1)
        shell.add_section("The answer", answer, extent=360, stretch=3)
        shell.link_style(panel.question_report)
        shell.add_section("How far to trust it", panel.question_report, extent=180, stretch=1)

        shell.add_action("Answer this question", panel.run_question, primary=True,
                         tip="Build the map from what the question may use, group it, and name the "
                             "unlabelled genes in groups that are mostly one label; the control is "
                             "scored on the same groups.")
        shell.add_action("Check what it may use", panel.check_question,
                         tip="Without building anything: which inputs the question may use, and "
                             "what was removed because it would restate the answer.")
        shell.add_action_widget(panel.stop_button())
        shell.add_action("Read the chosen gene", lambda: self._read_row(panel.question_table.currentRow()),
                         tip="Open the selected named gene in the Gene module.")
        shell.add_action("Put the answer on the map", self._to_map,
                         tip="Mark every named gene in the Map module.")
        shell.add_action("Keep the chosen answers", self.keep_selected,
                         tip="Add the selected named genes to your own layer of kept claims, with "
                             "the question, the purity and the control's verdict. Kept claims show "
                             "in the Gene module as kept, never as facts, and are never training "
                             "data.")
        shell.borrow("Console and AI chat", [("Console", getattr(window, "console_dock", None)),
                                             ("AI chat", getattr(window, "chat_dock", None))],
                     folded=True)
        shell.home_requested.connect(lambda: _go_home(window))
        panel.status.connect(shell.say)
        panel.status.connect(lambda *_: self._after_answer())
        _replace_tab(panel, window)
        panel.show_question()
        self.refresh_uses()

    # ------------------------------------------------------------------ what the answer may use
    def refresh_uses(self) -> None:
        """One tick box per kind of evidence the question's inputs hold; the label's own kind is
        shown held out, with the reason, and cannot be ticked (instruction 74, the data step)."""
        from .recipes import _categories, input_categories
        for box in self.uses_boxes.values():
            box.setParent(None)
            box.deleteLater()
        self.uses_boxes = {}
        self.panel.question_leave_out = ()
        q = self.panel.current_question()
        nodes = getattr(self.window, "nodes", None)
        if q is None or nodes is None:
            self.uses_note.setText("Choose a question to see what its answer may use.")
            return
        from .questions import as_recipe
        recipe = as_recipe(q)
        try:
            counts = input_categories(nodes, recipe)
            held = set(_categories(recipe.organism)(recipe.holdout)) - {"identity"}
        except Exception as error:
            self.uses_note.setText(f"Could not read the inputs: {error}")
            return
        own = ", ".join(sorted(c.replace("_", " ") for c in held)) or "none recorded"
        self.uses_note.setText("Kinds of evidence the question's inputs hold. Untick one to leave it "
                               f"out. Held out always: the label's own kind ({own}), so the answer "
                               "is not circular.")
        for category, n in counts.items():
            if category == "identity":
                continue
            box = QtWidgets.QCheckBox(f"{category.replace('_', ' ')} ({n} column{'s' if n != 1 else ''})")
            if category in held:
                box.setChecked(False)
                box.setEnabled(False)
                box.setToolTip(f"Held out: {category.replace('_', ' ')} is the kind of knowledge the "
                               "label itself is, so using it would predict the label from itself.")
            else:
                box.setChecked(True)
                box.setToolTip(f"Untick to build the answer without any {category.replace('_', ' ')} "
                               "column.")
                box.toggled.connect(lambda *_: self._uses_changed())
            self.uses_layout.addWidget(box)
            self.uses_boxes[category] = box

    def _uses_changed(self) -> None:
        left = tuple(c for c, b in self.uses_boxes.items() if b.isEnabled() and not b.isChecked())
        self.panel.question_leave_out = left
        self.shell.say("leaving out: " + (", ".join(c.replace("_", " ") for c in left) or "nothing"))

    # ------------------------------------------------------------------ your genes
    def _ids(self) -> np.ndarray:
        return self.window.nodes["gene_id"].astype(str).to_numpy()

    def _set_mine(self, genes: list, note: str) -> None:
        known = set(self._ids())
        self.mine = [g for g in genes if g in known]
        self.mine_note.setText(f"{len(self.mine)} genes of yours" + (f" ({note})" if note else "")
                               if self.mine else "No genes of yours: the whole answer is shown.")
        self._after_answer()

    def _typed(self) -> None:
        tokens = [t for t in self.field.text().replace(",", " ").split() if t]
        lower = {g.lower(): g for g in self._ids()}
        self._set_mine([lower[t.lower()] for t in tokens if t.lower() in lower], "typed")

    def load_list(self, path: str) -> None:
        """Read a file of genes to read against the answer."""
        from . import gene_drop as GD
        try:
            drop = GD.read_drop(path, known_ids=list(self._ids()))
        except Exception as error:
            self.shell.say(f"{os.path.basename(path)} could not be read: {error}")
            return
        self._set_mine(list(drop.genes), os.path.basename(path))

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
            self.shell.say("No gene list loaded: drop a file on ParaScape or choose one.")
            return
        self._set_mine(list(genes), "the gene list")

    # ------------------------------------------------------------------ the answer
    def answer(self):
        """The named-gene table of the last answer, or None."""
        result = getattr(self.panel, "_question_result", None)
        frame = getattr(result, "inference", None)
        return frame if frame is not None and len(frame) else None

    @staticmethod
    def _gene_column(frame) -> str:
        for name in ("gene_id", "gene", "id"):
            if name in frame.columns:
                return name
        return frame.columns[0]

    def _after_answer(self) -> None:
        """Say how many genes were named, and which of yours, for what."""
        frame = self.answer()
        if frame is None:
            self.answer_note.setText("Answer the question to see the genes it names.")
            return
        column = self._gene_column(frame)
        named = frame[column].astype(str)
        agrees = int(frame["control_agrees"].fillna(False).sum()) if "control_agrees" in frame else 0
        text = (f"<b>{len(frame):,}</b> genes named; the control agrees for <b>{agrees:,}</b>. "
                "Double-click a gene to read it.")
        if self.mine:
            hits = frame[named.isin(set(self.mine))]
            label = next((c for c in ("label", "predicted", "call", "holdout_label") if c in frame), None)
            if len(hits):
                what = ", ".join(f"{escape(str(r[column]))}" + (f" ({escape(str(r[label]))})" if label else "")
                                 for _i, r in hits.head(8).iterrows())
                text += f"<br>Of your {len(self.mine)} genes, <b>{len(hits)}</b> are named: {what}."
            else:
                text += f"<br>None of your {len(self.mine)} genes is named."
        text += self._artefact_warning(frame[column].astype(str).tolist())
        self.answer_note.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.answer_note.setText(text)

    def _artefact_warning(self, genes: list) -> str:
        """A warning when the named genes rest on filled-in values more than on measurements."""
        result = getattr(self.panel, "_question_result", None)
        columns = getattr(getattr(result, "closure", None), "columns", ())
        if not columns:
            return ""
        from .diagnostics import missingness
        check = missingness(self.window.nodes, columns, genes)
        return f"<br><b>Caution:</b> {escape(check['sentence'])}" if check.get("warn") else ""

    def _row_gene(self, row: int) -> str:
        frame = self.answer()
        if frame is None or not (0 <= row < len(frame)):
            item = self.panel.question_table.item(row, 0) if row >= 0 else None
            return item.text() if item is not None else ""
        return str(frame.iloc[row][self._gene_column(frame)])

    def _read_row(self, row: int) -> None:
        gene = self._row_gene(row)
        if gene:
            from .gene_module import open_in
            open_in(self.window, gene)

    def _to_map(self) -> None:
        frame = self.answer()
        if frame is None or getattr(self.window, "map_module", None) is None:
            return
        ids = self._ids()
        wanted = set(frame[self._gene_column(frame)].astype(str))
        from .map_module import open_in
        open_in(self.window)
        self.window.map_module.module.select([int(i) for i in np.flatnonzero(np.isin(ids, list(wanted)))],
                                             "the answer")

    def keep_selected(self) -> int:
        """Keep the selected named genes as claims in the user's own layer. Returns how many."""
        from . import kept_claims as KC
        frame = self.answer()
        if frame is None:
            self.shell.say("Answer the question first.")
            return 0
        rows = sorted({i.row() for i in self.panel.question_table.selectedIndexes()})
        if not rows:
            self.shell.say("Select the answers to keep in the table first.")
            return 0
        q = self.panel.current_question() or {}
        try:
            from .questions import as_recipe
            trait = as_recipe(q).holdout
        except Exception:
            trait = str(q.get("question", "answer"))
        column = self._gene_column(frame)
        label = next((c for c in ("label", "predicted", "call", "holdout_label") if c in frame), None)
        purity = next((c for c in ("purity", "precision", "cluster_purity") if c in frame), None)
        kept = 0
        for r in rows:
            if r >= len(frame):
                continue
            row = frame.iloc[r]
            p = float(row[purity]) if purity and row[purity] == row[purity] else float("nan")
            KC.keep(KC.KeptClaim(
                gene=str(row[column]), trait=str(trait), value=str(row[label]) if label else "named",
                probability=p, uncertainty=("the control agrees" if bool(row.get("control_agrees"))
                                            else "the control does not corroborate"),
                method=f"infer:{trait}", method_version=str(q.get("id", q.get("question", ""))[:60]),
                question=str(q.get("question", "")), accepted_score=p))
            kept += 1
        self.shell.say(f"Kept {kept} claims in your own layer.")
        return kept

    def _export_kept(self) -> None:
        from . import kept_claims as KC
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self.shell, "Export your kept claims",
                                                        "kept_claims.csv", "CSV (*.csv)")
        if path:
            self.shell.say(f"Exported {KC.export(path)} kept claims")

    def open_bundle(self, folder: str = "") -> dict:
        """Read an Infer bundle back (64.39): verify it, select its question, show its named genes."""
        from . import export_bundle
        if not folder:
            folder = export_bundle.ask_folder(self.shell)
            if not folder:
                return {}
        try:
            manifest = export_bundle.load(folder)
        except (OSError, ValueError) as error:
            self.shell.say(f"bundle refused: {error}")
            return {}
        query = manifest.get("query", {})
        texts = [str(q.get("question", "")) for q in getattr(self.panel, "_questions", [])]
        if query.get("question") in texts:
            self.panel.question_choice.setCurrentIndex(texts.index(query["question"]))
        left = tuple(query.get("leave_out_categories") or ())
        for category, box in self.uses_boxes.items():
            if box.isEnabled():
                box.setChecked(category not in left)
        self.panel.question_leave_out = left
        named = os.path.join(folder, "named_genes.csv")
        if os.path.exists(named):
            frame = pd.read_csv(named)
            self.answer_note.setText(f"From the bundle of {escape(str(manifest.get('written_at', '')))}: "
                                     f"<b>{len(frame):,}</b> genes named (the saved answer, not a re-run).")
            self.loaded_bundle = frame
        self.shell.say(f"bundle opened: {folder}")
        return manifest

    def save_bundle(self, folder: str = "") -> str:
        """Write the current answer as an export bundle (`export_bundle`); returns the folder."""
        result = getattr(self.panel, "_question_result", None)
        if result is None or getattr(result, "inference", None) is None:
            self.shell.say("answer a question first; the bundle holds its answer")
            return ""
        if not folder:
            folder = QtWidgets.QFileDialog.getExistingDirectory(self.shell, "Folder for the bundle")
            if not folder:
                return ""
        from . import export_bundle
        recipe = result.recipe
        tables, kinds = {"named genes": result.inference}, {"named genes": "inferred"}
        if len(getattr(result, "per_cluster", [])):
            tables["clusters"], kinds["clusters"] = result.per_cluster, "inferred"
        if len(getattr(result, "recovery", [])):
            tables["recovery of the held-out label"], kinds["recovery of the held-out label"] = \
                result.recovery, "measured"
        figures = {}
        view = getattr(getattr(getattr(self.window, "map_module", None), "module", None), "view", None)
        if view is not None and view.width() > 10:
            figures["map"] = view
        export_bundle.write(folder, organism=self.panel.organism,
                            query={"question": recipe.question, "holdout": recipe.holdout,
                                   "control": recipe.validation_holdout,
                                   "confidence": recipe.min_precision,
                                   "leave_out_categories": list(self.panel.question_leave_out)},
                            tables=tables, kinds=kinds, columns=result.closure.columns,
                            figures=figures, nodes=self.window.nodes,
                            strategy="map neighbourhoods and clusters",
                            settings={"confidence": recipe.min_precision},
                            n_genes=len(self.window.nodes), held_out=recipe.holdout)
        self.shell.say(f"bundle written: {folder}")
        return folder

    def _export(self) -> None:
        frame = self.answer()
        if frame is None:
            self.shell.say("Answer the question first.")
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self.shell, "Save the named genes",
                                                        "named_genes.csv", "CSV (*.csv)")
        if path:
            frame.to_csv(path, index=False)
            self.shell.say(f"Saved {len(frame)} genes to {os.path.basename(path)}")


def _replace_tab(panel, window) -> None:
    tabs = panel.findChild(QtWidgets.QTabWidget)
    if tabs is None:
        return
    for i in range(tabs.count()):
        if tabs.tabText(i) == QUESTIONS_TAB:
            page = QtWidgets.QWidget()
            v = QtWidgets.QVBoxLayout(page)
            note = QtWidgets.QLabel("Asking a question and reading its answer is now the Infer "
                                    "module.")
            note.setWordWrap(True)
            go = QtWidgets.QPushButton("Open the Infer module")
            go.setProperty("primary", True)
            go.setToolTip("Open the Infer module: choose a question, check what it may use, answer "
                          "it, and read the genes it names.")
            go.clicked.connect(lambda: open_in(window))
            v.addWidget(note)
            v.addWidget(go, 0, QtCore.Qt.AlignmentFlag.AlignLeft)
            v.addStretch(1)
            scroll = tabs.widget(i)
            if isinstance(scroll, QtWidgets.QScrollArea):
                scroll.setWidget(page)
            return


def build(window) -> ModuleShell | None:
    """The Infer module for `window`; None without the analysis panel."""
    panel = getattr(window, "panel", None)
    if panel is None or not hasattr(panel, "question_choice"):
        return None
    return InferModule(window).shell


def open_in(window) -> str:
    """Show the Infer module over the window."""
    from .landing_install import enter_landing
    landing = getattr(window, "landing", None)
    shell = getattr(window, "infer_module", None)
    if landing is None or shell is None:
        return "the Infer module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    return "Infer"


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)
