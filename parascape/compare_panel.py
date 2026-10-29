"""The Compare tab: two groups of samples from one study -- strains, a knockout and its parent, a drug and
its vehicle -- and which genes differ, by how much, and how sure (instruction 69).

Condensed first: one line saying how many genes went up and down, a volcano, and the ranked table. Every
convention is a visible choice with its reason in the tooltip -- quantifier, reference, normalisation --
because they are defended differently in the literature and the user asked to choose. Groups are guessed
from sample titles and can be reassigned, because titles are free text.
"""
from __future__ import annotations

from . import paths as _data_paths

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtWidgets

from . import expression_contrast as X
from . import theme as TH

TIPS = {
    "study": "A public study quantified from its raw reads. Its title, accession and sample count are "
             "shown; the Discoveries catalogue lists every study found and why some are not loaded.",
    "quantifier": "salmon: selective alignment with the genome as decoys, GC and sequence bias "
                  "corrected. kallisto: pseudoalignment. They usually agree; where they do not, the "
                  "gene is worth a second look.",
    "reference": "ME49+host: parasite and host transcripts compete for each read, so host reads from "
                 "conserved genes are not counted as parasite -- the default for intracellular "
                 "samples. ME49: parasite only. GT1 / VEG: the strain's own genome, for Type I / III "
                 "reads, which avoids penalising the genes that differ most between strains.",
    "normalisation": "TMM (edgeR) and median-of-ratios (DESeq2) both resist a few very abundant genes "
                     "shifting everything else; upper quartile is older; CPM scales by library size "
                     "only; TPM also corrects for gene length.",
    "groups": "Each sample's group. Guessed from its title (replicate numbers removed); change any "
              "that is wrong. A sample set to '—' is left out.",
    "compare": "Moderated comparison of B against A (limma-trend): positive log2 fold change = higher "
               "in B. Genes too lowly expressed to test are left out rather than called unchanged.",
    "fdr": "False discovery rate: of the genes called at this level, about this share are expected "
           "to be false.",
    "lfc": "Minimum absolute log2 fold change to call a gene (1 = a two-fold change).",
    "use": "Hide the called genes together and ask the fast strategies what each one is -- would "
           "the rest of the data have found them? The answer opens in the evidence panel, as Test on "
           "my genes does.",
    "save": "Save the whole table: every tested gene with its fold change, statistics and FDR.",
}

GROUPS = ("A", "B", "—")


class ComparePanel(QtWidgets.QWidget):
    """Choose a study and two groups; read which genes differ."""

    gene_chosen = QtCore.pyqtSignal(str)
    genes_chosen = QtCore.pyqtSignal(list)

    def __init__(self, root: str = X.NAS, titles: dict | None = None, parent=None):
        super().__init__(parent)
        self.root = root
        self.titles = titles or {}
        self.result: X.Contrast | None = None
        self.samples = pd.DataFrame()
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)

        self.available = X.available(root)
        row = QtWidgets.QHBoxLayout()
        self.study = QtWidgets.QComboBox()
        for study in sorted(set(self.available["study"])):
            title = self.titles.get(study, "")
            self.study.addItem(f"{study} — {title[:70]}" if title else study, study)
        self.study.setToolTip(TH.tip(TIPS["study"]))
        row.addWidget(self.study, 1)
        lay.addLayout(row)

        row = QtWidgets.QHBoxLayout()
        self.quantifier = QtWidgets.QComboBox()
        self.reference = QtWidgets.QComboBox()
        self.normalisation = QtWidgets.QComboBox()
        for box, key, items in ((self.quantifier, "quantifier", X.QUANTIFIERS),
                                (self.reference, "reference", X.REFERENCES),
                                (self.normalisation, "normalisation", X.NORMALISATIONS)):
            for i in items:
                box.addItem(i.replace("_", " "), i)
            box.setToolTip(TH.tip(TIPS[key]))
            row.addWidget(box)
        lay.addLayout(row)

        self.sample_table = QtWidgets.QTableWidget(0, 4)
        self.sample_table.setHorizontalHeaderLabels(["run", "sample", "host share", "group"])
        self.sample_table.verticalHeader().setVisible(False)
        self.sample_table.setToolTip(TH.tip(TIPS["groups"]))
        self.sample_table.setMinimumWidth(120)
        lay.addWidget(self.sample_table, 1)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("FDR ≤"))
        self.fdr = QtWidgets.QDoubleSpinBox()
        self.fdr.setRange(0.001, 0.5)
        self.fdr.setSingleStep(0.01)
        self.fdr.setDecimals(3)
        self.fdr.setValue(0.05)
        self.fdr.setToolTip(TH.tip(TIPS["fdr"]))
        row.addWidget(self.fdr)
        row.addWidget(QtWidgets.QLabel("|log2FC| ≥"))
        self.lfc = QtWidgets.QDoubleSpinBox()
        self.lfc.setRange(0.0, 10.0)
        self.lfc.setSingleStep(0.5)
        self.lfc.setValue(1.0)
        self.lfc.setToolTip(TH.tip(TIPS["lfc"]))
        row.addWidget(self.lfc)
        self.compare_btn = QtWidgets.QPushButton("Compare")
        self.compare_btn.setToolTip(TH.tip(TIPS["compare"]))
        row.addWidget(self.compare_btn)
        self.use_btn = QtWidgets.QPushButton("Test these genes")
        self.use_btn.setToolTip(TH.tip(TIPS["use"]))
        row.addWidget(self.use_btn)
        self.save_btn = QtWidgets.QPushButton("Save…")
        self.save_btn.setToolTip(TH.tip(TIPS["save"]))
        row.addWidget(self.save_btn)
        holder = QtWidgets.QWidget()
        holder.setLayout(row)
        strip = QtWidgets.QScrollArea()
        strip.setWidget(holder)
        strip.setWidgetResizable(True)
        strip.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        strip.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        strip.setFixedHeight(holder.sizeHint().height() + 2)
        lay.addWidget(strip)

        self.summary = QtWidgets.QLabel("")
        self.summary.setWordWrap(True)
        lay.addWidget(self.summary)
        # Render the scientific plot through Agg. A pyqtgraph graphics scene here crashed during
        # deferred Qt window destruction, including when switching between comparison windows.
        self.volcano = None
        self._volcano_slot = QtWidgets.QVBoxLayout()
        lay.addLayout(self._volcano_slot, 1)
        self.result_table = QtWidgets.QTableWidget(0, 5)
        self.result_table.setHorizontalHeaderLabels(["gene", "log2FC", "mean log2CPM", "p", "FDR"])
        self.result_table.verticalHeader().setVisible(False)
        self.result_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.result_table.setMinimumWidth(120)
        lay.addWidget(self.result_table, 1)

        self.study.currentIndexChanged.connect(self._study_changed)
        self.quantifier.currentIndexChanged.connect(self._study_changed)
        self.reference.currentIndexChanged.connect(self._study_changed)
        self.compare_btn.clicked.connect(lambda: self.compare())
        self.use_btn.clicked.connect(self._use)
        self.save_btn.clicked.connect(lambda: self.save())
        self.result_table.cellClicked.connect(self._gene_clicked)
        self.fdr.valueChanged.connect(self._show)
        self.lfc.valueChanged.connect(self._show)
        self._study_changed()

    # ------------------------------------------------------------------ choices
    def _conventions(self, study: str) -> pd.DataFrame:
        return self.available[self.available["study"] == study]

    def _study_changed(self):
        study = self.study.currentData()
        if not study:
            self.summary.setText("No studies have been quantified here yet. Set PARASCAPE_OMICS to the "
                                 "folder holding the expression pack.")
            return
        have = self._conventions(study)
        ref = self.reference.currentData()
        if ref not in set(have[have["quantifier"] == self.quantifier.currentData()]["reference"]):
            refs = list(have[have["quantifier"] == self.quantifier.currentData()]["reference"])
            if refs:
                self.reference.blockSignals(True)
                self.reference.setCurrentIndex(max(self.reference.findData(
                    "ME49_host" if "ME49_host" in refs else refs[0]), 0))
                self.reference.blockSignals(False)
        try:
            self.counts, self.tpm, _len, self.samples = X.load(
                study, self.quantifier.currentData(), self.reference.currentData(), self.root)
        except FileNotFoundError as e:
            self.summary.setText(str(e))
            self.samples = pd.DataFrame()
            self.sample_table.setRowCount(0)
            return
        conditions = list(dict.fromkeys(self.samples["condition"].astype(str)))
        self.sample_table.setRowCount(len(self.samples))
        for i, r in self.samples.reset_index(drop=True).iterrows():
            self.sample_table.setItem(i, 0, QtWidgets.QTableWidgetItem(str(r["run_accession"])))
            title = str(r.get("sample_title") or r.get("sample_alias") or "")
            self.sample_table.setItem(i, 1, QtWidgets.QTableWidgetItem(title))
            share = r.get("host_share")
            self.sample_table.setItem(i, 2, QtWidgets.QTableWidgetItem(
                "" if share is None or pd.isna(share) else f"{share:.0%}"))
            box = QtWidgets.QComboBox()
            box.setToolTip(TH.tip(TIPS['groups']))
            box.addItems(GROUPS)
            # First condition A, second B, the rest left out: a starting point, never a decision.
            k = conditions.index(str(r["condition"]))
            box.setCurrentIndex(0 if k == 0 else 1 if k == 1 else 2)
            self.sample_table.setCellWidget(i, 3, box)
        self.sample_table.resizeColumnsToContents()
        self.summary.setText(f"{len(self.samples)} samples in {len(conditions)} apparent conditions. "
                             f"Check the groups, then Compare.")

    def groups(self) -> tuple:
        """(group A runs, group B runs) as currently assigned in the sample table."""
        a, b = [], []
        for i in range(self.sample_table.rowCount()):
            run = self.sample_table.item(i, 0).text()
            g = self.sample_table.cellWidget(i, 3).currentText()
            (a if g == "A" else b if g == "B" else []).append(run)
        return a, b

    def set_groups(self, a: list, b: list):
        """Assign groups by run accession (everything else left out) -- for scripts and tests."""
        for i in range(self.sample_table.rowCount()):
            run = self.sample_table.item(i, 0).text()
            self.sample_table.cellWidget(i, 3).setCurrentIndex(0 if run in a else 1 if run in b else 2)

    # ------------------------------------------------------------------ actions
    def compare(self) -> X.Contrast | None:
        """Run the moderated comparison of B against A on the chosen conventions, and draw it."""
        a, b = self.groups()
        from .rna_design import independent_runs
        if not independent_runs(self.samples, a + b):
            self.summary.setText("Selected runs share a submitted biological sample. They cannot be counted as independent replicates; merge the repeated runs or select one per sample before comparison.")
            self.result = None
            return None
        try:
            self.result = X.contrast(self.counts, a, b, normalisation=self.normalisation.currentData(),
                                     tpm=self.tpm, quantifier=self.quantifier.currentData(),
                                     reference=self.reference.currentData())
        except ValueError as e:
            self.summary.setText(str(e))
            self.result = None
            return None
        self._show()
        return self.result

    def called(self) -> pd.DataFrame:
        """The genes called at the current FDR and fold-change thresholds."""
        if self.result is None:
            return pd.DataFrame()
        return self.result.significant(self.fdr.value(), self.lfc.value())

    def _show(self):
        if self.result is None:
            return
        t = self.result.table
        hit = self.called()
        up, down = int((hit["log2fc"] > 0).sum()), int((hit["log2fc"] < 0).sum())
        self.summary.setText(
            f"<b>{up:,} up, {down:,} down</b> in B against A at FDR ≤ {self.fdr.value():g} and "
            f"|log2FC| ≥ {self.lfc.value():g}, of {len(t):,} genes expressed enough to test "
            f"({len(self.result.group_a)} vs {len(self.result.group_b)} samples; "
            f"{self.result.quantifier} / {self.result.reference} / "
            f"{self.result.normalisation.replace('_', ' ')}; prior df {self.result.prior_df:.1f}).")
        import io
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from PyQt6 import QtGui
        if self.volcano is None:
            self.volcano = QtWidgets.QLabel(parent=self)
            self.volcano.setScaledContents(True)
            self.volcano.setMinimumHeight(160)
            self.volcano.setMinimumWidth(120)
            self._volcano_slot.addWidget(self.volcano)
        y = -np.log10(np.clip(t["p"].to_numpy(), 1e-300, 1))
        is_hit = t["gene"].isin(set(hit["gene"])).to_numpy()
        figure = Figure(figsize=(6, 2.8), tight_layout=True)
        canvas = FigureCanvasAgg(figure)
        axes = figure.add_subplot(111)
        axes.scatter(t["log2fc"].to_numpy()[~is_hit], y[~is_hit], s=5, c='#999999', alpha=.4)
        colours = np.where(t["log2fc"].to_numpy()[is_hit] > 0, "#c0392b", "#2e6db4")
        axes.scatter(t["log2fc"].to_numpy()[is_hit], y[is_hit], s=7, c=colours)
        axes.set(xlabel='log2 fold change (B over A)', ylabel='−log10 p')
        buffer = io.BytesIO()
        canvas.print_png(buffer)
        pixmap = QtGui.QPixmap()
        pixmap.loadFromData(buffer.getvalue(), 'PNG')
        self.volcano.setPixmap(pixmap)
        show = pd.concat([hit, t[~t["gene"].isin(set(hit["gene"]))]]).head(500)
        self.result_table.setRowCount(len(show))
        for i, r in enumerate(show.itertuples(index=False)):
            for j, v in enumerate((r.gene, f"{r.log2fc:+.2f}", f"{r.mean_expression:.1f}",
                                   f"{r.p:.2e}", f"{r.fdr:.2e}")):
                self.result_table.setItem(i, j, QtWidgets.QTableWidgetItem(str(v)))
        self.result_table.resizeColumnsToContents()

    def _gene_clicked(self, i: int, _j: int):
        item = self.result_table.item(i, 0)
        if item is not None:
            self.gene_chosen.emit(item.text())

    def _use(self):
        hit = self.called()
        if len(hit):
            self.genes_chosen.emit(hit["gene"].astype(str).tolist())

    def save(self, path: str = "") -> int:
        """Write every tested gene, with a `called` column, to CSV or TSV. Returns the rows written."""
        if self.result is None:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save comparison", "comparison.csv",
                                                            "CSV (*.csv);;TSV (*.tsv)")
        if not path:
            return 0
        t = self.result.table.assign(called=lambda d: d["gene"].isin(set(self.called()["gene"])))
        t.to_csv(path, index=False, sep="\t" if path.endswith(".tsv") else ",")
        return len(t)


def catalogue_titles() -> dict:
    """Study titles from the shipped omics catalogue, if present."""
    import os
    path = os.path.join(_data_paths.data_dir(), "omics_catalogue.tsv")
    if not os.path.exists(path):
        return {}
    c = pd.read_csv(path, sep="\t", usecols=["accession", "title"])
    return dict(zip(c["accession"].astype(str), c["title"].astype(str)))


def install(window) -> ComparePanel | None:
    """Add the Compare tab beside Discoveries and wire it to the gene card and the strategies."""
    panel = ComparePanel(titles=catalogue_titles())
    panel.gene_chosen.connect(window._workflow_gene)
    strategies = getattr(window, "strategy_panel", None)
    if strategies is not None and hasattr(window, "run_my_list"):
        panel.genes_chosen.connect(window.run_my_list)
    dock = QtWidgets.QDockWidget("compare", window)
    dock.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
    from .rna_panel import RNAProfilesPanel
    from .rna_design_panel import StudyDesignPanel
    from .published_expression_panel import PublishedExpressionPanel
    from .dataset_coverage_panel import DatasetCoveragePanel
    from .source_matrices_panel import SourceMatricesPanel
    from .deposited_values_panel import DepositedValuesPanel
    profiles = RNAProfilesPanel()
    profiles.gene_chosen.connect(window._workflow_gene)
    tabs = QtWidgets.QTabWidget()
    tabs.setMinimumWidth(120)
    for page, title in ((profiles, 'Gene expression'), (StudyDesignPanel(), 'Studies & samples'),
                        (PublishedExpressionPanel(), 'Published expression'),
                        (SourceMatricesPanel(), 'Deposited matrices'),
                        (DepositedValuesPanel(), 'Other measurements'),
                        (DatasetCoveragePanel(), 'Dataset coverage'),
                        (panel, 'Condition comparison')):
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setMinimumWidth(120)
        scroll.setWidget(page)
        tabs.addTab(scroll, title)
    dock.setWidget(tabs)
    window.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, dock)
    anchor = getattr(window, "discoveries_dock", None) or getattr(window, "strategies_dock", None)
    if anchor is not None:
        window.tabifyDockWidget(anchor, dock)
    if getattr(window, "right_dock", None) is not None:
        window.right_dock.raise_()
    window.compare, window.compare_dock = panel, dock
    window.rna_profiles = profiles
    return panel
