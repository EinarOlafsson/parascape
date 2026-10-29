"""Compare exact gene accessions across the bundled RNA-seq samples, offline.

Researchers need to see the same gene across samples without treating a missing measurement as zero or silently pooling different protocols. The browser keeps study identity and quantification choices visible, exposes the snapshot coverage, and exports the measurements with their sample metadata.
"""
from __future__ import annotations

from PyQt6 import QtCore, QtWidgets

from . import expression_contrast as X
from . import rna_profiles as R


class RNAProfilesPanel(QtWidgets.QWidget):
    """A searchable, exportable per-sample expression table with explicit coverage."""

    gene_chosen = QtCore.pyqtSignal(str)

    def __init__(self, root: str = X.NAS, parent=None):
        super().__init__(parent)
        self.root = root
        self.result = None
        layout = QtWidgets.QVBoxLayout(self)
        info = R.snapshot(root)
        self.coverage = QtWidgets.QLabel(
            f"RNA-seq snapshot: {info.get('verified_runs', 'external')} verified runs; "
            f"{info.get('genes', 'external')} measured genes. "
            "Partial studies may lack controls. TPM comparisons across studies are descriptive; "
            "protocol and batch effects remain. Use Condition comparison for within-study groups.")
        self.coverage.setWordWrap(True)
        layout.addWidget(self.coverage)
        self.genes = QtWidgets.QLineEdit()
        self.genes.setPlaceholderText('Toxoplasma gene accessions, separated by spaces or commas')
        self.genes.setToolTip('Enter one or more exact reference gene IDs. Missing genes are reported; zeros remain measured zeros.')
        layout.addWidget(self.genes)
        row = QtWidgets.QHBoxLayout()
        self.study = QtWidgets.QComboBox()
        self.study.setToolTip('Limit expression measurements to one study, or browse all completed studies. Cross-study TPM retains protocol and batch effects.')
        self.study.addItem('All studies', None)
        available = X.available(root)
        for study in sorted(set(available.study)):
            self.study.addItem(study, study)
        row.addWidget(self.study)
        self.quantifier = QtWidgets.QComboBox()
        self.quantifier.setToolTip('Choose Salmon or kallisto estimates for the same submitted reads and competitive reference.')
        self.quantifier.addItems(list(X.QUANTIFIERS))
        row.addWidget(self.quantifier)
        self.reference = QtWidgets.QComboBox()
        self.reference.setToolTip('Reference transcript set. ME49_host allows parasite and human/mouse host transcripts to compete for reads.')
        for ref in sorted(set(available.reference)):
            self.reference.addItem(ref)
        self.reference.setCurrentText('ME49_host')
        row.addWidget(self.reference)
        layout.addLayout(row)
        row = QtWidgets.QHBoxLayout()
        self.show_button = QtWidgets.QPushButton('Compare genes')
        self.show_button.setToolTip('Show per-sample TPM and estimated counts for the entered exact gene accessions; report unmeasured genes separately.')
        self.show_button.clicked.connect(self.compare)
        self.genes.returnPressed.connect(self.compare)
        row.addWidget(self.show_button)
        self.export_button = QtWidgets.QPushButton('Export CSV…')
        self.export_button.setToolTip('Save all displayed gene measurements with study, sample and run metadata to CSV.')
        self.export_button.clicked.connect(lambda: self.save())
        row.addWidget(self.export_button)
        self.qc_button = QtWidgets.QPushButton('Export sample technical QC…')
        self.qc_button.setToolTip('Export whole-matrix numerical validity, missingness and measured zeros for the selected studies, quantifier and reference, with source hashes. Runs remain separate.')
        self.qc_button.clicked.connect(lambda: self.save_qc())
        row.addWidget(self.qc_button)
        layout.addLayout(row)
        self.summary = QtWidgets.QLabel('Enter gene accessions to compare their measured expression across samples.')
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(['gene', 'study', 'run', 'sample', 'parasite TPM', 'estimated count'])
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        self.table.cellDoubleClicked.connect(self._gene)
        layout.addWidget(self.table, 1)

    def compare(self):
        """Show measurements without silently dropping unmeasured query genes."""
        import re
        requested = list(dict.fromkeys(g.upper() for g in re.split(r'[\s,;]+', self.genes.text().strip()) if g))
        try:
            self.result = R.profiles(requested, self.quantifier.currentText(),
                                     self.reference.currentText(), self.root, self.study.currentData())
        except (ValueError, OSError) as exc:
            self.result = None
            self.table.setRowCount(0)
            self.summary.setText(str(exc))
            return
        found = set(self.result.gene)
        missing = [g for g in requested if g not in found]
        self.summary.setText(f"{len(found)} genes, {len(self.result)} measurements. "
                             + ("Not measured in this selection: " + ', '.join(missing) if missing else ''))
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(self.result))
        for i, row in enumerate(self.result.to_dict('records')):
            values = [row['gene'], row['study'], row['run_accession'], row.get('sample_title', ''),
                      f"{row['tpm']:.4g}", f"{row['estimated_count']:.4g}"]
            for j, value in enumerate(values):
                self.table.setItem(i, j, QtWidgets.QTableWidgetItem(str(value)))
        self.table.resizeColumnsToContents()
        self.table.setSortingEnabled(True)

    def save(self, path: str = '') -> int:
        """Export all selected measurements and their sample metadata to CSV."""
        if self.result is None or self.result.empty:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export RNA profiles', 'rna_profiles.csv', 'CSV (*.csv)')
        if not path:
            return 0
        self.result.to_csv(path, index=False)
        return len(self.result)

    def save_qc(self, path: str = '') -> int:
        """Export QC for complete selected matrices independently of queried genes."""
        try:
            frame = R.technical_qc(self.quantifier.currentText(), self.reference.currentText(),
                                   self.root, self.study.currentData())
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return 0
        if frame.empty:
            self.summary.setText('No matrices available for technical QC in this selection.')
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export RNA sample technical QC',
                                                         'rna_sample_technical_qc.csv', 'CSV (*.csv)')
        if not path:
            return 0
        frame.to_csv(path, index=False)
        return len(frame)

    def _gene(self, row, column):
        item = self.table.item(row, 0)
        if item:
            self.gene_chosen.emit(item.text())
