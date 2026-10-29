"""Browse original deposited expression with explicit sample and annotation identity.

Published tables may use older annotations and different quantifiers from the raw-read
pipeline. This view retains those choices alongside source units and sample labels,
so alternative processing of the same samples cannot become extra replicates or an
unlabelled cross-study comparison. Export keeps the source metadata and checksums.
"""
from __future__ import annotations

from PyQt6 import QtWidgets

from . import published_expression as P


class PublishedExpressionPanel(QtWidgets.QWidget):
    """Query and export validated source measurements without pooling processing variants."""

    def __init__(self, root=None, parent=None):
        super().__init__(parent)
        self.result = None
        layout = QtWidgets.QVBoxLayout(self)
        self.description = QtWidgets.QLabel(
            'Original deposited expression. Annotation variants reuse the same samples; '
            'they are not extra replicates. Reported TPM, FPKM and expected counts retain '
            'their original units and normalization. Different pipelines remain separate.')
        self.description.setWordWrap(True)
        layout.addWidget(self.description)
        self.pack = QtWidgets.QComboBox()
        self.pack.setToolTip('Choose a validated study pack of original deposited measurements.')
        packs, errors = P.catalogue(root)
        for pack in packs:
            self.pack.addItem(f"{pack['study']} · {pack['samples']} submitted samples · {pack['path'].split('/')[-1]}", pack)
        layout.addWidget(self.pack)
        self.reference = QtWidgets.QComboBox()
        self.reference.setToolTip('Choose the source annotation. Annotation variants reuse the same samples.')
        layout.addWidget(self.reference)
        self.genes = QtWidgets.QLineEdit()
        self.genes.setPlaceholderText('Exact source gene identifiers, separated by spaces or commas')
        layout.addWidget(self.genes)
        self.compare_button = QtWidgets.QPushButton('Show deposited measurements')
        self.compare_button.setToolTip('Show original measurements for the entered gene identifiers in the selected study and annotation.')
        self.compare_button.clicked.connect(self.compare)
        self.genes.returnPressed.connect(self.compare)
        layout.addWidget(self.compare_button)
        self.export_button = QtWidgets.QPushButton('Export CSV…')
        self.export_button.setToolTip('Save the displayed measurements with their units, sample metadata and source provenance.')
        self.export_button.clicked.connect(lambda: self.save())
        layout.addWidget(self.export_button)
        self.qc_button = QtWidgets.QPushButton('Export sample technical QC…')
        self.qc_button.setToolTip('Export whole-matrix numerical validity, missingness and measured zeros with source hashes. Biological quality scores are not calibrated.')
        self.qc_button.clicked.connect(lambda: self.save_qc())
        layout.addWidget(self.qc_button)
        self.summary = QtWidgets.QLabel('Enter exact identifiers from the selected annotation.')
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.table = QtWidgets.QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(['gene', 'sample', 'sample title', 'reference', 'TPM', 'FPKM', 'expected count'])
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        layout.addWidget(self.table)
        self.pack.currentIndexChanged.connect(self._select_pack)
        self.reference.currentIndexChanged.connect(self._clear)
        self._select_pack()
        if errors:
            warning = QtWidgets.QLabel('Some source packs could not be verified: ' + '; '.join(errors))
            warning.setWordWrap(True)
            layout.addWidget(warning)
        if not packs:
            self.summary.setText('No validated deposited expression packs are installed. The Gene expression tab contains common-pipeline measurements.')
            self.compare_button.setEnabled(False)
        self.setMinimumWidth(120)

    def _clear(self):
        self.result = None
        self.table.setRowCount(0)
        self.summary.setText('Enter exact identifiers from the selected annotation.')

    def save_qc(self, path='') -> int:
        """Export full-matrix QC independently of the current gene query."""
        info = self.pack.currentData()
        if not info:
            return 0
        try:
            frame = P.technical_qc(info['path'], self.reference.currentText())
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export sample technical QC',
                                                         'sample_technical_qc.csv', 'CSV files (*.csv)')
        if not path:
            return 0
        frame.to_csv(path, index=False)
        return len(frame)

    def _select_pack(self):
        self.reference.clear()
        info = self.pack.currentData()
        if info:
            self.reference.addItems(sorted(info['variants']))
        self._clear()

    def compare(self):
        """Show one original source measurement per exact gene and submitted sample."""
        info = self.pack.currentData()
        if not info:
            return
        try:
            frame = P.profiles(info['path'], self.genes.text(), self.reference.currentText())
        except (OSError, ValueError, KeyError, TypeError) as error:
            self._clear()
            self.summary.setText(str(error))
            return
        self.result = frame
        self.table.setRowCount(len(frame))
        for i, row in enumerate(frame.to_dict('records')):
            for j, field in enumerate(['gene', 'sample', 'sample_title', 'reference', 'TPM', 'FPKM', 'expected_count']):
                item = QtWidgets.QTableWidgetItem(str(row[field]))
                item.setToolTip(row['sample_metadata'] + '\nSource: ' + row['source_url'])
                self.table.setItem(i, j, item)
        import re
        requested = set(re.split(r'[\s,;]+', self.genes.text().strip())) - {''}
        missing = requested - set(frame.gene)
        self.summary.setText(f'{len(frame)} measurements from {frame["sample"].nunique()} submitted samples. '
                             + ('Not measured: ' + ', '.join(sorted(missing)) if missing else ''))
        self.table.resizeColumnsToContents()

    def save(self, path='') -> int:
        """Export source values, software, metadata, URLs and checksums to CSV."""
        if self.result is None or self.result.empty:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export deposited expression', 'published_expression.csv', 'CSV (*.csv)')
        if not path:
            return 0
        self.result.to_csv(path, index=False)
        return len(self.result)
