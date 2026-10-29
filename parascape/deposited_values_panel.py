"""Browse independently reconciled deposited values with their original source context.

Other quantitative assays cannot safely inherit gene-expression units or sample
group assumptions. This view selects one source analysis and shows exact feature
labels, row numbers, literal tokens and submitted sample labels together. Repeated
feature labels retain separate rows. Exported values and whole-matrix numerical
QC carry provenance, while biological independence, compatible scales and feature
equivalence remain questions for source review rather than hidden UI defaults.
"""
from PyQt6 import QtWidgets

from . import deposited_values as D


class DepositedValuesPanel(QtWidgets.QWidget):
    """Inspect one reconciled source representation and export values or sample QC."""

    def __init__(self, root=None, parent=None):
        super().__init__(parent)
        self.root, self.result = root, None
        layout = QtWidgets.QVBoxLayout(self)
        label = QtWidgets.QLabel('Original deposited measurements. Feature rows, value tokens and source units are retained. '
            'Choose one analysis; assays, repeated names and submitted columns remain separate. '
            'Cross-study compatibility and biological independence require review.')
        label.setWordWrap(True)
        layout.addWidget(label)
        self.pack = QtWidgets.QComboBox()
        self.pack.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.pack.setMinimumContentsLength(12)
        self.pack.setToolTip('Choose one independently reconciled source analysis with its original quantity and units.')
        self.pack.currentIndexChanged.connect(self._clear)
        layout.addWidget(self.pack)
        self.refresh_button = QtWidgets.QPushButton('Refresh deposited measurements')
        self.refresh_button.setToolTip('Reload the validated measurement catalogue; incomplete preparations remain unavailable.')
        self.refresh_button.clicked.connect(self.refresh)
        layout.addWidget(self.refresh_button)
        self.features = QtWidgets.QPlainTextEdit()
        self.features.setMaximumHeight(80)
        self.features.setPlaceholderText('One exact original feature label per line')
        self.features.setToolTip('Match case-sensitive source labels, retaining every matching row. No aliases, feature collapse or gene mappings are inferred.')
        layout.addWidget(self.features)
        for attribute, caption, tooltip, callback in [
            ('show_button', 'Show original measurements', 'Verify source hashes and show original values and literal source tokens.', self.show_values),
            ('export_button', 'Export original measurements…', 'Export each original row/sample value with its source token, units and provenance.', lambda: self.save()),
            ('qc_button', 'Export column technical QC…', 'Check all source rows independently of the feature query. Biological quality scores remain unavailable.', lambda: self.save_qc())]:
            button = QtWidgets.QPushButton(caption)
            button.setToolTip(tooltip)
            button.clicked.connect(callback)
            setattr(self, attribute, button)
            layout.addWidget(button)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(['source row', 'source feature', 'source column', 'numeric value', 'original token', 'source units'])
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        layout.addWidget(self.table)
        self.refresh()
        self.setMinimumWidth(120)

    def _clear(self):
        self.result = None
        if hasattr(self, 'table'):
            self.table.setRowCount(0)
        if hasattr(self, 'summary'):
            self.summary.setText('Enter exact feature labels for the selected source analysis.')

    def refresh(self):
        """Load completed source analyses without mixing quantities or units."""
        self._clear()
        self.pack.clear()
        packs, errors = D.catalogue(self.root)
        for pack in packs:
            self.pack.addItem(f"{pack['study']} · {pack['analysis']} · {pack['assay']} · {pack['samples']} source columns · {pack['units']}", pack)
        self.summary.setText(f'{len(packs)} reconciled deposited measurement packs. ' +
                             (f'{len(errors)} catalogue entries require review.' if errors else 'Enter exact source feature labels.'))
        self.show_button.setEnabled(bool(packs))
        self.qc_button.setEnabled(bool(packs))

    def show_values(self):
        """Show exact source rows and original tokens for the selected feature labels."""
        info = self.pack.currentData()
        if not info:
            return
        self._clear()
        try:
            frame = D.profiles(info['path'], self.features.toPlainText())
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return
        self.result = frame
        self.table.setRowCount(len(frame))
        for i, row in enumerate(frame.to_dict('records')):
            for j, field in enumerate(('source_row', 'feature_label', 'source_sample', 'value', 'source_token', 'units')):
                item = QtWidgets.QTableWidgetItem('' if row[field] is None else str(row[field]))
                item.setToolTip(row['sample_metadata'] + '\n' + row.get('source_row_metadata', '') + '\nSource: ' + row['source_url'] + '\n' + row['comparison_scope'])
                self.table.setItem(i, j, item)
        missing = set(self.features.toPlainText().splitlines()) - set(frame.feature_label)
        self.summary.setText(f'{len(frame)} original row/column measurements. ' +
            ('Labels not matched: ' + ', '.join(sorted(missing)) if missing else 'Repeated labels retain separate source rows.'))
        self.table.resizeColumnsToContents()

    def save(self, path='') -> int:
        """Export original row/sample measurements and provenance without collapsing names."""
        if self.result is None or self.result.empty:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export deposited measurements', 'deposited_measurements.csv', 'CSV (*.csv)')
        if not path:
            return 0
        self.result.to_csv(path, index=False)
        return len(self.result)

    def save_qc(self, path='') -> int:
        """Export all-source-row numerical QC independently of the current feature query."""
        info = self.pack.currentData()
        if not info:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export sample technical QC', 'deposited_sample_qc.csv', 'CSV (*.csv)')
        if not path:
            return 0
        try:
            frame = D.technical_qc(info['path'])
            frame.to_csv(path, index=False)
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return 0
        self.summary.setText(f'Exported technical QC for {len(frame)} source columns. Biological quality scores remain unavailable.')
        return len(frame)
