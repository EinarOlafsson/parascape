"""Browse original-value matrices with source context visible beside measurements.

Public matrices can mix platforms, assay quantities and sample-channel designs.
Selecting one submitted representation at a time keeps those differences visible
instead of pooling values on incompatible scales. Exact identifier queries and
optional source annotation matches preserve individual probes and ambiguous
associations. CSV exports retain sample metadata and provenance so downstream
review can establish the appropriate units and independent experimental samples
before quantitative comparisons are made.
"""
from PyQt6 import QtWidgets

from . import source_matrices as S


class SourceMatricesPanel(QtWidgets.QWidget):
    """Inspect exact source IDs and export original values with sample and assay provenance."""

    def __init__(self, root=None, parent=None):
        super().__init__(parent)
        self.root, self.result = root, None
        layout = QtWidgets.QVBoxLayout(self)
        description = QtWidgets.QLabel('Original submitted VALUE matrices. Original probes and values are retained. Optional annotation matches reproduce source-declared associations. '
            'Units and scale require source review. Assays and platforms remain separate; missing values are unknown. '
            'Channels, technical repeats and shared samples are not extra biological replicates.')
        description.setWordWrap(True)
        layout.addWidget(description)
        self.pack = QtWidgets.QComboBox()
        self.pack.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.pack.setMinimumContentsLength(12)
        self.pack.setToolTip('Choose a source matrix, retaining its study, assay, platform-specific file and submitted samples.')
        self.pack.currentIndexChanged.connect(self._clear)
        layout.addWidget(self.pack)
        self.refresh_button = QtWidgets.QPushButton('Refresh source matrices')
        self.refresh_button.setToolTip('Find completed original-value packs as background conversion finishes.')
        self.refresh_button.clicked.connect(self.refresh)
        layout.addWidget(self.refresh_button)
        self.probes = QtWidgets.QPlainTextEdit()
        self.probes.setPlaceholderText('One exact original probe or submitted gene identifier per line')
        self.probes.setToolTip('Case-sensitive original identifiers; spaces within an identifier are retained. No aliases or probe-to-gene mappings are inferred.')
        self.probes.setMaximumHeight(80)
        layout.addWidget(self.probes)
        self.annotation = QtWidgets.QCheckBox('Match source platform identifiers')
        self.annotation.setToolTip('Match exact source-declared identifiers through independently reconciled platform annotations. Keep every probe and ambiguous source association; current gene equivalence is not inferred.')
        self.annotation.toggled.connect(self._clear)
        layout.addWidget(self.annotation)
        self.show_button = QtWidgets.QPushButton('Show original values')
        self.show_button.setToolTip('Verify the selected pack and show submitted values for the exact source identifiers.')
        self.show_button.clicked.connect(self.compare)
        layout.addWidget(self.show_button)
        self.export_button = QtWidgets.QPushButton('Export source values…')
        self.export_button.setToolTip('Save original values with source sample metadata, assay, quantity, scale limits and checksums.')
        self.export_button.clicked.connect(lambda: self.save())
        layout.addWidget(self.export_button)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(['query ID', 'source probe ID', 'study', 'sample', 'sample title', 'submitted VALUE'])
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        layout.addWidget(self.table)
        self.refresh()
        self.setMinimumWidth(120)

    def _clear(self):
        self.result = None
        if hasattr(self, 'table'):
            self.table.setRowCount(0)
        if hasattr(self, 'summary'):
            self.summary.setText('Enter exact source identifiers for the selected matrix.')

    def refresh(self):
        """Reload completed pack manifests while preserving assay and file identity."""
        self._clear()
        self.pack.clear()
        packs, errors = S.catalogue(self.root)
        for pack in packs:
            self.pack.addItem(f"{pack['study']} · {pack['readout']} · {pack['samples']} submitted samples · {pack['source_url'].split('/')[-1]}", pack)
        self.summary.setText(f'{len(packs)} original-value matrices available. ' +
                             (f'{len(errors)} pack manifests require review.' if errors else 'Enter exact source identifiers.'))
        self.show_button.setEnabled(bool(packs))

    def compare(self):
        """Show one selected source representation without mixing assay units."""
        info = self.pack.currentData()
        if not info:
            return
        self._clear()
        try:
            method = S.annotated_profiles if self.annotation.isChecked() else S.profiles
            frame = method(info['path'], self.probes.toPlainText())
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return
        self.result = frame
        self.table.setRowCount(len(frame))
        for i, row in enumerate(frame.to_dict('records')):
            for j, field in enumerate(('query_identifier', 'probe_id', 'study', 'geo_sample', 'sample_title', 'value')):
                item = QtWidgets.QTableWidgetItem(str(row.get(field, row['probe_id'])))
                item.setToolTip(row['sample_metadata'] + '\nSource: ' + row['source_url'] + '\n' +
                                row.get('annotation_evidence', '') + '\n' + row.get('identifier_interpretation', ''))
                self.table.setItem(i, j, item)
        matched = frame['query_identifier'] if 'query_identifier' in frame else frame.probe_id
        missing = set(self.probes.toPlainText().splitlines()) - set(matched)
        self.summary.setText(f'{len(frame)} original-value records. ' + ('Identifiers not matched: ' + ', '.join(sorted(missing)) if missing else ''))
        self.table.resizeColumnsToContents()

    def save(self, path='') -> int:
        """Export original values, source identities and comparison limits to CSV."""
        if self.result is None or self.result.empty:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export original source values', 'source_values.csv', 'CSV (*.csv)')
        if not path:
            return 0
        self.result.to_csv(path, index=False)
        return len(self.result)
