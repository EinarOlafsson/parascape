"""Browse one original-value matrix at a time without conflating probes, assays or scales."""
from PyQt6 import QtWidgets

from . import source_matrices as S


class SourceMatricesPanel(QtWidgets.QWidget):
    """Inspect exact source IDs and export original values with sample and assay provenance."""

    def __init__(self, root=None, parent=None):
        super().__init__(parent)
        self.root, self.result = root, None
        layout = QtWidgets.QVBoxLayout(self)
        description = QtWidgets.QLabel('Original submitted VALUE matrices. Probe/gene IDs retain source spelling; gene remapping is not applied. '
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
        self.table = QtWidgets.QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(['source ID', 'study', 'sample', 'sample title', 'submitted VALUE'])
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
            frame = S.profiles(info['path'], self.probes.toPlainText())
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return
        self.result = frame
        self.table.setRowCount(len(frame))
        for i, row in enumerate(frame.to_dict('records')):
            for j, field in enumerate(('probe_id', 'study', 'geo_sample', 'sample_title', 'value')):
                item = QtWidgets.QTableWidgetItem(str(row[field]))
                item.setToolTip(row['sample_metadata'] + '\nSource: ' + row['source_url'])
                self.table.setItem(i, j, item)
        missing = set(self.probes.toPlainText().splitlines()) - set(frame.probe_id)
        self.summary.setText(f'{len(frame)} original-value records. ' + ('Source IDs not measured: ' + ', '.join(sorted(missing)) if missing else ''))
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
