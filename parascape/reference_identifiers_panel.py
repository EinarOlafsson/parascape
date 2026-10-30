"""Browse source reference identifiers while keeping reference versions and evidence scope visible.

The original header index helps inspect explicitly declared identifier relationships
before comparing measurements. Each lookup stays inside one acquired source file,
retains repeated header records and exposes literal gene declarations with source
hashes. It does not choose a canonical identifier, infer which search database was
used, convert gene symbols or combine measurements. Pending sources and missing
declarations remain gaps rather than evidence of biological absence.
"""
from PyQt6 import QtWidgets

from . import reference_identifiers as R


class ReferenceIdentifiersPanel(QtWidgets.QWidget):
    """Inspect exact original header identifiers and source gene declarations with CSV provenance."""

    def __init__(self, root=None, parent=None):
        super().__init__(parent)
        self.root, self.result = root, None
        layout = QtWidgets.QVBoxLayout(self)
        label = QtWidgets.QLabel('Original reference header identifiers. Choose a source file and an exact identifier. '
            'Explicit gene declarations retain source meanings; current gene versions and measurement database use remain unreviewed.')
        label.setWordWrap(True)
        layout.addWidget(label)
        self.source = QtWidgets.QComboBox()
        self.source.setToolTip('Choose an acquired reference file with independently reconciled original header identifiers.')
        self.source.setMinimumContentsLength(12)
        self.source.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.source.currentIndexChanged.connect(self._clear)
        layout.addWidget(self.source)
        refresh = QtWidgets.QPushButton('Refresh reference sources')
        refresh.setToolTip('Reload reconciled reference sources and show files still awaiting acquisition or review.')
        refresh.clicked.connect(self.refresh)
        layout.addWidget(refresh)
        self.kind = QtWidgets.QComboBox()
        self.kind.setToolTip('Choose exact original primary identifiers or explicit source gene declarations; identifier versions remain distinct.')
        self.kind.addItem('Original primary identifier', 'primary')
        self.kind.addItem('Explicit source gene declaration', 'gene')
        self.kind.currentIndexChanged.connect(self._clear)
        layout.addWidget(self.kind)
        self.query = QtWidgets.QLineEdit()
        self.query.setPlaceholderText('One exact original identifier; case and version suffixes are retained')
        self.query.textChanged.connect(self._clear)
        layout.addWidget(self.query)
        self.show_button = QtWidgets.QPushButton('Show original identifier records')
        self.show_button.setToolTip('Verify the selected source and display every exact matching original header record without collapsing repeated identifiers.')
        self.show_button.clicked.connect(self.show_records)
        layout.addWidget(self.show_button)
        export = QtWidgets.QPushButton('Export identifier records…')
        export.setToolTip('Export displayed original headers, literal declarations and source hashes as CSV.')
        export.clicked.connect(lambda: self.save())
        layout.addWidget(export)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.table = QtWidgets.QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(['source row', 'original primary identifier', 'source identifier declarations', 'original header'])
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        layout.addWidget(self.table)
        self.refresh()
        self.setMinimumWidth(120)

    def _clear(self):
        self.result = None
        if hasattr(self, 'table'):
            self.table.setRowCount(0)
        if hasattr(self, 'summary'):
            self.summary.setText('Choose an original reference and exact source identifier.')

    def refresh(self):
        """Reload reconciled sources and disclose files still awaiting acquisition or review."""
        self._clear()
        self.source.clear()
        try:
            sources, _ = R.catalogue(self.root)
            for source in sources:
                if source['status'] == 'all_header_cells_reconciled':
                    caption = '; '.join(source['declared_projects']) + ' · ' + source['filename'] + f" · {source['counts']['headers']} headers"
                    self.source.addItem(caption, source)
            pending = sum(source['status'] != 'all_header_cells_reconciled' for source in sources)
            self.summary.setText(f'{self.source.count()} reconciled source files; {pending} awaiting acquisition or review. '
                'Identical file contents may occur in multiple projects; file occurrences are not independent reference versions.')
        except FileNotFoundError:
            self.summary.setText('Original reference identifier data is not installed.')
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
        self.show_button.setEnabled(bool(self.source.count()))

    def show_records(self):
        """Display exact source matches without canonicalisation or biological inference."""
        self._clear()
        source = self.source.currentData()
        if not source:
            return
        try:
            frame = R.lookup(source['source_url'], self.query.text(), kind=self.kind.currentData(), root=self.root)
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return
        self.result = frame
        self.table.setRowCount(len(frame))
        for i, row in enumerate(frame.to_dict('records')):
            for j, key in enumerate(('header_row', 'primary_identifier', 'declared_identifier_fields', 'original_header')):
                item = QtWidgets.QTableWidgetItem(str(row[key]))
                item.setToolTip(str(row[key]) + '\nSource: ' + row['source_url'] + '\n' + row['scope'])
                self.table.setItem(i, j, item)
        self.summary.setText(f'{len(frame)} original header records. ' +
            ('Repeated identifiers retain separate rows; declarations do not establish current gene equivalence.' if len(frame) else
             'No exact declaration in this source file; other references and unsupported fields were not matched.'))
        self.table.resizeColumnsToContents()

    def save(self, path='') -> int:
        """Export original header records and exact declaration/source provenance without collapsing rows."""
        if self.result is None or self.result.empty:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export original reference identifiers', 'reference_identifiers.csv', 'CSV (*.csv)')
        if not path:
            return 0
        self.result.to_csv(path, index=False)
        return len(self.result)
