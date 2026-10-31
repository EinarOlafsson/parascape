"""Browse original sparse measurements with background reads and complete exports."""
from __future__ import annotations

import json
from pathlib import Path
import threading

from PyQt6 import QtCore, QtWidgets

from . import sparse_values as S
from .sparse_query import prepare_query, QueryCancelled


class _Signals(QtCore.QObject):
    finished = QtCore.pyqtSignal(int, object, str)


class _Query(QtCore.QRunnable):
    def __init__(self, generation, descriptor, labels, coordinate_root, axis_root):
        super().__init__()
        self.signals = _Signals()
        self.generation = generation
        self.descriptor, self.labels = descriptor, labels
        self.coordinate_root, self.axis_root = coordinate_root, axis_root
        self.cancel = threading.Event()

    def run(self):
        result, error = None, ''
        try:
            result = prepare_query(self.descriptor, self.labels, coordinate_root=self.coordinate_root,
                                   axis_root=self.axis_root, cancel=self.cancel)
        except QueryCancelled:
            error = 'Query cancelled.'
        except Exception as exception:
            error = str(exception)
        self.signals.finished.emit(self.generation, result, error)


class SparseValuesPanel(QtWidgets.QWidget):
    """Select one original sparse source and inspect exact labels without pooling studies."""

    def __init__(self, root=None, parent=None, *, coordinate_root=None, axis_root=None):
        super().__init__(parent)
        self.root, self.coordinate_root, self.axis_root = root, coordinate_root, axis_root
        self.result = None
        self._worker = None
        self._generation = 0
        self._catalogue = []
        layout = QtWidgets.QVBoxLayout(self)
        description = QtWidgets.QLabel(
            'Original sparse measurements. Choose a source and enter exact feature labels. '
            'Repeated features and stored coordinates remain separate. Source units, cell filtering, '
            'gene equivalence and biological sample independence require review.')
        description.setWordWrap(True)
        layout.addWidget(description)
        self.source_filter = QtWidgets.QLineEdit()
        self.source_filter.setPlaceholderText('Filter by study, sample or source filename')
        self.source_filter.setToolTip('Filter source metadata without changing case-sensitive feature matching.')
        self.source_filter.textChanged.connect(self._filter_sources)
        layout.addWidget(self.source_filter)
        self.pack = QtWidgets.QComboBox()
        self.pack.setMinimumContentsLength(12)
        self.pack.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.pack.setToolTip('Choose one source matrix. Original source files and quantity families remain distinct.')
        self.pack.currentIndexChanged.connect(self._clear)
        layout.addWidget(self.pack)
        self.refresh_button = QtWidgets.QPushButton('Refresh sparse measurements')
        self.refresh_button.setToolTip('Include newly completed matrix and metadata bundles.')
        self.refresh_button.clicked.connect(self.refresh)
        layout.addWidget(self.refresh_button)
        self.catalogue_status = QtWidgets.QLabel()
        self.catalogue_status.setWordWrap(True)
        layout.addWidget(self.catalogue_status)
        self.features = QtWidgets.QPlainTextEdit()
        self.features.setMaximumHeight(80)
        self.features.setPlaceholderText('One exact original first-field feature label per line')
        self.features.setToolTip('Keep original case, spaces and suffixes. Every matching original row is retained; aliases and canonical genes are not inferred.')
        self.features.textChanged.connect(self._clear)
        layout.addWidget(self.features)
        self.show_button = QtWidgets.QPushButton('Show sparse measurements')
        self.show_button.setToolTip('Verify checksums in the background, stream all matching records and display the first 500. No zeros are filled and repeated coordinates are not summed.')
        self.show_button.clicked.connect(self.show_values)
        layout.addWidget(self.show_button)
        self.cancel_button = QtWidgets.QPushButton('Cancel query')
        self.cancel_button.setToolTip('Cancel between streamed batches. A checksum read already in progress finishes before cancellation takes effect.')
        self.cancel_button.clicked.connect(self._clear)
        layout.addWidget(self.cancel_button)
        self.export_button = QtWidgets.QPushButton('Export all matching records…')
        self.export_button.setToolTip('Save the complete CSV, including records beyond the preview, and a checksum/source/treatment sidecar.')
        self.export_button.clicked.connect(lambda: self.save())
        layout.addWidget(self.export_button)
        self.context_button = QtWidgets.QPushButton('Inspect selected measurement…')
        self.context_button.setToolTip('Show original coordinates, tokens and all submitted feature/cell metadata for the selected displayed record.')
        self.context_button.clicked.connect(self.show_context)
        layout.addWidget(self.context_button)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.table = QtWidgets.QTableWidget(0, 9)
        self.table.setHorizontalHeaderLabels(['entry index', 'source row', 'source column', 'feature position',
                                             'source feature', 'cell position', 'source cell label', 'numeric value', 'original token'])
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.currentCellChanged.connect(lambda *_: self.context_button.setEnabled(self.result is not None and self.table.currentRow() >= 0))
        self.table.cellDoubleClicked.connect(lambda *_: self.show_context())
        layout.addWidget(self.table)
        self.refresh()
        self.setMinimumWidth(120)

    def _buttons(self):
        self.show_button.setEnabled(bool(self.pack.currentData()) and self._worker is None)
        self.cancel_button.setEnabled(self._worker is not None)
        self.export_button.setEnabled(self.result is not None)
        self.context_button.setEnabled(self.result is not None and self.table.currentRow() >= 0)

    def _clear(self):
        self._generation += 1
        if self._worker is not None:
            self._worker.cancel.set()
        if self.result is not None:
            self.result.close()
            self.result = None
        if getattr(self, 'context_dialog', None) is not None:
            self.context_dialog.close()
            self.context_dialog.deleteLater()
            self.context_dialog = None
        if hasattr(self, 'table'):
            self.table.setRowCount(0)
            self.summary.setText('Query changed; waiting for the current read to finish.' if self._worker is not None else
                                 'Enter exact original feature labels for this source.')
            self._buttons()

    def refresh(self):
        self._catalogue, self._errors = S.catalogue(self.root)
        self._filter_sources()

    def _filter_sources(self):
        previous = (self.pack.currentData() or {}).get('path')
        text = self.source_filter.text().casefold().strip()
        sources = [item for item in self._catalogue if text in json.dumps(item['source_record'], ensure_ascii=False).casefold()]
        blocker = QtCore.QSignalBlocker(self.pack)
        self.pack.clear()
        for item in sources:
            record = item['source_record']
            self.pack.addItem(f"{record['studies']} · {record['source_entities']} · {Path(item['source_url']).name}", item)
        for index, item in enumerate(sources):
            if item['path'] == previous:
                self.pack.setCurrentIndex(index)
                break
        del blocker
        self._clear()
        self.catalogue_status.setText(f'{len(sources)} of {len(self._catalogue)} checked sparse sources shown. '
                                     + (f'{len(self._errors)} entries require review.' if self._errors else 'Choose one source; quantities are not pooled.'))
        self.catalogue_status.setToolTip('\n'.join(self._errors))
        if not sources:
            self.summary.setText('No sparse sources match this filter.' if self._catalogue else
                                 'No completed sparse measurement collection is installed.')

    def show_values(self):
        item = self.pack.currentData()
        if item is None or self._worker is not None:
            return
        labels = [line for line in self.features.toPlainText().splitlines() if line != '']
        self._clear()
        if not labels:
            self.summary.setText('Enter at least one exact original feature label.')
            return
        self._worker = _Query(self._generation, item['path'], labels, self.coordinate_root, self.axis_root)
        self._worker.signals.finished.connect(self._finished)
        self._buttons()
        self.summary.setText('Checking source files and reading matching records in the background…')
        QtCore.QThreadPool.globalInstance().start(self._worker)

    @QtCore.pyqtSlot(int, object, str)
    def _finished(self, generation, result, error):
        self._worker = None
        if generation != self._generation:
            if result is not None:
                result.close()
            self.summary.setText('Previous query discarded. Enter labels for the current source.')
            self._buttons()
            return
        if error:
            self.summary.setText(error)
            self._buttons()
            return
        self.result = result
        fields = ['entry_index', 'source_row', 'source_column', 'source_feature_axis_row', 'source_feature_label',
                  'source_cell_axis_row', 'source_cell_label', 'value', 'source_token']
        self.table.setRowCount(len(result.preview))
        for i, row in enumerate(result.preview):
            tooltip = json.dumps(dict(original_feature_fields=json.loads(row['original_feature_fields']),
                                      original_cell_fields=row['original_cell_fields'],
                                      source_url=result.record['source_url']), indent=2, ensure_ascii=False)
            for j, field in enumerate(fields):
                item = QtWidgets.QTableWidgetItem(str(row[field]))
                item.setToolTip(tooltip)
                self.table.setItem(i, j, item)
        self.table.resizeColumnsToContents()
        for column in range(self.table.columnCount()):
            self.table.setColumnWidth(column, self.table.columnWidth(column) + 12)
        missing = result.record['unmatched_labels']
        self.summary.setText(f"{result.record['exported_records']} original stored records; {len(result.preview)} displayed. "
                             'Export includes every matching record. Absent coordinates are not filled. '
                             + ('Labels not matched: ' + ', '.join(missing) if missing else ''))
        self._buttons()

    def save(self, path='') -> int:
        if self.result is None:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export complete sparse query', 'sparse_records.csv', 'CSV (*.csv)')
        if not path:
            return 0
        try:
            count = self.result.export(path)
        except (OSError, ValueError) as error:
            self.summary.setText(str(error))
            return 0
        self.summary.setText(f'Exported all {count} stored records and treatment record to {path}. '
                             f'The table shows {len(self.result.preview)} records.')
        return count

    def show_context(self):
        row = self.table.currentRow()
        if self.result is None or row < 0:
            return
        if getattr(self, 'context_dialog', None) is not None:
            self.context_dialog.close()
            self.context_dialog.deleteLater()
        self.context_dialog = QtWidgets.QDialog(self)
        self.context_dialog.setWindowTitle('Original sparse measurement context')
        layout = QtWidgets.QVBoxLayout(self.context_dialog)
        tree = QtWidgets.QTreeWidget()
        tree.setHeaderLabels(['Original field / context', 'Value'])
        tree.setToolTip('Submitted fields and their original positions are retained; biological roles and gene equivalence are not inferred.')
        measurement = self.result.preview[row]
        record = self.result.record

        def group(title, pairs):
            node = QtWidgets.QTreeWidgetItem([title, ''])
            tree.addTopLevelItem(node)
            for key, value in pairs:
                child = QtWidgets.QTreeWidgetItem([str(key), str(value)])
                child.setToolTip(1, str(value))
                node.addChild(child)
            node.setExpanded(True)

        group('Original measurement', [(key, measurement[key]) for key in
              ('entry_index', 'source_row', 'source_column', 'row_token', 'column_token',
               'source_feature_axis_row', 'source_feature_label', 'source_cell_axis_row',
               'source_cell_label', 'value', 'source_token')])
        for role, caption, fields in (
                ('features', 'Submitted feature fields', json.loads(measurement['original_feature_fields'])),
                ('barcodes', 'Submitted cell fields', measurement['original_cell_fields'])):
            headers = record['axis_roles'][role]['original_axis_artifacts'][0]['original_headers']
            group(caption, [(headers[index] if index < len(headers) else f'Field {index + 1}', value)
                            for index, value in enumerate(fields)])
        group('Source', [(key, value) for key, value in record['source_record'].items()]
              + [('Original source URL', record['source_url'])])
        group('Treatment', [('Feature coordinate', record['operations'][0]['source_coordinate']),
                            ('Original positions selected', record['operations'][0]['positions']),
                            ('Normalization / imputation / transposition', 'None applied'),
                            ('Repeated coordinates', 'Retained separately; no summation'),
                            ('Absent coordinates', 'No value filled or interpretation inferred'),
                            ('Biological comparability', 'Requires separate review')])
        group('Provenance', [(key, record[key]) for key in
                            ('descriptor_sha256', 'source_sha256', 'reader_sha256', 'query_materializer_sha256')])
        tree.resizeColumnToContents(0)
        layout.addWidget(tree)
        self.context_dialog.resize(700, 500)
        self.context_dialog.show()

    def closeEvent(self, event):
        self._clear()
        super().closeEvent(event)
