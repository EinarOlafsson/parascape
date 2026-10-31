"""Browse original sparse measurements with background reads and complete exports."""
from __future__ import annotations

import json
from pathlib import Path
import threading

from PyQt6 import QtCore, QtWidgets

from . import sparse_values as S
from .sparse_query import prepare_query, QueryCancelled
from . import optional_sparse_data as O
from . import sparse_record_qc as Q
from .component_download_dialog import ComponentDownloadDialog


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
        self.unique_coordinates = QtWidgets.QCheckBox('Require zero repeated coordinates')
        self.unique_coordinates.setToolTip('Only show sources whose full coordinate QC proves zero repeated row/column entries. Sources with missing or unknown global counts are excluded while checked. This does not establish biological quality or independent cells.')
        self.unique_coordinates.toggled.connect(self._filter_sources)
        layout.addWidget(self.unique_coordinates)
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
        self.download_button = QtWidgets.QPushButton('Download / verify selected dataset…')
        self.download_button.setToolTip('Fetch missing hosted files or verify and repair a previously downloaded dataset. Explicit local collection overrides are preserved.')
        self.download_button.clicked.connect(self.download_selected)
        layout.addWidget(self.download_button)
        self.export_button = QtWidgets.QPushButton('Export all matching records…')
        self.export_button.setToolTip('Save the complete CSV, including records beyond the preview, and a checksum/source/treatment sidecar.')
        self.export_button.clicked.connect(lambda: self.save())
        layout.addWidget(self.export_button)
        self.qc_button = QtWidgets.QPushButton('Inspect coordinate QC…')
        self.qc_button.setToolTip('Show source-bound ordering and duplicate-record facts. Unordered sources retain unknown global duplicate counts; these facts are not a biological quality score.')
        self.qc_button.clicked.connect(self.show_qc)
        layout.addWidget(self.qc_button)
        self.context_button = QtWidgets.QPushButton('Inspect selected measurement…')
        self.context_button.setToolTip('Show source coordinates, value representations and all submitted feature/cell metadata for the selected displayed record.')
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
        selected = self.pack.currentData()
        self.download_button.setEnabled(bool(selected) and self._worker is None and
                                        (Path(selected['path']).parent.parent / 'data_components.json').exists())
        self.export_button.setEnabled(self.result is not None)
        self.qc_button.setEnabled(bool(self.pack.currentData()))
        self.context_button.setEnabled(self.result is not None and self.table.currentRow() >= 0)

    def _clear(self):
        self._generation += 1
        if self._worker is not None:
            self._worker.cancel.set()
        if self.result is not None:
            self.result.close()
            self.result = None
        if getattr(self, 'qc_dialog', None) is not None:
            self.qc_dialog.close()
            self.qc_dialog.deleteLater()
            self.qc_dialog = None
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
        """Reload completed source descriptors and their coordinate QC."""
        self._catalogue, self._errors = S.catalogue(self.root)
        self._qc_reports, self._qc_errors = {}, []
        for source in self._catalogue:
            try:
                self._qc_reports[source['path']] = Q.report(source['path'])
            except (OSError, ValueError, KeyError, TypeError) as error:
                self._qc_reports[source['path']] = None
                self._qc_errors.append(str(error))
        self._filter_sources()

    def _filter_sources(self):
        previous = (self.pack.currentData() or {}).get('path')
        text = self.source_filter.text().casefold().strip()
        sources = [item for item in self._catalogue if text in json.dumps(
            dict(source=item['source_record'], url=item['source_url'], quantity=item['quantity'],
                 layer=item.get('source_layer')), ensure_ascii=False).casefold()]
        if self.unique_coordinates.isChecked():
            sources = [source for source in sources if self._qc_reports.get(source['path']) is not None
                       and self._qc_reports[source['path']]['facts']['global_duplicate_count_proven']
                       and self._qc_reports[source['path']]['facts']['global_repeated_records'] == 0]
        blocker = QtCore.QSignalBlocker(self.pack)
        self.pack.clear()
        for item in sources:
            record = item['source_record']
            label = f"{record['studies']} · {record['source_entities']} · {Path(item['source_url']).name}"
            if item.get('source_layer'):
                label += f" · layer {item['source_layer']}"
            self.pack.addItem(label, item)
            self.pack.setItemData(self.pack.count() - 1, item['quantity'], QtCore.Qt.ItemDataRole.ToolTipRole)
        for index, item in enumerate(sources):
            if item['path'] == previous:
                self.pack.setCurrentIndex(index)
                break
        del blocker
        self._clear()
        self.catalogue_status.setText(f'{len(sources)} of {len(self._catalogue)} checked sparse sources shown. '
                                     + (f'{len(self._errors)} entries require review.' if self._errors else 'Choose one source; quantities are not pooled.'))
        self.catalogue_status.setToolTip('\n'.join(self._errors + self._qc_errors))
        if not sources:
            self.summary.setText('No sparse sources match this filter.' if self._catalogue else
                                 'No completed sparse measurement collection is installed.')

    def show_values(self):
        """Read every matching original record in the background with a bounded preview."""
        item = self.pack.currentData()
        if item is None or self._worker is not None:
            return
        labels = [line for line in self.features.toPlainText().splitlines() if line != '']
        self._clear()
        if not labels:
            self.summary.setText('Enter at least one exact original feature label.')
            return
        try:
            plan = O.plan(item['path'], coordinate_root=self.coordinate_root, axis_root=self.axis_root)
            generation = self._generation
            if plan['missing']:
                dialog = ComponentDownloadDialog(plan['missing'], self)
                if dialog.exec() != QtWidgets.QDialog.DialogCode.Accepted:
                    self.summary.setText('Dataset download cancelled. Existing files are retained.')
                    return
                if generation != self._generation:
                    self.summary.setText('Source or query changed; choose Show sparse measurements again.')
                    return
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return
        self._worker = _Query(self._generation, item['path'], labels, plan['coordinate_root'], plan['axis_root'])
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
        try:
            selected = self.pack.currentData()
            result.record['coordinate_QC'] = Q.report(selected['path']) if selected else None
        except (OSError, ValueError, KeyError, TypeError) as error:
            result.record['coordinate_QC'] = None
            result.record['coordinate_QC_gap'] = str(error)
        fields = ['entry_index', 'source_row', 'source_column', 'source_feature_axis_row', 'source_feature_label',
                  'source_cell_axis_row', 'source_cell_label', 'value', 'source_token']
        derived_tokens = result.record.get('source_storage', {}).get('decimal_tokens_are_derived', False)
        self.table.horizontalHeaderItem(8).setText('derived decimal value' if derived_tokens else 'original token')
        self.table.setRowCount(len(result.preview))
        for i, row in enumerate(result.preview):
            tooltip = json.dumps(dict(original_feature_fields=json.loads(row['original_feature_fields']),
                                      original_cell_fields=row['original_cell_fields'],
                                      source_url=result.record['source_url']), indent=2, ensure_ascii=False)
            for j, field in enumerate(fields):
                item = QtWidgets.QTableWidgetItem(str(row[field]))
                item.setToolTip(tooltip)
                if j == 8 and derived_tokens:
                    item.setToolTip('Decimal rendering of the original HDF5 binary value; no original text token exists.\n' + tooltip)
                self.table.setItem(i, j, item)
        self.table.resizeColumnsToContents()
        for column in range(self.table.columnCount()):
            self.table.setColumnWidth(column, self.table.columnWidth(column) + 12)
        missing = result.record['unmatched_labels']
        self.summary.setText(f"{result.record['exported_records']} original stored records; {len(result.preview)} displayed. "
                             'Export includes every matching record. Absent coordinates are not filled. ' + self._qc_summary() + ' '
                             + ('Labels not matched: ' + ', '.join(missing) if missing else ''))
        if derived_tokens:
            self.summary.setText(self.summary.text() + ' Value strings are derived decimal representations of original HDF5 binary values.')
        self._buttons()

    def download_selected(self):
        """Offer checksum verification or retrieval of this dataset and its shared axes."""
        item = self.pack.currentData()
        if item is None or self._worker is not None:
            return
        self._clear()
        try:
            plan = O.plan(item['path'], coordinate_root=self.coordinate_root, axis_root=self.axis_root)
            if not plan['components']:
                self.summary.setText('These files are supplied by a local collection. Show measurements to verify them.')
                return
            dialog = ComponentDownloadDialog(plan['components'], self)
            if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
                self.summary.setText('Dataset files verified. Choose Show sparse measurements to browse them.')
            else:
                self.summary.setText('Dataset download cancelled. Existing files are retained.')
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))

    def save(self, path='') -> int:
        """Export all matching records and their source-treatment sidecar."""
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
                             f'The table shows {len(self.result.preview)} records. ' + self._qc_summary())
        return count

    def _qc_summary(self):
        report = self.result.record.get('coordinate_QC') if self.result is not None else None
        if report is None:
            return 'Coordinate QC unavailable; biological score not assigned.'
        facts = report['facts']
        if facts['global_duplicate_count_proven']:
            return f"Repeated coordinate records: {facts['global_repeated_records']} (complete-source check)."
        return 'Global repeated-coordinate count remains unknown.'

    def show_qc(self):
        """Show whole-source coordinate facts and their explicit biological limits."""
        selected = self.pack.currentData()
        if selected is None:
            return
        try:
            report = Q.report(selected['path'])
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return
        if report is None:
            self.summary.setText('Coordinate QC is not available for this source. No quality score is assigned.')
            return
        if getattr(self, 'qc_dialog', None) is not None:
            self.qc_dialog.close()
            self.qc_dialog.deleteLater()
        self.qc_dialog = QtWidgets.QDialog(self)
        self.qc_dialog.setWindowTitle('Original sparse coordinate QC')
        layout = QtWidgets.QVBoxLayout(self.qc_dialog)
        label = QtWidgets.QLabel('Technical source-record facts. Repeated entries remain separate. '
                                'A global count is proven only by complete monotone ordering; '
                                'otherwise it remains unknown. Biological quality and comparability require separate review.')
        label.setWordWrap(True)
        layout.addWidget(label)
        tree = QtWidgets.QTreeWidget()
        tree.setHeaderLabels(['QC fact / provenance', 'Value'])
        captions = dict(original_records='Stored source records',
                        global_duplicate_count_proven='Global duplicate count verified',
                        global_repeated_records='Repeated stored records',
                        global_repeated_coordinate_groups='Coordinates with repeated records',
                        unique_stored_coordinates='Distinct stored coordinates',
                        duplicate_record_fraction='Repeated-record fraction',
                        adjacent_repeated_records='Adjacent repeated records',
                        adjacent_duplicate_runs='Adjacent runs of repeated records',
                        source_row_then_column_nondecreasing='Source sorted by row, then column',
                        source_column_then_row_nondecreasing='Source sorted by column, then row',
                        global_duplicate_scope='Scope of global count',
                        coordinate_aggregation_applied='Coordinate aggregation applied',
                        biological_quality_score='Biological quality score')
        for key, value in report['facts'].items():
            item = QtWidgets.QTreeWidgetItem([captions.get(key, key), 'Unknown / not assigned' if value is None else str(value)])
            item.setToolTip(1, str(value))
            tree.addTopLevelItem(item)
        for key in ('receipt_sha256', 'processor_sha256', 'source_input_sha256', 'output_sha256'):
            item = QtWidgets.QTreeWidgetItem([key, report[key]])
            item.setToolTip(1, report[key])
            tree.addTopLevelItem(item)
        tree.resizeColumnToContents(0)
        layout.addWidget(tree)
        self.qc_dialog.resize(800, 560)
        self.qc_dialog.show()

    def show_context(self):
        """Inspect original coordinates, numeric tokens and submitted axis metadata."""
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
