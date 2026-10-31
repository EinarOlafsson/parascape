"""Browse complete original workbook rows with explicit verification and scope."""
from __future__ import annotations

import json
from pathlib import Path
import threading

from PyQt6 import QtCore, QtWidgets

from . import workbook_values as W
from . import optional_workbook_data as O
from . import workbook_qc as Q
from .component_download_dialog import ComponentDownloadDialog


class _Signals(QtCore.QObject):
    finished = QtCore.pyqtSignal(int, object, str)


class _Task(QtCore.QRunnable):
    def __init__(self, generation, function):
        super().__init__()
        self.generation, self.function = generation, function
        self.signals = _Signals()
        self.cancel = threading.Event()

    def run(self):
        result, error = None, ''
        try:
            result = self.function(self.cancel)
        except Exception as exception:
            error = str(exception)
        self.signals.finished.emit(self.generation, result, error)


class WorkbookValuesPanel(QtWidgets.QWidget):
    """Inspect source rows without pooling measurement, annotation and test columns."""

    COLUMNS = ('row', 'column', 'source_address', 'storage_type', 'source_value', 'formula_xml')

    def __init__(self, root=None, parent=None):
        super().__init__(parent)
        self.root, self.result = root, None
        self._worker = None
        self._generation = 0
        self._entries = []
        self.context_dialog = None
        self.context_record = None
        self.qc_dialog = None
        layout = QtWidgets.QVBoxLayout(self)
        description = QtWidgets.QLabel(
            'Original workbook rows. Labels, measurements, annotations and statistical results retain '
            'their submitted positions. Quantity units and sample roles require review before comparison.')
        description.setWordWrap(True)
        layout.addWidget(description)
        self.source_filter = QtWidgets.QLineEdit()
        self.source_filter.setPlaceholderText('Filter by study or workbook filename')
        self.source_filter.textChanged.connect(self._filter)
        layout.addWidget(self.source_filter)
        self.require_verified = QtWidgets.QCheckBox('Require complete source-cell verification')
        self.require_verified.setToolTip('Show sources with independent whole-source cell reconciliation. This checks preservation; biological quality and comparable measurement units remain unreviewed.')
        self.require_verified.toggled.connect(self._filter)
        layout.addWidget(self.require_verified)
        self.require_agreement = QtWidgets.QCheckBox('Require decoded-value reader agreement')
        self.require_agreement.setToolTip('Exclude unknown verification and reported reader discrepancies. Sources with primary-byte verification and a reader disagreement remain available when unchecked.')
        self.require_agreement.toggled.connect(self._filter)
        layout.addWidget(self.require_agreement)
        self.pack = QtWidgets.QComboBox()
        self.pack.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.pack.setMinimumContentsLength(12)
        self.pack.currentIndexChanged.connect(self._source_changed)
        layout.addWidget(self.pack)
        self.sheet = QtWidgets.QComboBox()
        self.sheet.currentIndexChanged.connect(self._clear)
        layout.addWidget(self.sheet)
        self.refresh_button = QtWidgets.QPushButton('Refresh workbooks')
        self.refresh_button.clicked.connect(self.refresh)
        layout.addWidget(self.refresh_button)
        self.catalogue_status = QtWidgets.QLabel()
        self.catalogue_status.setWordWrap(True)
        layout.addWidget(self.catalogue_status)
        self.labels = QtWidgets.QPlainTextEdit()
        self.labels.setMaximumHeight(75)
        self.labels.setPlaceholderText('Exact text labels, one per line. Leave empty to preview the first 20 rows with stored cells.')
        self.labels.setToolTip('Matching preserves case, spaces and duplicate source rows. Numeric cells are not interpreted as gene labels.')
        self.labels.textChanged.connect(self._clear)
        layout.addWidget(self.labels)
        column_row = QtWidgets.QHBoxLayout()
        column_row.addWidget(QtWidgets.QLabel('Label column:'))
        self.label_column = QtWidgets.QSpinBox()
        self.label_column.setRange(0, 16384)
        self.label_column.setSpecialValueText('Any text column')
        self.label_column.setToolTip('Use the original column number shown in the row preview to restrict label matching. Zero searches all text columns.')
        self.label_column.valueChanged.connect(self._clear)
        column_row.addWidget(self.label_column)
        layout.addLayout(column_row)
        self.show_button = QtWidgets.QPushButton('Show source rows')
        self.show_button.clicked.connect(self.show_values)
        layout.addWidget(self.show_button)
        self.cancel_button = QtWidgets.QPushButton('Cancel')
        self.cancel_button.clicked.connect(self._clear)
        layout.addWidget(self.cancel_button)
        self.download_button = QtWidgets.QPushButton('Download / verify selected workbook…')
        self.download_button.clicked.connect(self.download_selected)
        layout.addWidget(self.download_button)
        self.export_button = QtWidgets.QPushButton('Export all selected cells…')
        self.export_button.clicked.connect(lambda: self.save())
        layout.addWidget(self.export_button)
        self.context_button = QtWidgets.QPushButton('Inspect source and treatment…')
        self.context_button.clicked.connect(self.show_context)
        layout.addWidget(self.context_button)
        self.qc_button = QtWidgets.QPushButton('Inspect source QC…')
        self.qc_button.setToolTip('Read source-bound QC metadata before downloading the workbook. Original errors, formula bytes and reader differences remain separate from biological quality.')
        self.qc_button.clicked.connect(self.show_qc)
        layout.addWidget(self.qc_button)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.table = QtWidgets.QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(['source row', 'source column', 'cell address', 'cell kind',
                                              'original value', 'original formula'])
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.cellDoubleClicked.connect(lambda *_: self.show_context())
        layout.addWidget(self.table)
        self.setMinimumWidth(120)
        self.refresh()

    def _buttons(self):
        self.show_button.setEnabled(self.sheet.currentData() is not None and self._worker is None)
        self.cancel_button.setEnabled(self._worker is not None)
        self.export_button.setEnabled(self.result is not None and self._worker is None)
        self.context_button.setEnabled(self.result is not None)
        self.qc_button.setEnabled(bool(self.pack.currentData()))
        self.download_button.setEnabled(bool(self.pack.currentData()) and self._worker is None)

    def _clear(self, *_):
        self._generation += 1
        if self._worker is not None:
            self._worker.cancel.set()
        self.result = None
        self.table.setRowCount(0)
        self.summary.clear()
        self.context_record = None
        if self.context_dialog is not None:
            self.context_dialog.close()
            self.context_dialog = None
        if self.qc_dialog is not None:
            self.qc_dialog.close()
            self.qc_dialog = None
        self._buttons()

    def refresh(self):
        """Include newly completed workbooks and independent verification records."""
        self._clear()
        self._entries, errors = W.catalogue(self.root)
        self.catalogue_status.setText(f'{len(self._entries)} complete source workbooks; {len(errors)} source issues. '
                                      'Storage verification does not establish biological quality or comparable units.')
        self.catalogue_status.setToolTip('\n'.join(errors))
        self._filter()

    def _filter(self, *_):
        self._clear()
        needle = self.source_filter.text().casefold()
        self.pack.blockSignals(True)
        self.pack.clear()
        for row in self._entries:
            if self.require_verified.isChecked() and not row['independent_whole_source_reconciliation']:
                continue
            if self.require_agreement.isChecked() and row['decoded_value_checks_agree'] is not True:
                continue
            context = row['source_record']
            label = f"{context.get('studies', '')} · {Path(row['source_url']).name}"
            if needle in (label + json.dumps(context)).casefold():
                self.pack.addItem(label, row)
        self.pack.blockSignals(False)
        self._source_changed()

    def _source_changed(self, *_):
        self._clear()
        self.sheet.blockSignals(True)
        self.sheet.clear()
        selected = self.pack.currentData()
        try:
            if selected:
                info = W.manifest(selected['path'])
                for sheet in info['extraction']['sheets']:
                    if sheet.get('extraction_status') == 'non_worksheet_member_retained_in_original_only':
                        continue
                    attributes = sheet.get('source_attributes', {})
                    name = attributes.get('name', sheet.get('name', str(sheet['sheet_index'])))
                    state = attributes.get('state', 'hidden' if sheet.get('visibility') else 'visible')
                    self.sheet.addItem(f"{name} · {sheet['cells']:,} cells · {state}", sheet['sheet_index'])
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
        self.sheet.blockSignals(False)
        self._buttons()

    def _start(self, function, mode):
        self._worker = _Task(self._generation, function)
        self._worker.signals.finished.connect(lambda generation, result, error: self._finished(generation, result, error, mode))
        QtCore.QThreadPool.globalInstance().start(self._worker)
        self._buttons()

    def show_values(self):
        """Read original source rows, offering missing hosted files first."""
        if self._worker is not None or self.sheet.currentData() is None:
            return
        self._clear()
        folder = self.pack.currentData()['path']
        try:
            dependencies = O.plan(folder)
            if dependencies['missing']:
                dialog = ComponentDownloadDialog(dependencies['missing'], self)
                if dialog.exec() != QtWidgets.QDialog.DialogCode.Accepted:
                    self.summary.setText('Workbook download cancelled; no query or export is active.')
                    return
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return
        sheet, column = self.sheet.currentData(), self.label_column.value()
        labels = self.labels.toPlainText().splitlines()
        self.summary.setText('Verifying source files and reading original rows…')
        self._start(lambda cancel: W.prepare_query(folder, sheet, labels, label_column=column, cancel=cancel), 'query')

    def _finished(self, generation, result, error, mode):
        self._worker = None
        if generation != self._generation:
            self._buttons()
            return
        if error:
            if mode == 'query':
                self.result = None
            self.summary.setText(error)
        elif mode == 'export':
            self.summary.setText(f'Exported all {result:,} selected physical cells with a checksum and treatment sidecar.')
        else:
            self.result = result
            self.table.setRowCount(len(result.preview))
            supplementary = result.record['validation'].get('supplementary_biff_qc')
            for r, row in enumerate(result.preview):
                original_formula = Q.formula_for(supplementary, result.sheet, row)
                for c, key in enumerate(self.COLUMNS):
                    if key == 'source_value':
                        value = row['decoded_text'] if row['decoded_text'] is not None else row['source_value_token']
                    elif key == 'storage_type':
                        value = {'s':'Text', 'inlineStr':'Text', 'str':'Text', 'b':'Boolean', 'e':'Error',
                                 'd':'Date text', 'BIFF:text':'Text', 'BIFF:number':'Number',
                                 'BIFF:date_serial':'Date serial', 'BIFF:boolean':'Boolean',
                                 'BIFF:error_code':'Error code', 'BIFF:blank':'Blank'}.get(row[key], row[key])
                        if row[key] == 'n':
                            value = 'Formula cache' if row['formula_xml'] else ('Number' if row['source_value_token'] is not None else 'Blank')
                        if original_formula:
                            value = 'Formula cache / ' + original_formula['cached_result_kind']
                    elif key == 'formula_xml' and row[key] is not None:
                        import xml.etree.ElementTree as ET
                        value = ET.fromstring(row[key]).text or 'Shared formula; inspect source details'
                    elif key == 'formula_xml' and original_formula:
                        value = 'Original BIFF tokens; inspect source details'
                    else:
                        value = row[key]
                    item = QtWidgets.QTableWidgetItem('' if value is None else str(value))
                    hint = ('Missing source field' if value is None else 'Original source field; no numeric or gene-ID conversion')
                    if key == 'source_value':
                        hint = 'Preserved source value. Inspect source and treatment for storage tokens, legacy binary64 values and validation limits.'
                    item.setToolTip(hint)
                    self.table.setItem(r, c, item)
            proof = result.record['validation']
            grade = ('Original source cells independently reconciled' if proof['independent_whole_source_reconciliation']
                     else 'Independent source verification incomplete')
            if proof['reported_reader_disagreements']:
                grade += f"; {proof['reported_reader_disagreements']:,} reader text disagreements retained"
            self.summary.setText(f"{result.record['selected_physical_cells']:,} stored cells in selected row coordinates; "
                                 f"showing {len(result.preview):,}. {grade}. Units, sample roles and biological quality remain unreviewed. "
                                 f"Unmatched labels: {len(result.record['unmatched_labels'])}.")
        self._buttons()

    def show_context(self):
        """Expose exact source cells, original treatment and verification limitations."""
        if self.result is None:
            return
        selected = self.table.currentRow()
        self.context_record = dict(self.result.record, selected_source_cell=(
            self.result.preview[selected] if selected >= 0 else None))
        self.context_record['selected_source_formula'] = Q.formula_for(
            self.result.record['validation'].get('supplementary_biff_qc'), self.result.sheet,
            self.context_record['selected_source_cell'])
        if self.context_dialog is not None:
            self.context_dialog.close()
        self.context_dialog = QtWidgets.QDialog(self)
        self.context_dialog.setWindowTitle('Original workbook source and treatment')
        layout = QtWidgets.QVBoxLayout(self.context_dialog)
        text = QtWidgets.QPlainTextEdit(json.dumps(self.context_record, indent=2, ensure_ascii=False))
        text.setReadOnly(True)
        layout.addWidget(text)
        self.context_dialog.resize(800, 550)
        self.context_dialog.show()

    def show_qc(self):
        """Show bound source-preservation facts without reading optional payloads."""
        if not self.pack.currentData():
            return
        try:
            record = W.validation(self.pack.currentData()['path'])
            if self.qc_dialog is not None:
                self.qc_dialog.close()
            self.qc_dialog = QtWidgets.QDialog(self)
            self.qc_dialog.setWindowTitle('Original workbook source QC')
            layout = QtWidgets.QVBoxLayout(self.qc_dialog)
            label = QtWidgets.QLabel('Source preservation and reader agreement are separate checks. '
                                    'Biological quality, units and sample roles remain unreviewed.')
            label.setWordWrap(True)
            layout.addWidget(label)
            text = QtWidgets.QPlainTextEdit(json.dumps(record, indent=2, ensure_ascii=False))
            text.setReadOnly(True)
            layout.addWidget(text)
            self.qc_dialog.resize(800, 550)
            self.qc_dialog.show()
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))

    def download_selected(self):
        """Download or reverify only this hosted workbook's pinned source artifacts."""
        if self._worker is not None or not self.pack.currentData():
            return
        self._clear()
        try:
            dependencies = O.plan(self.pack.currentData()['path'])
            if not dependencies['components']:
                self.summary.setText('This workbook uses local source files. Show source rows to verify their checksums.')
                return
            dialog = ComponentDownloadDialog(dependencies['components'], self)
            if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
                self.summary.setText('Workbook downloaded and verified. Show source rows to browse it offline.')
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))

    def save(self, destination=None):
        """Export complete records; interactive exports run in the background."""
        if self.result is None or self._worker is not None:
            return 0
        result = self.result
        if destination is not None:
            return W.export_csv(result, destination)
        destination, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export original workbook cells', '', 'CSV files (*.csv)')
        if not destination:
            return 0
        self.summary.setText('Exporting every selected source cell…')
        self._start(lambda cancel: W.export_csv(result, destination, cancel=cancel), 'export')
        return 0
