"""RNA/protein coexpression with explicit source, column and transformation choices."""
from __future__ import annotations

import json
import threading

import numpy as np
from PyQt6 import QtCore, QtWidgets

from . import coexpression as C
from . import optional_deposited_data as O
from .component_download_dialog import ComponentDownloadDialog
from .comparison_provenance import export


class _Worker(QtCore.QThread):
    ready = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, action, parent):
        super().__init__(parent)
        self.action = action
        self.cancelled = threading.Event()

    def run(self):
        try:
            result = self.action(self.cancelled.is_set)
            self.ready.emit(result)
        except (OSError, ValueError, KeyError, TypeError, InterruptedError) as error:
            self.failed.emit(str(error))


class CoexpressionPanel(QtWidgets.QWidget):
    """Explore associations across selected columns of a single measurement source."""
    gene_chosen = QtCore.pyqtSignal(str)

    def __init__(self, rna_root=None, deposited_root=None, parent=None):
        super().__init__(parent)
        self.rna_root = rna_root
        self.matrix = self.result = self.record = None
        self.worker = None
        self._generation = 0
        self.sources, errors = C.catalogue(rna_root, deposited_root)
        layout = QtWidgets.QVBoxLayout(self)
        intro = QtWidgets.QLabel(
            'Find features that vary together across samples within one RNA or protein source. '
            'Choose a source, load it, review sample columns, then enter one exact feature identifier. '
            'Protein groups and transcripts remain separate source rows. '
            'Sources can include host and contaminant features; taxon filtering is not inferred. '
            'RNA–protein agreement across genes is a separate analysis.')
        intro.setWordWrap(True)
        layout.addWidget(intro)
        self.level = QtWidgets.QComboBox()
        self.level.addItems(['RNA', 'Protein'])
        self.filter = QtWidgets.QLineEdit()
        self.filter.setPlaceholderText('Filter by study, analysis, quantity or reference')
        self.source = QtWidgets.QComboBox()
        self.source.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.source.setMinimumContentsLength(20)
        self.load_button = QtWidgets.QPushButton('Load selected source')
        for widget in (self.level, self.filter, self.source, self.load_button):
            layout.addWidget(widget)
        self.context = QtWidgets.QLabel()
        self.context.setWordWrap(True)
        layout.addWidget(self.context)
        layout.addWidget(QtWidgets.QLabel('Select columns after reviewing their roles; columns can represent channels or repeated measurements.'))
        self.samples = QtWidgets.QTableWidget(0, 3)
        self.samples.setHorizontalHeaderLabels(['Use', 'Original column', 'Submitted context (hover for full record)'])
        self.samples.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.samples.setMaximumHeight(190)
        layout.addWidget(self.samples)
        buttons = QtWidgets.QHBoxLayout()
        self.all_button = QtWidgets.QPushButton('Select all columns')
        self.none_button = QtWidgets.QPushButton('Clear selection')
        buttons.addWidget(self.all_button); buttons.addWidget(self.none_button)
        layout.addLayout(buttons)
        form = QtWidgets.QFormLayout()
        self.seed = QtWidgets.QLineEdit()
        self.seed.setPlaceholderText('Exact source feature label; case-sensitive')
        self.seed_row = QtWidgets.QLineEdit()
        self.seed_row.setPlaceholderText('Optional original row identity, required for repeated labels')
        self.method = QtWidgets.QComboBox()
        self.method.addItems(['Spearman', 'Pearson'])
        self.transform = QtWidgets.QComboBox()
        self.transform.addItem('Original source scale', 'none')
        self.transform.addItem('log2(x + 1); nonnegative values only', 'log2p1')
        self.minimum = QtWidgets.QSpinBox()
        self.minimum.setRange(4, 100000); self.minimum.setValue(6)
        form.addRow('Seed feature', self.seed)
        form.addRow('Original row identity', self.seed_row)
        form.addRow('Correlation', self.method)
        form.addRow('Transformation', self.transform)
        form.addRow('Minimum paired columns', self.minimum)
        layout.addLayout(form)
        actions = QtWidgets.QHBoxLayout()
        self.run_button = QtWidgets.QPushButton('Find coexpressed features')
        self.cancel_button = QtWidgets.QPushButton('Cancel'); self.cancel_button.setEnabled(False)
        self.export_button = QtWidgets.QPushButton('Export all results + treatment record')
        self.export_button.setEnabled(False)
        for widget in (self.run_button, self.cancel_button, self.export_button):
            actions.addWidget(widget)
        layout.addLayout(actions)
        self.summary = QtWidgets.QLabel('Load a source to begin.' + (f' {len(errors)} catalogue errors: ' + '; '.join(errors[:3]) if errors else ''))
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.table = QtWidgets.QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(['Original row', 'Feature', 'Correlation', 'Paired columns', 'Status'])
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setMinimumHeight(180)
        layout.addWidget(self.table)
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
        self.figure = Figure(figsize=(8, 3.5), constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumHeight(260)
        layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        layout.addWidget(self.canvas)
        self.level.currentTextChanged.connect(self._filter_sources)
        self.filter.textChanged.connect(self._filter_sources)
        self.source.currentIndexChanged.connect(self._source_changed)
        self.load_button.clicked.connect(self.load_selected)
        self.all_button.clicked.connect(lambda: self._select_all(True))
        self.none_button.clicked.connect(lambda: self._select_all(False))
        self.samples.itemChanged.connect(self._invalidate)
        for widget in (self.seed, self.seed_row):
            widget.textChanged.connect(self._invalidate)
        self.method.currentIndexChanged.connect(self._invalidate)
        self.transform.currentIndexChanged.connect(self._invalidate)
        self.minimum.valueChanged.connect(self._invalidate)
        self.run_button.clicked.connect(self.calculate)
        self.cancel_button.clicked.connect(self._cancel)
        self.export_button.clicked.connect(lambda: self.save())
        self.table.cellClicked.connect(self._show_pair)
        self._filter_sources()

    def _invalidate(self, *_):
        self._generation += 1
        if self.worker is not None:
            self.worker.cancelled.set()
        self.result = self.record = None
        self.table.setRowCount(0)
        self.export_button.setEnabled(False)
        self.figure.clear(); self.canvas.draw_idle()
        self.summary.setText('Selection changed. Calculate again before exporting.')

    def _source_changed(self, *_):
        self._invalidate()
        self.matrix = None
        self.samples.setRowCount(0)
        self.context.setText(self.source.currentData().get('title', '') if self.source.currentData() else 'No matching sources available.')

    def _filter_sources(self, *_):
        query = self.filter.text().casefold()
        self.source.blockSignals(True)
        self.source.clear()
        for source in self.sources:
            if source['level'] == self.level.currentText() and query in json.dumps(source).casefold():
                self.source.addItem(source['title'], source)
        self.source.blockSignals(False)
        self._source_changed()

    def _start(self, action, callback):
        if self.worker is not None:
            self.summary.setText('Wait for the current operation to finish or cancel it first.')
            return
        generation = self._generation
        worker = _Worker(action, self)
        self.worker = worker
        worker.ready.connect(lambda value: callback(value) if generation == self._generation else None)
        worker.failed.connect(lambda message: self.summary.setText(message) if generation == self._generation else None)
        worker.finished.connect(self._finished)
        self.load_button.setEnabled(False); self.run_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        worker.start()

    def _finished(self):
        worker = self.worker
        self.worker = None
        self.load_button.setEnabled(True); self.run_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        if worker is not None:
            worker.deleteLater()

    def _cancel(self):
        self._invalidate()
        self.summary.setText('Cancelled; verified source files are retained. Any active file verification will finish safely.')

    def load_selected(self):
        source = self.source.currentData()
        if not source or self.worker is not None:
            return
        self._source_changed()
        try:
            if source['kind'] == 'deposited':
                plan = O.plan(source['path'])
                generation = self._generation
                if plan['missing']:
                    dialog = ComponentDownloadDialog(plan['missing'], self)
                    if dialog.exec() != QtWidgets.QDialog.DialogCode.Accepted or generation != self._generation:
                        self.summary.setText('Source download cancelled or selection changed.')
                        return
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error)); return
        self.summary.setText('Loading and verifying the selected source…')
        self._start(lambda cancelled: C.load(source, self.rna_root), self._loaded)

    def _loaded(self, matrix):
        self.matrix = matrix
        self.samples.blockSignals(True)
        self.samples.setRowCount(len(matrix.samples))
        for i, row in enumerate(matrix.samples):
            check = QtWidgets.QTableWidgetItem()
            check.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            check.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.samples.setItem(i, 0, check)
            self.samples.setItem(i, 1, QtWidgets.QTableWidgetItem(str(row['column'])))
            text = json.dumps(row, ensure_ascii=False, default=str)
            item = QtWidgets.QTableWidgetItem(text[:250]); item.setToolTip(text)
            self.samples.setItem(i, 2, item)
        self.samples.blockSignals(False)
        self.samples.resizeColumnsToContents()
        self.context.setText(f"{matrix.record['quantity']} · {matrix.record['units']} · {len(matrix.values):,} source rows. "
                             + matrix.record['source_scope'])
        self.summary.setText('Review the original column roles and select samples. Conditions and replicate independence are not inferred.')

    def _select_all(self, selected):
        self.samples.blockSignals(True)
        for i in range(self.samples.rowCount()):
            self.samples.item(i, 0).setCheckState(QtCore.Qt.CheckState.Checked if selected else QtCore.Qt.CheckState.Unchecked)
        self.samples.blockSignals(False)
        self._invalidate()

    def selected_columns(self):
        return [self.matrix.samples[i]['column'] for i in range(self.samples.rowCount())
                if self.samples.item(i, 0).checkState() == QtCore.Qt.CheckState.Checked]

    def calculate(self):
        if self.matrix is None:
            self.summary.setText('Load and review a source first.'); return
        if self.worker is not None:
            return
        self._invalidate()
        matrix, seed = self.matrix, self.seed.text()
        columns = self.selected_columns()
        row = self.seed_row.text()
        if row:
            matches = [key for key in matrix.values.index if str(key) == row]
            if len(matches) != 1:
                self.summary.setText('Original row identity is absent or ambiguous.'); return
            row = matches[0]
        else:
            row = None
        method, transform, minimum = self.method.currentText().lower(), self.transform.currentData(), self.minimum.value()
        self.summary.setText('Calculating within-source descriptive correlations…')
        self._start(lambda cancelled: C.correlate(matrix, seed, columns, seed_row=row,
                    method=method, transform=transform, min_pairs=minimum, cancelled=cancelled), self._calculated)

    def _calculated(self, result):
        self.result, self.record = result
        self.export_button.setEnabled(True)
        shown = self.result.head(500)
        self.table.setRowCount(len(shown))
        for i, row in enumerate(shown.itertuples(index=False)):
            for j, value in enumerate((row.source_row, row.feature_label,
                    f'{row.correlation:.3f}' if np.isfinite(row.correlation) else '', row.paired_columns, row.status)):
                self.table.setItem(i, j, QtWidgets.QTableWidgetItem(str(value)))
        self.table.resizeColumnsToContents()
        count = int(self.result.correlation.notna().sum())
        small = ' Small sample set: interpret correlations cautiously.' if len(self.record['selected_columns']) < 15 else ''
        self.summary.setText(f'{count:,} descriptive correlations across {len(self.record["selected_columns"])} selected columns. '
            f'Showing {len(shown):,} of {len(self.result):,} rows; export retains all rows and exclusion reasons. '
            'No significance, biological quality, interaction or causal claim.' + small)
        if count:
            self._show_pair(0, 0)

    def _show_pair(self, index, _column):
        if self.result is None or index >= len(self.result):
            return
        row = self.result.iloc[index]
        columns = self.record['selected_columns']
        keys = [self.record['seed_row'], row.source_row]
        values = self.matrix.values.loc[keys, columns].to_numpy(dtype=float)
        if self.record['transform'] == 'log2p1':
            values = np.log2(values + 1)
        keep = np.isfinite(values).all(axis=0)
        self.figure.clear()
        scatter = self.figure.add_subplot(121)
        scatter.scatter(values[0, keep], values[1, keep], s=22)
        scatter.set(xlabel=self.record['seed_label'], ylabel=str(row.feature_label),
                    title=f'{self.record["method"]}: {row.correlation:.3f} · n={int(keep.sum())}')
        scatter.text(.02, .98, self.record['transform'] + ' · source units', transform=scatter.transAxes, va='top', fontsize=8)
        # A display-only within-row z score makes profile shapes visible, without
        # changing the source values used for correlations or exported results.
        top = self.result[self.result.correlation.notna()].head(12)
        heat_keys = [self.record['seed_row'], *top.source_row.tolist()]
        display = self.matrix.values.loc[heat_keys, columns].to_numpy(dtype=float, copy=True)
        if self.record['transform'] == 'log2p1':
            display = np.log2(display + 1)
        for i, vector in enumerate(display):
            finite = np.isfinite(vector)
            if finite.any():
                scaled = vector[finite] / max(float(np.max(np.abs(vector[finite]))), 1.)
                deviation = scaled.std()
                display[i, finite] = (scaled - scaled.mean()) / deviation if deviation else 0.
            display[i, ~finite] = np.nan
        heat = self.figure.add_subplot(122)
        image = heat.imshow(np.ma.masked_invalid(display), aspect='auto', cmap='coolwarm', vmin=-2, vmax=2)
        heat.set_yticks(range(len(heat_keys)), [str(self.matrix.labels.loc[key])[:25] for key in heat_keys], fontsize=7)
        heat.set_xticks(range(len(columns)), list(map(str, columns)), rotation=90, fontsize=6)
        heat.set_title('Top profiles · display-only row z score', fontsize=9)
        self.figure.colorbar(image, ax=heat, fraction=.04)
        self.canvas.draw_idle()

    def save(self, path=''):
        if self.result is None:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save coexpression', 'coexpression.csv', 'CSV (*.csv);;TSV (*.tsv)')
        if not path:
            return 0
        try:
            return export(self.result, path, dict(self.record,
                visualization='Scatter: selected transformed source values; heatmap: display-only within-row z score, top 12 valid correlations'))
        except (OSError, ValueError) as error:
            self.summary.setText(str(error)); return 0

    def closeEvent(self, event):
        if self.worker is not None:
            self.worker.cancelled.set()
            self.summary.setText('Finishing source verification or cancellation before closing.')
            event.ignore()
            return
        super().closeEvent(event)
