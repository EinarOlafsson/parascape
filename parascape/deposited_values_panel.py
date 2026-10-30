"""Browse independently reconciled deposited values with their original source context.

Other quantitative assays cannot safely inherit gene-expression units or sample
group assumptions. This view selects one source analysis and shows exact feature
labels, row numbers, literal tokens and submitted sample labels together. Repeated
feature labels retain separate rows. Exported values and whole-matrix numerical
QC carry provenance, while biological independence, compatible scales and feature
equivalence remain questions for source review rather than hidden UI defaults.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import uuid

from PyQt6 import QtCore, QtWidgets

from . import deposited_values as D
from . import protein_identifier_evidence as I
from . import source_flags as F


class DepositedValuesPanel(QtWidgets.QWidget):
    """Inspect one reconciled source representation and export values or sample QC."""

    def __init__(self, root=None, parent=None):
        super().__init__(parent)
        self.root, self.result, self.treatment_record = root, None, None
        self._catalogue, self._catalogue_errors = [], []
        layout = QtWidgets.QVBoxLayout(self)
        label = QtWidgets.QLabel('Original deposited measurements. Feature rows, value tokens and source units are retained. '
            'Choose one analysis; assays, repeated names and submitted columns remain separate. '
            'Cross-study compatibility and biological independence require review.')
        label.setWordWrap(True)
        layout.addWidget(label)
        self.source_filter = QtWidgets.QLineEdit()
        self.source_filter.setPlaceholderText('Filter by study, file, assay or quantity label')
        self.source_filter.setToolTip('Find source analyses by case-insensitive text in study, source analysis, assay or units. Original feature identifiers remain case-sensitive.')
        self.source_filter.textChanged.connect(self._filter_sources)
        layout.addWidget(self.source_filter)
        self.pack = QtWidgets.QComboBox()
        self.pack.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.pack.setMinimumContentsLength(12)
        self.pack.setToolTip('Choose one independently reconciled source analysis with its original quantity and units.')
        self.pack.currentIndexChanged.connect(self._clear)
        layout.addWidget(self.pack)
        self.catalogue_status = QtWidgets.QLabel()
        self.catalogue_status.setWordWrap(True)
        layout.addWidget(self.catalogue_status)
        self.refresh_button = QtWidgets.QPushButton('Refresh deposited measurements')
        self.refresh_button.setToolTip('Reload the validated measurement catalogue; incomplete preparations remain unavailable.')
        self.refresh_button.clicked.connect(self.refresh)
        layout.addWidget(self.refresh_button)
        self.identifier_kind = QtWidgets.QComboBox()
        self.identifier_kind.addItem('Exact original feature label', 'exact')
        for namespace, caption in I.NAMESPACES.items():
            self.identifier_kind.addItem(caption, namespace)
        self.identifier_kind.currentIndexChanged.connect(self._clear)
        self.identifier_kind.setToolTip('Keep original labels, submitted names and reference candidates distinct. Whole protein-group quantities are retained.')
        layout.addWidget(self.identifier_kind)
        self.features = QtWidgets.QPlainTextEdit()
        self.features.setMaximumHeight(80)
        self.features.setPlaceholderText('One exact original feature label per line')
        self.features.setToolTip('Match case-sensitive source labels, retaining every matching row. No aliases, feature collapse or gene mappings are inferred.')
        self.features.textChanged.connect(self._clear)
        layout.addWidget(self.features)
        self.exclude_source_marked = QtWidgets.QCheckBox('Exclude rows with source "+" flags')
        self.exclude_source_marked.setToolTip('Opt in to excluding protein rows with literal + in Reverse, Contaminant, Potential contaminant or Only identified by site. Missing or unfamiliar flags remain unreviewed and retained. Whole-matrix technical QC stays unfiltered.')
        self.exclude_source_marked.toggled.connect(self._clear)
        layout.addWidget(self.exclude_source_marked)
        for attribute, caption, tooltip, callback in [
            ('show_button', 'Show original measurements', 'Verify source hashes and show original values and literal source tokens.', self.show_values),
            ('measurement_button', 'Inspect selected measurement…', 'Inspect the selected original row, source column, condition and replicate declarations. Repeated labels remain separate; sample independence is not inferred.', self.show_measurement_context),
            ('export_button', 'Export displayed measurements…', 'Export retained original values and source tokens with a separate treatment record containing source hashes and any flag exclusions.', lambda: self.save()),
            ('treatment_button', 'Export view treatment record…', 'Record the source hashes, query, flag policy and excluded original rows, including when every match was excluded.', lambda: self.save_treatment()),
            ('context_button', 'Source columns and submitted factors…', 'Inspect original column labels and exact submitted filename associations without a feature query.', self.show_column_context),
            ('qc_button', 'Export column technical QC…', 'Check all source rows independently of the feature query. Biological quality scores remain unavailable.', lambda: self.save_qc())]:
            button = QtWidgets.QPushButton(caption)
            button.setToolTip(tooltip)
            button.clicked.connect(callback)
            setattr(self, attribute, button)
            layout.addWidget(button)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.table = QtWidgets.QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(['source row', 'source feature', 'source column', 'numeric value', 'original token', 'source units', 'original source flags'])
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.currentCellChanged.connect(self._measurement_selection_changed)
        self.table.cellDoubleClicked.connect(lambda *_: self.show_measurement_context())
        layout.addWidget(self.table)
        self.refresh()
        self.setMinimumWidth(120)

    def _clear(self):
        self.result = None
        self.treatment_record = None
        self._close_measurement_context()
        if hasattr(self, 'measurement_button'):
            self.measurement_button.setEnabled(False)
        if getattr(self, 'context_dialog', None) is not None:
            self.context_dialog.close()
            self.context_dialog = None
        if hasattr(self, 'identifier_kind'):
            source = self.pack.currentData() or {}
            protein = source.get('assay') == 'deposited protein-group quantities'
            for index in range(1, self.identifier_kind.count()):
                self.identifier_kind.model().item(index).setEnabled(protein)
            if not protein and self.identifier_kind.currentData() != 'exact':
                self.identifier_kind.setCurrentIndex(0)
        if hasattr(self, 'table'):
            self.table.setRowCount(0)
        if hasattr(self, 'export_button'):
            self.export_button.setEnabled(False)
        if hasattr(self, 'treatment_button'):
            self.treatment_button.setEnabled(False)
        if hasattr(self, 'summary'):
            self.summary.setText('Enter exact feature labels for the selected source analysis.')

    def refresh(self):
        """Load completed source analyses without mixing quantities or units."""
        self._catalogue, self._catalogue_errors = D.catalogue(self.root)
        self._filter_sources()

    def _filter_sources(self):
        """Filter source choices while retaining a visible selection and invalidating all earlier query results."""
        previous = (self.pack.currentData() or {}).get('path')
        query = self.source_filter.text().casefold().strip()
        packs = [pack for pack in self._catalogue if query in ' '.join(
            str(pack[key]) for key in ('study', 'analysis', 'assay', 'units')).casefold()]
        blocked = self.pack.blockSignals(True)
        try:
            self.pack.clear()
            for pack in packs:
                self.pack.addItem(f"{pack['study']} · {pack['analysis']} · {pack['assay']} · {pack['samples']} source columns · {pack['units']}", pack)
            for index, pack in enumerate(packs):
                if pack['path'] == previous:
                    self.pack.setCurrentIndex(index)
                    break
        finally:
            self.pack.blockSignals(blocked)
        self._clear()
        self.catalogue_status.setText(f'{len(packs)} of {len(self._catalogue)} reconciled source analyses shown. ' +
            (f'{len(self._catalogue_errors)} catalogue entries require review.' if self._catalogue_errors else
             'Quantity families and source units remain separate.'))
        self.catalogue_status.setToolTip('\n'.join(self._catalogue_errors))
        if not packs:
            self.summary.setText('No source analyses match this filter.' if self._catalogue else
                                 'No reconciled deposited measurement data is installed.')
        for button in (self.show_button, self.qc_button, self.context_button):
            button.setEnabled(bool(packs))

    def _close_measurement_context(self):
        if getattr(self, 'measurement_dialog', None) is not None:
            self.measurement_dialog.close()
            self.measurement_dialog.deleteLater()
            self.measurement_dialog = None
            self.measurement_tree = None

    def _selected_measurement(self):
        """Use the item's result position so sorting cannot substitute another repeated label."""
        if self.result is None:
            return None
        item = self.table.item(self.table.currentRow(), 0)
        position = item.data(QtCore.Qt.ItemDataRole.UserRole) if item is not None else None
        if not isinstance(position, int) or isinstance(position, bool) or not 0 <= position < len(self.result):
            return None
        return self.result.iloc[position].to_dict()

    def _measurement_selection_changed(self, *_):
        self._close_measurement_context()
        self.measurement_button.setEnabled(self._selected_measurement() is not None)

    def show_measurement_context(self):
        """Inspect one original row/column measurement without inferring roles from its label."""
        self._close_measurement_context()
        row = self._selected_measurement()
        if row is None:
            return
        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle('Original measurement · ' + str(row['feature_label']) + ' · row ' + str(row['source_row']))
        layout = QtWidgets.QVBoxLayout(dialog)
        label = QtWidgets.QLabel('Original annotations belong to this source row and quantity column. '
            'Condition, filename and replicate declarations are preserved as submitted; they do not establish '
            'biological independence, control pairing, gene equivalence or compatible scales. '
            'Numeric missingness and the original token remain separate.')
        label.setWordWrap(True)
        layout.addWidget(label)
        tree = QtWidgets.QTreeWidget()
        tree.setHeaderLabels(['Original field or context', 'Value or declaration'])
        tree.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        layout.addWidget(tree)

        def add(parent, name, value):
            text = '' if isinstance(value, (dict, list)) else value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)
            item = QtWidgets.QTreeWidgetItem(parent, [str(name), text])
            item.setToolTip(0, str(name))
            item.setToolTip(1, text)
            if isinstance(value, dict):
                for key, child in value.items():
                    add(item, key, child)
            elif isinstance(value, list):
                for position, child in enumerate(value, 1):
                    add(item, 'Occurrence ' + str(position), child)
            return item

        add(tree, 'Selected original measurement', {caption: row[field] for caption, field in (
            ('Source row ordinal', 'source_row'), ('Original feature label', 'feature_label'),
            ('Original quantity column', 'source_sample'), ('Numeric value', 'value'),
            ('Original source token', 'source_token'), ('Source units', 'units'))}).setExpanded(True)
        for caption, field in (('Original row annotations', 'source_row_metadata'),
                               ('Original column metadata', 'sample_metadata'),
                               ('Qualified identifier evidence', 'identifier_evidence')):
            encoded = row.get(field)
            add(tree, caption, json.loads(encoded) if isinstance(encoded, str) and encoded else
                'No declaration supplied for this measurement').setExpanded(True)
        add(tree, 'Source flag view', dict(original_source_flags=json.loads(row['source_flags']),
            interpretation=row['source_flag_status'], view_policy=json.loads(row['source_flag_view_policy'])))
        add(tree, 'Source and comparison scope', {field: row[field] for field in (
            'study', 'source_analysis', 'assay', 'quantity', 'units', 'source_url', 'input_sha256',
            'manifest_sha256', 'matrix_sha256', 'literal_tokens_sha256', 'comparison_scope')})
        tree.setColumnWidth(0, 370)
        dialog.resize(1100, 700)
        self.measurement_dialog, self.measurement_tree = dialog, tree
        dialog.show()


    def show_column_context(self):
        """Inspect source columns and every exact submitted factor occurrence without inferring sample roles."""
        info = self.pack.currentData()
        if getattr(self, 'context_dialog', None) is not None:
            self.context_dialog.close()
            self.context_dialog = None
        if not info:
            return
        try:
            contexts = D.column_context(info['path'])
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return
        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle(info['study'] + ' · ' + info['analysis'] + ' · Source columns')
        layout = QtWidgets.QVBoxLayout(dialog)
        label = QtWidgets.QLabel('Original source-column metadata and submitted filename associations. '
            'All repeated rows and duplicate fields remain visible. Biological independence and control pairing are not inferred. '
            'Experiment-label candidates require review of input-table linkage.')
        label.setWordWrap(True)
        layout.addWidget(label)
        tree = QtWidgets.QTreeWidget()
        tree.setHeaderLabels(['Source column / field', 'Original value or association status'])
        layout.addWidget(tree)
        for context in contexts:
            association = context['source_design_metadata']
            status = association.get('primary_file_status', association.get('association_status', ''))
            caption = {'exact_unique_primary_file': 'Exact primary file match',
                       'ambiguous_exact_primary_file': 'Multiple primary file matches',
                       'unresolved_exact_primary_file': 'No exact primary file match',
                       'filename_role_unreviewed': 'Filename role unreviewed'}.get(status, status)
            if status == 'filename_role_unreviewed' and association.get('submitted_experiment_rows'):
                caption = 'Exact experiment-label candidates'
            parent = QtWidgets.QTreeWidgetItem(tree, [context['source_sample'], caption])
            if association.get('source_filename_literal') is not None:
                QtWidgets.QTreeWidgetItem(parent, ['Submitted filename', association['source_filename_literal']])
                QtWidgets.QTreeWidgetItem(parent, ['Primary file identifiers', '; '.join(association['exact_primary_file_ids'])])
                QtWidgets.QTreeWidgetItem(parent, ['Submitted design row occurrences', str(len(association['submitted_factor_rows']))])
            if association.get('source_experiment_label_literal') is not None:
                QtWidgets.QTreeWidgetItem(parent, ['Submitted experiment label', association['source_experiment_label_literal']])
                QtWidgets.QTreeWidgetItem(parent, ['Candidate design row occurrences', str(len(association.get('submitted_experiment_rows', [])))])
            original = QtWidgets.QTreeWidgetItem(parent, ['Original column metadata', 'Unchanged submitted record'])
            for key, value in context.items():
                if key != 'source_design_metadata':
                    QtWidgets.QTreeWidgetItem(original, [key, value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)])
            audit = QtWidgets.QTreeWidgetItem(parent, ['Association provenance', 'Exact filename scope and source receipts'])
            for key, value in association.items():
                if key not in ('submitted_factor_rows', 'submitted_experiment_rows'):
                    QtWidgets.QTreeWidgetItem(audit, [key, value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)])
            occurrences = [(row, 'Submitted row ', 'matching_file_column_indices') for row in association.get('submitted_factor_rows', [])]
            occurrences += [(row, 'Experiment row candidate ', 'matching_experiment_column_indices') for row in association.get('submitted_experiment_rows', [])]
            for row, title, match_key in occurrences:
                submitted = QtWidgets.QTreeWidgetItem(parent, [title + str(row['source_row']), row['source_url']])
                for position, (header, cell) in enumerate(zip(row['headers'], row['cells']), 1):
                    QtWidgets.QTreeWidgetItem(submitted, [str(position) + ': ' + header, cell])
                for key in ('source_sha256', 'source_row_id', 'metadata_origin', 'biological_truth_status', match_key):
                    QtWidgets.QTreeWidgetItem(submitted, [key, json.dumps(row[key], ensure_ascii=False)])
        iterator = QtWidgets.QTreeWidgetItemIterator(tree)
        while iterator.value() is not None:
            item = iterator.value()
            for column in range(2):
                item.setToolTip(column, item.text(column))
            iterator += 1
        tree.setColumnWidth(0, 400)
        dialog.resize(1000, 650)
        self.context_dialog, self.context_tree = dialog, tree
        dialog.show()

    def show_values(self):
        """Show exact source rows and original tokens for the selected feature labels."""
        info = self.pack.currentData()
        if not info:
            return
        self._clear()
        try:
            namespace = self.identifier_kind.currentData()
            frame = (D.profiles(info['path'], self.features.toPlainText()) if namespace == 'exact' else
                     I.profiles(info['path'], self.features.toPlainText(), namespace=namespace))
            original = frame
            manifest = D.manifest(info['path'], verify_outputs=False)
            manifest_sha256 = hashlib.sha256((Path(info['path']) / 'manifest.json').read_bytes()).hexdigest()
            if not frame.empty and not frame.manifest_sha256.eq(manifest_sha256).all():
                raise ValueError('Source manifest changed after the original query')
            frame, treatment = F.apply(frame, exclude_marked=self.exclude_source_marked.isChecked())
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.summary.setText(str(error))
            return
        self.result = frame
        identifier_contexts = ([dict(source_row=int(row.source_row), evidence=json.loads(row.identifier_evidence))
            for row in original.drop_duplicates('source_row').itertuples()]
            if not original.empty and 'identifier_evidence' in original else [])
        self.treatment_record = dict(schema_version=1, recorded_at=datetime.now(timezone.utc).isoformat(),
            manifest_sha256=manifest_sha256, matrix_sha256=manifest['outputs']['values.parquet'],
            literal_tokens_sha256=manifest['outputs']['source_tokens.parquet'],
            source_url=manifest['source_url'], input_sha256=manifest['input_sha256'],
            source_study=manifest['study'], source_analysis=manifest['source_analysis'],
            quantity=manifest['quantity'], units=manifest['units'], identifier_namespace=namespace,
            exact_query_literals=self.features.toPlainText().splitlines(), source_selection_filter=self.source_filter.text(),
            original_identifier_evidence=identifier_contexts,
            source_flag_view=treatment, original_files_modified=False, numerical_values_transformed=False,
            biological_quality_score_assigned=False)
        self.export_button.setEnabled(not frame.empty)
        self.treatment_button.setEnabled(True)
        self.table.setRowCount(len(frame))
        sorting = self.table.isSortingEnabled()
        self.table.setSortingEnabled(False)
        for i, row in enumerate(frame.to_dict('records')):
            for j, field in enumerate(('source_row', 'feature_label', 'source_sample', 'value', 'source_token', 'units', 'source_flags')):
                item = QtWidgets.QTableWidgetItem('' if row[field] is None else str(row[field]))
                if j == 0:
                    item.setData(QtCore.Qt.ItemDataRole.UserRole, i)
                item.setToolTip(row['sample_metadata'] + '\n' + row.get('source_row_metadata', '') + '\n' + row.get('identifier_evidence', '') + '\nSource flags: ' + row['source_flag_status'] + '\nSource: ' + row['source_url'] + '\n' + row['comparison_scope'])
                self.table.setItem(i, j, item)
        self.table.setSortingEnabled(sorting)
        self._measurement_selection_changed()
        if namespace == 'exact':
            matched = set(original.feature_label)
        else:
            matched = {row['query_literal'] for text in original.identifier_evidence for row in json.loads(text)}
        missing = set(self.features.toPlainText().splitlines()) - matched
        self.summary.setText(f'{len(frame)} original row/column measurements. ' +
            ('Labels not matched: ' + ', '.join(sorted(missing)) if missing else 'Repeated labels retain separate source rows.'))
        if namespace != 'exact':
            self.summary.setText(self.summary.text() + ' Original group quantities retained; unique-gene abundance and current gene equivalence remain unreviewed.')
        self.summary.setText(self.summary.text() + f" Source-flag view excluded {treatment['excluded_source_rows']} of {treatment['input_source_rows']} queried original rows ({treatment['excluded_measurements']} measurements). Missing or unfamiliar flags are retained; this is not a biological quality score.")
        if treatment['input_source_rows'] and not treatment['retained_source_rows']:
            self.summary.setText('All matching source rows were excluded by the selected source-flag rule. ' + self.summary.text())
        self.table.resizeColumnsToContents()
        self.table.setColumnWidth(6, 240)

    def save(self, path='') -> int:
        """Export original row/sample measurements and provenance without collapsing names."""
        if self.result is None or self.result.empty:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export deposited measurements', 'deposited_measurements.csv', 'CSV (*.csv)')
        if not path:
            return 0
        self.result.to_csv(path, index=False)
        csv_path = Path(path)
        record = dict(self.treatment_record, exported_csv=str(csv_path.resolve()),
                      exported_csv_sha256=hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                      exported_measurements=len(self.result))
        sidecar = csv_path.with_name(csv_path.name + '.treatment-' + uuid.uuid4().hex[:12] + '.json')
        sidecar.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
        self.summary.setText(f'Exported {len(self.result)} displayed original measurements. Treatment record: {sidecar.name}')
        return len(self.result)

    def save_treatment(self, path='') -> int:
        """Export the exact view operation even when source flags exclude every matched measurement."""
        if self.treatment_record is None:
            return 0
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export view treatment record', 'measurement_view_treatment.json', 'JSON (*.json)')
        if not path:
            return 0
        Path(path).write_text(json.dumps(self.treatment_record, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
        return 1

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
        self.summary.setText(f'Exported unfiltered whole-source technical QC for {len(frame)} source columns. The view flag filter is not applied; biological quality scores remain unavailable.')
        return len(frame)
