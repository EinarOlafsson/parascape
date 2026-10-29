"""Browse catalogue coverage and original sample labels without overstating processing.

A downloaded archive is not an admitted expression matrix. This view keeps the
whole catalogue visible, distinguishes source acquisition from local measurements,
and exposes exact submitted identities and factors. Shared samples remain linked;
their presence in multiple studies does not establish independent replication.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
import re
from pathlib import Path

from PyQt6 import QtWidgets

from . import paths


def audit_root() -> Path:
    """Locate the separately installed audit data or an explicit local override."""
    return Path(os.environ.get('PARASCAPE_AUDIT', str(Path(paths.data_dir()) / 'processing_audit')))


def coverage(root: Path) -> tuple[list[dict], dict]:
    """Read a coverage table only when its bytes match the accompanying ledger."""
    data = (root / 'dataset_coverage.tsv').read_bytes()
    info = json.loads((root / 'dataset_coverage_status.json').read_text())
    if hashlib.sha256(data).hexdigest() != info['coverage_sha256']:
        raise ValueError('Coverage ledger is updating or its checksum differs; refresh to retry')
    rows = list(csv.DictReader(io.StringIO(data.decode()), delimiter='\t'))
    if len(rows) != info['catalogue_records'] or len({(r['source'], r['accession']) for r in rows}) != len(rows):
        raise ValueError('Catalogue coverage identities differ from the ledger')
    return rows, info


def source_samples(root: Path, accession: str) -> dict:
    """Load original submitted labels only for a reconciled primary GEO source."""
    if not accession.startswith('GSE') or not accession[3:].isdigit():
        raise ValueError('No primary GEO sample metadata for this catalogue source')
    folder = root / 'all_geo_samples'
    entries = json.loads((folder / 'coverage.json').read_text())
    matches = [r for r in entries if r['accession'] == accession]
    if len(matches) != 1 or matches[0]['status'] != 'validated_sample_metadata':
        raise ValueError('Primary sample metadata has not been validated for this source')
    data = (folder / accession / 'geo-samples.json').read_bytes()
    if hashlib.sha256(data).hexdigest() != matches[0]['sample_metadata_sha256']:
        raise ValueError('Primary sample metadata checksum differs')
    report = json.loads(data)
    if report['accession'] != accession or report['source_samples'] != len(report['samples']):
        raise ValueError('Primary sample identity or count differs')
    return report


def shared_samples(root: Path) -> dict:
    """Read audited cross-study links for identical submitted GEO sample identifiers."""
    folder = root / 'all_geo_samples'
    info = json.loads((folder / 'sample_identity_manifest.json').read_text())
    data = (folder / 'sample_identity_links.tsv').read_bytes()
    if hashlib.sha256(data).hexdigest() != info['outputs']['sample_identity_links.tsv']:
        raise ValueError('Shared-sample identity index checksum differs')
    rows = csv.DictReader(io.StringIO(data.decode()), delimiter='\t')
    return {row['geo_sample']: row['same_geo_sample_in_studies'].split(';') for row in rows}


def workbench_samples(root: Path, accession: str) -> dict:
    """Read independently reconciled submitted Workbench factors and analysis metadata."""
    if not re.fullmatch(r'ST\d{6}', accession):
        raise ValueError('Invalid Workbench study identity')
    collection = root / 'metabolomics_workbench'
    pointer = json.loads((collection / 'collection.json').read_text())
    directory = pointer['directory']
    if pointer.get('schema_version') != 1 or not re.fullmatch(r'[0-9a-f]{64}', directory):
        raise ValueError('Unsupported Workbench collection identity')
    folder = collection / directory
    if folder.resolve().parent != collection.resolve():
        raise ValueError('Workbench collection lies outside its installed root')
    receipt = (folder / 'validation.json').read_bytes()
    if hashlib.sha256(receipt).hexdigest() != pointer['validation_sha256']:
        raise ValueError('Workbench reconciliation receipt checksum differs')
    outputs = json.loads(receipt)['outputs']

    def read(name):
        data = (folder / name).read_bytes()
        if hashlib.sha256(data).hexdigest() != outputs[name]:
            raise ValueError('Workbench metadata checksum differs: ' + name)
        return json.loads(data)

    entries = [row for row in read('coverage.json') if row['accession'] == accession]
    if len(entries) != 1 or entries[0]['status'] != 'independently_reconciled_primary_metadata':
        raise ValueError('Primary Workbench metadata has not been reconciled')
    report = read(accession + '/sample_design.json')
    analyses = read(accession + '/analyses.json')
    if (report['study'] != accession or analyses['study'] != accession
            or len(report['submitted_rows']) != entries[0]['submitted_factor_records']
            or len(analyses['analyses']) != entries[0]['analyses']):
        raise ValueError('Workbench study identity or metadata coverage differs')
    return dict(**report, analysis_inventory=analyses['analyses'], source_coverage=entries[0])


def pride_samples(root: Path, accession: str) -> dict:
    """Read original PRIDE project metadata and independently reconciled submitted SDRF rows."""
    if not re.fullmatch(r'PXD\d+', accession):
        raise ValueError('Invalid PRIDE project identity')
    projects = root / 'pride_projects'
    coverage_bytes = (projects / 'coverage.json').read_bytes()
    records = json.loads(coverage_bytes)
    matches = [row for row in records if row['accession'] == accession]
    if len(matches) != 1 or matches[0]['status'] != 'validated_primary_metadata':
        raise ValueError('Primary PRIDE project metadata has not been reconciled')
    record = matches[0]
    summary_bytes = (projects / accession / 'project_summary.json').read_bytes()
    if hashlib.sha256(summary_bytes).hexdigest() != record['summary_sha256']:
        raise ValueError('Primary PRIDE project summary checksum differs')
    summary = json.loads(summary_bytes)
    if summary['accession'] != accession:
        raise ValueError('Primary PRIDE project identity differs')
    result = dict(project=summary, designs=[], design_status='Submitted design collection is not installed.')
    collection = root / 'pride_sample_designs'
    pointer_path = collection / 'collection.json'
    if not pointer_path.exists():
        return result
    pointer = json.loads(pointer_path.read_text())
    directory = pointer['directory']
    if pointer.get('schema_version') != 1 or not re.fullmatch(r'[0-9a-f]{64}', directory):
        raise ValueError('Unsupported PRIDE design collection identity')
    folder = collection / directory
    if folder.resolve().parent != collection.resolve():
        raise ValueError('PRIDE design collection lies outside its installed root')
    receipt_bytes = (folder / 'validation.json').read_bytes()
    if hashlib.sha256(receipt_bytes).hexdigest() != pointer['validation_sha256']:
        raise ValueError('PRIDE design reconciliation receipt checksum differs')
    receipt = json.loads(receipt_bytes)
    if (receipt['status'] != 'source_designs_reconciled' or receipt['catalogue_projects'] != len(records)
            or receipt['bindings']['project_coverage_sha256'] != hashlib.sha256(coverage_bytes).hexdigest()):
        raise ValueError('PRIDE design collection project context differs')

    def read(name):
        if name != 'index.json' and not re.fullmatch(r'[0-9a-f]{64}\.json', name):
            raise ValueError('Invalid PRIDE design output identity')
        path = folder / name
        if path.resolve().parent != folder.resolve():
            raise ValueError('PRIDE design output lies outside its collection')
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != receipt['outputs'][name]:
            raise ValueError('PRIDE design output checksum differs')
        return json.loads(data)

    index = read('index.json')
    names = index['projects'][accession]
    if (index['schema_version'] != 1 or len(names) != len(set(names))
            or set(index['projects']) != {row['accession'] for row in records}):
        raise ValueError('PRIDE design index project identities differ')
    for name in names:
        design = read(name)
        headers, rows = design['headers'], design['rows']
        if (design['status'] != 'structurally_validated_design' or len(rows) != design['source_rows']
                or design['primary_file_metadata_sha256'].get(accession) != record['files_sha256']
                or any(not isinstance(header, str) for header in headers)):
            raise ValueError('Submitted PRIDE design coverage or project binding differs')
        for position, row in enumerate(rows, 2):
            if (row['source_row'] != position or len(row['cells']) != len(headers)
                    or any(not isinstance(cell, str) for cell in row['cells'])
                    or row['source_row_id'] != hashlib.sha256((design['source_sha256'] + ':' + str(position)).encode()).hexdigest()):
                raise ValueError('Submitted PRIDE design row identities differ')
        result['designs'].append(design)
    result['design_status'] = ('Submitted design rows retained; biological independence and control pairing are not inferred.'
                               if names else 'No submitted SDRF design is represented for this project in the installed collection.')
    return result


class DatasetCoveragePanel(QtWidgets.QWidget):
    """Search all catalogue records and inspect source sample labels and acquisition gaps."""

    def __init__(self, root=None, parent=None):
        super().__init__(parent)
        self.root = Path(root) if root is not None else audit_root()
        self.rows = []
        self.visible_rows = []
        self.links = None
        layout = QtWidgets.QVBoxLayout(self)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.search = QtWidgets.QLineEdit()
        self.search.setPlaceholderText('Find an accession, study title, source or assay')
        self.search.textChanged.connect(self._filter)
        layout.addWidget(self.search)
        button = QtWidgets.QPushButton('Refresh coverage')
        button.setToolTip('Reload the local catalogue coverage ledger and its acquisition and measurement status.')
        button.clicked.connect(self.refresh)
        layout.addWidget(button)
        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(['accession', 'source', 'study', 'assay', 'source files', 'local measurements'])
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.itemSelectionChanged.connect(self._select)
        layout.addWidget(self.table)
        self.details = QtWidgets.QLabel('Select a study to inspect its submitted sample labels.')
        self.details.setWordWrap(True)
        layout.addWidget(self.details)
        self.samples = QtWidgets.QTreeWidget()
        self.samples.setHeaderLabels(['Submitted sample / source field', 'Original label or value'])
        layout.addWidget(self.samples)
        self.refresh()

    def refresh(self):
        """Reload an audited local coverage ledger without implying public release coverage."""
        self.links = None
        try:
            self.rows, info = coverage(self.root)
            observed = info.get('observed_epoch')
            date = datetime.fromtimestamp(observed, timezone.utc).strftime('%Y-%m-%d %H:%M UTC') if observed else 'time unspecified'
            self.summary.setText(f"{len(self.rows)} catalogue records; "
                                 f"{info['local_runtime_catalogue_records']} with local measurements. "
                                 'Downloaded files and schema previews still require assay and comparison review. '
                                 f'Coverage snapshot: {date}. Source-file counts refer to the processing archive.')
        except FileNotFoundError:
            self.rows = []
            self.summary.setText('Dataset coverage data is not installed.')
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.rows = []
            self.summary.setText('Local acquisition coverage is not available: ' + str(error))
        self._filter()

    def _filter(self):
        query = self.search.text().casefold().strip()
        self.visible_rows = [r for r in self.rows if query in ' '.join(
            r.get(k, '') for k in ('accession', 'source', 'title', 'readout')).casefold()]
        self.table.blockSignals(True)
        self.table.clearSelection()
        self.table.setRowCount(len(self.visible_rows))
        for i, row in enumerate(self.visible_rows):
            total = row['supplementary_urls_in_plan']
            files = f"{row['source_files_with_acquisition_receipts']}/{total}" if int(total) else (
                'no files declared' if row['primary_sample_metadata']=='validated_sample_metadata' else 'not catalogued')
            if int(row.get('series_matrix_files_declared', 0)):
                matrices = f"{row.get('series_matrix_files_acquired', 0)}/{row['series_matrix_files_declared']} Series Matrices"
                files = (files + ' supplementary; ' + matrices) if int(total) else matrices
            if int(row.get('workbench_acquired_api_responses', 0)):
                files = row['workbench_acquired_api_responses'] + ' API responses; raw payloads pending'
            measurement = 'Available locally' if row['measurement_status']=='local runtime measurements present' else 'Pending measurements'
            values = [row['accession'], row['source'], row['title'], row['readout'], files, measurement]
            for j, value in enumerate(values):
                self.table.setItem(i, j, QtWidgets.QTableWidgetItem(value))
        self.table.blockSignals(False)
        self.samples.clear()
        self.details.setText('Select a study to inspect its submitted sample labels. Shared source identities and multi-channel samples remain explicit.')
        self.table.resizeColumnsToContents()
        self.table.setColumnWidth(2, 300)

    def _select(self):
        self.samples.clear()
        index = self.table.currentRow()
        if index < 0 or index >= len(self.visible_rows):
            return
        row = self.visible_rows[index]
        measurement = 'Measurements available locally' if row['measurement_status']=='local runtime measurements present' else 'Measurements pending'
        metadata = 'verified' if row['primary_sample_metadata'] in (
            'validated_sample_metadata', 'independently_reconciled_primary_metadata') else 'not collected'
        self.details.setText(f"{row['accession']}: {measurement}. "
                             f"Primary sample metadata: {metadata}. "
                             'Biological independence, technical replicate type and control pairing are not inferred. '
                             'Multiple channels and annotation variants are not extra replicates.')
        if int(row.get('series_matrix_files_declared', 0)):
            self.details.setText(self.details.text() +
                f" Series Matrices: {row.get('series_matrix_files_acquired', 0)}/{row['series_matrix_files_declared']} acquired; "
                f"{row.get('series_matrices_numeric_checked', 0)} with checked numeric structure; "
                f"{row.get('series_matrices_without_embedded_values', 0)} without embedded measurements; "
                f"{row.get('series_matrices_review_required', 0)} require review. "
                'Original units, scale and platform/probe mappings still require review before comparison.')
        try:
            if row['source'] == 'PRIDE':
                report = pride_samples(self.root, row['accession'])
                self.details.setText(f"{row['accession']}: {measurement}. Primary project/file metadata: verified. "
                                     'Biological independence, technical replicate type and control pairing are not inferred.')
                project = QtWidgets.QTreeWidgetItem(self.samples, [row['accession'], 'Submitted project metadata'])
                for key, value in report['project'].items():
                    text = value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)
                    QtWidgets.QTreeWidgetItem(project, [key, text])
                for design in report['designs']:
                    parent = QtWidgets.QTreeWidgetItem(self.samples, ['Submitted SDRF', design['source_url']])
                    QtWidgets.QTreeWidgetItem(parent, ['Metadata origin', design['metadata_origin']])
                    QtWidgets.QTreeWidgetItem(parent, ['Biological truth status', design['biological_truth_status']])
                    for submitted in design['rows']:
                        item = QtWidgets.QTreeWidgetItem(parent, ['Source row ' + str(submitted['source_row']),
                            '; '.join(submitted['submitted_source_names'])])
                        for position, (header, cell) in enumerate(zip(design['headers'], submitted['cells']), 1):
                            QtWidgets.QTreeWidgetItem(item, [str(position) + ': ' + header, cell])
                        for link in submitted['file_links']:
                            if link['project'] == row['accession']:
                                child = QtWidgets.QTreeWidgetItem(item, ['Submitted file', link['submitted_file']])
                                QtWidgets.QTreeWidgetItem(child, ['Exact filename association', link['status']])
                                QtWidgets.QTreeWidgetItem(child, ['Primary file identifiers', json.dumps(link['exact_primary_file_ids'])])
                iterator = QtWidgets.QTreeWidgetItemIterator(self.samples)
                while iterator.value() is not None:
                    item = iterator.value()
                    for column in range(2):
                        item.setToolTip(column, item.text(column))
                    iterator += 1
                self.details.setText(self.details.text() + ' ' + report['design_status'] +
                    ' Community annotations remain distinct from independently verified biological truth. '
                    'Repeated filenames, fractions and channels do not establish independent replicates; '
                    'file identifiers do not prove downloaded raw payloads or measurement-column pairing.')
                self.samples.resizeColumnToContents(0)
                return
            if row['source'] == 'Metabolomics Workbench':
                report = workbench_samples(self.root, row['accession'])
                for position, sample in enumerate(report['submitted_rows'], 1):
                    parent = QtWidgets.QTreeWidgetItem(self.samples, [sample['local_sample_id'], sample.get('sample_source', '')])
                    QtWidgets.QTreeWidgetItem(parent, ['Submitted factor row', str(position)])
                    for key, value in sample.items():
                        if key != 'local_sample_id':
                            text = value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)
                            QtWidgets.QTreeWidgetItem(parent, [key, text])
                for analysis in report['analysis_inventory']:
                    parent = QtWidgets.QTreeWidgetItem(self.samples, [analysis['analysis'], 'Submitted analysis metadata'])
                    for key, value in analysis.items():
                        if key != 'analysis':
                            text = value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)
                            QtWidgets.QTreeWidgetItem(parent, [key, text])
                self.details.setText(self.details.text() +
                    ' Workbench factor rows and analysis representations remain separate. Source units and row-count discrepancies are preserved; '
                    'raw filenames are declarations, not downloaded payloads. Cross-study comparisons remain pending.')
                self.samples.resizeColumnToContents(0)
                return
            report = source_samples(self.root, row['accession'])
            if self.links is None:
                try:
                    self.links = shared_samples(self.root)
                except (OSError, ValueError, KeyError, TypeError):
                    self.links = {}
            for sample in report['samples']:
                parent = QtWidgets.QTreeWidgetItem(self.samples, [sample['geo_accession'], sample['title']])
                QtWidgets.QTreeWidgetItem(parent, ['Source channels', str(len(sample['channels']))])
                linked = self.links.get(sample['geo_accession'])
                QtWidgets.QTreeWidgetItem(parent, ['Same submitted sample in', '; '.join(linked) if linked else 'Cross-study links unavailable'])
                for channel in sample['channels']:
                    child = QtWidgets.QTreeWidgetItem(parent, ['Channel ' + str(channel['position']), channel['source']])
                    QtWidgets.QTreeWidgetItem(child, ['Organisms', '; '.join(channel['organisms'])])
                    for item in channel['characteristics']:
                        QtWidgets.QTreeWidgetItem(child, [item['tag'], item['value']])
                for relation in sample['relations']:
                    QtWidgets.QTreeWidgetItem(parent, [relation['type'], relation['target']])
            self.samples.resizeColumnToContents(0)
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.details.setText(self.details.text() + ' Sample metadata unavailable: ' + str(error))
