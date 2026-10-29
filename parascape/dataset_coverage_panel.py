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
        metadata = 'verified' if row['primary_sample_metadata']=='validated_sample_metadata' else 'not collected'
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
