"""Stream a checked sparse query into a complete local CSV and bounded preview."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import tempfile
import threading

from .sparse_values import SparseDataset, _sha

FIELDS = ['entry_index', 'source_row', 'source_column', 'row_token', 'column_token',
          'source_token', 'value', 'source_feature_axis_row', 'source_feature_label',
          'original_feature_fields', 'source_cell_axis_row']


class QueryCancelled(Exception):
    """The source/query changed before the complete result was ready."""


@dataclass
class SparseQueryResult:
    """Keep a bounded preview and complete source-record export for one query."""
    temporary: tempfile.TemporaryDirectory
    csv_path: Path
    record: dict
    preview: list[dict]
    feature_matches: list[dict]

    def close(self):
        """Release the temporary complete-query spool after browsing or export."""
        self.temporary.cleanup()

    def export(self, destination: str | Path) -> int:
        """Save all matching stored records and the exact treatment sidecar."""
        destination = Path(destination)
        sidecar = destination.with_suffix(destination.suffix + '.json')
        if destination.exists() or sidecar.exists():
            raise FileExistsError('Choose a new filename for the CSV and its treatment record')
        if _sha(self.csv_path) != self.record['csv_sha256']:
            raise ValueError('Prepared sparse query CSV changed')
        # Exclusive creation protects an unrelated file if another process writes
        # between the dialog selection and the actual save.
        created = []
        try:
            with destination.open('xb') as target, self.csv_path.open('rb') as source:
                created.append(destination)
                shutil.copyfileobj(source, target)
            with sidecar.open('x') as target:
                created.append(sidecar)
                target.write(json.dumps(self.record, indent=2, allow_nan=False) + '\n')
        except Exception:
            for path in created:
                path.unlink(missing_ok=True)
            raise
        return self.record['exported_records']


def prepare_query(descriptor, labels: list[str], *, coordinate_root=None, axis_root=None,
                  cancel: threading.Event | None = None, preview_limit: int = 500,
                  temporary_root=None) -> SparseQueryResult:
    """Read all matching records in source order; preview limits do not limit export."""
    def checkpoint():
        if cancel is not None and cancel.is_set():
            raise QueryCancelled()

    if type(preview_limit) is not int or preview_limit < 0:
        raise ValueError('Preview limit must be a nonnegative integer')
    if not labels or any(not isinstance(label, str) or label == '' for label in labels):
        raise ValueError('Enter at least one exact original feature label')
    checkpoint()
    dataset = SparseDataset(descriptor, coordinate_root=coordinate_root, axis_root=axis_root)
    checkpoint()
    labels = list(dict.fromkeys(labels))
    features = {}
    unmatched = []
    for label in labels:
        matched = False
        for batch in dataset.iter_feature_matches(label):
            checkpoint()
            for row in batch.to_pylist():
                features[row['source_axis_row']] = row
                matched = True
        if not matched:
            unmatched.append(label)
    positions = sorted(features)
    temporary = tempfile.TemporaryDirectory(prefix='parascape-sparse-query-', dir=temporary_root)
    try:
        path = Path(temporary.name) / 'original_sparse_records.csv'
        preview, count = [], 0
        cell_coordinate = dataset.info['axis_roles']['barcodes']['source_coordinate']
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            for batch in dataset.iter_entries(positions):
                checkpoint()
                for row in batch.to_pylist():
                    feature_position = row[dataset.feature_coordinate]
                    feature = features[feature_position]
                    row.update(source_feature_axis_row=feature_position,
                               source_feature_label=feature['first_field'],
                               original_feature_fields=json.dumps(feature['original_fields'], ensure_ascii=False),
                               source_cell_axis_row=row[cell_coordinate])
                    writer.writerow(row)
                    count += 1
                    if len(preview) < preview_limit:
                        preview.append(row)
        checkpoint()
        # Read original cell metadata only for displayed records. The complete
        # original cell-axis artifact remains pinned in the export sidecar.
        cells = {}
        needed = {row['source_cell_axis_row'] for row in preview}
        for batch in dataset.iter_axis_rows('barcodes', sorted(needed)):
            checkpoint()
            cells.update({row['source_axis_row']: row for row in batch.to_pylist()})
        for row in preview:
            cell = cells[row['source_cell_axis_row']]
            row.update(source_cell_label=cell['first_field'], original_cell_fields=cell['original_fields'])
        record = dataset.query_record(positions)
        record.update(created_utc=datetime.now(timezone.utc).isoformat(), query_labels=labels,
                      unmatched_labels=unmatched, matched_feature_rows=[features[p] for p in positions],
                      exported_records=count, csv_sha256=_sha(path), csv_fields=FIELDS,
                      displayed_records=len(preview), preview_limit=preview_limit,
                      cell_metadata_scope='CSV retains original cell-axis positions; all original cell metadata remain in the checksum-pinned barcodes artifact. Cell labels and fields in the browser preview are exact original rows.',
                      query_materializer_sha256=_sha(Path(__file__)))
        return SparseQueryResult(temporary, path, record, preview, [features[p] for p in positions])
    except Exception:
        temporary.cleanup()
        raise
