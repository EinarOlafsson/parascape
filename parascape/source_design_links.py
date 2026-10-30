"""Read source-column metadata associations without guessing experimental sample roles.

Some quantity headers explicitly contain a filename. An independently reconciled
collection can associate that literal name with primary file records and submitted
SDRF rows, preserving duplicate headings and repeated occurrences. Exact experiment
labels can also identify source design candidates; input-table linkage remains
a separate question. These links do
not establish biological independence, control pairing, downloaded raw payloads or
compatible scales. Associations require the exact measurement manifest and source
column identities; missing collections and unsupported header roles stay visible.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re

from . import paths


def column_links(pack: str | Path, *, audit_root: str | Path | None = None) -> tuple[dict, str]:
    """Load pinned associations for one exact measurement manifest, retaining explicit gaps."""
    root = Path(audit_root) if audit_root is not None else Path(os.environ.get(
        'PARASCAPE_AUDIT', str(Path(paths.data_dir()) / 'processing_audit')))
    collection = root / 'protein_column_design_links'
    pointer_path = collection / 'collection.json'
    if not pointer_path.exists():
        return {}, 'Source-column association collection is not installed.'
    pointer = json.loads(pointer_path.read_text())
    directory = pointer['directory']
    if pointer.get('schema_version') != 1 or not re.fullmatch(r'[0-9a-f]{64}', directory):
        raise ValueError('Unsupported source-column association collection identity')
    folder = collection / directory
    if folder.resolve().parent != collection.resolve():
        raise ValueError('Source-column association collection lies outside its installed root')
    receipt_bytes = (folder / 'validation.json').read_bytes()
    if hashlib.sha256(receipt_bytes).hexdigest() != pointer['validation_sha256']:
        raise ValueError('Source-column reconciliation receipt checksum differs')
    receipt = json.loads(receipt_bytes)
    if receipt.get('schema_version') != 1 or receipt['status'] != 'exact_source_filename_associations_reconciled':
        raise ValueError('Source-column filename associations have not been reconciled')

    def read(name):
        if name != 'index.json' and not re.fullmatch(r'[0-9a-f]{64}\.json', name):
            raise ValueError('Invalid source-column association output identity')
        path = folder / name
        if path.resolve().parent != folder.resolve():
            raise ValueError('Source-column association output lies outside its collection')
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != receipt['outputs'][name]:
            raise ValueError('Source-column association output checksum differs')
        return json.loads(data)

    index = read('index.json')
    if index.get('schema_version') != 1 or len(index['manifests']) != receipt['measurement_packs']:
        raise ValueError('Source-column association index coverage differs')
    manifest_bytes = (Path(pack) / 'manifest.json').read_bytes()
    digest = hashlib.sha256(manifest_bytes).hexdigest()
    if digest not in index['manifests']:
        return {}, 'This exact measurement manifest is not represented in the installed association collection.'
    report = read(index['manifests'][digest])
    info = json.loads(manifest_bytes)
    columns = report['columns']
    if (report.get('schema_version') != 1 or report['measurement_manifest_sha256'] != digest
            or report['input_sha256'] != info['input_sha256'] or report['study'] != info['study']
            or report['source_url'] != info['source_url'] or report['source_analysis'] != info['source_analysis']
            or [row['column'] for row in columns] != info['sample_columns']
            or [row['source_column'] for row in columns] != info['source_samples']):
        raise ValueError('Source-column association measurement context differs')
    for column in columns:
        if column['biological_independence'] != 'not inferred':
            raise ValueError('Unsupported biological independence inferred from source filenames')
        for row in [*column['submitted_factor_rows'], *column.get('submitted_experiment_rows', [])]:
            if (len(row['headers']) != len(row['cells']) or any(not isinstance(s, str) for s in [*row['headers'], *row['cells']])
                    or row['source_row_id'] != hashlib.sha256((row['source_sha256'] + ':' + str(row['source_row'])).encode()).hexdigest()):
                raise ValueError('Submitted source-column factor row identity differs')
    context = dict(collection=pointer, measurement_manifest_sha256=digest,
                   primary_file_metadata_sha256=report['primary_file_metadata_sha256'])
    return {row['column']: dict(row, association_provenance=context) for row in columns}, 'Exact filename associations only; biological sample and control pairing remain unreviewed.'
