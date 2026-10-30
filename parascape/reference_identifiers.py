"""Inspect original reference header declarations without converting them to canonical gene mappings.

A deposited FASTA may be a complete database, a contaminant list or a small
sequence set. This reader preserves repeated header rows and explicit source
identifier declarations, requiring independent full-header reconciliation and
pinned input/output hashes. Matching a source gene field does not establish the
current gene version, protein-group equivalence or the database used for a
measurement table. Gene symbols and unrecognised fields are never silently
substituted for source gene identifiers.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from . import paths


def data_root() -> Path:
    """Locate original reference identifiers in the external audit data or an explicit override."""
    audit = Path(os.environ.get('PARASCAPE_AUDIT', str(Path(paths.data_dir()) / 'processing_audit')))
    return Path(os.environ.get('PARASCAPE_REFERENCE_IDENTIFIERS', str(audit / 'reference_identifiers')))


def _sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def catalogue(root: str | Path | None = None) -> tuple[list[dict], dict]:
    """List all planned reference sources, keeping pending acquisition and source review explicit."""
    root = Path(root) if root is not None else data_root()
    data = (root / 'catalogue.json').read_bytes()
    status = json.loads((root / 'status.json').read_text())
    if hashlib.sha256(data).hexdigest() != status['catalogue_sha256']:
        raise ValueError('Reference identifier catalogue is updating or its checksum differs')
    info = json.loads(data)
    sources = info['sources']
    if (info.get('schema_version') != 1 or len(sources) != status['planned']
            or len({row['source_url'] for row in sources}) != len(sources)):
        raise ValueError('Reference identifier source identities or coverage differ')
    return sources, status


def _index(root, source):
    if source['status'] != 'all_header_cells_reconciled' or not re.fullmatch(r'[0-9a-f]{64}', source['pack']):
        raise ValueError('Original reference headers have not been reconciled')
    folder = root / source['pack']
    if folder.resolve().parent != root.resolve():
        raise ValueError('Reference identifier pack lies outside its collection')
    receipt_bytes = (folder / 'validation.json').read_bytes()
    if hashlib.sha256(receipt_bytes).hexdigest() != source['validation_sha256']:
        raise ValueError('Reference header reconciliation receipt checksum differs')
    receipt = json.loads(receipt_bytes)
    if (receipt['status'] != 'all_header_cells_reconciled'
            or receipt['bindings']['input_sha256'] != source['input_sha256']
            or set(receipt['outputs']) != {'manifest.json', 'headers.parquet'}):
        raise ValueError('Reference header reconciliation input context differs')
    for name, digest in receipt['outputs'].items():
        if _sha(folder / name) != digest:
            raise ValueError('Reference identifier output checksum differs')
    info = json.loads((folder / 'manifest.json').read_text())
    table = pq.ParquetFile(folder / 'headers.parquet')
    schema = table.schema_arrow
    expected = ['header_row', 'primary_identifier', 'original_header', 'declared_identifier_fields']
    if (info.get('schema_version') != 1 or info['input_sha256'] != source['input_sha256']
            or info['fields'] != expected or schema.names != expected
            or table.metadata.num_rows != info['header_rows'] or info['header_rows'] != receipt['counts']['headers']
            or info['header_rows'] != source['counts']['headers'] or not pa.types.is_int64(schema.field(0).type)
            or any(not (pa.types.is_string(schema.field(name).type) or pa.types.is_large_string(schema.field(name).type)) for name in expected[1:])):
        raise ValueError('Reference header row/schema coverage differs')
    return folder, info


def lookup(source_url: str, value: str, *, kind: str = 'primary', root: str | Path | None = None) -> pd.DataFrame:
    """Match an exact original primary identifier or explicit gene declaration, preserving repeated rows."""
    if kind not in ('primary', 'gene') or not isinstance(value, str) or not value:
        raise ValueError('Choose an exact primary identifier or explicit source gene declaration')
    root = Path(root) if root is not None else data_root()
    sources, _ = catalogue(root)
    matches = [row for row in sources if row['source_url'] == source_url]
    if len(matches) != 1 or matches[0]['status'] != 'all_header_cells_reconciled':
        raise ValueError('Selected reference source is not independently reconciled')
    source = matches[0]
    folder, info = _index(root, source)
    filters = [('primary_identifier', '=', value)] if kind == 'primary' else None
    frame = pq.read_table(folder / 'headers.parquet', filters=filters).to_pandas()
    fields = frame.declared_identifier_fields.map(json.loads)
    selected = fields.map(lambda rows: [row for row in rows if row['key'] in ('gene', 'gene_id') and row['literal_value'] == value])
    if kind == 'gene':
        frame = frame[selected.map(bool)].copy()
        selected = selected.loc[frame.index]
    if frame.header_row.isna().any() or frame.header_row.duplicated().any() or (frame.header_row < 0).any():
        raise ValueError('Original reference header row identities differ')
    frame['matching_gene_declarations'] = selected.map(lambda rows: json.dumps(rows, ensure_ascii=False))
    frame['match_kind'] = 'original_primary_identifier' if kind == 'primary' else 'explicit_source_gene_declaration'
    frame['matched_literal'] = value
    for key in ('source_url', 'filename', 'input_sha256'):
        frame[key] = source[key]
    frame['declared_projects'] = json.dumps(source['declared_projects'])
    frame['header_index_sha256'] = info['outputs']['headers.parquet']
    frame['scope'] = info['scope']
    return frame.reset_index(drop=True)
