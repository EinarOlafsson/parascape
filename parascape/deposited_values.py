"""Read original deposited measurements without erasing their experimental context.

Assays outside gene-expression matrices still need exact feature rows, original
sample labels, source units and auditable numerical conversion. These packs retain
literal source tokens beside binary64 values and use source row numbers to keep
repeated feature names distinct. Readers require an independent full-cell
reconciliation receipt before exposing values. Matching a label does not establish
gene equivalence, compatible scales or independent biological replicates, and no
normalization or missing-value imputation is performed during browsing.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from . import paths
from . import source_design_links


def data_root() -> Path:
    """Locate separately installed original measurements or a local collection override."""
    return Path(os.environ.get('PARASCAPE_DEPOSITED_VALUES', str(Path(paths.data_dir()) / 'deposited_values')))


def _sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def manifest(pack: str | Path, *, verify_outputs: bool = True) -> dict:
    """Require complete source-cell reconciliation and pin every runtime output."""
    pack = Path(pack)
    receipt = json.loads((pack / 'validation.json').read_text())
    if receipt['status'] != 'source_cells_reconciled' or _sha(pack / 'manifest.json') != receipt['manifest_sha256']:
        raise ValueError('Deposited values lack matching independent source-cell reconciliation')
    info = json.loads((pack / 'manifest.json').read_text())
    columns, samples = info['sample_columns'], info['source_samples']
    if (info.get('schema_version') != 1 or info['features'] <= 0 or not columns
            or len(columns) != len(set(columns)) or len(samples) != len(columns)
            or len(samples) != len(set(samples)) or any(not isinstance(c, str) or not c for c in [*columns, *samples])
            or set(columns) & {'source_row', 'feature_label'}
            or set(info['outputs']) != {'values.parquet', 'source_tokens.parquet', 'samples.json', 'source_metadata.json'}):
        raise ValueError('Unsupported deposited value schema, sample identities or dimensions')
    if (receipt['source_sha256'] != info['input_sha256'] or receipt['verified_rows'] != info['features']
            or receipt['verified_cells'] != info['features'] * len(columns)
            or not receipt['literal_tokens_verified']):
        raise ValueError('Deposited source-cell validation coverage differs')
    if verify_outputs:
        for name, expected in info['outputs'].items():
            if _sha(pack / name) != expected:
                raise ValueError('Deposited value output checksum differs: ' + name)
    return info


def catalogue(root: str | Path | None = None) -> tuple[list[dict], list[str]]:
    """List independently reconciled packs while retaining assay and source-analysis identity."""
    root = Path(root) if root is not None else data_root()
    packs, errors = [], []
    indexes = [root / 'catalogue.json', *sorted(root.glob('*/catalogue.json'))]
    for index in indexes:
        namespace = index.parent
        if namespace.name.startswith('.') or not index.exists():
            continue
        try:
            if namespace != root and namespace.resolve().parent != root.resolve():
                raise ValueError('Deposited namespace lies outside its collection')
            names = json.loads(index.read_text())['packs']
            if len(names) != len(set(names)):
                raise ValueError('Duplicated pack identities in deposited catalogue')
            for name in names:
                try:
                    if not isinstance(name, str) or Path(name).name != name or name.startswith('.'):
                        raise ValueError('Invalid deposited pack path')
                    folder = namespace / name
                    if folder.resolve().parent != namespace.resolve():
                        raise ValueError('Deposited pack lies outside its collection')
                    info = manifest(folder, verify_outputs=False)
                    packs.append(dict(path=str(folder), study=info['study'], analysis=info['source_analysis'],
                        assay=info['assay'], units=info['units'], features=info['features'], samples=len(info['sample_columns'])))
                except (OSError, ValueError, KeyError, TypeError) as error:
                    errors.append(str(name) + ': ' + str(error))
        except (OSError, ValueError, KeyError, TypeError) as error:
            errors.append(str(index) + ': ' + str(error))
    return packs, errors


def _context(pack):
    pack = Path(pack)
    info = manifest(pack)
    values, tokens = pq.ParquetFile(pack / 'values.parquet'), pq.ParquetFile(pack / 'source_tokens.parquet')
    expected = ['source_row', 'feature_label', *info['sample_columns']]
    for source in (values, tokens):
        schema = source.schema_arrow
        if source.metadata.num_rows != info['features'] or schema.names != expected:
            raise ValueError('Deposited feature/sample axes differ')
        if (not pa.types.is_int64(schema.field('source_row').type) or not
                (pa.types.is_string(schema.field('feature_label').type) or pa.types.is_large_string(schema.field('feature_label').type))):
            raise ValueError('Unsupported deposited source identity types')
    if (any(not pa.types.is_float64(values.schema_arrow.field(c).type) for c in info['sample_columns'])
            or any(not (pa.types.is_string(tokens.schema_arrow.field(c).type) or
                        pa.types.is_large_string(tokens.schema_arrow.field(c).type)) for c in info['sample_columns'])):
        raise ValueError('Unsupported deposited numeric or literal token types')
    metadata = json.loads((pack / 'samples.json').read_text())
    samples = metadata['samples']
    if (metadata['study'] != info['study'] or metadata['source_analysis'] != info['source_analysis']
            or [s['column'] for s in samples] != info['sample_columns']
            or [s['source_sample'] for s in samples] != info['source_samples']):
        raise ValueError('Deposited sample metadata differs from matrix identities')
    return info, values, samples


def _column_context(pack, samples, audit_root=None):
    links, status = source_design_links.column_links(pack, audit_root=audit_root)
    return [dict(sample, source_design_metadata=links.get(sample['column'], {'association_status': status}))
            for sample in samples]


def column_context(pack: str | Path, *, audit_root: str | Path | None = None) -> list[dict]:
    """Inspect original source columns and separately pinned submitted metadata without a feature query."""
    _, _, samples = _context(pack)
    return _column_context(pack, samples, audit_root)


def _reported_parent_rows(pack, info, labels):
    if info['assay'] != 'deposited original transcript-table quantities':
        raise ValueError('Reported parent gene matching requires an original transcript-table source')
    metadata = json.loads((pack / 'source_metadata.json').read_text())
    if metadata.get('original_feature_column') != 'transcript_id' or 'gene_id' not in metadata.get('original_annotation_columns', []):
        raise ValueError('This transcript source has no checked reported gene_id annotation')
    annotations = metadata.get('row_annotations')
    if not isinstance(annotations, list) or len(annotations) != info['features']:
        raise ValueError('Reported parent gene annotation coverage differs')
    axes = pq.read_table(pack / 'values.parquet', columns=['source_row', 'feature_label']).to_pylist()
    rows = []; found = set(); wanted = set(labels)
    for position, (annotation, axis) in enumerate(zip(annotations, axes, strict=True)):
        if (type(annotation.get('source_row')) is not int or annotation['source_row'] != position
                or axis['source_row'] != position or annotation.get('feature_label') != axis['feature_label']
                or not isinstance(annotation.get('source_fields', {}).get('gene_id'), str)):
            raise ValueError('Reported parent gene annotations differ from original transcript rows')
        parent = annotation['source_fields']['gene_id']
        if parent in wanted:
            rows.append(position); found.add(parent)
    return rows, dict(schema_version=1, match='reported_gene_id', exact_query_literals=labels,
        unmatched_literals=[label for label in labels if label not in found], selected_source_rows=rows,
        original_feature_column='transcript_id', original_annotation_column='gene_id',
        source_metadata_sha256=info['outputs']['source_metadata.json'],
        scope='Exact reported parent gene_id within this source only; every original transcript row remains separate',
        canonical_gene_equivalence_established=False, gene_aggregation_performed=False,
        transcript_suffix_stripping_performed=False, biological_independence_inferred=False)


def profiles(pack: str | Path, features, *, audit_root: str | Path | None = None,
             match: str = 'feature_label') -> pd.DataFrame:
    """Return exact labels or reported parent gene IDs, preserving original transcript rows.

    ``match='reported_gene_id'`` uses only the checked gene_id annotation in an
    original transcript table. It neither sums transcripts nor maps gene versions.
    The source-specific query record is returned in ``frame.attrs['identifier_query']``.
    """
    pack = Path(pack)
    labels = features.splitlines() if isinstance(features, str) else list(features)
    if not labels or any(not isinstance(label, str) or not label for label in labels):
        raise ValueError('Enter exact deposited feature labels, one per line')
    info, _, samples = _context(pack)
    labels = list(dict.fromkeys(labels))
    query_record = None
    if match == 'reported_gene_id':
        rows, query_record = _reported_parent_rows(pack, info, labels)
        filters = [('source_row', 'in', rows)]
    elif match == 'feature_label':
        filters = [('feature_label', 'in', labels)]
    else:
        raise ValueError('Unknown deposited label matching mode')
    matrix = pq.read_table(pack / 'values.parquet', filters=filters).to_pandas()
    tokens = pq.read_table(pack / 'source_tokens.parquet', filters=[('source_row', 'in', matrix.source_row.tolist())]).to_pandas()
    if (matrix.source_row.isna().any() or matrix.source_row.duplicated().any()
            or not matrix[['source_row', 'feature_label']].equals(tokens[['source_row', 'feature_label']])):
        raise ValueError('Deposited values and literal-token source rows differ')
    frame = matrix.melt(id_vars=['source_row', 'feature_label'], var_name='sample_column', value_name='value')
    literal = tokens.melt(id_vars=['source_row', 'feature_label'], var_name='sample_column', value_name='source_token')
    frame['source_token'] = literal['source_token']
    lookup = {s['column']: s for s in _column_context(pack, samples, audit_root)}
    frame['source_sample'] = frame.sample_column.map(lambda c: lookup[c]['source_sample'])
    frame['sample_metadata'] = frame.sample_column.map(lambda c: json.dumps(lookup[c], ensure_ascii=False))
    source_metadata = json.loads((pack / 'source_metadata.json').read_text())
    annotations = source_metadata.get('row_annotations')
    if annotations is not None:
        if not isinstance(annotations, list) or len(annotations) != info['features'] or any(
                not isinstance(row, dict) or row.get('source_row') != i for i, row in enumerate(annotations)):
            raise ValueError('Deposited row-annotation identities differ')
        selected = {int(row.source_row): row.feature_label for row in matrix.itertuples()}
        if any(i < 0 or i >= len(annotations) or annotations[i].get('feature_label') != label for i, label in selected.items()):
            raise ValueError('Deposited row-annotation labels differ')
        frame['source_row_metadata'] = frame.source_row.map(lambda i: json.dumps(annotations[int(i)], ensure_ascii=False))
    for field in ('study', 'source_analysis', 'assay', 'quantity', 'units', 'quantity_domain', 'source_url', 'input_sha256', 'comparison_scope'):
        frame[field] = info[field]
    frame['matrix_sha256'] = info['outputs']['values.parquet']
    frame['literal_tokens_sha256'] = info['outputs']['source_tokens.parquet']
    frame['manifest_sha256'] = _sha(pack / 'manifest.json')
    if query_record is not None:
        frame.attrs['identifier_query'] = query_record
    return frame


def technical_qc(pack: str | Path) -> pd.DataFrame:
    """Describe every original sample across all source rows without changing values or units."""
    from . import measurement_qc as Q

    pack = Path(pack)
    info, source, samples = _context(pack)
    matrix = source.read().to_pandas()
    if matrix.source_row.tolist() != list(range(info['features'])):
        raise ValueError('Deposited source-row ordinals differ from the full matrix')
    result = Q.summarize(matrix.set_index('source_row')[info['sample_columns']], domain=info['quantity_domain'])
    lookup = {s['column']: s for s in samples}
    result['source_sample'] = result['sample'].map(lambda c: lookup[c]['source_sample'])
    result['sample_metadata'] = result['sample'].map(lambda c: json.dumps(lookup[c], ensure_ascii=False))
    for field in ('study', 'source_analysis', 'assay', 'quantity', 'units', 'source_url', 'input_sha256', 'comparison_scope'):
        result[field] = info[field]
    result['matrix_sha256'] = info['outputs']['values.parquet']
    result['manifest_sha256'] = _sha(pack / 'manifest.json')
    return result
