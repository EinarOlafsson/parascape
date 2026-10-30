"""Query original protein identifiers and qualified source declarations without gene-level abundance inference.

Protein groups, submitted names and reference gene fields provide different kinds
of evidence. This reader preserves exact components and every ambiguous occurrence,
requires the original measurement manifest and independently reconciled association
bytes, and returns whole original group measurements. It never converts symbols,
collapses groups, infers a search database or establishes current gene-version
identity. Multiple matching declarations annotate one original measurement rather
than creating extra measurements or independent biological samples.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re

import pandas as pd
import pyarrow.parquet as pq

from . import deposited_values as D
from . import paths

NAMESPACES = {
    'original_protein_group_component': 'Original protein-group component',
    'submitted_gene_label': 'Submitted gene label',
    'reference_declared_gene_identifier': 'Reference gene declaration · candidate',
    'reference_GN_gene_name_candidate': 'Reference GN name · candidate',
}


def data_root() -> Path:
    """Locate separately installed identifier evidence or an explicit local override."""
    audit = Path(os.environ.get('PARASCAPE_AUDIT', str(Path(paths.data_dir()) / 'processing_audit')))
    return Path(os.environ.get('PARASCAPE_PROTEIN_IDENTIFIERS', str(audit / 'protein_identifier_evidence')))


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(pack: str | Path, queries, *, namespace: str, root: str | Path | None = None) -> pd.DataFrame:
    """Return exact, manifest-bound identifier evidence while preserving ambiguous source occurrences."""
    if namespace not in NAMESPACES:
        raise ValueError('Choose an original identifier or a declared source evidence class')
    labels = queries.splitlines() if isinstance(queries, str) else list(queries)
    if not labels or any(not isinstance(label, str) or not label for label in labels):
        raise ValueError('Enter exact source identifiers, one per line')
    pack = Path(pack)
    info = D.manifest(pack)
    root = Path(root) if root is not None else data_root()
    if not (root / 'collection.json').exists():
        raise ValueError('Protein identifier evidence is not installed; exact original feature-label lookup remains available')
    pointer = json.loads((root / 'collection.json').read_text())
    directory = pointer['directory']
    if pointer.get('schema_version') != 1 or not re.fullmatch(r'[0-9a-f]{64}', directory):
        raise ValueError('Unsupported protein identifier collection identity')
    folder = root / directory
    if folder.resolve().parent != root.resolve():
        raise ValueError('Protein identifier collection lies outside its installed root')
    receipt_bytes = (folder / 'validation.json').read_bytes()
    if hashlib.sha256(receipt_bytes).hexdigest() != pointer['validation_sha256']:
        raise ValueError('Protein identifier reconciliation receipt checksum differs')
    receipt = json.loads(receipt_bytes)
    if receipt['status'] != 'all_identifier_evidence_reconciled' or _sha(folder / 'index.json') != receipt['outputs']['index.json']:
        raise ValueError('Protein identifier evidence has not been reconciled')
    index = json.loads((folder / 'index.json').read_text())
    if index.get('schema_version') != 1 or len(index['measurements']) != receipt['measurement_packs']:
        raise ValueError('Protein identifier collection coverage differs')
    digest = _sha(pack / 'manifest.json')
    matches = [row for row in index['measurements'] if row['measurement_manifest_sha256'] == digest]
    if len(matches) != 1:
        raise ValueError('No reconciled identifier evidence for this exact measurement manifest')
    source = matches[0]
    if (source['input_sha256'] != info['input_sha256'] or source['source_url'] != info['source_url']
            or source['features'] != info['features'] or source['file'] != digest + '.parquet'):
        raise ValueError('Protein identifier measurement context differs')
    path = folder / source['file']
    if path.resolve().parent != folder.resolve() or _sha(path) != receipt['outputs'][source['file']]:
        raise ValueError('Protein identifier output checksum or location differs')
    frame = pq.read_table(path, filters=[('query_namespace', '=', namespace), ('query_literal', 'in', list(dict.fromkeys(labels)))]).to_pandas()
    if ((frame.source_row < 0).any() or (frame.source_row >= info['features']).any()
            or not frame.query_namespace.eq(namespace).all() or not frame.query_literal.isin(labels).all()):
        raise ValueError('Protein identifier source-row or query scope differs')
    axes = pq.read_table(pack / 'values.parquet', columns=['source_row', 'feature_label'],
                         filters=[('source_row', 'in', frame.source_row.unique().tolist())]).to_pandas().set_index('source_row')
    if any(axes.loc[int(row.source_row), 'feature_label'] != row.feature_label for row in frame.itertuples()):
        raise ValueError('Protein identifier feature labels differ from original measurement rows')
    frame['measurement_manifest_sha256'] = digest
    frame['identifier_collection_sha256'] = pointer['validation_sha256']
    frame['identifier_scope'] = index['scope']
    return frame


def profiles(pack: str | Path, queries, *, namespace: str, root: str | Path | None = None) -> pd.DataFrame:
    """Return whole original protein-group values with all matching identifier evidence attached."""
    evidence = records(pack, queries, namespace=namespace, root=root)
    if evidence.empty:
        return pd.DataFrame(columns=['source_row', 'feature_label', 'source_sample', 'value', 'source_token',
                                     'units', 'sample_metadata', 'source_url', 'comparison_scope', 'identifier_evidence'])
    frame = D.profiles(pack, evidence.feature_label.unique().tolist())
    frame = frame[frame.source_row.isin(evidence.source_row.unique())].copy()
    grouped = {int(row): group.to_dict('records') for row, group in evidence.groupby('source_row', sort=False)}
    frame['identifier_evidence'] = frame.source_row.map(lambda row: json.dumps(grouped[int(row)], ensure_ascii=False))
    frame['identifier_scope'] = evidence.identifier_scope.iloc[0]
    return frame.reset_index(drop=True)
