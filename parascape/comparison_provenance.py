"""Check exact comparison inputs and export reproducible treatment records.

Submitted sample identity prevents known repeated sequencing from becoming
extra replicates. Distinct accessions alone never prove biological independence.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd


def selection(counts, group_a, group_b, metadata=None):
    """Validate selected run identities and retain known/unknown sample relationships."""
    a, b = list(group_a), list(group_b)
    if not counts.columns.is_unique or not counts.index.is_unique:
        raise ValueError('Comparison requires unique run columns and gene identifiers')
    if len(set(a)) != len(a) or len(set(b)) != len(b):
        raise ValueError('Repeated run selections cannot count as extra replicates')
    if set(a) & set(b):
        raise ValueError('Groups A and B must not share a run')
    missing = [run for run in a + b if run not in counts.columns]
    if missing:
        raise ValueError('Requested runs are absent from the matrix: ' + ', '.join(map(str, missing)))
    if len(a) < 2 or len(b) < 2:
        raise ValueError('each group needs at least two replicates for a variance')
    rows = []
    if metadata is not None:
        if 'run_accession' not in metadata:
            raise ValueError('Comparison metadata must identify its runs')
        selected = metadata[metadata.run_accession.isin(a + b)]
        if selected.run_accession.duplicated().any() or set(selected.run_accession) != set(a + b):
            raise ValueError('Comparison metadata must contain exactly one row per selected run')
        rows = selected.astype(object).where(selected.notna(), None).to_dict('records')
        if 'study' in selected and selected.study.dropna().astype(str).str.strip().replace('', np.nan).dropna().nunique() > 1:
            raise ValueError('Compare runs within one study; cross-study batch correction is not provided')
    by_run = {row['run_accession']: row for row in rows}
    samples = {}; unknown = []
    for run in a + b:
        value = by_run.get(run, {}).get('sample_accession')
        sample = str(value).strip() if value is not None else ''
        if not sample or sample.lower() in {'missing', 'unknown', 'not provided', 'not applicable', 'nan', 'n/a'}:
            unknown.append(run)
        else:
            samples.setdefault(sample, []).append(run)
    repeats = {sample: runs for sample, runs in samples.items() if len(runs) > 1}
    if repeats:
        raise ValueError('Selected runs share a submitted sample; repeated sequencing cannot count as independent replicates: '
                         + json.dumps(repeats, sort_keys=True))
    return dict(schema_version=1, group_a=a, group_b=b, selected_metadata=rows,
                known_submitted_samples=samples, runs_without_sample_identity=unknown,
                biological_independence='Not established by run or submitted-sample accessions; review the study design',
                group_roles='User-selected A and B; control pairing and biological replicate roles are not inferred',
                sample_identity_matching='Trim surrounding whitespace for relationship checks; empty/missing/unknown/not provided/not applicable/nan/n/a are unknown; original metadata retained',
                repeated_run_aggregation=False, shared_sample_replicates=False)


def matrix_record(matrix):
    """Hash the ordered numerical matrix actually supplied to the model, with typed labels."""
    values = np.asarray(matrix, dtype='<f8', order='C')
    labels = lambda items: [dict(type=type(x).__module__ + '.' + type(x).__qualname__, value=repr(x)) for x in items]
    description = dict(shape=list(values.shape), index=labels(matrix.index), columns=labels(matrix.columns),
                       original_column_dtypes=[str(dtype) for dtype in matrix.dtypes], index_dtype=str(matrix.index.dtype),
                       representation='IEEE754 float64 little-endian C order; typed ordered labels encoded as UTF8 JSON')
    digest = hashlib.sha256(json.dumps(description, sort_keys=True, ensure_ascii=False).encode())
    digest.update(values.tobytes(order='C'))
    return dict(**description, sha256=digest.hexdigest(), scope='Numerical model input, not original source-file bytes')


def export(table, destination, record):
    """Atomically publish a fresh CSV/TSV and checksum-bound treatment sidecar without overwrites."""
    destination = Path(destination)
    sidecar = destination.with_suffix(destination.suffix + '.json')
    if destination.exists() or sidecar.exists():
        raise FileExistsError('Choose a new comparison export filename')
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix='.' + destination.name, suffix='.partial', dir=destination.parent)
    temporary = Path(name); side_temp = temporary.with_suffix('.json.partial')
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8', newline='') as stream:
            table.to_csv(stream, index=False, sep='\t' if destination.suffix == '.tsv' else ',')
        digest = hashlib.sha256(temporary.read_bytes()).hexdigest()
        side_temp.write_text(json.dumps(dict(record, exported_rows=len(table), export_sha256=digest),
                                      indent=2, ensure_ascii=False, allow_nan=False, default=str) + '\n')
        os.link(temporary, destination)
        try:
            os.link(side_temp, sidecar)
        except Exception:
            destination.unlink()
            raise
        return len(table)
    finally:
        temporary.unlink(missing_ok=True); side_temp.unlink(missing_ok=True)
