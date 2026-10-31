"""Descriptive seed-to-feature coexpression within one source and quantity.

Rows remain source features, including protein groups and transcripts. Correlation
does not establish interaction, regulation, gene equivalence or replicate identity.

Why it works this way: each search stays inside one deposited table, because values from
different studies, units or quantity families are not on a common scale, and mixing them
would make study effects look like coexpression. Only columns where both features are
finite are paired, so missing measurements never become zeros, while measured zeros are
kept. Runs that share one submitted sample are refused, because repeated sequencing of a
sample would otherwise count as extra evidence.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from . import deposited_values as D
from . import expression_contrast as X
from .comparison_provenance import matrix_record


RNA_ASSAYS = {
    'deposited original gene-table quantities',
    'deposited original transcript-table quantities',
    'deposited original target-table quantities',
}


@dataclass
class ExpressionMatrix:
    """One source table as stored: original row identities, exact labels, sample context and lineage."""

    values: pd.DataFrame  # original row identity x original sample column
    labels: pd.Series
    samples: list[dict]
    record: dict


def catalogue(rna_root=None, deposited_root=None):
    """List source-specific RNA and protein quantities, without pooling conventions."""
    sources = []
    for row in X.available(X.NAS if rna_root is None else str(rna_root)).to_dict('records'):
        sources.append(dict(kind='rna', level='RNA', **row,
            title=f"{row['study']} · {row['quantifier']} / {row['reference']} · TPM"))
    packs, errors = D.catalogue(deposited_root)
    for row in packs:
        if row['assay'] in RNA_ASSAYS:
            level = 'RNA'
        elif 'protein' in row['assay']:
            level = 'Protein'
        else:
            continue
        info = D.manifest(row['path'], verify_outputs=False)
        # Ratios / isoform proportions answer different questions from abundance.
        if any(term in info['quantity'].lower() for term in ('ratio', 'isopct')):
            continue
        sources.append(dict(kind='deposited', level=level, **row,
            title=f"{row['study']} · {row['analysis']} · {info['quantity']} ({row['samples']} columns)"))
    return sources, errors


def load(source, rna_root=None):
    """Load checked original quantities and pin source bytes and sample metadata."""
    if source['kind'] == 'rna':
        root = str(Path(X.NAS if rna_root is None else rna_root).resolve())
        base = Path(X.matrix_dir(root)) / source['study'] / f"{source['quantifier']}_{source['reference']}"
        files = {str(path): D._sha(path) for path in
                 (base.with_suffix('.parquet'), Path(str(base) + '_samples.tsv'))}
        _, values, _, samples = X.load(source['study'], source['quantifier'], source['reference'], root)
        if any(D._sha(Path(path)) != digest for path, digest in files.items()):
            raise ValueError('RNA source files changed during loading')
        if 'run_accession' not in samples or samples.run_accession.duplicated().any():
            raise ValueError('RNA sample metadata must identify each run exactly once')
        if set(samples.run_accession) != set(values.columns):
            raise ValueError('RNA sample metadata and measurement columns differ')
        metadata = samples.astype(object).where(samples.notna(), None).to_dict('records')
        by_run = {row['run_accession']: row for row in metadata}
        metadata = [dict(by_run[column], column=column) for column in values.columns]
        labels = pd.Series(values.index.astype(str), index=values.index)
        record = dict(source=source, quantity='TPM', units='TPM', source_files_sha256=files,
            source_scope='One study, quantifier and reference; gene-level matrix supplied by RNA pipeline')
    elif source['kind'] == 'deposited':
        pack = D.runtime_pack(source['path'])
        info, parquet, metadata = D._context(pack)
        frame = parquet.read().to_pandas()
        if frame.source_row.tolist() != list(range(info['features'])):
            raise ValueError('Deposited source-row identities differ')
        values = frame.set_index('source_row')[info['sample_columns']]
        labels = pd.Series(frame.feature_label.to_numpy(), index=values.index)
        if any(D._sha(pack / name) != info['outputs'][name] for name in ('values.parquet', 'samples.json')):
            raise ValueError('Deposited source files changed during loading')
        record = dict(source=source, quantity=info['quantity'], units=info['units'],
            manifest=info, manifest_sha256=D._sha(pack / 'manifest.json'),
            source_scope='Original source rows; protein groups, transcripts and targets are not collapsed to genes')
    else:
        raise ValueError('Unsupported coexpression source')
    if (not values.index.is_unique or not values.columns.is_unique
            or not labels.index.equals(values.index)):
        raise ValueError('Coexpression requires distinct source row and column identities')
    return ExpressionMatrix(values, labels, metadata, record)


def correlate(matrix, seed, columns, *, seed_row=None, method='spearman',
              transform='none', min_pairs=6, cancelled=lambda: False):
    """Compare a seed with every source row using pairwise finite observations.

    Spearman ranks are recomputed on each pair's overlap. No imputation, zero
    replacement, cross-study pooling, p-values or significance claims are made.
    Duplicate feature labels require an explicit original row identity.
    """
    columns = list(columns)
    if method not in {'spearman', 'pearson'} or transform not in {'none', 'log2p1'}:
        raise ValueError('Choose Spearman/Pearson and an explicit supported transformation')
    if not isinstance(min_pairs, int) or isinstance(min_pairs, bool) or min_pairs < 4:
        raise ValueError('At least four paired columns are required for exploratory correlation')
    if len(set(columns)) != len(columns) or not set(columns) <= set(matrix.values.columns):
        raise ValueError('Choose distinct columns from the selected source')
    if len(columns) < min_pairs:
        raise ValueError(f'Select at least {min_pairs} columns; review their biological relationships')
    matches = matrix.labels[matrix.labels.eq(seed)]
    if seed_row is not None:
        matches = matches.loc[matches.index == seed_row]
    if len(matches) != 1:
        raise ValueError(f'Exact seed matches {len(matches)} rows. Specify its original row identity: '
                         + ', '.join(map(str, matches.index[:20])))
    key = matches.index[0]
    selected = matrix.values.loc[:, columns]
    metadata = [row for row in matrix.samples if row['column'] in columns]
    known_samples = {}
    for row in metadata:
        sample = row.get('sample_accession')
        if sample is not None and str(sample).strip().lower() not in {
                '', 'nan', 'n/a', 'unknown', 'missing', 'not provided', 'not applicable'}:
            known_samples.setdefault(str(sample).strip(), []).append(row['column'])
    if any(len(runs) > 1 for runs in known_samples.values()):
        raise ValueError('Selected runs share a submitted sample; repeated sequencing must not count as extra samples')
    data = selected.to_numpy(dtype=float, copy=True)
    if transform == 'log2p1':
        if (data[np.isfinite(data)] < 0).any():
            raise ValueError('log2(x+1) requires nonnegative finite source values; use the original scale')
        data = np.log2(data + 1)
    seed_position = matrix.values.index.get_loc(key)
    x = data[seed_position]
    rows = []
    for i, (row_key, label) in enumerate(matrix.labels.items()):
        if cancelled():
            raise InterruptedError('Coexpression cancelled')
        if i == seed_position:
            continue
        y = data[i]
        keep = np.isfinite(x) & np.isfinite(y)
        n = int(keep.sum())
        coefficient, status = np.nan, 'insufficient_pairs'
        if n >= min_pairs:
            a, b = x[keep], y[keep]
            if method == 'spearman':
                a, b = rankdata(a), rankdata(b)
            # Scale before centering to avoid overflow with extreme original units.
            a = a / (float(np.max(np.abs(a))) or 1.)
            b = b / (float(np.max(np.abs(b))) or 1.)
            a, b = a - a.mean(), b - b.mean()
            denominator = np.linalg.norm(a) * np.linalg.norm(b)
            if denominator == 0:
                status = 'constant_on_overlap'
            else:
                coefficient = float(np.clip(np.dot(a, b) / denominator, -1., 1.))
                status = 'descriptive'
        rows.append(dict(source_row=row_key, feature_label=label, correlation=coefficient,
                         paired_columns=n, excluded_columns=len(columns) - n,
                         seed_nonzero_pairs=int(np.count_nonzero(x[keep])),
                         feature_nonzero_pairs=int(np.count_nonzero(y[keep])), status=status))
    table = pd.DataFrame(rows, columns=['source_row', 'feature_label', 'correlation',
                                      'paired_columns', 'excluded_columns', 'seed_nonzero_pairs',
                                      'feature_nonzero_pairs', 'status'])
    if len(table):
        table = table.sort_values('correlation', key=lambda s: s.abs(), ascending=False,
                                  kind='stable', na_position='last').reset_index(drop=True)
    record = dict(schema_version=1, analysis='within-source descriptive seed coexpression',
        source_record=matrix.record, selected_columns=columns, selected_sample_metadata=metadata,
        input_matrix=matrix_record(selected), method=method, transform=transform,
        seed_label=seed, seed_row=key.item() if isinstance(key, np.generic) else key,
        minimum_paired_columns=min_pairs,
        operations=['Select explicit source columns', 'Preserve original row identities and measured zeros',
                    'Apply ' + transform, 'Exclude nonfinite observations separately for each pair',
                    'Recompute ranks on pairwise overlap' if method == 'spearman' else 'Center each overlapping pair',
                    'Calculate correlation; retain ineligible rows with reason'],
        missing_policy='Pairwise finite observations; no imputation or zero replacement',
        cross_study_pooling=False, gene_aggregation=False, biological_quality_score=None,
        biological_independence='Not established; review submitted sample identity and study design',
        p_values_computed=False, interpretation='Descriptive association; not physical interaction, regulation or causality',
        implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    return table, record
