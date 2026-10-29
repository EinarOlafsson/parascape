"""Read original submitted VALUE matrices with exact source identities and explicit scale limits."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from . import paths


def data_root() -> Path:
    """Locate separately installed source-matrix packs or a local collection override."""
    return Path(os.environ.get('PARASCAPE_SOURCE_MATRICES', str(Path(paths.data_dir()) / 'source_matrices')))


def _sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def manifest(pack: str | Path, *, verify_outputs: bool = True) -> dict:
    """Verify a conversion receipt and, by default, every pack output before reading values."""
    pack = Path(pack)
    validation = json.loads((pack / 'validation.json').read_text())
    if _sha(pack / 'manifest.json') != validation['manifest_sha256']:
        raise ValueError('Source-matrix manifest differs from its conversion validation')
    info = json.loads((pack / 'manifest.json').read_text())
    samples = info['sample_order']
    required = {'values.parquet', 'samples.json', 'original_scale_and_platform.json'}
    if (info.get('schema_version') != 1 or set(info['outputs']) != required or not samples
            or len(samples) != len(set(samples)) or info['features'] <= 0):
        raise ValueError('Unsupported source-matrix schema, samples or dimensions')
    if (validation['original_source_sha256'] != info['input_sha256']
            or validation['verified_rows'] != info['features']
            or validation['verified_numeric_or_missing_cells'] != info['features'] * len(samples)):
        raise ValueError('Source-cell validation coverage differs from the matrix dimensions')
    if verify_outputs:
        for name, digest in info['outputs'].items():
            if _sha(pack / name) != digest:
                raise ValueError('Source-matrix output checksum differs: ' + name)
    return info


def catalogue(root: str | Path | None = None) -> tuple[list[dict], list[str]]:
    """List validated manifest metadata; selected value queries verify complete output hashes."""
    root = Path(root) if root is not None else data_root()
    packs, errors = [], []
    for source in sorted(root.glob('*/manifest.json')):
        if source.parent.name.startswith('.'):
            continue
        try:
            info = manifest(source.parent, verify_outputs=False)
            packs.append(dict(path=str(source.parent), study=info['study'],
                              readout=info['catalogue_readout'], features=info['features'],
                              samples=len(info['sample_order']), source_url=info['source_url']))
        except (OSError, ValueError, KeyError, TypeError) as error:
            errors.append(source.parent.name + ': ' + str(error))
    return packs, errors


def profiles(pack: str | Path, probes) -> pd.DataFrame:
    """Return original values for exact, case-sensitive source IDs without gene remapping."""
    pack = Path(pack)
    info = manifest(pack)
    probes = probes.splitlines() if isinstance(probes, str) else list(probes)
    if not probes or any(not isinstance(p, str) or not p for p in probes):
        raise ValueError('Enter exact source probe or gene identifiers, one per line')
    probes = list(dict.fromkeys(probes))
    source = pq.ParquetFile(pack / 'values.parquet')
    if source.metadata.num_rows != info['features'] or source.schema_arrow.names != ['probe_id', *info['sample_order']]:
        raise ValueError('Source-matrix row or sample-column identities differ')
    probe_type = source.schema_arrow.field('probe_id').type
    if (not (pa.types.is_string(probe_type) or pa.types.is_large_string(probe_type)) or
            any(not pa.types.is_float64(source.schema_arrow.field(s).type) for s in info['sample_order'])):
        raise ValueError('Unsupported source identifier or numeric value types')
    metadata = json.loads((pack / 'samples.json').read_text())
    samples = metadata['samples']
    if metadata['study'] != info['study'] or [s['geo_accession'] for s in samples] != info['sample_order']:
        raise ValueError('Source sample metadata differs from matrix columns')
    matrix = pq.read_table(pack / 'values.parquet', filters=[('probe_id', 'in', probes)]).to_pandas()
    if matrix.probe_id.duplicated().any() or matrix.probe_id.isna().any():
        raise ValueError('Ambiguous source feature identity')
    frame = matrix.melt(id_vars='probe_id', var_name='geo_sample', value_name='value')
    if np.isinf(frame['value'].to_numpy(dtype=float)).any():
        raise ValueError('Nonfinite observed source values require review')
    lookup = {s['geo_accession']: s for s in samples}
    frame['sample_title'] = frame.geo_sample.map(lambda s: lookup[s]['title'])
    frame['sample_metadata'] = frame.geo_sample.map(lambda s: json.dumps(lookup[s], ensure_ascii=False))
    for key in ('study', 'quantity', 'units', 'catalogue_readout', 'source_url', 'input_sha256', 'comparison_scope'):
        frame[key] = info[key]
    frame['matrix_sha256'] = info['outputs']['values.parquet']
    frame['scale_and_platform_metadata_sha256'] = info['outputs']['original_scale_and_platform.json']
    frame['feature_identity'] = info['feature_identity']
    return frame
