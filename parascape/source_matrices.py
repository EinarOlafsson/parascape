"""Read submitted matrices while preserving source identities and scale limits.

Deposited VALUE fields may describe different assays, platforms or normalization
conventions. Readers therefore verify conversion receipts and output hashes,
retain original probes and sample metadata, and expose comparison limits beside
measurements. Optional identifier matching reproduces explicit source annotation
associations without averaging probes or claiming current gene equivalence.
Preserving these distinctions prevents source acquisition from silently becoming
an unsupported gene-level or cross-study comparison.
"""
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


def _context(pack):
    """Verify complete source identities, schema and submitted sample metadata."""
    pack = Path(pack)
    info = manifest(pack)
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
    return info, source, samples


def technical_qc(pack: str | Path) -> pd.DataFrame:
    """Describe every submitted sample using all original probes and an unknown scale.

    Negative values remain valid real measurements because deposited VALUE scales
    may be signed. Neither biological quality nor source units are inferred.
    """
    from . import measurement_qc as Q

    pack = Path(pack)
    info, source, samples = _context(pack)
    matrix = source.read().to_pandas().set_index('probe_id')
    result = Q.summarize(matrix, domain='unspecified_real')
    lookup = {sample['geo_accession']: sample for sample in samples}
    result['geo_sample'] = result['sample']
    result['sample_title'] = result['sample'].map(lambda sample: lookup[sample]['title'])
    result['sample_metadata'] = result['sample'].map(lambda sample: json.dumps(lookup[sample], ensure_ascii=False))
    for key in ('study', 'quantity', 'units', 'catalogue_readout', 'source_url', 'input_sha256',
                'comparison_scope', 'feature_identity'):
        result[key] = info[key]
    result['matrix_sha256'] = info['outputs']['values.parquet']
    result['sample_metadata_sha256'] = info['outputs']['samples.json']
    result['scale_and_platform_metadata_sha256'] = info['outputs']['original_scale_and_platform.json']
    result['manifest_sha256'] = _sha(pack / 'manifest.json')
    result['qc_scope'] = ('All original source probes; no identifier expansion or probe collapse. '
                          'Channels, shared samples and repeated measurements are not independent biological replicates.')
    return result


def profiles(pack: str | Path, probes) -> pd.DataFrame:
    """Return original values for exact, case-sensitive source IDs without gene remapping."""
    pack = Path(pack)
    probes = probes.splitlines() if isinstance(probes, str) else list(probes)
    if not probes or any(not isinstance(p, str) or not p for p in probes):
        raise ValueError('Enter exact source probe or gene identifiers, one per line')
    probes = list(dict.fromkeys(probes))
    info, source, samples = _context(pack)
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


def annotated_profiles(pack: str | Path, identifiers, links_root: str | Path | None = None) -> pd.DataFrame:
    """Match source-declared platform identifiers while retaining every original probe.

    Associations reproduce the source annotation. They do not establish independently
    validated current gene equivalence, and measurements are never averaged by identifier.
    """
    pack = Path(pack)
    info = manifest(pack)
    identifiers = identifiers.splitlines() if isinstance(identifiers, str) else list(identifiers)
    if not identifiers or any(not isinstance(value, str) or not value for value in identifiers):
        raise ValueError('Enter exact source annotation identifiers, one per line')
    metadata = json.loads((pack / 'original_scale_and_platform.json').read_text())
    platforms = set()
    for field in ('!Sample_platform_id', '!Series_platform_id'):
        platforms = {value for row in metadata.get(field, []) for value in row if value}
        if platforms:
            break
    if len(platforms) != 1:
        raise ValueError('Identifier lookup requires one explicitly declared platform for this matrix')
    platform = next(iter(platforms))
    root = Path(links_root) if links_root is not None else Path(os.environ.get(
        'PARASCAPE_PLATFORM_IDENTIFIERS', str(Path(paths.data_dir()) / 'platform_identifiers')))
    coverage = json.loads((root / 'independent_validation_coverage.json').read_text())
    candidates = [r for r in coverage if r['platform'] == platform and r['status'] == 'source_identifier_links_reconciled']
    if len(candidates) != 1:
        raise ValueError('No uniquely reconciled source identifier index for this platform')
    candidate = candidates[0]
    folder = root / Path(candidate['path']).name
    if folder.resolve().parent != root.resolve() or not folder.name.startswith(platform + '-'):
        raise ValueError('Source identifier index lies outside its collection')
    if _sha(folder / 'independent_validation.json') != candidate['independent_validation_sha256']:
        raise ValueError('Source identifier validation receipt checksum differs')
    validation = json.loads((folder / 'independent_validation.json').read_text())
    if _sha(folder / 'manifest.json') != validation['manifest_sha256']:
        raise ValueError('Source identifier manifest differs from its independent reconciliation')
    annotation = json.loads((folder / 'manifest.json').read_text())
    if (annotation['platform'] != platform or validation['verified_links'] != annotation['identifier_link_records']
            or annotation['source_sha256'] != validation['source_sha256']):
        raise ValueError('Source identifier platform or validation coverage differs')
    if _sha(pack / 'manifest.json') not in annotation['observed_matrix_pack_manifests'].values():
        raise ValueError('Identifier index was not prepared for this matrix snapshot')
    source = folder / 'links.parquet'
    if _sha(source) != annotation['outputs']['links.parquet']:
        raise ValueError('Source identifier links checksum differs')
    links = pq.read_table(source, filters=[('identifier', 'in', list(dict.fromkeys(identifiers)))]).to_pandas()
    if links.empty:
        return pd.DataFrame(columns=['query_identifier', 'probe_id', 'geo_sample', 'value', 'sample_title'])
    associations = []
    for (identifier, probe), evidence in links.groupby(['identifier', 'probe_id'], sort=False):
        associations.append(dict(query_identifier=identifier, probe_id=probe,
            annotation_evidence=json.dumps(evidence.to_dict('records'), ensure_ascii=False)))
    association = pd.DataFrame(associations)
    values = profiles(pack, association.probe_id.unique())
    result = association.merge(values, on='probe_id', validate='many_to_many')
    result['identifier_platform'] = platform
    result['identifier_source_url'] = annotation['source_url']
    result['identifier_source_sha256'] = annotation['source_sha256']
    result['identifier_manifest_sha256'] = validation['manifest_sha256']
    result['identifier_interpretation'] = ('Source-declared association only; current gene equivalence not independently validated. '
        'Multiple probes and query aliases retain original measurements and are not extra biological replicates.')
    return result
