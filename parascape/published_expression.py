"""Read audited deposited expression packs without mixing annotation variants.

Source RSEM measurements retain their reported units and exact gene identifiers.
Alternative annotations of the same samples never create extra replicates. Packs
are external data, independent of the common-pipeline RNA snapshots and models.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re

import numpy as np
import pandas as pd

from . import paths


def data_root() -> Path:
    """Locate installed deposited expression packs or an explicit local override."""
    return Path(os.environ.get('PARASCAPE_PUBLISHED', str(Path(paths.data_dir()) / 'published_expression')))


def _sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _file(pack, name):
    candidate = (pack / name).resolve()
    if not candidate.is_relative_to(pack.resolve()):
        raise ValueError('Pack file lies outside its data directory')
    return candidate


def manifest(pack: str | Path) -> dict:
    """Verify a pack's validation receipt, output hashes and source sample coverage."""
    pack = Path(pack)
    source = pack / 'manifest.json'
    validation = json.loads((pack / 'independent_validation.json').read_text())
    if validation['manifest_sha256'] != _sha(source):
        raise ValueError('Deposited-expression manifest differs from its validation receipt')
    info = json.loads(source.read_text())
    if info.get('schema_version') != 1 or info.get('software', {}).get('RSEM') is None:
        raise ValueError('Unsupported deposited-expression pack schema or measurement software')
    outputs = info['outputs']
    fields = {'expected_count', 'TPM', 'FPKM', 'effective_length', 'length'}
    for variant in info['variants'].values():
        if (len(variant['samples']) != len(set(variant['samples'])) or
                len(variant['fields']) != len(fields) or set(variant['fields']) != fields):
            raise ValueError('Duplicated samples or unsupported measurement fields')
    required = {'samples.json'} | {ref + '.parquet' for ref in info['variants']}
    if not required.issubset(outputs):
        raise ValueError('Pack does not declare every matrix and sample metadata hash')
    for name, digest in outputs.items():
        if _sha(_file(pack, name)) != digest:
            raise ValueError('Deposited-expression output checksum differs: ' + name)
    samples = json.loads((pack / 'samples.json').read_text())
    ids = [s['geo_accession'] for s in samples['samples']]
    if len(set(ids)) != len(ids) or samples['primary_archive_sha256'] != info['primary_archive_sha256']:
        raise ValueError('Primary sample identity or source archive differs')
    keys = [(r['reference'], r['sample']) for r in info['inputs']]
    expected = {(ref, sample) for ref, v in info['variants'].items() for sample in v['samples']}
    if len(keys) != len(set(keys)) or set(keys) != expected:
        raise ValueError('Variant/sample source identities are missing or duplicated')
    if {sample for _, sample in keys} != set(ids):
        raise ValueError('Matrix samples differ from source sample metadata')
    return info


def catalogue(root: str | Path | None = None) -> tuple[list[dict], list[str]]:
    """List validated source packs and visible errors; unfinished packs stay unlisted."""
    root = Path(root) if root is not None else data_root()
    packs, errors = [], []
    for source in sorted(root.glob('*/manifest.json')):
        if not (source.parent / 'independent_validation.json').exists():
            continue
        try:
            info = manifest(source.parent)
            packs.append(dict(path=str(source.parent), study=info['study'], variants=info['variants'],
                              samples=len({r['sample'] for r in info['inputs']})))
        except (OSError, ValueError, KeyError, TypeError) as error:
            errors.append(f'{source.parent.name}: {error}')
    return packs, errors


def _matrix(pack: str | Path, reference: str):
    """Load the complete checked matrix before selecting genes or reporting QC."""
    pack = Path(pack)
    info = manifest(pack)
    if reference not in info['variants']:
        raise ValueError('Annotation variant is not present in this pack')
    variant = info['variants'][reference]
    matrix = pd.read_parquet(_file(pack, reference + '.parquet'))
    if matrix.index.has_duplicates or matrix.index.isna().any() or matrix.columns.has_duplicates:
        raise ValueError('Duplicate or missing matrix identifiers')
    expected = pd.MultiIndex.from_product([variant['fields'], variant['samples']])
    if matrix.shape[0] != variant['genes'] or set(matrix.columns) != set(expected):
        raise ValueError('Matrix dimensions or sample/field coverage differ from the manifest')
    values = matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Invalid deposited numeric measurements')
    return pack, info, matrix


def technical_qc(pack: str | Path, reference: str) -> pd.DataFrame:
    """Export whole-matrix sample QC with original measurement and source provenance."""
    from . import measurement_qc as Q
    pack, info, matrix = _matrix(pack, reference)
    inputs = {r['sample']: r for r in info['inputs'] if r['reference'] == reference}
    rows = []
    for field in ('TPM', 'FPKM', 'expected_count'):
        report = Q.summarize(matrix[field])
        report['field'] = field
        report['study'] = info['study']
        report['reference'] = reference
        report['normalization'] = info['normalization']
        report['source_url'] = report['sample'].map(lambda sample: inputs[sample]['url'])
        report['source_sha256'] = report['sample'].map(lambda sample: inputs[sample]['sha256'])
        report['manifest_sha256'] = _sha(pack / 'manifest.json')
        report['matrix_sha256'] = info['outputs'][reference + '.parquet']
        rows.append(report)
    return pd.concat(rows, ignore_index=True)


def profiles(pack: str | Path, genes, reference: str) -> pd.DataFrame:
    """Return original source measurements and metadata for exact gene identifiers."""
    if isinstance(genes, str):
        genes = re.split(r'[\s,;]+', genes.strip())
    genes = list(dict.fromkeys(str(g).strip() for g in genes if str(g).strip()))
    if not genes:
        raise ValueError('Enter at least one exact source gene identifier')
    pack, info, matrix = _matrix(pack, reference)
    variant = info['variants'][reference]
    samples = json.loads((pack / 'samples.json').read_text())['samples']
    source_samples = {s['geo_accession']: s for s in samples}
    source_inputs = {r['sample']: r for r in info['inputs'] if r['reference'] == reference}
    rows = []
    for gene in genes:
        if gene not in matrix.index:
            continue
        for sample in variant['samples']:
            original = source_samples[sample]
            row = dict(gene=gene, study=info['study'], sample=sample,
                       sample_title=original['title'], reference=reference,
                       quantifier='RSEM', quantifier_version=info['software']['RSEM'],
                       source_url=source_inputs[sample]['url'],
                       source_sha256=source_inputs[sample]['sha256'],
                       sample_metadata=json.dumps(original, ensure_ascii=False),
                       normalization=info['normalization'], pack=pack.name)
            row.update({field: matrix.loc[gene, (field, sample)] for field in variant['fields']})
            rows.append(row)
    return pd.DataFrame(rows, columns=['gene', 'study', 'sample', 'sample_title', 'reference',
                                     'quantifier', 'quantifier_version', *variant['fields'],
                                     'source_url', 'source_sha256', 'sample_metadata', 'normalization', 'pack'])


def main(argv=None):
    """Export an exact-gene query from one source pack and annotation variant."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('genes', nargs='+')
    parser.add_argument('--pack', type=Path, required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    frame = profiles(args.pack, args.genes, args.reference)
    if frame.empty:
        parser.error('No measured exact identifiers in this annotation variant')
    if args.output:
        frame.to_csv(args.output, index=False)
    else:
        print(frame.to_csv(index=False), end='')


if __name__ == '__main__':
    main()
