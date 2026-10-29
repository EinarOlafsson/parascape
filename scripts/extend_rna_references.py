"""Add all already-verified alternative reference matrices to a frozen RNA snapshot."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from parascape import expression_contrast as X
from scripts.build_rna_snapshot import sha256


def main():
    """Include alternative conventions without admitting newly completed or unverified runs."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--snapshot', type=Path, required=True)
    p.add_argument('--verified', type=Path, required=True)
    a = p.parse_args()
    manifest = json.loads((a.snapshot / 'snapshot.json').read_text())
    evidence = json.loads(a.verified.read_text())
    identities = {tuple(key.split('/')[:2]) for key in manifest['inputs']}
    variants = {}
    for reference in ('ME49', 'GT1', 'VEG'):
        mapping = a.root / 'references' / reference
        variants[reference] = {'reference': json.loads((mapping / 'reference.json').read_text()),
                               'tx2gene_sha256': sha256(mapping / 'tx2gene.tsv')}
        for study in sorted({study for study, run in identities}):
            allowed = []
            for s, run in sorted(identities):
                if s != study:
                    continue
                receipt = a.root / 'rnaseq' / study / run / 'done.json'
                meta = json.loads(receipt.read_text())
                if reference in meta['references']:
                    allowed.append(run)
                    for q in X.QUANTIFIERS:
                        table = receipt.parent / f'{q}_{reference}' / ('quant.sf.gz' if q == 'salmon' else 'abundance.tsv.gz')
                        expected = evidence[f'{study}/{run}']['sha256'][str(table)]
                        if sha256(table) != expected:
                            raise ValueError(f'Changed verified input: {table}')
                        manifest['inputs'][f'{study}/{run}/{reference}/{q}'] = {
                            'quantification_sha256': expected, 'receipt_sha256': sha256(receipt), 'tool_version': meta[q]}
            if not allowed:
                continue
            for q in X.QUANTIFIERS:
                c, t, length, samples = X.matrices(study, q, reference, str(a.root))
                samples = samples[samples.run_accession.isin(allowed)]
                base = a.snapshot / 'matrices' / study / f'{q}_{reference}'
                values = pd.concat({'count': c[allowed], 'tpm': t[allowed], 'length': length[allowed]}, axis=1).astype('float32')
                values.to_parquet(str(base) + '.parquet', compression='zstd', compression_level=19)
                samples.to_csv(str(base) + '_samples.tsv', sep='\t', index=False)
            print(f'{reference} / {study}: {len(allowed)} verified runs', flush=True)
    manifest['references'] = list(X.REFERENCES)
    manifest['reference_variants'] = variants
    manifest['scope'] = f"Partial queue snapshot; all already-verified quantifier/reference conventions for its frozen {manifest['verified_runs']} study/run entries"
    manifest['limitations'] = 'Descriptive expression. Partial studies may lack controls. Cross-study TPM is not batch-corrected. Alternative strain references retain their own exact gene identifiers.'
    for path in (a.snapshot / 'matrices').rglob('*'):
        if path.is_file():
            manifest['outputs'][str(path.relative_to(a.snapshot))] = sha256(path)
    (a.snapshot / 'snapshot.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('All verified alternative references added', flush=True)


if __name__ == '__main__':
    main()
