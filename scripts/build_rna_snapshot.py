"""Freeze verified competitive-reference gene matrices for an offline PyPI release.

The full quantifications and alternative reference builds stay on the NAS. Only
audited study/run entries are admitted. A manifest retains input and output hashes,
tool versions, reference lineage, partial-study coverage and numeric precision.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from parascape import expression_contrast as X


def sha256(path):
    """Hash a file without loading it into memory."""
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def build(root: Path, verified: Path, plan: Path, out: Path, compact: bool = False):
    """Publish a complete snapshot only after all selected input/output checks pass."""
    evidence = json.loads(verified.read_text())
    planned = [json.loads(line) for line in plan.read_text().splitlines() if line.strip()]
    by_study = defaultdict(list)
    for row in planned:
        key = row['study'] + '/' + row['run_accession']
        if key in evidence:
            by_study[row['study']].append(row)
    if not by_study or len(evidence) != sum(map(len, by_study.values())):
        raise ValueError('Verified entries do not reconcile to the frozen plan')
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'snapshot.json').exists():
        raise ValueError('Use an empty output directory; snapshots are immutable')
    ref = 'ME49_host'
    mapping = X.tx2gene(ref, str(root))
    ref_dir = root / 'references' / ref
    manifest = {'schema_version': 1, 'created_utc': datetime.now(timezone.utc).isoformat(),
                'scope': 'Partial queue snapshot; competitive ME49 + human/mouse reference',
                'planned_runs': len(planned), 'verified_runs': len(evidence),
                'studies': [], 'quantifiers': list(X.QUANTIFIERS), 'references': [ref],
                'plan_sha256': sha256(plan), 'inputs': {}, 'outputs': {},
                'reference': json.loads((ref_dir / 'reference.json').read_text()),
                'tx2gene_sha256': sha256(ref_dir / 'tx2gene.tsv'),
                'precision': 'Gene count, TPM and effective length stored as float32; full precision inputs retained externally',
                'normalization': 'Sum transcript counts and TPM by exact tx2gene mapping; drop HOST: genes; reclose parasite TPM to 1e6',
                'missingness': 'Missing genes remain NaN; measured zeros remain zero',
                'matrix_fields': ['count', 'tpm'] if compact else ['count', 'tpm', 'length'],
                'effective_lengths': 'Omitted from compact pack; load returns NaN lengths' if compact else 'Included',
                'limitations': 'Descriptive expression only. Partial studies may lack controls. Cross-study TPM is not a batch-corrected effect or differential expression. Alternative references remain external.'}
    gene_ids = set()
    for study, runs in sorted(by_study.items()):
        entry = {'study': study, 'completed_runs': len(runs),
                 'planned_runs': sum(r['study'] == study for r in planned)}
        manifest['studies'].append(entry)
        destination = out / 'matrices' / study
        destination.mkdir(parents=True, exist_ok=True)
        for q in X.QUANTIFIERS:
            counts, tpms, lengths, samples = {}, {}, {}, []
            for row in runs:
                run = row['run_accession']
                key = study + '/' + run
                folder = root / 'rnaseq' / study / run
                receipt = folder / 'done.json'
                before = receipt.read_bytes()
                meta = json.loads(before)
                if meta['study'] != study or meta['run_accession'] != run or ref not in meta['references']:
                    raise ValueError(f'Receipt mismatch: {key}')
                table = folder / f'{q}_{ref}' / ('quant.sf.gz' if q == 'salmon' else 'abundance.tsv.gz')
                expected = evidence[key]['sha256'][str(table)]
                if sha256(table) != expected:
                    raise ValueError(f'Changed quantification: {table}')
                f = X.read_run(str(table.parent), q)
                f['gene'] = f.tx.map(mapping)
                if f.gene.isna().any():
                    raise ValueError(f'Unmapped transcripts: {table}')
                g = f.groupby('gene')
                c, t = g['count'].sum(), g['tpm'].sum()
                length = (f.assign(w=f.tpm * f.length).groupby('gene').w.sum()
                          / t.replace(0, np.nan)).fillna(g['length'].mean())
                host = c.index.str.startswith(X.HOST_PREFIX)
                share = float(c[host].sum() / c.sum()) if c.sum() else None
                c, t, length = c[~host], t[~host], length[~host]
                if t.sum() > 0:
                    t = t / t.sum() * 1e6
                else:
                    t[:] = np.nan
                counts[run], tpms[run], lengths[run] = c, t, length
                samples.append({**meta, 'references': ','.join(meta['references']),
                                'condition': X._condition(meta), 'host_share': share})
                if receipt.read_bytes() != before:
                    raise ValueError(f'Receipt changed during build: {key}')
                manifest['inputs'][f'{key}/{q}'] = {
                    'quantification_sha256': expected,
                    'receipt_sha256': hashlib.sha256(before).hexdigest(),
                    'tool_version': meta[q]}
                gene_ids.update(c.index)
            base = destination / f'{q}_{ref}'
            fields = {'count': pd.DataFrame(counts), 'tpm': pd.DataFrame(tpms)}
            if not compact:
                fields['length'] = pd.DataFrame(lengths)
            pd.concat(fields, axis=1).astype('float32').to_parquet(
                str(base) + '.parquet', compression='zstd', compression_level=19,
                use_dictionary=['gene'], use_byte_stream_split=True)
            pd.DataFrame(samples).to_csv(str(base) + '_samples.tsv', sep='\t', index=False)
        print(f'{study}: {len(runs)} runs built', flush=True)
    manifest['genes'] = len(gene_ids)
    for path in sorted((out / 'matrices').rglob('*')):
        if path.is_file():
            manifest['outputs'][str(path.relative_to(out))] = sha256(path)
    (out / 'snapshot.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Complete: {len(evidence)} runs, {len(by_study)} studies, {len(gene_ids)} genes', flush=True)


def main():
    """Build a release snapshot from explicit, immutable queue evidence."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--verified', type=Path, required=True)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--compact', action='store_true', help='Omit effective length matrices from the compact PyPI pack')
    a = p.parse_args()
    build(a.root, a.verified, a.plan, a.out, a.compact)


if __name__ == '__main__':
    main()
