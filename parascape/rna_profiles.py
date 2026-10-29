"""Offline gene expression profiles from a versioned snapshot of verified RNA-seq runs.

These are descriptive measurements, not batch-corrected cross-study effects. Compare
conditions within the same study; different protocols and host fractions affect TPM.
Missing genes remain missing, while measured zero TPM remains zero. No expression
snapshot is added to the existing prediction features or their calibration.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import pandas as pd

from . import expression_contrast as X


def snapshot(root: str = X.NAS) -> dict:
    """Return snapshot coverage and provenance, or an explicit external-root status."""
    p = Path(root) / 'snapshot.json'
    return json.loads(p.read_text()) if p.exists() else {'scope': 'external expression matrices'}


def profiles(genes, quantifier: str = 'salmon', reference: str = 'ME49_host',
             root: str = X.NAS, study: str | None = None) -> pd.DataFrame:
    """Return one row per requested gene and measured sample, retaining study identity.

    Accept exact gene accessions, including genes absent from the map's node table.
    An unmeasured gene is not silently replaced by a similarly named gene.
    """
    if isinstance(genes, str):
        genes = re.split(r'[\s,;]+', genes.strip())
    genes = list(dict.fromkeys(str(g).upper() for g in genes if str(g).strip()))
    if not genes:
        raise ValueError('Enter at least one gene accession')
    if quantifier not in X.QUANTIFIERS or reference not in X.REFERENCES:
        raise ValueError('Unknown quantifier or reference')
    available = X.available(root)
    have = available[(available.quantifier == quantifier) & (available.reference == reference)]
    if study:
        have = have[have.study == study]
    rows = []
    for accession in have.study:
        counts, tpm, _length, samples = X.load(accession, quantifier, reference, root)
        found = tpm.index.intersection(genes)
        if not len(found):
            continue
        measured = tpm.loc[found].rename_axis('gene').reset_index().melt(
            id_vars='gene', var_name='run_accession', value_name='tpm')
        c = counts.loc[found].rename_axis('gene').reset_index().melt(
            id_vars='gene', var_name='run_accession', value_name='estimated_count')
        measured = measured.merge(c, on=['gene', 'run_accession'], validate='one_to_one')
        measured = measured.merge(samples.drop(columns=['study'], errors='ignore'),
                                  on='run_accession', validate='many_to_one')
        measured['study'] = accession
        measured['quantifier'], measured['reference'] = quantifier, reference
        rows.append(measured)
    if not rows:
        return pd.DataFrame(columns=['gene', 'study', 'run_accession', 'tpm', 'estimated_count'])
    return pd.concat(rows, ignore_index=True).sort_values(['study', 'run_accession', 'gene'])


def main(argv=None):
    """Export gene profiles from the installed snapshot without a NAS or network."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('genes', nargs='+', help='Exact Toxoplasma gene accessions')
    parser.add_argument('--study')
    parser.add_argument('--quantifier', choices=X.QUANTIFIERS, default='salmon')
    parser.add_argument('--reference', choices=X.REFERENCES, default='ME49_host')
    parser.add_argument('--root', default=X.NAS)
    parser.add_argument('--output', help='CSV destination (default: standard output)')
    args = parser.parse_args(argv)
    frame = profiles(args.genes, args.quantifier, args.reference, args.root, args.study)
    if frame.empty:
        parser.error('No measured values for these genes and conventions in this snapshot')
    if args.output:
        frame.to_csv(args.output, index=False)
    else:
        print(frame.to_csv(index=False), end='')


if __name__ == '__main__':
    main()
