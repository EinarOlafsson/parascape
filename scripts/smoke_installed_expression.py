"""Check the installed, downloaded expression pack without network or NAS access.

Run with the wheel environment's isolated Python after installing its pinned data
archive. Source identifiers are chosen mechanically to exercise storage and query
routes, without interpreting gene function or treating annotations as replicates.
"""
import json
import os
from pathlib import Path
import sys


def _offline(event, args):
    if event in ('open', 'os.listdir', 'os.scandir') and args:
        path = args[0]
        if isinstance(path, (str, bytes, os.PathLike)) and os.fsdecode(path).startswith('/nas_mnt/'):
            raise PermissionError('Installed-pack smoke test forbids NAS reads')
    if event in ('socket.connect', 'socket.getaddrinfo'):
        raise PermissionError('Installed-pack smoke test forbids network access')


def main():
    """Exercise the installed priority study, deposited variants and primary sample metadata."""
    sys.addaudithook(_offline)
    from parascape import paths
    from parascape import expression_contrast as X
    from parascape import rna_profiles as R
    from parascape import published_expression as P
    from parascape import dataset_coverage_panel as C
    import pandas as pd

    root = Path(paths.data_dir())
    snapshot = R.snapshot(str(root / 'rnaseq'))
    priority = next(s for s in snapshot['studies'] if s['study']=='PRJNA1008424')
    assert priority['completed_runs']==priority['planned_runs']==8
    for tool in X.QUANTIFIERS:
        gene = X.load('PRJNA1008424', tool, 'ME49_host', str(root / 'rnaseq'))[1].index[0]
        values = R.profiles([gene], tool, 'ME49_host', str(root / 'rnaseq'), 'PRJNA1008424')
        assert len(values)==values.run_accession.nunique()==8
    packs, errors = P.catalogue(root / 'published_expression')
    assert packs and not errors
    deposited = 0
    for pack in packs:
        for reference, variant in pack['variants'].items():
            gene = pd.read_parquet(Path(pack['path']) / (reference + '.parquet')).index[0]
            values = P.profiles(pack['path'], [gene], reference)
            assert values['sample'].nunique()==len(values)==len(variant['samples'])
            deposited += len(values)
    rows, status = C.coverage(root / 'processing_audit')
    assert len(rows)==status['catalogue_records']
    assert len(C.source_samples(root / 'processing_audit', 'GSE134099')['samples'])==6
    assert C.shared_samples(root / 'processing_audit')
    print(json.dumps(dict(verified_RNA_entries=snapshot['verified_runs'], priority_study_runs=8,
                          deposited_query_rows=deposited, catalogue_records=len(rows),
                          network_disabled=True, NAS_reads_disabled=True)))


if __name__ == '__main__':
    main()
