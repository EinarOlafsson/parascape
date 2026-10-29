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
    from parascape import source_matrices as S
    import pandas as pd
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq

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
    sources, errors = S.catalogue(root / 'source_matrices')
    assert len(sources) >= 147 and not errors
    for source in sources:
        path = Path(source['path'])
        probe = pq.ParquetFile(path / 'values.parquet').read_row_group(0, columns=['probe_id']).column(0)[0].as_py()
        values = S.profiles(path, [probe])
        assert len(values) == source['samples'] and values.geo_sample.nunique() == source['samples']
    identifiers = root / 'platform_identifiers'
    records = json.loads((identifiers / 'independent_validation_coverage.json').read_text())
    assert sum(r['status']=='source_identifier_links_reconciled' for r in records) >= 54
    annotated_rows = 0
    for record in records:
        if record['status'] != 'source_identifier_links_reconciled':
            continue
        folder = identifiers / Path(record['path']).name
        annotation = json.loads((folder / 'manifest.json').read_text())
        table = pq.ParquetFile(folder / 'links.parquet')
        if not table.metadata.num_rows:
            continue
        query_candidates = table.read_row_group(0, columns=['identifier']).column(0).to_pylist()
        dataset = ds.dataset(folder / 'links.parquet')
        for identifier in list(dict.fromkeys(query_candidates))[:10]:
            if dataset.count_rows(filter=ds.field('identifier') == identifier) > 4:
                continue
            probes = pq.read_table(folder / 'links.parquet', columns=['probe_id'],
                                  filters=[('identifier', '=', identifier)]).column(0).to_pylist()
            for name in annotation['observed_matrix_pack_manifests']:
                pack = root / 'source_matrices' / name
                info = S.manifest(pack, verify_outputs=False)
                if len(info['sample_order']) > 200:
                    continue
                present = pq.read_table(pack / 'values.parquet', columns=['probe_id'],
                                        filters=[('probe_id', 'in', probes)])
                if not present.num_rows:
                    continue
                frame = S.annotated_profiles(pack, [identifier], identifiers)
                assert len(frame) == present.num_rows * len(info['sample_order'])
                annotated_rows = len(frame)
                break
            if annotated_rows:
                break
        if annotated_rows:
            break
    assert annotated_rows > 0
    qc = R.technical_qc(root=str(root / 'rnaseq'), study='PRJNA1008424')
    assert len(qc) == 16 and qc.quality_score.isna().all()
    print(json.dumps(dict(verified_RNA_entries=snapshot['verified_runs'], priority_study_runs=8,
                          deposited_query_rows=deposited, catalogue_records=len(rows),
                          original_value_matrices=len(sources), annotation_platforms=len(records),
                          annotated_query_rows=annotated_rows, whole_matrix_QC_rows=len(qc),
                          network_disabled=True, NAS_reads_disabled=True)))


if __name__ == '__main__':
    main()
