"""Check the installed, downloaded expression pack without network or NAS access.

Run with the wheel environment's isolated Python after installing its pinned data
archive. Source identifiers are chosen mechanically to exercise storage and query
routes, without interpreting gene function or treating annotations as replicates.
"""
import hashlib
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
    from parascape import deposited_values as D
    from parascape import reference_identifiers as F
    from parascape import protein_identifier_evidence as I
    import parascape
    assert Path(parascape.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()), 'Package was not imported from the wheel environment'
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
    original_values, errors = D.catalogue(root / 'deposited_values')
    assert len(original_values) == 136 and not errors
    protein_packs = workbench_packs = quantity_columns = cells = missing = 0
    filename_matches = filename_gaps = factor_occurrences = experiment_candidates = 0
    chosen_pack = reference_choice = None
    identifier_root = root / 'processing_audit/protein_identifier_evidence'
    identifier_pointer = json.loads((identifier_root / 'collection.json').read_text())
    identifier_folder = identifier_root / identifier_pointer['directory']
    identifier_index = json.loads((identifier_folder / 'index.json').read_text())
    identifier_sources = {entry['measurement_manifest_sha256']: entry for entry in identifier_index['measurements']}
    assert len(identifier_sources) == 129
    identifier_occurrences = sum(entry['evidence_records'] for entry in identifier_sources.values())
    assert identifier_occurrences == 601605
    identifier_queries = {}
    for pack in original_values:
        folder = Path(pack['path'])
        info = D.manifest(folder)
        feature = pq.read_table(folder / 'values.parquet', columns=['feature_label']).column(0)[0].as_py()
        frame = D.profiles(folder, [feature])
        assert not frame.empty and frame.feature_label.eq(feature).all()
        sample_qc = D.technical_qc(folder)
        assert len(sample_qc) == len(info['sample_columns']) and sample_qc.quality_score.isna().all()
        quantity_columns += len(sample_qc)
        cells += info['features'] * len(sample_qc)
        missing += int(sample_qc.missing.sum())
        if info['assay'] == 'deposited protein-group quantities':
            protein_packs += 1
            digest = hashlib.sha256((folder / 'manifest.json').read_bytes()).hexdigest()
            identifier_source = identifier_sources[digest]
            evidence = pq.read_table(identifier_folder / identifier_source['file'],
                                     filters=[('query_namespace', '=', 'original_protein_group_component')])
            query = evidence.column('query_literal')[0].as_py()
            matches = I.records(folder, [query], namespace='original_protein_group_component')
            assert not matches.empty and matches.query_literal.eq(query).all()
            identifier_queries[str(folder)] = query

            for context in D.column_context(folder):
                association = context['source_design_metadata']
                filename_matches += association.get('primary_file_status') == 'exact_unique_primary_file'
                filename_gaps += association.get('primary_file_status') == 'unresolved_exact_primary_file'
                factor_occurrences += len(association.get('submitted_factor_rows', []))
                experiment_candidates += bool(association.get('submitted_experiment_rows'))
                if association.get('submitted_factor_rows'):
                    chosen_pack = folder, feature
        else:
            workbench_packs += 1
    assert (protein_packs, workbench_packs, quantity_columns, cells, missing) == (129, 7, 1085, 1266353, 38130)
    assert (filename_matches, filename_gaps, factor_occurrences, experiment_candidates) == (11, 12, 5, 357)
    reference_sources, reference_status = F.catalogue(root / 'processing_audit/reference_identifiers')
    acquired = [row for row in reference_sources if row['status'] == 'all_header_cells_reconciled']
    assert len(acquired) >= 47
    for source in acquired:
        table = pq.read_table(root / 'processing_audit/reference_identifiers' / source['pack'] / 'headers.parquet')
        primary = table.column('primary_identifier')[0].as_py()
        assert not F.lookup(source['source_url'], primary).empty
        for literal in table.column('declared_identifier_fields').to_pylist():
            fields = [field for field in json.loads(literal) if field['key'] in ('gene', 'gene_id')]
            if fields:
                value = fields[0]['literal_value']
                assert not F.lookup(source['source_url'], value, kind='gene').empty
                reference_choice = source['source_url'], value
                break
    projects = json.loads((root / 'processing_audit/pride_projects/coverage.json').read_text())
    designs = {}
    for project in projects:
        report = C.pride_samples(root / 'processing_audit', project['accession'])
        assert report['project']['accession'] == project['accession']
        for design in report['designs']:
            designs[design['source_url']] = design
    assert len(projects) == 209 and len(designs) == 52 and sum(d['source_rows'] for d in designs.values()) == 1392
    for accession in ('ST000817', 'ST001304', 'ST002051', 'ST002053', 'ST002544', 'ST002545', 'ST002549'):
        assert C.workbench_samples(root / 'processing_audit', accession)['submitted_rows']
    from PyQt6 import QtWidgets
    from parascape import compare_panel
    from parascape import deposited_values_panel as DP
    from parascape import reference_identifiers_panel as FP
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = QtWidgets.QMainWindow()
    window._workflow_gene = lambda *_args: None
    compare_panel.install(window)
    panels = window.findChildren(DP.DepositedValuesPanel)
    reference_panels = window.findChildren(FP.ReferenceIdentifiersPanel)
    assert len(panels) == len(reference_panels) == 1
    panel = panels[0]
    for index in range(panel.pack.count()):
        if Path(panel.pack.itemData(index)['path']) == chosen_pack[0]:
            panel.pack.setCurrentIndex(index)
            break
    panel.show_column_context()
    assert panel.context_dialog is not None and panel.context_tree.topLevelItemCount() > 0
    panel.features.setPlainText(chosen_pack[1]); panel.show_values()
    output = Path(os.environ['PARASCAPE_SMOKE_OUTPUT'])
    output.mkdir(exist_ok=True)
    assert panel.save(str(output / 'original-measurements.csv')) > 0
    assert panel.save_qc(str(output / 'column-qc.csv')) > 0
    panel.source_filter.setText('Ratio H/L normalized')
    assert panel.pack.count() == 4 and panel.result is None
    panel.source_filter.clear()
    for index in range(panel.pack.count()):
        if Path(panel.pack.itemData(index)['path']) == chosen_pack[0]:
            panel.pack.setCurrentIndex(index)
            break
    panel.identifier_kind.setCurrentIndex(panel.identifier_kind.findData('original_protein_group_component'))
    panel.features.setPlainText(identifier_queries[str(chosen_pack[0])]); panel.show_values()
    assert panel.result is not None and not panel.result.empty
    assert panel.save(str(output / 'qualified-protein-components.csv')) > 0
    exported = pd.read_csv(output / 'qualified-protein-components.csv', dtype=str, keep_default_na=False)
    assert exported.identifier_evidence.tolist() == panel.result.identifier_evidence.tolist()
    panel.features.setPlainText('a missing literal source identifier')
    assert panel.result is None and panel.save(str(output / 'stale-query.csv')) == 0

    reference_panel = reference_panels[0]
    for index in range(reference_panel.source.count()):
        if reference_panel.source.itemData(index)['source_url'] == reference_choice[0]:
            reference_panel.source.setCurrentIndex(index)
            break
    reference_panel.kind.setCurrentIndex(1); reference_panel.query.setText(reference_choice[1]); reference_panel.show_records()
    assert reference_panel.save(str(output / 'original-reference-identifiers.csv')) > 0
    window.close()
    print(json.dumps(dict(verified_RNA_entries=snapshot['verified_runs'], priority_study_runs=8,
                          deposited_query_rows=deposited, catalogue_records=len(rows),
                          original_value_matrices=len(sources), annotation_platforms=len(records),
                          annotated_query_rows=annotated_rows, whole_matrix_QC_rows=len(qc),
                          network_disabled=True, NAS_reads_disabled=True, original_measurement_packs=len(original_values),
                          protein_quantity_packs=protein_packs, workbench_quantity_packs=workbench_packs,
                          original_quantity_cells=cells, missing_original_cells=missing, source_quantity_columns=quantity_columns,
                          exact_primary_filename_matches=filename_matches, submitted_factor_occurrences=factor_occurrences,
                          experiment_label_candidate_columns=experiment_candidates, reference_source_occurrences=len(acquired),
                          reference_header_records=reference_status['unique_content_header_rows'],
                          PRIDE_projects=len(projects), submitted_PRIDE_designs=len(designs),
                          installed_actual_Compare_tabs_and_CSV_exports=True,
                          protein_identifier_packs=len(identifier_sources), protein_identifier_evidence_records=identifier_occurrences,
                          installed_source_search_and_stale_export_refusal=True)))


if __name__ == '__main__':
    main()
