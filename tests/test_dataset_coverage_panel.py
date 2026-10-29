"""Coverage and sample browsing preserve identities and expose incomplete measurement admission."""
import csv
import hashlib
import json

import pytest

from parascape import dataset_coverage_panel as D


@pytest.fixture
def audit(tmp_path):
    """Two synthetic source studies share one explicitly multi-channel submitted sample."""
    rows = []
    for accession, title, state in [('GSE1', 'original one', 'local runtime measurements present'),
                                    ('GSE2', 'original two', 'not integrated into local comparison runtime')]:
        rows.append(dict(accession=accession, source='GEO', title=title, readout='synthetic assay',
                         supplementary_urls_in_plan='2', source_files_with_acquisition_receipts='1',
                         measurement_status=state, primary_sample_metadata='validated_sample_metadata'))
    path = tmp_path / 'dataset_coverage.tsv'
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    (tmp_path / 'dataset_coverage_status.json').write_text(json.dumps(dict(
        coverage_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), catalogue_records=2,
        local_runtime_catalogue_records=1)))
    folder = tmp_path / 'all_geo_samples'
    (folder / 'GSE1').mkdir(parents=True)
    source = dict(accession='GSE1', source_samples=1, samples=[dict(
        geo_accession='GSM1', title='Original sample label', relations=[dict(type='BioSample', target='SAMN1')],
        channels=[dict(position='1', source='Original channel', organisms=['synthetic organism'],
                       characteristics=[dict(tag='replicate', value='source label 1')]),
                  dict(position='2', source='Second source channel', organisms=[], characteristics=[])])])
    source_path = folder / 'GSE1/geo-samples.json';source_path.write_text(json.dumps(source))
    (folder / 'coverage.json').write_text(json.dumps([dict(accession='GSE1', status='validated_sample_metadata',
        sample_metadata_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest())]))
    links = folder / 'sample_identity_links.tsv'
    links.write_text('geo_sample\tsame_geo_sample_in_studies\nGSM1\tGSE1;GSE2\n')
    (folder / 'sample_identity_manifest.json').write_text(json.dumps(dict(outputs={links.name:hashlib.sha256(links.read_bytes()).hexdigest()})))
    return tmp_path


def test_full_catalogue_and_local_admission_remain_separate(audit):
    rows, status = D.coverage(audit)
    assert len(rows)==2 and status['local_runtime_catalogue_records']==1
    assert rows[1]['measurement_status'].startswith('not integrated')
    assert D.shared_samples(audit)['GSM1']==['GSE1','GSE2']
    assert len(D.source_samples(audit, 'GSE1')['samples'])==1


def test_changed_coverage_and_sample_labels_are_rejected(audit):
    path = audit / 'dataset_coverage.tsv';path.write_bytes(path.read_bytes()+b'\n')
    with pytest.raises(ValueError, match='checksum'):
        D.coverage(audit)
    path = audit / 'all_geo_samples/GSE1/geo-samples.json';path.write_text('{}')
    with pytest.raises(ValueError, match='checksum'):
        D.source_samples(audit,'GSE1')


def test_browse_original_labels_shared_identity_and_channels(qtbot,audit):
    panel = D.DatasetCoveragePanel(audit);qtbot.addWidget(panel)
    assert panel.table.rowCount()==2 and '1 with local measurements' in panel.summary.text()
    panel.search.setText('original one');assert panel.table.rowCount()==1
    panel.table.selectRow(0)
    assert panel.samples.topLevelItemCount()==1
    sample = panel.samples.topLevelItem(0)
    assert sample.text(0)=='GSM1' and sample.text(1)=='Original sample label'
    assert sample.child(0).text(1)=='2'
    assert sample.child(1).text(1)=='GSE1; GSE2'
    assert 'not inferred' in panel.details.text()
    panel.search.setText('no-such-study');assert panel.table.rowCount()==0
    assert panel.samples.topLevelItemCount()==0


def test_missing_optional_audit_pack_is_visible(qtbot,tmp_path):
    panel = D.DatasetCoveragePanel(tmp_path);qtbot.addWidget(panel)
    assert panel.table.rowCount()==0 and 'not installed' in panel.summary.text()


def test_series_matrix_acquisition_does_not_claim_measurement_admission(qtbot, audit):
    rows, _ = D.coverage(audit)
    for row in rows:
        row.update(series_matrix_files_declared='2', series_matrix_files_acquired='1',
                   series_matrices_numeric_checked='1', series_matrices_without_embedded_values='0',
                   series_matrices_review_required='0')
    panel = D.DatasetCoveragePanel(audit)
    qtbot.addWidget(panel)
    panel.rows = rows
    panel._filter()
    assert '1/2 Series Matrices' in panel.table.item(1, 4).text()
    assert panel.table.item(1, 5).text() == 'Pending measurements'
    panel.table.selectRow(1)
    assert 'Original units, scale' in panel.details.text()


@pytest.fixture
def workbench_audit(audit):
    """Synthetic factor rows retain repeated labels and distinct analysis representations."""
    root = audit / 'metabolomics_workbench'
    folder = root / ('1' * 64)
    study = folder / 'ST000001'
    study.mkdir(parents=True)
    design = dict(study='ST000001', submitted_rows=[
        dict(local_sample_id='001', sample_source='Original source', factors={'batch': 'source A'}, raw_data='source.raw'),
        dict(local_sample_id='001', sample_source='Original source', factors={'batch': 'source B'}, raw_data='source-repeat.raw')])
    analyses = dict(study='ST000001', analyses=[dict(analysis='AN000001', source_units=['source units'],
        mwtab_deposited_rows=3, study_data_endpoint_rows=2, source_representations_have_equal_row_counts=False)])
    (study / 'sample_design.json').write_text(json.dumps(design))
    (study / 'analyses.json').write_text(json.dumps(analyses))
    (folder / 'coverage.json').write_text(json.dumps([dict(accession='ST000001',
        status='independently_reconciled_primary_metadata', submitted_factor_records=2, analyses=1)]))
    outputs = {p.relative_to(folder).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in folder.rglob('*.json')}
    (folder / 'validation.json').write_text(json.dumps(dict(outputs=outputs)))
    (root / 'collection.json').write_text(json.dumps(dict(schema_version=1, directory=folder.name,
        validation_sha256=hashlib.sha256((folder / 'validation.json').read_bytes()).hexdigest())))
    return audit


def test_workbench_metadata_retains_repeated_labels_units_and_source_discrepancies(workbench_audit):
    report = D.workbench_samples(workbench_audit, 'ST000001')
    assert [row['local_sample_id'] for row in report['submitted_rows']] == ['001', '001']
    assert report['submitted_rows'][0]['factors'] == {'batch': 'source A'}
    assert report['analysis_inventory'][0]['source_units'] == ['source units']
    assert not report['analysis_inventory'][0]['source_representations_have_equal_row_counts']
    pointer = json.loads((workbench_audit / 'metabolomics_workbench/collection.json').read_text())
    path = workbench_audit / 'metabolomics_workbench' / pointer['directory'] / 'ST000001/sample_design.json'
    path.write_text('{}')
    with pytest.raises(ValueError, match='checksum'):
        D.workbench_samples(workbench_audit, 'ST000001')


def test_workbench_browser_shows_factors_without_admitting_measurements(qtbot, workbench_audit):
    panel = D.DatasetCoveragePanel(workbench_audit)
    qtbot.addWidget(panel)
    panel.rows = [dict(accession='ST000001', source='Metabolomics Workbench', title='Original source study',
        readout='synthetic assay', supplementary_urls_in_plan='0', source_files_with_acquisition_receipts='0',
        workbench_acquired_api_responses='12', measurement_status='not integrated into local comparison runtime',
        primary_sample_metadata='independently_reconciled_primary_metadata')]
    panel._filter()
    panel.table.selectRow(0)
    assert '12 API responses' in panel.table.item(0, 4).text()
    assert panel.table.item(0, 5).text() == 'Pending measurements'
    assert panel.samples.topLevelItemCount() == 3
    assert panel.samples.topLevelItem(0).text(0) == panel.samples.topLevelItem(1).text(0) == '001'
    assert 'not inferred' in panel.details.text() and 'raw filenames are declarations' in panel.details.text()


def test_workbench_collection_path_and_receipt_must_be_valid(workbench_audit):
    path = workbench_audit / 'metabolomics_workbench/collection.json'
    pointer = json.loads(path.read_text())
    pointer['validation_sha256'] = '0' * 64
    path.write_text(json.dumps(pointer))
    with pytest.raises(ValueError, match='receipt checksum'):
        D.workbench_samples(workbench_audit, 'ST000001')
    pointer['directory'] = '../outside'
    path.write_text(json.dumps(pointer))
    with pytest.raises(ValueError, match='collection identity'):
        D.workbench_samples(workbench_audit, 'ST000001')


@pytest.fixture
def pride_audit(audit):
    """Submitted design rows repeat identities and column headings without creating new replicates."""
    projects = audit / 'pride_projects'
    records = []
    for accession in ['PXD000001', 'PXD000002']:
        folder = projects / accession
        folder.mkdir(parents=True)
        summary = folder / 'project_summary.json'
        summary.write_text(json.dumps(dict(accession=accession, title='Original project title',
            sampleAttributes=[{'value': 'Submitted source label'}])))
        records.append(dict(accession=accession, status='validated_primary_metadata',
            summary_sha256=hashlib.sha256(summary.read_bytes()).hexdigest(), files_sha256='a'*64))
    coverage = projects / 'coverage.json'
    coverage.write_text(json.dumps(records))
    collection = audit / 'pride_sample_designs'
    folder = collection / ('1'*64)
    folder.mkdir(parents=True)
    source_hash = 'b'*64
    rows = [dict(source_row=position, cells=['001', batch, '', 'original.raw'], submitted_source_names=['001'],
        source_row_id=hashlib.sha256((source_hash+':'+str(position)).encode()).hexdigest(),
        file_links=[dict(column_index=3,submitted_file='original.raw',project='PXD000001',
            exact_primary_file_ids=['file-001'], status='exact_unique_file')])
        for position,batch in [(2,'not available'),(3,'second submitted label')]]
    name = 'c'*64+'.json'
    (folder/name).write_text(json.dumps(dict(source_url='https://example.org/original.sdrf.tsv', source_sha256=source_hash,
        headers=['source name','characteristics[batch]','characteristics[batch]','comment[data file]'],
        rows=rows,source_rows=2,status='structurally_validated_design',
        primary_file_metadata_sha256={'PXD000001':'a'*64},metadata_origin='community_annotated',
        biological_truth_status='Submitted metadata; not independently verified biological truth')))
    (folder/'index.json').write_text(json.dumps(dict(schema_version=1, projects={'PXD000001':[name],'PXD000002':[]})))
    outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.iterdir()}
    receipt=folder/'validation.json'
    receipt.write_text(json.dumps(dict(status='source_designs_reconciled',catalogue_projects=2,
        bindings={'project_coverage_sha256':hashlib.sha256(coverage.read_bytes()).hexdigest()},outputs=outputs)))
    (collection/'collection.json').write_text(json.dumps(dict(schema_version=1,directory=folder.name,
        validation_sha256=hashlib.sha256(receipt.read_bytes()).hexdigest())))
    return audit


def test_pride_submitted_rows_duplicate_headers_and_exact_file_links_are_preserved(pride_audit):
    report=D.pride_samples(pride_audit,'PXD000001')
    design=report['designs'][0]
    assert design['headers'][1:3]==['characteristics[batch]','characteristics[batch]']
    assert [row['submitted_source_names'] for row in design['rows']]==[['001'],['001']]
    assert design['rows'][0]['cells'][1:3]==['not available','']
    assert design['rows'][0]['file_links'][0]['exact_primary_file_ids']==['file-001']
    assert 'not inferred' in report['design_status']
    empty=D.pride_samples(pride_audit,'PXD000002')
    assert empty['project']['accession']=='PXD000002' and not empty['designs']
    assert 'No submitted SDRF' in empty['design_status']


def test_pride_receipts_outputs_and_project_context_are_required(pride_audit):
    pointer_path=pride_audit/'pride_sample_designs/collection.json'
    pointer=json.loads(pointer_path.read_text())
    folder=pointer_path.parent/pointer['directory']
    design=folder/('c'*64+'.json')
    original=design.read_bytes()
    design.write_text('{}')
    with pytest.raises(ValueError,match='checksum'):
        D.pride_samples(pride_audit,'PXD000001')
    design.write_bytes(original)
    project=pride_audit/'pride_projects/PXD000001/project_summary.json'
    project.write_text('{}')
    with pytest.raises(ValueError,match='summary checksum'):
        D.pride_samples(pride_audit,'PXD000001')
    pointer['directory']='../outside'
    pointer_path.write_text(json.dumps(pointer))
    with pytest.raises(ValueError,match='collection identity'):
        D.pride_samples(pride_audit,'PXD000002')


def test_pride_browser_keeps_submitted_columns_and_repeated_source_rows(qtbot,pride_audit):
    panel=D.DatasetCoveragePanel(pride_audit)
    qtbot.addWidget(panel)
    panel.rows=[dict(accession='PXD000001',source='PRIDE',title='Original project',readout='synthetic assay',
        supplementary_urls_in_plan='2',source_files_with_acquisition_receipts='1',
        measurement_status='not integrated into local comparison runtime',primary_sample_metadata='project metadata only')]
    panel._filter()
    panel.table.selectRow(0)
    assert panel.samples.topLevelItemCount()==2
    design=panel.samples.topLevelItem(1)
    assert design.childCount()==4
    first,second=design.child(2),design.child(3)
    assert first.text(1)==second.text(1)=='001'
    assert first.child(1).text(0)=='2: characteristics[batch]'
    assert first.child(1).toolTip(0)=='2: characteristics[batch]'
    assert first.child(2).text(0)=='3: characteristics[batch]'
    assert first.child(1).text(1)=='not available' and first.child(2).text(1)==''
    assert first.child(4).child(0).text(1)=='exact_unique_file'
    assert 'Primary project/file metadata: verified' in panel.details.text()
    assert 'not inferred' in panel.details.text() and 'do not prove downloaded' in panel.details.text()
    assert panel.table.item(0,5).text()=='Pending measurements'
