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
