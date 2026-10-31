"""QC reports bind to exact source artifacts and keep unknown duplication unknown."""
import hashlib
import json

import pytest

from parascape import sparse_record_qc as Q
from parascape.sparse_values_panel import SparseValuesPanel
from test_sparse_values import fixture_collection


def write_report(root, descriptor, *, unordered=False, mutation=None):
    info=json.loads(descriptor.read_text());reference=info['coordinate_artifact'];processor='f'*64
    facts=dict(original_records=4,source_row_then_column_nondecreasing=not unordered,
               source_column_then_row_nondecreasing=False,adjacent_repeated_records=1,adjacent_duplicate_runs=1,
               global_duplicate_count_proven=not unordered,global_repeated_records=None if unordered else 1,
               global_repeated_coordinate_groups=None if unordered else 1,unique_stored_coordinates=None if unordered else 3,
               duplicate_record_fraction=None if unordered else .25,coordinate_aggregation_applied=False,biological_quality_score=None)
    record=dict(bindings={**{key:reference[key] for key in ('directory','manifest_sha256','validation_sha256','output_sha256')},'processor_sha256':processor},
                source_url=info['source_url'],source_input_sha256=info['input_sha256'],source_record=info['source_record'],
                original_values_changed=False,biological_or_comparative_admission=False,result=facts)
    if mutation:mutation(record)
    root.mkdir(exist_ok=True);filename='a'*64+'.json';body=json.dumps(record,sort_keys=True)+'\n';(root/filename).write_text(body)
    collection=dict(schema_version=1,processor_sha256=processor,packs=[dict(source_url=info['source_url'],file=filename,
                    sha256=hashlib.sha256(body.encode()).hexdigest(),global_duplicate_count_proven=facts['global_duplicate_count_proven'],global_repeated_records=facts['global_repeated_records'])])
    body=json.dumps(collection,sort_keys=True)+'\n';digest=hashlib.sha256(body.encode()).hexdigest();(root/(digest+'.json')).write_text(body);(root/'current.json').write_text(json.dumps(dict(file=digest+'.json',sha256=digest)))


def test_source_bound_complete_and_unknown_facts(tmp_path):
    roots,descriptor=fixture_collection(tmp_path);root=tmp_path/'qc';write_report(root,descriptor)
    report=Q.report(descriptor,root);assert report['facts']['global_repeated_records']==1 and report['facts']['duplicate_record_fraction']==.25
    write_report(root,descriptor,unordered=True);report=Q.report(descriptor,root);assert report['facts']['global_repeated_records'] is None
    assert report['facts']['biological_quality_score'] is None
    assert Q.report(descriptor,tmp_path/'missing') is None


@pytest.mark.parametrize('change',['output','source','fraction','unknown_as_zero','biology'])
def test_rehashed_incompatible_reports_are_refused(tmp_path,change):
    roots,descriptor=fixture_collection(tmp_path);root=tmp_path/'qc'
    def mutation(record):
        if change=='output':record['bindings']['output_sha256']='0'*64
        elif change=='source':record['source_record']['studies']='different'
        elif change=='fraction':record['result']['duplicate_record_fraction']=.5
        elif change=='unknown_as_zero':record['result']['global_repeated_records']=0
        else:record['result']['biological_quality_score']=1.
    write_report(root,descriptor,unordered=change=='unknown_as_zero',mutation=mutation)
    with pytest.raises(ValueError):Q.report(descriptor,root)


def test_qt_context_and_export_keep_qc_source_and_scope(qtbot,tmp_path,monkeypatch):
    roots,descriptor=fixture_collection(tmp_path);audit=tmp_path/'audit';audit.mkdir();write_report(audit/'sparse_coordinate_qc',descriptor)
    monkeypatch.setenv('PARASCAPE_AUDIT',str(audit));panel=SparseValuesPanel(roots['sparse_datasets'],coordinate_root=roots['sparse_coordinates'],axis_root=roots['sparse_axes']);qtbot.addWidget(panel)
    panel.show_qc();assert panel.qc_dialog.isVisible()
    panel.features.setPlainText('001');assert panel.qc_dialog is None
    panel.show_values();qtbot.waitUntil(lambda:panel._worker is None,timeout=10000);assert panel.result is not None
    target=tmp_path/'query.csv';assert panel.save(str(target))==4
    sidecar=json.loads(target.with_suffix('.csv.json').read_text());assert sidecar['coordinate_QC']['facts']['global_repeated_records']==1
    assert sidecar['coordinate_QC']['facts']['biological_quality_score'] is None


def test_coordinate_filter_is_opt_in_and_excludes_repeats_or_unknowns(qtbot,tmp_path,monkeypatch):
    roots,descriptor=fixture_collection(tmp_path);audit=tmp_path/'audit';audit.mkdir();write_report(audit/'sparse_coordinate_qc',descriptor)
    monkeypatch.setenv('PARASCAPE_AUDIT',str(audit));panel=SparseValuesPanel(roots['sparse_datasets'],coordinate_root=roots['sparse_coordinates'],axis_root=roots['sparse_axes']);qtbot.addWidget(panel)
    assert panel.pack.count()==1 and not panel.unique_coordinates.isChecked()
    panel.unique_coordinates.setChecked(True);assert panel.pack.count()==0 and not panel.show_button.isEnabled()
    panel.unique_coordinates.setChecked(False);assert panel.pack.count()==1
    write_report(audit/'sparse_coordinate_qc',descriptor,unordered=True);panel.refresh();panel.unique_coordinates.setChecked(True)
    assert panel.pack.count()==0
