"""Source-bound storage QC is distinct from biological quality and reader agreement."""
import copy
import json
import struct

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from parascape import workbook_qc as Q
from parascape import workbook_values as W
from test_workbook_values import fixture_collection


def biff_fixture(root):
    folder,_=fixture_collection(root,legacy=True)
    path=folder/'cells.parquet';table=pq.read_table(path);rows=table.to_pylist()
    row=next(r for r in rows if r['row']==2 and r['column']==3)
    row.update(storage_type='BIFF:error_code',decoded_text='29',numeric_binary64_hex=None)
    pq.write_table(pa.Table.from_pylist(rows,schema=table.schema),path,row_group_size=2)
    info=json.loads((folder/'manifest.json').read_text())
    info['outputs']['cells.parquet']=dict(sha256=W._sha(path),size=path.stat().st_size)
    (folder/'manifest.json').write_text(json.dumps(info))
    original=json.loads((folder/'independent_validation.json').read_text())
    original['measurement_manifest_sha256']=W._sha(folder/'manifest.json')
    (folder/'independent_validation.json').write_text(json.dumps(original))
    return folder


def wrapper(folder,*,differences=0,formula=True):
    info=W.manifest(folder);count=info['extraction']['cells'];digest=W._sha(folder/'manifest.json')
    tokens=b'\x43\x01\x00\x00\x00\x13';cached=b'\x02\x00\x1d\x00\xa3\x00\xff\xff'
    payload=struct.pack('<HHH',1,2,0)+cached+struct.pack('<HIH',0,0,len(tokens))+tokens
    preview=[dict(sheet_index=0,row=2,column=3,record_offset=100,cached_result_bytes_hex=cached.hex(),
                  formula_tokens_hex=tokens.hex(),complete_record_payload_hex=payload.hex())] if formula else []
    raw=dict(schema_version=1,status='all_BIFF_physical_cells_and_decoded_values_reconciled',
        measurement_manifest_sha256=digest,input_sha256=info['input_sha256'],all_physical_cells_and_decoded_values_verified=True,
        source_and_existing_pack_bytes_changed=False,formula_bytes_preserved=True,original_formulas_recalculated=False,
        measurement_admission=False,biological_quality_score=None,verifier_sha256='9'*64,
        physical_cells=count,decoded_cells_checked=count,global_repeated_cell_coordinates=0,
        original_workbook_stream_sha256='d'*64,
        outputs={name:dict(sha256=sha*64,size=4096) for name,sha in [('workbook_stream.bin','d'),('physical_cells.parquet','e'),('formula_records.parquet','f')]},
        physical_cell_kinds=dict(text=6,numeric=3,blank=1,error=1),original_error_code_counts={'29':1},
        text_confirmation_paths={'raw_LABEL_declared_codepage':6},independent_Calamine_text_disagreements=differences,
        independent_Calamine_whole_text_agreement=differences==0,
        sheets=[dict(sheet_index=0,name='Source quantities',cells=count)],formula_records=len(preview),
        mixed_BIFF_version_declarations=True,calamine_version='fixture-only',olefile_version='fixture-only')
    return dict(schema_version=1,source_directory=folder.name,measurement_manifest_sha256=digest,input_sha256=info['input_sha256'],
        raw_audit_directory='c'*64,raw_report_sha256='8'*64,raw_report=raw,assembler_sha256='7'*64,
        formula_preview_limit=100,formula_records_preview=preview,formula_preview_is_complete=True,
        numerical_values_changed=False,original_receipts_replaced=False,biological_or_comparative_admission=False,biological_quality_score=None)


def install_qc(folder,value):
    root=folder.parent/'source_qc';root.mkdir(exist_ok=True)
    body=(json.dumps(value,sort_keys=True,indent=2)+'\n').encode()
    import hashlib
    digest=hashlib.sha256(body).hexdigest();(root/(digest+'.json')).write_bytes(body)
    index=dict(schema_version=1,measurement_admission=False,biological_quality_score=None,sources={folder.name:dict(
        file=digest+'.json',sha256=digest,measurement_manifest_sha256=value['measurement_manifest_sha256'],input_sha256=value['input_sha256'])})
    (root/'catalogue.json').write_text(json.dumps(index))
    return root,index


def test_supplement_preserves_original_receipt_csv_and_separate_agreement(tmp_path):
    folder=biff_fixture(tmp_path/'collection')
    original=(folder/'independent_validation.json').read_bytes()
    before=W.prepare_query(folder,0,['001'],label_column=1)
    assert not before.record['validation']['independent_whole_source_reconciliation']
    W.export_csv(before,tmp_path/'before.csv')
    install_qc(folder,wrapper(folder,differences=1))
    proof=W.validation(folder)
    assert proof['independent_whole_source_reconciliation']
    assert proof['decoded_value_checks_agree'] is False and proof['reported_reader_disagreements']==1
    assert proof['biological_quality_score'] is None and not proof['original_validation']['independent_whole_source_reconciliation']
    assert (folder/'independent_validation.json').read_bytes()==original
    query=W.prepare_query(folder,0,['001'],label_column=1,preview_cells=1)
    assert W.export_csv(query,tmp_path/'after.csv')==7
    assert (tmp_path/'after.csv').read_bytes()==(tmp_path/'before.csv').read_bytes()
    record=json.loads((tmp_path/'after.csv.json').read_text())
    qc=record['validation']['supplementary_biff_qc']
    formula=Q.formula_for(qc,0,dict(row=2,column=3))
    assert formula['cached_result_kind']=='error' and formula['cached_result_code']==29
    assert formula['formula_evaluated'] is False
    assert Q.formula_for(qc,0,dict(row=2,column=2)) is None


@pytest.mark.parametrize('change,expected',[
    (lambda x:x.update(input_sha256='f'*64),'binding'),
    (lambda x:x.update(source_directory='0'*64),'identity'),
    (lambda x:x.update(numerical_values_changed=True),'transformation'),
    (lambda x:x['raw_report'].update(decoded_cells_checked=999),'coverage'),
    (lambda x:x['raw_report'].update(biological_quality_score=0.99),'claim'),
    (lambda x:x['raw_report'].update(independent_Calamine_whole_text_agreement=True),'agreement'),
    (lambda x:x['raw_report'].update(original_error_code_counts={}), 'error coverage'),
    (lambda x:x['raw_report'].update(text_confirmation_paths={}), 'text coverage'),
    (lambda x:x['raw_report'].update(global_repeated_cell_coordinates=1),'repeated'),
    (lambda x:x['formula_records_preview'][0].update(formula_tokens_hex='00'),'token'),
    (lambda x:x['formula_records_preview'][0].update(row=99),'coordinates'),
    (lambda x:x.update(formula_preview_is_complete=False),'preview'),
])
def test_internally_rehashed_but_false_qc_is_refused(tmp_path,change,expected):
    folder=biff_fixture(tmp_path/'collection');value=wrapper(folder,differences=1)
    change(value);install_qc(folder,value)
    with pytest.raises(ValueError,match=expected):W.validation(folder)


def test_unknown_qc_stays_unknown_and_unsafe_qc_path_is_refused(tmp_path):
    folder=biff_fixture(tmp_path/'collection')
    assert Q.report(folder) is None and W.validation(folder)['decoded_value_checks_agree'] is None
    root,index=install_qc(folder,wrapper(folder))
    index['sources'][folder.name]['file']='../other.json';(root/'catalogue.json').write_text(json.dumps(index))
    with pytest.raises(ValueError,match='filename'):Q.report(folder)


def test_portable_query_keeps_core_qc_origin_after_payload_install(tmp_path,monkeypatch):
    from test_optional_workbook_data import portable_collection
    from parascape import data_components as C
    source,metadata,component,_=portable_collection(tmp_path,monkeypatch,legacy=True)
    # The standard fixture has four numeric cells and no original cached error.
    value=wrapper(source,differences=1,formula=False)
    value['raw_report']['physical_cell_kinds']=dict(text=6,numeric=4,blank=1)
    value['raw_report']['original_error_code_counts']={}
    install_qc(metadata,value)
    C.download(component)
    query=W.prepare_query(metadata,0,['001'],label_column=1)
    assert query.folder!=metadata
    assert query.record['validation']['independent_whole_source_reconciliation']
    assert query.record['validation']['decoded_value_checks_agree'] is False
    assert query.record['validation']['supplementary_biff_qc']['metadata_file_sha256']
    with pytest.raises(ValueError,match='immutable'):
        W.export_csv(query,metadata/'not-a-source.csv')
