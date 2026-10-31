"""Software source-storage fixtures do not establish biological quality."""
import csv
import json
from pathlib import Path
import threading

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from parascape import workbook_values as W


def fixture_collection(root, *, verified=True, legacy=False):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    key = 'a' * 64
    folder = root / key
    folder.mkdir()
    schema = pa.schema([('sheet_index',pa.int32()),('cell_order',pa.int64()),('row',pa.int32()),
                        ('column',pa.int32()),('source_address',pa.string()),('coordinate_inferred',pa.bool_()),
                        *[(name,pa.string()) for name in ('storage_type','source_value_token','decoded_text',
                          'numeric_binary64_hex','formula_xml','source_cell_xml','source_attributes_json')]])
    values = []
    def cell(row,column,kind,token=None,text=None,formula=None):
        values.append(dict(sheet_index=0,cell_order=len(values),row=row,column=column,
                           source_address=f'{chr(64+column)}{row}',coordinate_inferred=False,
                           storage_type=kind,source_value_token=token,decoded_text=text,
                           numeric_binary64_hex=None,formula_xml=formula,source_cell_xml='<c/>',
                           source_attributes_json='{}'))
    cell(1,1,'inlineStr',text='feature');cell(1,2,'inlineStr',text='source quantity')
    cell(2,1,'inlineStr',text='001');cell(2,2,'n',token='9007199254740993')
    cell(2,2,'n',token='-2.5000')  # Repeated coordinates remain physical records.
    cell(2,3,'n',token='0.3333333333333333333',formula='<f>A2/3</f>')
    cell(4,1,'inlineStr',text='001');cell(4,2,'n');cell(4,3,'inlineStr',text='')
    cell(5,1,'inlineStr',text='Other');cell(5,2,'n',token='1')
    if legacy:
        for row in values:
            row['source_address']=None;row['coordinate_inferred']=True;row['source_cell_xml']=None
            row['formula_xml']=None
            if row['storage_type']=='n':
                token=row['source_value_token'];row['source_value_token']=None
                row['numeric_binary64_hex']=float(token).hex() if token is not None else None
                row['decoded_text']=str(float(token)) if token is not None else ''
                row['storage_type']='BIFF:number' if token is not None else 'BIFF:blank'
            else:row['storage_type']='BIFF:text'
    pq.write_table(pa.Table.from_pylist(values,schema=schema),folder/'cells.parquet',row_group_size=2)
    pq.write_table(pa.table({'fixture':pa.array([],type=pa.string())}),folder/'strings.parquet')
    pq.write_table(pa.table({'fixture':pa.array([],type=pa.string())}),folder/'rows.parquet')
    info=dict(schema_version=1,operation_version='source-cells-v1',original_bytes_changed=False,
              input_sha256='b'*64,input_path='/nas_mnt/unavailable-original',source_url='https://ftp.ncbi.nlm.nih.gov/geo/fixture.xlsx',
              source_record=dict(studies='GSE_FIXTURE',filename='fixture.xlsx',source_entities='Series:fixture'),
              treatment=dict(normalization=False,aggregation=False,feature_mapping=False,header_inference=False,
                             formula_evaluation=False,biological_quality_score=None,stage_score=None,strain_score=None),
              extraction=dict(format='BIFF_decoded' if legacy else 'OOXML',cells=len(values),rows=4,shared_strings=0,
                              numeric_representation='Binary64 hex' if legacy else 'Original numeric text',
                              facts=dict(formula_cells=0 if legacy else 1),
                              sheets=[dict(sheet_index=0,source_attributes=dict(name='Source quantities',state='hidden'),
                                           cells=len(values),rows=4,extraction_status='source_cells_retained')]),
              outputs={p.name:dict(sha256=W._sha(p),size=p.stat().st_size) for p in folder.glob('*.parquet')})
    (folder/'manifest.json').write_text(json.dumps(info))
    if verified:
        proof=dict(measurement_manifest_sha256=W._sha(folder/'manifest.json'),input_sha256=info['input_sha256'],
                   measurement_admission=False,biological_quality_score=None,
                   independent_whole_source_reconciliation=not legacy,
                   status='BIFF_checksums_verified_independent_reader_review_required' if legacy else 'all_OOXML_source_cells_independently_reconciled',
                   cells=len(values),rows=4,shared_strings=0)
        (folder/'independent_validation.json').write_text(json.dumps(proof))
    (root/'catalogue.json').write_text(json.dumps(dict(workbooks=[dict(directory=key,status='complete_source_cell_extraction',
                       source_url=info['source_url'],input_sha256=info['input_sha256'],source_occurrences=[dict(studies='GSE_FIXTURE')])])) )
    return folder, values


def test_complete_rows_preserve_numeric_tokens_formulas_and_duplicate_coordinates(tmp_path):
    folder,values=fixture_collection(tmp_path)
    query=W.prepare_query(folder,0,['001'],label_column=1,preview_cells=1)
    assert len(query.preview)==1 and query.record['selected_physical_cells']==7
    target=tmp_path/'all.csv'
    assert W.export_csv(query,target)==7
    with target.open(newline='') as stream:exported=list(csv.DictReader(stream))
    assert [row['cell_order'] for row in exported]==[str(x['cell_order']) for x in values if x['row'] in (2,4)]
    assert [row['source_value_token'] for row in exported[:4]]==['','9007199254740993','-2.5000','0.3333333333333333333']
    assert exported[3]['formula_xml']=='<f>A2/3</f>'
    assert exported[5]['decoded_text_present']=='False' and exported[6]['decoded_text_present']=='True'
    sidecar=json.loads((tmp_path/'all.csv.json').read_text())
    assert sidecar['export_sha256']==W._sha(target) and sidecar['exported_physical_cells']==7
    assert sidecar['biological_or_comparative_admission'] is False
    assert sidecar['original_treatment']['biological_quality_score'] is None


def test_exact_text_does_not_coerce_numeric_cells_or_strip_labels(tmp_path):
    folder,_=fixture_collection(tmp_path)
    for label in ['1','001 ','9007199254740993','OTHER']:
        query=W.prepare_query(folder,0,[label])
        assert query.record['selected_physical_cells']==0
        assert query.record['unmatched_labels']==[label]
    assert W.prepare_query(folder,0,['001'],label_column=2).record['selected_physical_cells']==0
    with pytest.raises(ValueError):W.prepare_query(folder,0,[1])


@pytest.mark.parametrize('verified,legacy,grade',[(False,False,False),(True,False,True),(True,True,False)])
def test_unknown_verification_and_legacy_scope_remain_visible(tmp_path,verified,legacy,grade):
    folder,_=fixture_collection(tmp_path,verified=verified,legacy=legacy)
    entries,errors=W.catalogue(tmp_path)
    assert len(entries)==1 and not errors
    assert entries[0]['independent_whole_source_reconciliation'] is grade
    query=W.prepare_query(folder,0,['001'])
    assert query.record['selected_physical_cells']==7
    assert query.record['validation']['biological_quality_score'] is None
    assert query.record['source_record']['studies']=='GSE_FIXTURE'
    if legacy:
        numeric=next(row for block in query.blocks() for row in block.to_pylist() if row['numeric_binary64_hex'])
        assert numeric['source_value_token'] is None


def test_invalid_proof_scope_is_not_a_quality_pass(tmp_path):
    folder,_=fixture_collection(tmp_path)
    proof=folder/'independent_validation.json';record=json.loads(proof.read_text());record['cells']=999
    proof.write_text(json.dumps(record))
    entries,errors=W.catalogue(tmp_path)
    assert not entries and 'coverage differs' in errors[0]
    with pytest.raises(ValueError):W.prepare_query(folder,0,['001'])


def test_changed_file_refuses_export_and_leaves_no_partial_result(tmp_path):
    folder,_=fixture_collection(tmp_path)
    query=W.prepare_query(folder,0,['001'])
    path=folder/'cells.parquet';path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError):W.export_csv(query,tmp_path/'changed.csv')
    assert not (tmp_path/'changed.csv').exists() and not (tmp_path/'changed.csv.json').exists()
    assert not list(tmp_path.glob('*.partial'))


def test_cancellation_refuses_query_and_atomic_export(tmp_path):
    folder,_=fixture_collection(tmp_path)
    cancelled=threading.Event();cancelled.set()
    with pytest.raises(W.QueryCancelled):W.prepare_query(folder,0,['001'],cancel=cancelled)
    query=W.prepare_query(folder,0,['001'])
    with pytest.raises(W.QueryCancelled):W.export_csv(query,tmp_path/'cancelled.csv',cancel=cancelled)
    assert not (tmp_path/'cancelled.csv').exists()


def test_exports_do_not_replace_user_files_or_change_source_artifacts(tmp_path):
    folder,_=fixture_collection(tmp_path)
    query=W.prepare_query(folder,0,['001'])
    existing=tmp_path/'existing.csv';existing.write_text('user work')
    with pytest.raises(FileExistsError):W.export_csv(query,existing)
    assert existing.read_text()=='user work'
    with pytest.raises(ValueError):W.export_csv(query,folder/'new.csv')


def test_unsafe_catalogue_paths_and_missing_collections(tmp_path):
    assert W.catalogue(tmp_path/'absent')==([],[])
    fixture_collection(tmp_path)
    path=tmp_path/'catalogue.json';index=json.loads(path.read_text());index['workbooks'][0]['directory']='../outside'
    path.write_text(json.dumps(index))
    entries,errors=W.catalogue(tmp_path)
    assert not entries and 'directory' in errors[0]
