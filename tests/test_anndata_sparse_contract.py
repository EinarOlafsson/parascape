"""AnnData layers must retain observation/variable orientation and exact binary proof."""
import json

import pytest
import pyarrow as pa
import pyarrow.parquet as pq

from parascape import sparse_values as S
from scripts.build_data_components import build_component
from test_sparse_values import fixture_collection, save_descriptor, open_dataset


def fixture(tmp_path, layout='csr'):
    roots,path=fixture_collection(tmp_path,True,True)
    descriptor=json.loads(path.read_text());ref=descriptor['coordinate_artifact']
    folder=roots['sparse_coordinates']/ref['directory']
    matrix=json.loads((folder/'manifest.json').read_text())
    if layout=='dense':
        table=pq.read_table(folder/'entries.parquet')
        records=table.to_pylist()
        for position,record in enumerate(records):
            record.update(source_row=position//3+1,source_column=position%3+1,
                row_token=str(position//3+1),column_token=str(position%3+1))
        for position in [4,5]:
            records.append(dict(entry_index=position,source_row=2,source_column=position%3+1,
                row_token='2',column_token=str(position%3+1),source_token='0',value=0.))
        pq.write_table(pa.Table.from_pylist(records,schema=table.schema),folder/'entries.parquet')
        matrix['declared_dimensions'][2]=6;descriptor['declared_dimensions'][2]=6
        matrix['outputs']['entries.parquet']=S._sha(folder/'entries.parquet')
        ref['output_sha256']=matrix['outputs']['entries.parquet']
    storage=dict(format='AnnData_HDF5',layout=layout,matrix_path='layers/spliced',
        source_axis_order=['observations','variables'],source_index_base=0,output_index_base=1,
        decimal_tokens_are_derived=True,coordinate_aggregation_applied=False,
        original_dtype_roundtrip_verified=True,layer_aggregation_applied=False)
    matrix.update(source_format='AnnData_HDF5',source_storage=storage,
        original_numeric_array_logical_bytes_sha256='a'*64)
    (folder/'manifest.json').write_text(json.dumps(matrix));ref['manifest_sha256']=S._sha(folder/'manifest.json')
    receipt=json.loads((folder/'validation.json').read_text())
    receipt.update(manifest_sha256=ref['manifest_sha256'],state='all_original_AnnData_layer_records_reconciled',
        all_original_binary_values_and_order_checked=True,derived_coordinate_and_value_tokens_checked=True,
        original_text_tokens_present=False,independent_original_numeric_array_bytes_sha256='a'*64,
        all_original_dtype_bits_recoverable=True,every_dense_cell_retained=layout=='dense')
    receipt.update(original_records=matrix['declared_dimensions'][2],matrix_dimensions=matrix['declared_dimensions'])
    (folder/'validation.json').write_text(json.dumps(receipt));ref['validation_sha256']=S._sha(folder/'validation.json')
    descriptor['source_storage']=storage
    return roots,save_descriptor(roots,descriptor),folder


@pytest.mark.parametrize('layout',['csr','csc','dense'])
def test_original_ann_layer_is_queryable_and_packaged(tmp_path,layout):
    roots,path,folder=fixture(tmp_path,layout)
    dataset=open_dataset(roots,path)
    assert dataset.feature_coordinate=='source_column'
    rows=[r for batch in dataset.iter_entries([1]) for r in batch.to_pylist()]
    assert len(rows)==(2 if layout=='dense' else 3) and rows[0]['value']==0
    assert dataset.query_record([1])['source_storage']['matrix_path']=='layers/spliced'
    result=build_component(folder,'sparse_coordinates/'+folder.name,tmp_path/'parts')
    assert result['component']['files']['entries.parquet']['sha256']==S._sha(folder/'entries.parquet')


@pytest.mark.parametrize('field,value',[
    ('layout','unknown'),('matrix_path',''),('source_axis_order',['variables','observations']),
    ('original_dtype_roundtrip_verified',False),('layer_aggregation_applied',True)])
def test_incomplete_ann_layer_storage_refused_by_reader_and_builder(tmp_path,field,value):
    roots,path,folder=fixture(tmp_path)
    descriptor=json.loads(path.read_text());matrix=json.loads((folder/'manifest.json').read_text())
    matrix['source_storage'][field]=value;descriptor['source_storage'][field]=value
    (folder/'manifest.json').write_text(json.dumps(matrix));ref=descriptor['coordinate_artifact']
    ref['manifest_sha256']=S._sha(folder/'manifest.json')
    receipt=json.loads((folder/'validation.json').read_text());receipt['manifest_sha256']=ref['manifest_sha256']
    (folder/'validation.json').write_text(json.dumps(receipt));ref['validation_sha256']=S._sha(folder/'validation.json')
    with pytest.raises(ValueError,match='AnnData'):open_dataset(roots,save_descriptor(roots,descriptor))
    with pytest.raises(ValueError,match='AnnData'):build_component(folder,'sparse_coordinates/'+folder.name,tmp_path/'parts')


def test_ann_orientation_cannot_be_changed_even_when_axis_lengths_match(tmp_path):
    roots,path,_=fixture(tmp_path)
    descriptor=json.loads(path.read_text())
    descriptor['axis_roles']['features']['source_coordinate']='source_row'
    descriptor['axis_roles']['barcodes']['source_coordinate']='source_column'
    with pytest.raises(ValueError,match='without transposing'):open_dataset(roots,save_descriptor(roots,descriptor))
