"""Binary source proof and derived tokens remain distinct from text coordinates."""
import json

import pytest

from parascape import sparse_values as S
from scripts.build_data_components import build_component
from test_sparse_values import fixture_collection, save_descriptor, open_dataset


def binary_fixture(tmp_path):
    roots,path=fixture_collection(tmp_path,False,False)
    descriptor=json.loads(path.read_text());reference=descriptor['coordinate_artifact']
    folder=roots['sparse_coordinates']/reference['directory']
    matrix=json.loads((folder/'manifest.json').read_text());receipt=json.loads((folder/'validation.json').read_text())
    storage=dict(format='HDF5_CSC',source_index_base=0,output_index_base=1,decimal_tokens_are_derived=True,
                 coordinate_aggregation_applied=False,token_scope='Derived decimal rendering of original binary arrays')
    matrix.update(source_format='HDF5_CSC',source_storage=storage,original_numeric_array_logical_bytes_sha256='a'*64)
    (folder/'manifest.json').write_text(json.dumps(matrix));reference['manifest_sha256']=S._sha(folder/'manifest.json')
    receipt.update(manifest_sha256=reference['manifest_sha256'],state='all_original_HDF5_CSC_records_reconciled',
                   all_original_binary_values_and_order_checked=True,derived_coordinate_and_value_tokens_checked=True,original_text_tokens_present=False,
                   independent_original_numeric_array_bytes_sha256='a'*64)
    receipt.pop('all_original_coordinate_tokens_and_order_checked')
    (folder/'validation.json').write_text(json.dumps(receipt));reference['validation_sha256']=S._sha(folder/'validation.json')
    descriptor['source_storage']=storage
    return roots,save_descriptor(roots,descriptor),folder


def test_binary_source_receipt_reads_unchanged_values_and_records_derived_token_origin(tmp_path):
    roots,path,folder=binary_fixture(tmp_path)
    dataset=open_dataset(roots,path)
    assert dataset.query_record([1])['source_storage']['decimal_tokens_are_derived'] is True
    assert sum(batch.num_rows for batch in dataset.iter_entries([1]))>0
    component=build_component(folder,'sparse_coordinates/'+folder.name,tmp_path/'parts',tag='v0.0.9.0')
    assert component['component']['files']['entries.parquet']['sha256']==S._sha(folder/'entries.parquet')


@pytest.mark.parametrize('field,value',[('all_original_binary_values_and_order_checked',False),
                                        ('derived_coordinate_and_value_tokens_checked',False),
                                        ('original_text_tokens_present',True)])
def test_binary_source_cannot_inherit_or_exaggerate_text_validation(tmp_path,field,value):
    roots,path,folder=binary_fixture(tmp_path)
    descriptor=json.loads(path.read_text());receipt=json.loads((folder/'validation.json').read_text());receipt[field]=value
    (folder/'validation.json').write_text(json.dumps(receipt));descriptor['coordinate_artifact']['validation_sha256']=S._sha(folder/'validation.json')
    with pytest.raises(ValueError,match='source proof'):open_dataset(roots,save_descriptor(roots,descriptor))
    with pytest.raises(ValueError,match='source proof'):build_component(folder,'sparse_coordinates/'+folder.name,tmp_path/'parts')


def test_binary_descriptor_cannot_hide_original_token_scope(tmp_path):
    roots,path,_=binary_fixture(tmp_path);descriptor=json.loads(path.read_text());descriptor.pop('source_storage')
    with pytest.raises(ValueError,match='storage context'):open_dataset(roots,save_descriptor(roots,descriptor))
