"""Portable workbook discovery, on-demand install and offline source queries."""
import json
import shutil

import pytest

from parascape import data_components as C
from parascape import optional_workbook_data as O
from parascape import workbook_values as W
from scripts.build_data_components import build_component
from test_data_components import transport
from test_workbook_values import fixture_collection


def portable_collection(tmp_path,monkeypatch,*,legacy=False):
    source,_=fixture_collection(tmp_path/'source',legacy=legacy)
    result=build_component(source,'workbook_cells/'+source.name,tmp_path/'parts',tag='v0.0.9.0',part_limit=4096)
    component=result['component']
    payloads={p['url']:(tmp_path/'parts'/('workbook_cells-'+source.name)/p['url'].rsplit('/',1)[-1]).read_bytes()
              for p in component['parts']}
    core=tmp_path/'core';metadata=core/'workbook_cells'/source.name;metadata.mkdir(parents=True)
    for name in ['manifest.json','independent_validation.json']:shutil.copyfile(source/name,metadata/name)
    shutil.copyfile(source.parent/'catalogue.json',metadata.parent/'catalogue.json')
    # A software-fixture publication marker exercises local transport only.
    (core/'data_components.json').write_text(json.dumps(dict(schema_version=1,publication_status='published_verified',
         components={component['directory']:component})))
    monkeypatch.setenv('PARASCAPE_COMPONENTS',str(tmp_path/'shared-cache'))
    monkeypatch.delenv('PARASCAPE_WORKBOOK_VALUES',raising=False)
    calls=transport(monkeypatch,payloads)
    return source,metadata,component,calls


@pytest.mark.parametrize('legacy',[False,True])
def test_complete_component_roundtrip_literal_query_and_offline_reuse(tmp_path,monkeypatch,legacy):
    source,metadata,component,calls=portable_collection(tmp_path,monkeypatch,legacy=legacy)
    entries,errors=W.catalogue(metadata.parent)
    assert len(entries)==1 and not errors
    dependencies=O.plan(metadata)
    assert dependencies['missing']==[component]
    with pytest.raises(FileNotFoundError):W.prepare_query(metadata,0,['001'])
    installed=C.download(component)
    for name in component['files']:assert (source/name).read_bytes()==(installed/name).read_bytes()
    assert not O.plan(metadata)['missing']
    calls_before=len(calls)
    query=W.prepare_query(metadata,0,['001'],label_column=1,preview_cells=1)
    assert query.record['selected_physical_cells']==7
    assert W.export_csv(query,tmp_path/'portable.csv')==7
    assert len(calls)==calls_before
    assert query.record['validation']['independent_whole_source_reconciliation'] is not legacy
    assert not list(metadata.glob('*.parquet'))


def test_provisional_or_mismatched_components_cannot_download(tmp_path,monkeypatch):
    _,metadata,component,_=portable_collection(tmp_path,monkeypatch)
    path=metadata.parent.parent/'data_components.json';index=json.loads(path.read_text())
    index['publication_status']='provisional';path.write_text(json.dumps(index))
    with pytest.raises(ValueError,match='not yet available'):O.plan(metadata)
    index['publication_status']='published_verified'
    index['components'][component['directory']]['files']['cells.parquet']['sha256']='f'*64
    path.write_text(json.dumps(index))
    with pytest.raises(ValueError,match='does not match'):O.plan(metadata)


def test_explicit_local_override_is_preserved(tmp_path,monkeypatch):
    _,metadata,_,_=portable_collection(tmp_path,monkeypatch)
    monkeypatch.setenv('PARASCAPE_WORKBOOK_VALUES',str(metadata.parent))
    with pytest.raises(ValueError,match='explicit'):O.plan(metadata)


def test_core_builder_keeps_workbook_payload_outside_small_metadata_bundle(tmp_path):
    from scripts.build_data_release import build
    source,_=fixture_collection(tmp_path/'data/workbook_cells')
    (tmp_path/'data/nodes.parquet').write_bytes(b'core-node fixture')
    record=build(tmp_path/'data',tmp_path/'output',tmp_path/'descriptor.json')
    assert 'workbook_cells/catalogue.json' in record['files']
    assert f'workbook_cells/{source.name}/manifest.json' in record['files']
    assert f'workbook_cells/{source.name}/independent_validation.json' in record['files']
    assert not any(name.startswith('workbook_cells/') and name.endswith('.parquet') for name in record['files'])
