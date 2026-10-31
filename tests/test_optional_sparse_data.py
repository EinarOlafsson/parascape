"""Dataset downloads must preserve artifact identity and explicit local collections."""
import json

import pytest

from parascape import data_components as C, optional_sparse_data as O
from parascape.sparse_values import SparseDataset
from test_sparse_values import fixture_collection
from test_data_component_builder import B
from test_data_components import transport


@pytest.fixture
def hosted(tmp_path, monkeypatch):
    roots, descriptor = fixture_collection(tmp_path)
    info = json.loads(descriptor.read_text())
    components, payloads = {}, {}
    references = [('sparse_coordinates', info['coordinate_artifact'])]
    references.extend(('sparse_axes', ref) for role in ('features', 'barcodes') for ref in info['axis_roles'][role]['original_axis_artifacts'])
    for kind, ref in references:
        identity = kind + '/' + ref['directory']
        built = B.build_component(roots[kind]/ref['directory'], identity, tmp_path/'parts', part_limit=4096)
        components[identity] = built['component']
        payloads.update({part['url']:(tmp_path/'parts'/identity.replace('/','-')/part['url'].rsplit('/',1)[1]).read_bytes()
                         for part in built['component']['parts']})
    index = dict(schema_version=1, publication_status='published_verified', components=components,
                 datasets={descriptor.stem:dict(source_url=info['source_url'],components=list(components))})
    index_path = tmp_path/'data_components.json'
    index_path.write_text(json.dumps(index))
    original_roots = {}
    for kind in ('sparse_coordinates', 'sparse_axes'):
        original_roots[kind] = roots[kind].rename(tmp_path/('original-'+kind))
    monkeypatch.setenv('PARASCAPE_COMPONENTS', str(tmp_path/'shared'))
    monkeypatch.delenv('PARASCAPE_SPARSE_COORDINATES', raising=False)
    monkeypatch.delenv('PARASCAPE_SPARSE_AXES', raising=False)
    return dict(roots=roots,original_roots=original_roots,descriptor=descriptor,index=index,index_path=index_path,payloads=payloads)


def test_shared_download_and_offline_reader_routes(hosted, monkeypatch):
    item = hosted
    p = O.plan(item['descriptor'])
    assert p['registered'] and len(p['missing']) == 3
    assert p['coordinate_root'] == C.storage_root()/'sparse_coordinates'
    calls = transport(monkeypatch,item['payloads'])
    for component in p['missing']:
        C.download(component)
    downloaded = len(calls)
    ready = O.plan(item['descriptor'])
    assert not ready['missing']
    dataset = SparseDataset(item['descriptor'], coordinate_root=ready['coordinate_root'], axis_root=ready['axis_root'])
    records = [r for batch in dataset.iter_entries([1,2]) for r in batch.to_pylist()]
    assert len(records)==4 and records[1]['value']==2**53+1
    for component in ready['components']:
        C.download(component)
    assert len(calls)==downloaded


def test_explicit_local_paths_are_used_and_never_overwritten(hosted):
    item = hosted
    p = O.plan(item['descriptor'], coordinate_root=item['original_roots']['sparse_coordinates'],axis_root=item['original_roots']['sparse_axes'])
    assert p['components']==p['missing']==[]
    assert p['coordinate_root']==item['original_roots']['sparse_coordinates']
    with pytest.raises(ValueError,match='explicit'):
        O.plan(item['descriptor'],coordinate_root=item['roots']['sparse_coordinates'])
    assert not C.storage_root().exists()


@pytest.mark.parametrize('change',['dependency','artifact','unpublished'])
def test_incompatible_or_provisional_index_refused(hosted,change):
    item=hosted;index=item['index']
    if change=='dependency':index['datasets'][item['descriptor'].stem]['components'].pop()
    elif change=='artifact':next(iter(index['components'].values()))['files']['manifest.json']['sha256']='0'*64
    else:index['publication_status']='prepared locally; not uploaded'
    item['index_path'].write_text(json.dumps(index))
    with pytest.raises(ValueError):O.plan(item['descriptor'])
    assert not C.storage_root().exists()


def test_legacy_index_absence_keeps_existing_roots(hosted):
    item=hosted;item['index_path'].unlink()
    p=O.plan(item['descriptor'],coordinate_root=item['original_roots']['sparse_coordinates'])
    assert not p['registered'] and not p['missing']
    assert p['coordinate_root']==item['original_roots']['sparse_coordinates']
