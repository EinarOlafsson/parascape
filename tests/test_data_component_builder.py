"""Deterministic real artifact packaging round-trips through the independent installer."""
import importlib.util
import json
from pathlib import Path

import pytest

from parascape import data_components as C
from test_sparse_values import fixture_collection
from test_data_components import transport

spec = importlib.util.spec_from_file_location('component_builder', Path(__file__).resolve().parents[1]/'scripts/build_data_components.py')
B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)


def test_archive_parts_preserve_checked_axis_bytes_and_reuse(tmp_path,monkeypatch):
    roots,path=fixture_collection(tmp_path)
    descriptor=json.loads(path.read_text());ref=descriptor['axis_roles']['features']['original_axis_artifacts'][0]
    source=roots['sparse_axes']/ref['directory'];directory='sparse_axes/'+ref['directory']
    first=B.build_component(source,directory,tmp_path/'parts',part_limit=4096)
    assert not first['reused'] and len(first['component']['parts'])>1
    assert B.build_component(source,directory,tmp_path/'parts',part_limit=4096)['reused']
    second=B.build_component(source,directory,tmp_path/'other-parts',part_limit=4096)
    assert first['component']==second['component']
    folder=Path(first['path']);payloads={part['url']:(folder/part['url'].rsplit('/',1)[1]).read_bytes() for part in first['component']['parts']}
    transport(monkeypatch,payloads)
    installed=C.download(first['component'],tmp_path/'installed')
    for name in ('manifest.json','validation.json','axis.parquet'):assert (installed/name).read_bytes()==(source/name).read_bytes()
    damaged=folder/first['component']['parts'][0]['url'].rsplit('/',1)[1];damaged.write_bytes(b'changed')
    with pytest.raises(ValueError,match='changed'):B.build_component(source,directory,tmp_path/'parts',part_limit=4096)


def test_source_mutation_refused_without_published_parts(tmp_path):
    roots,path=fixture_collection(tmp_path);descriptor=json.loads(path.read_text());ref=descriptor['coordinate_artifact'];source=roots['sparse_coordinates']/ref['directory'];output=source/'entries.parquet';output.write_bytes(output.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='changed'):B.build_component(source,'sparse_coordinates/'+ref['directory'],tmp_path/'parts')
    assert not list((tmp_path/'parts').iterdir())
