"""Release assets include portable runtime evidence while excluding raw inputs and install markers."""
import hashlib
import json
from pathlib import Path
import tarfile

import pytest

from scripts.build_data_release import build


def test_portable_expression_and_audit_data_are_checksum_pinned(tmp_path):
    source = tmp_path / 'source'
    files = {'nodes.parquet': b'core nodes', 'rnaseq/snapshot.json': b'{}',
             'published_expression/source_pack/ToxoDB_release38.parquet': b'source measurements',
             'published_expression/source_pack/prepare_published_rsem.py': b'# original adapter source',
             'processing_audit/dataset_coverage.tsv': b'original coverage',
             'processing_audit/all_geo_samples/GSE1/geo-samples.json': b'{}'}
    for name, body in files.items():
        path=source/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(body)
    for name in ['data_release.json','processing_audit/source.body','processing_audit/download.part']:
        path=source/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('must not enter the release')
    descriptor = tmp_path / 'descriptor.json'
    info = build(source,tmp_path/'output',descriptor)
    assert set(info['files'])==set(files)
    asset=tmp_path/'output'/info['asset']
    assert hashlib.sha256(asset.read_bytes()).hexdigest()==info['sha256']
    assert json.loads(descriptor.read_text())==info
    with tarfile.open(asset) as archive:
        assert set(archive.getnames())==set(files)
        for name,body in files.items():
            assert archive.extractfile(name).read()==body
            assert info['files'][name]['sha256']==hashlib.sha256(body).hexdigest()


def test_release_builder_rejects_a_linked_runtime_data_directory(tmp_path):
    source=tmp_path/'source';source.mkdir();(source/'nodes.parquet').write_bytes(b'core')
    elsewhere=tmp_path/'elsewhere';elsewhere.mkdir();(elsewhere/'measurement.json').write_text('{}')
    (source/'published_expression').symlink_to(elsewhere,target_is_directory=True)
    with pytest.raises(ValueError,match='symlinks'):
        build(source,tmp_path/'output',tmp_path/'descriptor.json')
