"""Optional archives are independent, resumable and invisible until fully verified."""
import hashlib
import io
import json
import tarfile
import urllib.error

import pytest

from parascape import data_components as C


def fixture_component(members=None):
    files = {'manifest.json': b'{"source":"fixture"}', 'values.parquet': b'literal fixture payload'}
    content = io.BytesIO()
    with tarfile.open(fileobj=content, mode='w') as bundle:
        for name, data in (members or files).items():
            entry = tarfile.TarInfo(name)
            entry.size = len(data)
            bundle.addfile(entry, io.BytesIO(data))
    archive = content.getvalue()
    blocks = [archive[i:i+4096] for i in range(0, len(archive), 4096)]
    parts = [dict(url=f'https://github.com/EinarOlafsson/parascape/releases/download/v1/part-{i}',
                  size=len(block), sha256=hashlib.sha256(block).hexdigest()) for i, block in enumerate(blocks)]
    component = dict(directory='sparse_axes/' + 'a'*64, files={name: dict(size=len(data), sha256=hashlib.sha256(data).hexdigest())
                                                             for name, data in files.items()},
                     archive_size=len(archive), archive_sha256=hashlib.sha256(archive).hexdigest(), parts=parts)
    return component, dict(zip((part['url'] for part in parts), blocks)), files


class Response(io.BytesIO):
    def __init__(self, data, status=200, headers=None):
        super().__init__(data)
        self.status, self.headers = status, headers or {}


def transport(monkeypatch, payloads, *, calls=None):
    calls = [] if calls is None else calls
    def open_request(request, timeout):
        data = payloads[request.full_url]
        offset = int(request.headers.get('Range', 'bytes=0-')[6:-1])
        calls.append((request.full_url, offset))
        if offset:
            return Response(data[offset:], 206, {'Content-Range': f'bytes {offset}-{len(data)-1}/{len(data)}'})
        return Response(data)
    monkeypatch.setattr(C.urllib.request, 'urlopen', open_request)
    return calls


def test_all_parts_files_atomic_visibility_and_idempotent_reuse(tmp_path, monkeypatch):
    component, payloads, files = fixture_component()
    calls = transport(monkeypatch, payloads)
    root = tmp_path / 'components'
    other = root / 'user.txt'
    root.mkdir()
    other.write_text('preserve')
    assert not C.status(component, root)['ready']
    installed = C.download(component, root)
    assert C.status(component, root, verify=True)['ready']
    assert {name: (installed/name).read_bytes() for name in files} == files
    assert other.read_text() == 'preserve'
    downloads = len(calls)
    assert C.download(component, root) == installed and len(calls) == downloads
    assert not list(root.glob('.component-stage-*'))


def test_cancel_resume_verified_completed_parts_and_partial_range(tmp_path, monkeypatch):
    component, payloads, _ = fixture_component()
    calls = transport(monkeypatch, payloads)
    root = tmp_path/'components'
    # First cancellation occurs after a complete first response was written,
    # retaining its bytes for checksum validation on retry.
    with pytest.raises(InterruptedError):
        C.download(component, root, cancelled=lambda: len(calls) >= 2)
    assert not C.status(component, root)['ready']
    assert not (root/component['directory']).exists()
    C.download(component, root)
    assert C.status(component, root, verify=True)['ready']
    assert sum(url == component['parts'][0]['url'] for url, _ in calls) == 1


def test_existing_partial_uses_exact_range_and_final_sha(tmp_path, monkeypatch):
    component, payloads, _ = fixture_component()
    calls = transport(monkeypatch, payloads)
    cache = tmp_path/'.download-parts'
    cache.mkdir()
    part = component['parts'][0]
    (cache/(part['sha256']+'.part')).write_bytes(payloads[part['url']][:17])
    C.download(component, tmp_path)
    assert calls[0] == (part['url'],17)
    assert C.status(component, tmp_path, verify=True)['ready']


@pytest.mark.parametrize('bad', ['part', 'archive', 'file', 'extra', 'traversal', 'link'])
def test_corruption_or_unsafe_archive_never_becomes_visible(tmp_path, monkeypatch, bad):
    component, payloads, files = fixture_component({'../escape': b'bad'} if bad=='traversal' else
                                                 {'extra': b'bad'} if bad=='extra' else None)
    if bad=='part':
        payloads[component['parts'][0]['url']] = b'x'*component['parts'][0]['size']
    elif bad=='archive':
        component['archive_sha256']='0'*64
    elif bad=='file':
        component['files']['manifest.json']['sha256']='0'*64
    elif bad=='link':
        content=io.BytesIO()
        with tarfile.open(fileobj=content,mode='w') as bundle:
            member=tarfile.TarInfo('manifest.json');member.type=tarfile.SYMTYPE;member.linkname='/outside';bundle.addfile(member)
        body=content.getvalue();component['parts']=[dict(url=component['parts'][0]['url'],size=len(body),sha256=hashlib.sha256(body).hexdigest())];component['archive_size']=len(body);component['archive_sha256']=hashlib.sha256(body).hexdigest();payloads={component['parts'][0]['url']:body}
    transport(monkeypatch,payloads)
    with pytest.raises(ValueError):C.download(component,tmp_path/'components')
    assert not (tmp_path/'components'/component['directory']).exists()
    assert not (tmp_path/'escape').exists()


def test_unmanaged_data_and_bad_range_are_preserved(tmp_path,monkeypatch):
    component,payloads,_=fixture_component();target=tmp_path/component['directory'];target.mkdir(parents=True);(target/'user').write_text('preserve')
    with pytest.raises(ValueError,match='unmanaged'):C.download(component,tmp_path)
    assert (target/'user').read_text()=='preserve'
    other=tmp_path/'empty';cache=other/'.download-parts';cache.mkdir(parents=True);part=component['parts'][0];(cache/(part['sha256']+'.part')).write_bytes(b'prefix')
    monkeypatch.setattr(C.urllib.request,'urlopen',lambda *a,**k:Response(payloads[part['url']],206,{'Content-Range':'bytes 0-1/2'}))
    with pytest.raises(ValueError,match='range'):C.download(component,other)
    assert not (other/component['directory']).exists()


@pytest.mark.parametrize('change', ['directory','file','url','part_size'])
def test_contract_rejected_before_network_or_writes(tmp_path,monkeypatch,change):
    component,_,_=fixture_component()
    if change=='directory':component['directory']='../escape'
    elif change=='file':component['files']['../escape']=component['files']['manifest.json']
    elif change=='url':component['parts'][0]['url']='https://example.com/data'
    else:component['parts'][0]['size']=1024**3+1;component['archive_size']=sum(x['size'] for x in component['parts'])
    monkeypatch.setattr(C.urllib.request,'urlopen',lambda *a,**k:pytest.fail('network attempted'))
    with pytest.raises(ValueError):C.download(component,tmp_path/'new')
    assert not (tmp_path/'new').exists()


def test_corrupt_managed_component_is_replaced_with_recovery_copy(tmp_path, monkeypatch):
    component, payloads, files = fixture_component()
    transport(monkeypatch, payloads)
    installed = C.download(component, tmp_path)
    (installed/'manifest.json').write_bytes(b'x' * len(files['manifest.json']))
    assert not C.status(component, tmp_path, verify=True)['ready']
    C.download(component, tmp_path)
    assert C.status(component, tmp_path, verify=True)['ready']
    previous = list(installed.parent.glob(installed.name + '.previous-*'))
    assert len(previous) == 1 and (previous[0]/'manifest.json').read_bytes() == b'x' * len(files['manifest.json'])


def test_cached_part_link_is_refused_without_touching_target(tmp_path, monkeypatch):
    component, payloads, _ = fixture_component()
    cache = tmp_path/'.download-parts'
    cache.mkdir()
    user_file = tmp_path/'keep.bin'
    user_file.write_bytes(payloads[component['parts'][0]['url']])
    (cache/component['parts'][0]['sha256']).symlink_to(user_file)
    transport(monkeypatch, payloads)
    with pytest.raises(ValueError, match='Unsafe cached'):
        C.download(component, tmp_path)
    assert user_file.read_bytes() == payloads[component['parts'][0]['url']]
