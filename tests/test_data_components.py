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


def transport(monkeypatch, payloads, *, calls=None, attempts=None):
    """Serve the addresses in `payloads`; answer 404 to any other, as the host that has not got it would.

    `calls` records the addresses that answered, `attempts` every address asked, including the ones
    that answered 404 -- which is what a part still only on GitHub looks like on Hugging Face.
    """
    calls = [] if calls is None else calls
    def open_request(request, timeout):
        if attempts is not None:
            attempts.append(request.full_url)
        if request.full_url not in payloads:
            raise urllib.error.HTTPError(request.full_url, 404, 'Not Found', {}, None)
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


def test_a_stalled_part_is_retried_and_a_missing_one_is_not(monkeypatch, tmp_path):
    """Hundreds of small release parts: one stalled read must not end the install; a 404 must end it at once."""
    import urllib.error
    from parascape import data_components as DC
    monkeypatch.setattr(DC.time if hasattr(DC, "time") else __import__("time"), "sleep", lambda s: None)
    calls = []

    address = 'https://github.com/EinarOlafsson/parascape/releases/download/v1/part-0'

    def flaky(part, cache, progress, cancelled, addresses):
        calls.append(1)
        if len(calls) < 3:
            raise TimeoutError("The read operation timed out")
        return "done"
    monkeypatch.setattr(DC, "_download_part_once", flaky)
    assert DC._download_part({'url': address}, tmp_path, None, None) == "done" and len(calls) == 3

    def missing(part, cache, progress, cancelled, addresses):
        calls.append(1)
        raise urllib.error.HTTPError("u", 404, "Not Found", {}, None)
    calls.clear()
    monkeypatch.setattr(DC, "_download_part_once", missing)
    import pytest
    with pytest.raises(urllib.error.HTTPError):
        DC._download_part({'url': address}, tmp_path, None, None)
    assert len(calls) == 1


# ------------------------------------------- the data move to Hugging Face (2026-11-13)
"""Every component released so far has its parts on GitHub and keeps them, so both bases are
accepted: the one configured for the move and the GitHub releases the published indexes name."""


def test_both_the_configured_host_and_the_published_github_release_are_accepted(monkeypatch):
    monkeypatch.delenv(C.ENV_HOST, raising=False)
    assert C.hosts() == (C.COMPONENT_HOST + '/', C.GITHUB_RELEASES + '/')
    component, _, _ = fixture_component()
    assert C.validate(component)['download_bytes'] == component['archive_size']
    moved = dict(component, parts=[dict(part, url=f'{C.COMPONENT_HOST}/v1/part-{i}')
                                   for i, part in enumerate(component['parts'])])
    assert C.validate(moved) == C.validate(component), 'the same component, served from the new host'


def test_a_component_host_setting_moves_the_accepted_base(monkeypatch):
    monkeypatch.setenv(C.ENV_HOST, 'https://data.example.org/components')
    assert 'https://data.example.org/components/' in C.hosts()
    component, _, _ = fixture_component()
    component['parts'][0] = dict(component['parts'][0], url='https://data.example.org/components/v1/part-0')
    assert C.validate(component)


@pytest.mark.parametrize('bad', ['http://github.com/EinarOlafsson/parascape/releases/download/v1/p',
                                 'https://elsewhere.example/v1/p',
                                 'https://github.com/EinarOlafsson/parascape/releases/download/v1/../p',
                                 'https://github.com/EinarOlafsson/parascape/releases/download/v1/p?t=1'])
def test_an_address_off_the_accepted_bases_is_still_refused(monkeypatch, bad):
    monkeypatch.delenv(C.ENV_HOST, raising=False)
    component, _, _ = fixture_component()
    component['parts'][0] = dict(component['parts'][0], url=bad)
    with pytest.raises(ValueError, match='Unexpected optional-data download host'):
        C.validate(component)


# ------------------------------------- the parts themselves move to Hugging Face (2026-11-13)
"""The index inside an already-published core archive names a GitHub address for every part, and
that archive cannot be rewritten. So the first address is built from the configured host, the kind
and the part name, and the recorded address is the fallback. No new core release is needed, and a
part that has not moved yet is still served."""


def moved(component):
    """Where each part of this component is addressed on the configured host."""
    kind = component['directory'].split('/')[0]
    return [f"{C.COMPONENT_HOST}/{kind}/{part['url'].rsplit('/', 1)[1]}" for part in component['parts']]


def test_the_configured_host_is_tried_first_and_the_recorded_address_second(monkeypatch):
    monkeypatch.delenv(C.ENV_HOST, raising=False)
    component, _, _ = fixture_component()
    first, second = C.part_sources(component, component['parts'][0])
    assert first == (C.COMPONENT_HOST + '/sparse_axes/'
                     + component['parts'][0]['url'].rsplit('/', 1)[1])
    assert second == component['parts'][0]['url']


def test_a_part_on_hugging_face_is_taken_from_there(tmp_path, monkeypatch):
    monkeypatch.delenv(C.ENV_HOST, raising=False)
    component, payloads, files = fixture_component()
    on_hugging_face = dict(zip(moved(component), payloads.values()))
    attempts = []
    transport(monkeypatch, on_hugging_face, attempts=attempts)
    installed = C.download(component, tmp_path / 'components')
    assert C.status(component, tmp_path / 'components', verify=True)['ready']
    assert {name: (installed / name).read_bytes() for name in files} == files
    assert attempts and all('huggingface.co' in url for url in attempts), 'GitHub was never asked'


def test_a_part_still_only_on_github_is_taken_from_there(tmp_path, monkeypatch):
    monkeypatch.delenv(C.ENV_HOST, raising=False)
    component, payloads, files = fixture_component()
    attempts = []
    transport(monkeypatch, payloads, attempts=attempts)
    installed = C.download(component, tmp_path / 'components')
    assert C.status(component, tmp_path / 'components', verify=True)['ready']
    assert {name: (installed / name).read_bytes() for name in files} == files
    assert [url for url in attempts if 'huggingface.co' in url], 'the new host was tried first'
    assert [url for url in attempts if 'github.com' in url], 'and the published asset served it'


def test_a_failure_on_every_host_names_every_host(tmp_path, monkeypatch):
    monkeypatch.delenv(C.ENV_HOST, raising=False)
    monkeypatch.setattr(C, 'RETRY_PAUSE_SECONDS', (0, 0, 0))
    component, _, _ = fixture_component()
    attempts = []
    transport(monkeypatch, {}, attempts=attempts)
    with pytest.raises(C.PartDownloadError) as failure:
        C.download(component, tmp_path / 'components')
    assert 'huggingface.co' in str(failure.value) and 'github.com' in str(failure.value)
    assert not (tmp_path / 'components' / component['directory']).exists()


def test_the_wrong_bytes_from_the_new_host_are_refused_not_replaced(tmp_path, monkeypatch):
    """A checksum is checked once, after the hosts have been tried, so wrong bytes end the install."""
    monkeypatch.delenv(C.ENV_HOST, raising=False)
    component, payloads, _ = fixture_component()
    wrong = dict(zip(moved(component), (b'x' * len(block) for block in payloads.values())))
    attempts = []
    transport(monkeypatch, dict(payloads, **wrong), attempts=attempts)
    with pytest.raises(ValueError, match='checksum or size mismatch'):
        C.download(component, tmp_path / 'components')
    assert not (tmp_path / 'components' / component['directory']).exists()
    assert not [url for url in attempts if 'github.com' in url], 'wrong bytes are not retried elsewhere'


def test_a_component_host_setting_moves_where_the_parts_are_fetched_from(tmp_path, monkeypatch):
    monkeypatch.setenv(C.ENV_HOST, 'https://data.example.org/components')
    component, payloads, files = fixture_component()
    elsewhere = {f"https://data.example.org/components/sparse_axes/{url.rsplit('/', 1)[1]}": block
                 for url, block in payloads.items()}
    attempts = []
    transport(monkeypatch, elsewhere, attempts=attempts)
    installed = C.download(component, tmp_path / 'components')
    assert {name: (installed / name).read_bytes() for name in files} == files
    assert all('data.example.org' in url for url in attempts)
