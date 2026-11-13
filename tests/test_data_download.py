"""First-run installation verifies data and preserves working files on failures."""
import hashlib
import io
import json
import tarfile

import pytest

from PySide6 import QtWidgets

from parascape import data_manager as D, paths


def archive(tmp_path, members):
    """Create a small synthetic asset, with a descriptor independent of tar headers."""
    p = tmp_path / 'data.tar.gz'
    with tarfile.open(p, 'w:gz') as t:
        for name, payload in members.items():
            info = tarfile.TarInfo(name); info.size = len(payload)
            t.addfile(info, io.BytesIO(payload))
    descriptor = {'version': 'test', 'size': p.stat().st_size,
                  'sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
                  'files': {name: {'size': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
                            for name, data in members.items()}}
    return p, descriptor


def test_verified_install_and_missing_detection(tmp_path):
    asset, descriptor = archive(tmp_path, {'nodes.parquet': b'data', 'rnaseq/sample.json': b'{}'})
    target = tmp_path / 'installed'
    assert not D.status(target, descriptor)['ready']
    D.install_archive(asset, target, descriptor)
    assert D.status(target, descriptor)['ready']
    (target / 'nodes.parquet').unlink()
    assert D.status(target, descriptor)['missing'] == ['nodes.parquet']
    D.install_archive(asset, target, descriptor)
    assert D.status(target, descriptor)['ready']
    assert list(tmp_path.glob('installed.previous-*'))


def test_bad_checksum_never_replaces_existing_data(tmp_path):
    asset, descriptor = archive(tmp_path, {'nodes.parquet': b'data'})
    target = tmp_path / 'installed';target.mkdir();(target / 'keep').write_text('precious')
    with pytest.raises(D.DataLocationError, match='did not install'):
        D.install_archive(asset, target, descriptor)
    descriptor['sha256'] = '0' * 64
    with pytest.raises(ValueError, match='checksum'):
        D.install_archive(asset, target, descriptor)
    assert (target / 'keep').read_text() == 'precious'


@pytest.mark.parametrize('name', ['../escaped', '/absolute'])
def test_unsafe_tar_paths_are_rejected(tmp_path, name):
    asset, descriptor = archive(tmp_path, {name: b'bad'})
    with pytest.raises(ValueError, match='unsafe'):
        D.install_archive(asset, tmp_path / 'installed', descriptor)
    assert not (tmp_path / 'escaped').exists()


def test_member_checksum_and_cancellation_preserve_destination(tmp_path):
    asset, descriptor = archive(tmp_path, {'nodes.parquet': b'data'})
    target = tmp_path / 'installed'
    descriptor['files']['nodes.parquet']['sha256'] = '0' * 64
    with pytest.raises(ValueError, match='file checksum'):
        D.install_archive(asset, target, descriptor)
    assert not target.exists()
    with pytest.raises(InterruptedError):
        D.install_archive(asset, target, descriptor, cancelled=lambda: True)


def test_missing_data_offers_download_and_cancel_without_opening_window(qtbot, monkeypatch):
    from parascape import data_download_dialog as G
    monkeypatch.setattr(D, 'status', lambda: {'ready': False, 'version': 'test', 'download_bytes': 1024, 'directory': '/cache/data'})
    monkeypatch.setattr(G.DataDownloadDialog, 'exec', lambda self: G.QtWidgets.QDialog.DialogCode.Rejected)
    assert not G.ensure_data()
    dialog = G.DataDownloadDialog();qtbot.addWidget(dialog)
    assert dialog.start.text() == 'Download data'
    assert 'Hugging Face' in dialog.findChildren(G.QtWidgets.QLabel)[0].text()


def test_present_data_does_not_offer_download(monkeypatch):
    from parascape import data_download_dialog as G
    monkeypatch.setattr(D, 'status', lambda: {'ready': True})
    monkeypatch.setattr(G, 'DataDownloadDialog', lambda: pytest.fail('unexpected popup'))
    assert G.ensure_data()


def test_dialog_success_and_retry_are_responsive(qtbot, monkeypatch):
    from parascape import data_download_dialog as G
    monkeypatch.setattr(D, 'status', lambda: {'ready': False, 'version': 'test', 'download_bytes': 1024, 'directory': '/cache/data'})
    attempts = []

    def fetch(**kwargs):
        attempts.append(1)
        kwargs['progress']('Downloading data', 1, 2)
        if len(attempts) == 1:
            raise OSError('network unavailable')

    monkeypatch.setattr(D, 'download', fetch)
    dialog = G.DataDownloadDialog();qtbot.addWidget(dialog)
    dialog._start()
    qtbot.waitUntil(lambda: dialog.start.isEnabled(), timeout=5000)
    assert dialog.start.text() == 'Retry download'
    assert 'network unavailable' in dialog.message.text()
    dialog._start()
    qtbot.waitUntil(lambda: dialog.result() == G.QtWidgets.QDialog.DialogCode.Accepted, timeout=5000)


def test_dialog_cancels_worker_before_closing(qtbot, monkeypatch):
    import time
    from parascape import data_download_dialog as G
    monkeypatch.setattr(D, 'status', lambda: {'ready': False, 'version': 'test', 'download_bytes': 1024, 'directory': '/cache/data'})

    def fetch(**kwargs):
        while not kwargs['cancelled']():
            time.sleep(.01)
        raise InterruptedError('cancelled')

    monkeypatch.setattr(D, 'download', fetch)
    dialog = G.DataDownloadDialog();qtbot.addWidget(dialog)
    dialog._start()
    qtbot.waitUntil(dialog.worker.isRunning, timeout=5000)
    dialog.reject();dialog.reject()
    qtbot.waitUntil(lambda: not dialog.worker.isRunning(), timeout=5000)
    assert dialog.result() == G.QtWidgets.QDialog.DialogCode.Rejected


def test_missing_package_data_resolves_to_writable_versioned_cache(tmp_path, monkeypatch):
    monkeypatch.delenv(paths.ENV_CACHE, raising=False)
    monkeypatch.setenv(paths.ENV_STATE, str(tmp_path))
    monkeypatch.setattr(paths, '_cache_candidates', lambda: iter([]))
    assert paths.data_dir().startswith(str(tmp_path / 'data-'))


# ------------------------------------------------------- the three defects Einar hit on 2026-11-05
"""He ran it, the bar reached 100 %, nothing happened; he retried and it sat at 0 %. All three causes
are below. The download itself was never at fault -- the archive on disk was complete and its sha256
matched -- so every fix here is about when the refusal happens and whether anyone can see it."""


def test_a_populated_unmanaged_destination_is_refused_before_anything_is_downloaded(tmp_path):
    """The defect that cost 898 MB. `install_archive` always refused to overwrite somebody else's
    files, correctly -- but it refused AFTER the whole archive had been fetched and hashed."""
    target = tmp_path / 'checkout_data'
    target.mkdir()
    (target / 'nodes.parquet').write_bytes(b'committed in the repository')
    with pytest.raises(D.DataLocationError) as caught:
        D.check_destination(target)
    assert 'did not install' in str(caught.value)
    assert 'PARASCAPE_CACHE' in str(caught.value), 'the message must say what to do about it'


def test_a_managed_destination_is_accepted(tmp_path):
    target = tmp_path / 'managed'
    target.mkdir()
    (target / 'data_release.json').write_text('{}')
    (target / 'nodes.parquet').write_bytes(b'x')
    assert D.check_destination(target) == target.absolute()


def test_an_empty_or_absent_destination_is_accepted(tmp_path):
    assert D.check_destination(tmp_path / 'absent')
    (tmp_path / 'empty').mkdir()
    assert D.check_destination(tmp_path / 'empty')


def test_the_location_refusal_is_a_value_error_so_old_callers_still_catch_it():
    assert issubclass(D.DataLocationError, ValueError)


def test_hashing_reports_progress(tmp_path):
    """Why a finished download read 0 %: hashing emitted nothing, and on a complete part file the
    download loop -- the only thing that reported progress -- was skipped entirely."""
    blob = tmp_path / 'blob'
    blob.write_bytes(b'x' * (3 * 1024 * 1024))
    seen = []
    D._hash(blob, lambda phase, done, total: seen.append((phase, done, total)))
    assert seen, 'hashing 898 MB in silence is what looked like a hang'
    assert seen[-1][1] == seen[-1][2] == blob.stat().st_size
    assert all(p for p, _, _ in seen), 'every report names its phase'


def test_the_dialog_does_not_offer_a_retry_that_cannot_work(qtbot, monkeypatch):
    """The confusing half: the refusal appeared in a small label under a button saying
    "Retry download", so the only available action was one that could never succeed."""
    from parascape import data_download_dialog as G
    monkeypatch.setattr(D, 'status', lambda: {'ready': False, 'version': 'test',
                                              'download_bytes': 1024, 'directory': '/cache/data'})
    monkeypatch.setattr(QtWidgets.QMessageBox, 'warning', staticmethod(lambda *a, **k: None))

    def refuse(**kwargs):
        raise D.DataLocationError('/repo/parascape/data already holds files that ParaScape '
                                  'did not install. Point PARASCAPE_CACHE at an empty folder.')

    monkeypatch.setattr(D, 'download', refuse)
    dialog = G.DataDownloadDialog()
    qtbot.addWidget(dialog)
    dialog._start()
    qtbot.waitUntil(lambda: dialog.start.text() != 'Download data', timeout=5000)
    assert not dialog.start.isEnabled(), 'a non-retryable failure must not offer a retry'
    assert dialog.start.text() == 'Cannot download here'
    assert 'PARASCAPE_CACHE' in dialog.message.text()


def test_an_ordinary_network_failure_still_offers_a_retry(qtbot, monkeypatch):
    """The distinction has to cut both ways, or a transient failure becomes a dead end."""
    from parascape import data_download_dialog as G
    monkeypatch.setattr(D, 'status', lambda: {'ready': False, 'version': 'test',
                                              'download_bytes': 1024, 'directory': '/cache/data'})

    monkeypatch.setattr(D, 'download', lambda **k: (_ for _ in ()).throw(OSError('network unavailable')))
    dialog = G.DataDownloadDialog()
    qtbot.addWidget(dialog)
    dialog._start()
    qtbot.waitUntil(lambda: dialog.start.isEnabled(), timeout=5000)
    assert dialog.start.text() == 'Retry download'


# ------------------------------------------------- the data moved to Hugging Face (2026-11-13)
"""Einar: "all of the data that parascape uses should be hosted on huggingface in my account.
otherwise the space limits on github will fill quickly." The archive is the same archive and the
checksum contract is unchanged; what changed is where it is fetched from, that a second host is
tried when the first does not answer, and that the host is a setting. Nothing here touches the
network: the fetch is replaced by a stub that records the addresses it was asked for."""


def descriptor_for(tmp_path):
    asset, descriptor = archive(tmp_path, {'nodes.parquet': b'data'})
    descriptor.update(version='9.9.9', asset='parascape-data-9.9.9.tar.gz',
                      url=f'{D.HOST}/9.9.9/parascape-data-9.9.9.tar.gz',
                      fallback_url=f'{D.GITHUB_RELEASES}/v9.9.9/parascape-data-9.9.9.tar.gz')
    return asset, descriptor


def test_hugging_face_is_tried_first_and_the_github_release_second(tmp_path, monkeypatch):
    monkeypatch.delenv(D.ENV_HOST, raising=False)
    monkeypatch.delenv(D.ENV_URL, raising=False)
    _, descriptor = descriptor_for(tmp_path)
    assert D.sources(descriptor) == [descriptor['url'], descriptor['fallback_url']]
    assert D.sources(descriptor)[0].startswith('https://huggingface.co/datasets/einarolafsson/')
    assert D.host_name(descriptor) == 'Hugging Face'


def test_the_shipped_descriptor_points_at_hugging_face_with_github_behind_it(monkeypatch):
    monkeypatch.delenv(D.ENV_HOST, raising=False)
    monkeypatch.delenv(D.ENV_URL, raising=False)
    description = D.manifest()
    first, second = D.sources(description)
    assert first == f'{D.HOST}/{description["version"]}/{description["asset"]}'
    assert second == f'{D.GITHUB_RELEASES}/v{description["version"]}/{description["asset"]}'


def test_the_host_is_a_setting_not_an_edit(tmp_path, monkeypatch):
    _, descriptor = descriptor_for(tmp_path)
    monkeypatch.setenv(D.ENV_HOST, 'https://data.example.org/parascape/')
    assert D.sources(descriptor)[0] == 'https://data.example.org/parascape/9.9.9/parascape-data-9.9.9.tar.gz'
    assert D.host_name(descriptor) == 'data.example.org'
    monkeypatch.setenv(D.ENV_URL, 'https://data.example.org/exactly/this.tar.gz')
    assert D.sources(descriptor) == ['https://data.example.org/exactly/this.tar.gz']


@pytest.mark.parametrize('bad', ['http://huggingface.co/x', 'https://elsewhere.example/x',
                                 'https://huggingface.co/a/../b', 'https://huggingface.co/x?token=1'])
def test_an_address_off_the_known_hosts_is_refused_before_any_request(tmp_path, monkeypatch, bad):
    monkeypatch.delenv(D.ENV_HOST, raising=False)
    monkeypatch.setenv(D.ENV_URL, bad)
    _, descriptor = descriptor_for(tmp_path)
    with pytest.raises(ValueError, match='Unexpected data download host'):
        D.sources(descriptor)


def test_an_unanswering_first_host_falls_through_to_the_second_and_still_verifies(tmp_path, monkeypatch):
    """The reason the GitHub release stays published. The stub serves the real archive only from the
    fallback address, and the install must succeed with the pinned checksums unchanged."""
    import urllib.error
    asset, descriptor = descriptor_for(tmp_path)
    monkeypatch.delenv(D.ENV_HOST, raising=False)
    monkeypatch.delenv(D.ENV_URL, raising=False)
    asked = []

    def stub(url, part, description, progress=None, cancelled=None):
        asked.append(url)
        if url != description['fallback_url']:
            raise urllib.error.HTTPError(url, 404, 'Not Found', {}, None)
        part.write_bytes(asset.read_bytes())

    monkeypatch.setattr(D, '_fetch', stub)
    target = tmp_path / 'installed'
    assert D.download(target, descriptor) == target.absolute()
    assert asked == [descriptor['url'], descriptor['fallback_url']]
    assert D.status(target, descriptor)['ready']


def test_a_failure_on_every_host_names_every_host(tmp_path, monkeypatch):
    asset, descriptor = descriptor_for(tmp_path)
    monkeypatch.delenv(D.ENV_HOST, raising=False)
    monkeypatch.delenv(D.ENV_URL, raising=False)
    monkeypatch.setattr(D, '_fetch', lambda *a, **k: (_ for _ in ()).throw(OSError('refused')))
    with pytest.raises(OSError) as caught:
        D.download(tmp_path / 'installed', descriptor)
    assert 'huggingface.co' in str(caught.value) and 'github.com' in str(caught.value)


def test_a_hosted_archive_that_is_the_wrong_bytes_is_still_refused(tmp_path, monkeypatch):
    """Changing host changes where the bytes come from and nothing about what is accepted."""
    _, descriptor = descriptor_for(tmp_path)
    monkeypatch.delenv(D.ENV_HOST, raising=False)
    monkeypatch.delenv(D.ENV_URL, raising=False)
    monkeypatch.setattr(D, '_fetch', lambda url, part, *a, **k: part.write_bytes(b'x' * descriptor['size']))
    with pytest.raises(ValueError, match='checksum mismatch'):
        D.download(tmp_path / 'installed', descriptor)
    assert not (tmp_path / 'installed').exists()


def test_an_early_end_is_a_failure_rather_than_a_truncated_install(tmp_path, monkeypatch):
    _, descriptor = descriptor_for(tmp_path)

    class Response:
        status, headers = 200, {}

        def read(self, size):
            return b''

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(D.urllib.request, 'urlopen', lambda *a, **k: Response())
    with pytest.raises(OSError, match='ended early'):
        D._fetch(descriptor['url'], tmp_path / 'part', descriptor)


def test_the_download_dialog_names_the_host_it_will_use(qtbot, monkeypatch):
    from parascape import data_download_dialog as G
    monkeypatch.setattr(D, 'status', lambda: {'ready': False, 'version': 'test', 'host': 'Hugging Face',
                                              'download_bytes': 1024, 'directory': '/cache/data'})
    dialog = G.DataDownloadDialog();qtbot.addWidget(dialog)
    assert 'Hugging Face' in dialog.findChildren(G.QtWidgets.QLabel)[0].text()
