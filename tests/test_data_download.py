"""First-run installation verifies data and preserves working files on failures."""
import hashlib
import io
import json
import tarfile

import pytest

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
    with pytest.raises(ValueError, match='unmanaged'):
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
    assert 'GitHub' in dialog.findChildren(G.QtWidgets.QLabel)[0].text()


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
