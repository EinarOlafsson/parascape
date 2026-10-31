"""Download prompts remain responsive and cannot expose cancelled datasets."""
import threading

from PyQt6 import QtWidgets
import pytest

from parascape import component_download_dialog as D, sparse_values_panel as P
from parascape import data_components as C
from test_optional_sparse_data import hosted  # noqa: F401
from test_data_components import transport


def test_dialog_real_installer_success_and_retry(qtbot,hosted,monkeypatch):
    calls=transport(monkeypatch,hosted['payloads']);components=list(hosted['index']['components'].values())
    original=C.download;attempts=[]
    def fail_once(*args,**kwargs):
        attempts.append(True)
        if len(attempts)==1:raise OSError('temporary failure')
        return original(*args,**kwargs)
    monkeypatch.setattr(C,'download',fail_once)
    dialog=D.ComponentDownloadDialog(components);qtbot.addWidget(dialog);dialog.show();dialog._start()
    qtbot.waitUntil(lambda:dialog.start.isEnabled(),timeout=10000)
    assert 'temporary failure' in dialog.message.text() and dialog.start.text()=='Retry download'
    dialog._start()
    qtbot.waitUntil(lambda:dialog.result()==QtWidgets.QDialog.DialogCode.Accepted,timeout=10000)
    assert calls and all(C.status(component)['ready'] for component in components)
    assert not dialog.worker.isRunning()


def test_dialog_cancel_waits_for_worker_and_keeps_prompt_closed(qtbot,hosted,monkeypatch):
    entered=threading.Event()
    def cancellable(component,progress,cancelled):
        entered.set()
        assert entered.wait(1)
        while not cancelled():
            threading.Event().wait(.01)
        raise InterruptedError('cancelled safely')
    monkeypatch.setattr(C,'download',cancellable)
    dialog=D.ComponentDownloadDialog(list(hosted['index']['components'].values()));qtbot.addWidget(dialog);dialog.show();dialog._start()
    qtbot.waitUntil(entered.is_set,timeout=5000);dialog.reject()
    qtbot.waitUntil(lambda:not dialog.worker.isRunning() and not dialog.isVisible(),timeout=5000)
    assert dialog.result()==QtWidgets.QDialog.DialogCode.Rejected
    assert not C.storage_root().exists()


def test_show_downloads_missing_dependencies_then_reopens_offline(qtbot,hosted,monkeypatch):
    calls=transport(monkeypatch,hosted['payloads']);prompts=[]
    class Accepted:
        def __init__(self,components,parent):self.components=components;prompts.append(components)
        def exec(self):
            for component in self.components:C.download(component)
            return QtWidgets.QDialog.DialogCode.Accepted
    monkeypatch.setattr(P,'ComponentDownloadDialog',Accepted)
    panel=P.SparseValuesPanel(hosted['roots']['sparse_datasets']);qtbot.addWidget(panel);panel.features.setPlainText('001');panel.show_values()
    qtbot.waitUntil(lambda:panel._worker is None,timeout=10000)
    assert panel.result is not None,panel.summary.text()
    assert len(prompts)==1 and len(prompts[0])==3 and panel.table.rowCount()==4
    assert panel.download_button.isEnabled()
    downloads=len(calls);panel.show_values();qtbot.waitUntil(lambda:panel._worker is None,timeout=10000)
    assert panel.result is not None and len(prompts)==1 and len(calls)==downloads
    # Verify action repairs an owned corrupted component and invalidates the prior view.
    component=prompts[0][0];file=C.storage_root()/component['directory']/'entries.parquet';file.write_bytes(file.read_bytes()+b'changed')
    panel.download_selected()
    assert panel.result is None and C.status(component,verify=True)['ready']


def test_declining_download_does_not_start_query_or_write_data(qtbot,hosted,monkeypatch):
    class Rejected:
        def __init__(self,*args):pass
        def exec(self):return QtWidgets.QDialog.DialogCode.Rejected
    monkeypatch.setattr(P,'ComponentDownloadDialog',Rejected)
    panel=P.SparseValuesPanel(hosted['roots']['sparse_datasets']);qtbot.addWidget(panel);panel.features.setPlainText('001');panel.show_values()
    assert panel._worker is None and panel.result is None and not C.storage_root().exists()
    assert not panel.export_button.isEnabled()
