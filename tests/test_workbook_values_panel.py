"""Actual background workbook browsing, treatment and complete export routes."""
import json

from parascape.workbook_values_panel import WorkbookValuesPanel
from test_workbook_values import fixture_collection


def complete(qtbot,panel):
    qtbot.waitUntil(lambda:panel._worker is None,timeout=10000)


def test_actual_rows_context_export_and_stale_result_refusal(qtbot,tmp_path):
    fixture_collection(tmp_path/'collection')
    panel=WorkbookValuesPanel(tmp_path/'collection');qtbot.addWidget(panel)
    assert panel.pack.count()==1 and 'hidden' in panel.sheet.currentText()
    panel.labels.setPlainText('001');panel.show_values();complete(qtbot,panel)
    assert panel.result is not None,panel.summary.text()
    assert panel.table.rowCount()==7
    assert panel.table.item(1,4).text()=='9007199254740993'
    panel.table.setCurrentCell(3,0);panel.show_context()
    assert panel.context_dialog.isVisible()
    assert panel.context_record['selected_source_cell']['formula_xml']=='<f>A2/3</f>'
    assert panel.save(str(tmp_path/'gui.csv'))==7
    assert json.loads((tmp_path/'gui.csv.json').read_text())['exported_physical_cells']==7
    panel.labels.setPlainText('absent')
    assert panel.result is None and panel.context_dialog is None and panel.table.rowCount()==0
    assert not panel.export_button.isEnabled() and panel.save(str(tmp_path/'stale.csv'))==0


def test_source_filter_and_cancelled_generation_discard_late_results(qtbot,tmp_path):
    fixture_collection(tmp_path/'collection',verified=False)
    panel=WorkbookValuesPanel(tmp_path/'collection');qtbot.addWidget(panel)
    panel.show_values();assert panel._worker is not None
    panel.labels.setPlainText('Other');complete(qtbot,panel)
    assert panel.result is None and not panel.export_button.isEnabled()
    panel.show_values();complete(qtbot,panel)
    assert panel.result is not None and 'verification incomplete' in panel.summary.text()
    assert panel.table.rowCount()==2
    panel.source_filter.setText('nonexistent')
    assert panel.pack.count()==0 and panel.sheet.count()==0 and panel.result is None
    assert not panel.show_button.isEnabled()
    panel.source_filter.setText('GSE_FIXTURE')
    assert panel.pack.count()==1


def test_actual_missing_workbook_dialog_installs_and_then_opens_offline(qtbot,tmp_path,monkeypatch):
    from PyQt6 import QtCore
    from parascape.component_download_dialog import ComponentDownloadDialog
    from test_optional_workbook_data import portable_collection
    _,metadata,_,calls=portable_collection(tmp_path,monkeypatch)
    panel=WorkbookValuesPanel(metadata.parent);qtbot.addWidget(panel)
    panel.labels.setPlainText('001')
    def start_dialog():
        dialogs=panel.findChildren(ComponentDownloadDialog)
        assert len(dialogs)==1
        dialogs[0].start.click()
    QtCore.QTimer.singleShot(50,start_dialog)
    panel.show_values();complete(qtbot,panel)
    assert panel.result is not None,panel.summary.text()
    assert panel.table.item(1,4).text()=='9007199254740993'
    transferred=len(calls)
    panel.labels.setPlainText('Other');panel.show_values();complete(qtbot,panel)
    assert panel.table.rowCount()==2 and len(calls)==transferred
