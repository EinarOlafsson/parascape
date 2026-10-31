"""QC before payload download, exact filters and original formula inspection."""
from parascape.workbook_values_panel import WorkbookValuesPanel
from test_workbook_qc import biff_fixture,wrapper,install_qc
from test_workbook_values_panel import complete


def test_verified_source_with_reader_difference_is_not_hidden_by_default(qtbot,tmp_path):
    folder=biff_fixture(tmp_path/'collection');install_qc(folder,wrapper(folder,differences=1))
    panel=WorkbookValuesPanel(folder.parent);qtbot.addWidget(panel)
    assert panel.pack.count()==1
    panel.show_qc();assert panel.qc_dialog.isVisible()
    panel.require_verified.setChecked(True)
    assert panel.pack.count()==1 and panel.qc_dialog is None
    panel.require_agreement.setChecked(True)
    assert panel.pack.count()==0 and not panel.show_button.isEnabled()
    panel.require_agreement.setChecked(False)
    panel.labels.setPlainText('001');panel.show_values();complete(qtbot,panel)
    assert panel.result is not None,panel.summary.text()
    assert '1 reader text disagreements' in panel.summary.text()
    assert panel.table.item(3,3).text()=='Formula cache / error'
    assert 'BIFF tokens' in panel.table.item(3,5).text()
    panel.table.setCurrentCell(3,0);panel.show_context()
    assert panel.context_record['selected_source_formula']['cached_result_code']==29
    assert panel.context_record['selected_source_formula']['formula_evaluated'] is False


def test_qc_can_be_inspected_before_optional_values_are_present(qtbot,tmp_path,monkeypatch):
    from test_optional_workbook_data import portable_collection
    source,metadata,_,calls=portable_collection(tmp_path,monkeypatch,legacy=True)
    value=wrapper(source,formula=False)
    value['raw_report']['physical_cell_kinds']=dict(text=6,numeric=4,blank=1)
    value['raw_report']['original_error_code_counts']={}
    install_qc(metadata,value)
    panel=WorkbookValuesPanel(metadata.parent);qtbot.addWidget(panel)
    panel.show_qc()
    assert panel.qc_dialog.isVisible() and not calls
    assert not list(metadata.glob('*.parquet'))


def test_unknown_sources_are_visible_until_user_requests_verification(qtbot,tmp_path):
    folder=biff_fixture(tmp_path/'collection')
    panel=WorkbookValuesPanel(folder.parent);qtbot.addWidget(panel)
    assert panel.pack.count()==1
    panel.require_verified.setChecked(True);assert panel.pack.count()==0
    panel.require_verified.setChecked(False);assert panel.pack.count()==1
