"""Background UI queries cannot leave stale sources or truncated exports."""
import csv
import json

import pytest

from parascape.sparse_values_panel import SparseValuesPanel
from test_sparse_values import fixture_collection


def panel_for(qtbot, roots):
    panel = SparseValuesPanel(roots['sparse_datasets'], coordinate_root=roots['sparse_coordinates'], axis_root=roots['sparse_axes'])
    qtbot.addWidget(panel)
    return panel


def complete(qtbot, panel):
    qtbot.waitUntil(lambda: panel._worker is None, timeout=10000)


@pytest.mark.parametrize('inverse', [False, True])
def test_actual_background_query_context_and_complete_export(qtbot, tmp_path, inverse):
    roots, _ = fixture_collection(tmp_path, inverse=inverse)
    panel = panel_for(qtbot, roots)
    panel.features.setPlainText('001')
    panel.show_values()
    complete(qtbot, panel)
    assert panel.result is not None, panel.summary.text()
    assert panel.table.rowCount() == 4
    assert panel.table.item(0, 4).text() == '001'
    assert panel.table.item(0, 6).text() == 'bc-1'
    assert panel.table.item(1, 7).text() == str(2**53 + 1)
    assert 'submitted-genome' in panel.table.item(0, 4).toolTip()
    panel.table.setCurrentCell(2, 0)
    assert panel.context_button.isEnabled()
    panel.show_context()
    assert panel.context_dialog.isVisible()
    destination = tmp_path / 'gui.csv'
    assert panel.save(str(destination)) == 4
    with destination.open() as stream:
        rows = list(csv.DictReader(stream))
    assert rows[1]['source_token'] == str(2**53 + 1)
    sidecar = json.loads(destination.with_suffix('.csv.json').read_text())
    assert sidecar['operations'][0]['source_coordinate'] == ('source_column' if inverse else 'source_row')
    panel.features.setPlainText('absent')
    assert panel.table.rowCount() == 0 and panel.result is None
    assert panel.context_dialog is None and not panel.export_button.isEnabled()
    assert panel.save(str(tmp_path / 'stale.csv')) == 0


def test_editing_pending_query_discards_any_late_result(qtbot, tmp_path):
    roots, _ = fixture_collection(tmp_path)
    panel = panel_for(qtbot, roots)
    panel.features.setPlainText('001')
    panel.show_values()
    assert panel._worker is not None
    panel.features.setPlainText('ABC')
    complete(qtbot, panel)
    assert panel.result is None and panel.table.rowCount() == 0
    assert not panel.export_button.isEnabled()
    panel.show_values()
    complete(qtbot, panel)
    assert panel.result is not None and panel.result.record['exported_records'] == 0
    assert len(panel.result.feature_matches) == 1
    assert panel.export_button.isEnabled()  # Empty result still has a treatment record.


def test_filter_and_failed_integrity_clear_exports(qtbot, tmp_path):
    roots, _ = fixture_collection(tmp_path)
    panel = panel_for(qtbot, roots)
    panel.features.setPlainText('001')
    panel.show_values()
    complete(qtbot, panel)
    assert panel.result is not None
    panel.source_filter.setText('missing')
    assert panel.pack.count() == 0 and panel.result is None
    assert not panel.show_button.isEnabled() and not panel.export_button.isEnabled()
    panel.source_filter.clear()
    output = roots['sparse_coordinates'] / ('c' * 64) / 'entries.parquet'
    output.write_bytes(output.read_bytes() + b'changed')
    panel.show_values()
    complete(qtbot, panel)
    assert panel.result is None and 'changed' in panel.summary.text()
    assert not panel.export_button.isEnabled()


def test_compare_install_exposes_sparse_tab(qtbot, monkeypatch, tmp_path):
    from PyQt6 import QtCore, QtWidgets
    from parascape import compare_panel, rna_panel, rna_design_panel, published_expression_panel
    from parascape import dataset_coverage_panel, source_matrices_panel, deposited_values_panel, reference_identifiers_panel

    class Profiles(QtWidgets.QWidget):
        gene_chosen = QtCore.pyqtSignal(str)

    class Window(QtWidgets.QMainWindow):
        def _workflow_gene(self, gene):
            pass

    monkeypatch.setattr(compare_panel, 'catalogue_titles', lambda: {})
    monkeypatch.setattr(rna_panel, 'RNAProfilesPanel', Profiles)
    for module, name in [(rna_design_panel, 'StudyDesignPanel'), (published_expression_panel, 'PublishedExpressionPanel'),
                         (dataset_coverage_panel, 'DatasetCoveragePanel'), (source_matrices_panel, 'SourceMatricesPanel'),
                         (deposited_values_panel, 'DepositedValuesPanel'), (reference_identifiers_panel, 'ReferenceIdentifiersPanel')]:
        monkeypatch.setattr(module, name, QtWidgets.QWidget)
    monkeypatch.setenv('PARASCAPE_SPARSE_DATASETS', str(tmp_path / 'not-installed'))
    window = Window()
    qtbot.addWidget(window)
    compare_panel.install(window)
    tabs = window.compare_dock.widget()
    index = [tabs.tabText(i) for i in range(tabs.count())].index('Sparse measurements')
    assert isinstance(tabs.widget(index).widget(), SparseValuesPanel)
    assert not tabs.widget(index).widget().show_button.isEnabled()
