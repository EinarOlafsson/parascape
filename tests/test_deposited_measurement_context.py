"""Selected measurements retain row/column identities and original source declarations through sorting and edits."""
import hashlib
import json

import pytest
from PyQt6 import QtCore, QtWidgets

from parascape.deposited_values_panel import DepositedValuesPanel
from test_deposited_values import pack  # noqa: F401
from test_protein_identifier_evidence import evidence  # noqa: F401


@pytest.fixture
def annotated_pack(pack):
    metadata = dict(row_annotations=[dict(source_row=i, feature_label=label,
        source_fields={'R.FileName': 'Same.raw', 'R.Condition': condition, 'R.Replicate': '001',
                       'PG.Qvalue': qvalue}) for i, (label, condition, qvalue) in enumerate([
            ('001', 'Original condition A', '0.001'), ('001', 'Original condition B', '0.6'),
            ('CaseID', '', '')])])
    (pack / 'source_metadata.json').write_text(json.dumps(metadata))
    info = json.loads((pack / 'manifest.json').read_text())
    info['outputs']['source_metadata.json'] = hashlib.sha256((pack / 'source_metadata.json').read_bytes()).hexdigest()
    (pack / 'manifest.json').write_text(json.dumps(info))
    receipt = json.loads((pack / 'validation.json').read_text())
    receipt['manifest_sha256'] = hashlib.sha256((pack / 'manifest.json').read_bytes()).hexdigest()
    (pack / 'validation.json').write_text(json.dumps(receipt))
    return pack


def fields(tree):
    result = []
    iterator = QtWidgets.QTreeWidgetItemIterator(tree)
    while iterator.value() is not None:
        item = iterator.value()
        result.append((item.text(0), item.text(1)))
        iterator += 1
    return result


def test_repeated_label_details_select_original_row_and_column_after_sorting(qtbot, annotated_pack):
    panel = DepositedValuesPanel(annotated_pack.parent)
    qtbot.addWidget(panel)
    panel.features.setPlainText('001')
    panel.show_values()
    panel.table.setSortingEnabled(True)
    panel.table.sortItems(4, QtCore.Qt.SortOrder.DescendingOrder)
    wanted = None
    for visible in range(panel.table.rowCount()):
        panel.table.setCurrentCell(visible, 0)
        selected = panel._selected_measurement()
        if selected['source_row'] == 1 and selected['source_sample'] == 'Original 02':
            wanted = selected
            break
    assert wanted is not None and panel.measurement_button.isEnabled()
    panel.show_measurement_context()
    displayed = fields(panel.measurement_tree)
    assert ('R.Condition', 'Original condition B') in displayed
    assert ('R.Condition', 'Original condition A') not in displayed
    assert ('R.FileName', 'Same.raw') in displayed
    assert ('R.Replicate', '001') in displayed
    assert ('PG.Qvalue', '0.6') in displayed
    assert ('Original quantity column', 'Original 02') in displayed
    assert ('Original source token', '-2.500') in displayed
    assert ('Source row ordinal', '1') in displayed
    assert ('Numeric value', '-2.5') in displayed
    source_fields = panel.measurement_tree.findItems('source_fields',
        QtCore.Qt.MatchFlag.MatchExactly | QtCore.Qt.MatchFlag.MatchRecursive)
    assert len(source_fields) == 1 and source_fields[0].isExpanded()
    assert not panel.measurement_tree.editTriggers() & QtWidgets.QAbstractItemView.EditTrigger.DoubleClicked


def test_missing_numeric_value_and_original_empty_token_remain_distinct(qtbot, annotated_pack):
    panel = DepositedValuesPanel(annotated_pack.parent)
    qtbot.addWidget(panel)
    panel.features.setPlainText('CaseID')
    panel.show_values()
    panel.table.setCurrentCell(0, 0)
    panel.show_measurement_context()
    displayed = fields(panel.measurement_tree)
    assert ('Original source token', '') in displayed
    assert ('Numeric value', 'NaN') in displayed
    assert ('R.Condition', '') in displayed
    # All repeated submitted column-factor occurrences remain inspectable.
    assert ('batch', 'A') in displayed and ('batch', 'B') in displayed


@pytest.mark.parametrize('edit', ['query', 'source_filter', 'namespace', 'flag_rule', 'refresh', 'selection'])
def test_details_close_on_context_changes_and_do_not_modify_exports(qtbot, annotated_pack, tmp_path, edit):
    panel = DepositedValuesPanel(annotated_pack.parent)
    qtbot.addWidget(panel)
    panel.features.setPlainText('001')
    panel.show_values()
    original = panel.result.copy(deep=True)
    treatment = json.dumps(panel.treatment_record, sort_keys=True)
    panel.table.setCurrentCell(0, 0)
    panel.show_measurement_context()
    assert panel.measurement_dialog is not None
    assert panel.result.equals(original) and json.dumps(panel.treatment_record, sort_keys=True) == treatment
    if edit == 'query':
        panel.features.setPlainText('CaseID')
    elif edit == 'source_filter':
        panel.source_filter.setText('no matching source')
    elif edit == 'namespace':
        panel.identifier_kind.currentIndexChanged.emit(0)
    elif edit == 'flag_rule':
        panel.exclude_source_marked.setChecked(True)
    elif edit == 'refresh':
        panel.refresh()
    else:
        panel.table.setCurrentCell(1, 0)
    assert panel.measurement_dialog is None and panel.measurement_tree is None
    if edit != 'selection':
        assert panel.result is None and not panel.measurement_button.isEnabled()
        assert panel.save(str(tmp_path / 'stale.csv')) == 0
    else:
        assert panel.result.equals(original) and panel.measurement_button.isEnabled()


def test_no_selection_and_absent_row_annotations_do_not_create_inferred_details(qtbot, pack):
    panel = DepositedValuesPanel(pack.parent)
    qtbot.addWidget(panel)
    panel.show_measurement_context()
    assert getattr(panel, 'measurement_dialog', None) is None and not panel.measurement_button.isEnabled()
    panel.features.setPlainText('001')
    panel.show_values()
    panel.table.setCurrentCell(0, 0)
    panel.show_measurement_context()
    assert ('Original row annotations', 'No declaration supplied for this measurement') in fields(panel.measurement_tree)
    panel.table.clearSelection()
    panel.table.setCurrentCell(-1, -1)
    assert panel.measurement_dialog is None and not panel.measurement_button.isEnabled()


def test_qualified_identifier_details_retain_all_evidence_occurrences_for_the_selected_row(qtbot, evidence):
    pack, _, _ = evidence
    panel = DepositedValuesPanel(pack.parent)
    qtbot.addWidget(panel)
    panel.identifier_kind.setCurrentIndex(panel.identifier_kind.findData('submitted_gene_label'))
    panel.features.setPlainText('G001')
    panel.show_values()
    panel.table.setCurrentCell(0, 0)
    panel.show_measurement_context()
    displayed = fields(panel.measurement_tree)
    assert displayed.count(('query_literal', 'G001')) == 2
    assert ('query_namespace', 'submitted_gene_label') in displayed
    assert ('Original feature label', '001') in displayed
    assert ('Source row ordinal', '0') in displayed
