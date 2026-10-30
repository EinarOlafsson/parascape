"""Literal source flags remain distinct from biological quality and every exclusion has a source-bound record."""
import hashlib
import json

import pandas as pd
import pytest

from parascape import deposited_values as D
from parascape import source_flags as F
from parascape.deposited_values_panel import DepositedValuesPanel
from test_deposited_values import pack  # noqa: F401
from test_protein_identifier_evidence import evidence  # noqa: F401


def digest(path):
    """Pin fixture bytes independently of the view implementation."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def annotations(pack, *, mark_both=False):
    """A repeated group label has different row flags; unfamiliar whitespace is not silently interpreted."""
    metadata_path = pack / 'source_metadata.json'
    metadata = json.loads(metadata_path.read_text())
    metadata['row_annotations'] = [
        dict(source_row=0, feature_label='001', source_fields={'Reverse': '+', 'Potential contaminant': ''}),
        dict(source_row=1, feature_label='001', source_fields={'Reverse': '+' if mark_both else '', 'Potential contaminant': ''}),
        dict(source_row=2, feature_label='CaseID', source_fields={'Reverse': ' + ', 'Contaminant': ''})]
    metadata_path.write_text(json.dumps(metadata))
    info = json.loads((pack / 'manifest.json').read_text())
    info['assay'] = 'deposited protein-group quantities'
    info['outputs']['source_metadata.json'] = digest(metadata_path)
    (pack / 'manifest.json').write_text(json.dumps(info))
    receipt = json.loads((pack / 'validation.json').read_text())
    receipt['manifest_sha256'] = digest(pack / 'manifest.json')
    (pack / 'validation.json').write_text(json.dumps(receipt))


def test_flags_annotate_without_changing_values_and_exclude_by_original_row_not_group_label(pack):
    annotations(pack)
    raw = D.profiles(pack, ['001', 'CaseID'])
    before = digest(pack / 'values.parquet')
    shown, record = F.apply(raw)
    pd.testing.assert_frame_equal(shown[raw.columns], raw, check_exact=True)
    assert shown.value.to_numpy().tobytes() == raw.value.to_numpy().tobytes()
    assert record['excluded_source_rows'] == 0 and record['input_source_rows'] == 3
    selected, record = F.apply(raw, exclude_marked=True)
    assert set(selected.source_row) == {1, 2} and len(selected) == 4
    assert record['excluded_source_rows'] == 1 and record['excluded_measurements'] == 2
    assert record['uninterpreted_source_rows'] == 1
    assert selected.loc[selected.source_row.eq(2), 'value'].isna().sum() == 1
    assert selected.source_token.tolist() == ['1.2500', '', '-2.500', '2.5']
    assert record['excluded_rows'][0]['feature_label'] == '001'
    assert record['excluded_rows'][0]['source_marked_fields'] == ['Reverse']
    assert digest(pack / 'values.parquet') == before


def test_missing_flags_other_assays_and_unknown_tokens_are_not_called_clean_or_excluded(pack):
    raw = D.profiles(pack, '001')
    shown, record = F.apply(raw, exclude_marked=True)
    assert len(shown) == len(raw) and shown.source_flag_status.eq('flags_not_declared').all()
    annotations(pack)
    raw = D.profiles(pack, '001')
    raw['assay'] = 'other assay with unrelated fields'
    shown, record = F.apply(raw, exclude_marked=True)
    assert len(shown) == len(raw) and shown.source_flag_status.eq('not_interpreted_for_assay').all()
    assert record['excluded_source_rows'] == 0


def test_manifest_identity_prevents_flags_bleeding_between_sources_and_inconsistent_rows_are_refused(pack):
    annotations(pack)
    raw = D.profiles(pack, '001')
    other = raw[raw.source_row.eq(0)].copy()
    other['manifest_sha256'] = 'b' * 64
    metadata = json.loads(other.source_row_metadata.iloc[0])
    metadata['source_fields']['Reverse'] = ''
    other['source_row_metadata'] = json.dumps(metadata)
    combined = pd.concat([raw, other], ignore_index=True)
    selected, record = F.apply(combined, exclude_marked=True)
    assert len(selected) == 4 and set(selected.manifest_sha256) == {raw.manifest_sha256.iloc[0], 'b' * 64}
    assert record['input_source_rows'] == 3 and record['excluded_source_rows'] == 1
    broken = raw.copy()
    broken.loc[broken.index[-2], 'source_row_metadata'] = json.dumps(metadata)
    with pytest.raises(ValueError, match='Repeated source-row flag declarations'):
        F.apply(broken)


def test_panel_exports_sidecar_and_invalidates_treatment_when_filter_changes(qtbot, pack, tmp_path):
    annotations(pack)
    before = {name: digest(pack / name) for name in ['values.parquet', 'source_tokens.parquet', 'source_metadata.json']}
    panel = DepositedValuesPanel(pack.parent)
    qtbot.addWidget(panel)
    panel.features.setPlainText('001')
    panel.show_values()
    assert panel.table.columnCount() == 7 and panel.table.rowCount() == 4
    assert not panel.exclude_source_marked.isChecked()
    panel.exclude_source_marked.setChecked(True)
    assert panel.result is None and panel.treatment_record is None
    assert panel.save_treatment(str(tmp_path / 'stale.json')) == 0
    panel.show_values()
    assert panel.table.rowCount() == 2 and set(panel.result.source_row) == {1}
    assert 'Labels not matched' not in panel.summary.text()
    output = tmp_path / 'selected.csv'
    assert panel.save(str(output)) == 2
    sidecars = list(tmp_path.glob('selected.csv.treatment-*.json'))
    assert len(sidecars) == 1
    saved = json.loads(sidecars[0].read_text())
    assert saved['exported_csv_sha256'] == digest(output)
    assert saved['manifest_sha256'] == digest(pack / 'manifest.json')
    assert saved['source_flag_view']['excluded_source_rows'] == 1
    assert not saved['original_files_modified'] and not saved['numerical_values_transformed']
    assert not saved['biological_quality_score_assigned']
    assert panel.save(str(output)) == 2 and len(list(tmp_path.glob('selected.csv.treatment-*.json'))) == 2
    assert before == {name: digest(pack / name) for name in before}


def test_all_excluded_matches_can_export_a_treatment_record_without_becoming_missing_expression(qtbot, pack, tmp_path):
    annotations(pack, mark_both=True)
    panel = DepositedValuesPanel(pack.parent)
    qtbot.addWidget(panel)
    panel.features.setPlainText('001')
    panel.exclude_source_marked.setChecked(True)
    panel.show_values()
    assert panel.result.empty and panel.table.rowCount() == 0
    assert 'Labels not matched' not in panel.summary.text()
    assert panel.treatment_button.isEnabled() and not panel.export_button.isEnabled()
    assert panel.save(str(tmp_path / 'empty.csv')) == 0
    record = tmp_path / 'all-excluded.json'
    assert panel.save_treatment(str(record)) == 1
    report = json.loads(record.read_text())['source_flag_view']
    assert report['input_measurements'] == report['excluded_measurements'] == 4
    assert [row['source_row'] for row in report['excluded_rows']] == [0, 1]


def test_qualified_identifier_context_is_retained_in_the_view_treatment_record(qtbot, evidence, tmp_path):
    pack, _, _ = evidence
    panel = DepositedValuesPanel(pack.parent)
    qtbot.addWidget(panel)
    panel.identifier_kind.setCurrentIndex(panel.identifier_kind.findData('submitted_gene_label'))
    panel.features.setPlainText('G001')
    panel.exclude_source_marked.setChecked(True)
    panel.show_values()
    path = tmp_path / 'qualified-treatment.json'
    assert panel.save_treatment(str(path)) == 1
    record = json.loads(path.read_text())
    assert record['identifier_namespace'] == 'submitted_gene_label'
    assert record['exact_query_literals'] == ['G001']
    assert len(record['original_identifier_evidence']) == 1
    context = record['original_identifier_evidence'][0]
    assert context['source_row'] == 0 and len(context['evidence']) == 2
    assert all(row['query_literal'] == 'G001' and row['identifier_collection_sha256'] for row in context['evidence'])


def test_loaded_implementation_identity_does_not_depend_on_runtime_source_file_availability(pack, monkeypatch, tmp_path):
    raw = D.profiles(pack, '001')
    before = F.apply(raw)[1]['implementation_sha256']
    monkeypatch.setattr(F, '__file__', str(tmp_path / 'source-not-shipped-in-executable.py'))
    after = F.apply(raw)[1]['implementation_sha256']
    assert before == after == F.IMPLEMENTATION_SHA256
