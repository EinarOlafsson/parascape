"""Row-specific filename candidates preserve source ordinals, ambiguity and duplicate design fields."""
import hashlib
import json

import pandas as pd
import pytest
from PyQt6 import QtCore, QtWidgets

from parascape import deposited_values as D
from parascape import row_design_links as R
from parascape.deposited_values_panel import DepositedValuesPanel
from test_deposited_values import pack  # noqa: F401


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repin(root, folder, report_path, report):
    report_path.write_text(json.dumps(report))
    receipt = json.loads((folder / 'validation.json').read_text())
    receipt['outputs'][report_path.name] = digest(report_path)
    (folder / 'validation.json').write_text(json.dumps(receipt))
    pointer = json.loads((root / 'collection.json').read_text())
    pointer['validation_sha256'] = digest(folder / 'validation.json')
    (root / 'collection.json').write_text(json.dumps(pointer))


@pytest.fixture
def linked_rows(pack, tmp_path, monkeypatch):
    metadata = dict(row_layout=dict(layout='original_rowwise_quantities'), row_annotations=[
        dict(source_row=i, feature_label=label, source_fields={'R.FileName': filename,
             'R.Condition': condition, 'R.Replicate': '001'}) for i, (label, filename, condition) in enumerate([
                 ('001', 'Sample1', 'condition A'), ('001', 'Sample1', 'condition B'),
                 ('CaseID', 'Sample2', '')])])
    (pack / 'source_metadata.json').write_text(json.dumps(metadata))
    info = json.loads((pack / 'manifest.json').read_text())
    info['assay'] = 'deposited protein-accession quantities'
    info['outputs']['source_metadata.json'] = digest(pack / 'source_metadata.json')
    (pack / 'manifest.json').write_text(json.dumps(info))
    receipt = json.loads((pack / 'validation.json').read_text())
    receipt['manifest_sha256'] = digest(pack / 'manifest.json')
    (pack / 'validation.json').write_text(json.dumps(receipt))
    manifest = digest(pack / 'manifest.json')
    root = tmp_path / 'audit/protein_row_design_links'
    folder = root / ('1' * 64)
    folder.mkdir(parents=True)
    factors = [dict(source_row=position, source_sha256='b' * 64,
        source_row_id=hashlib.sha256(('b' * 64 + ':' + str(position)).encode()).hexdigest(),
        headers=['comment[data file]', 'comment[data file]', 'characteristics[biological replicate]', 'characteristics[biological replicate]'],
        cells=['Sample1.raw', 'Sample1.raw', '001', 'another declaration'],
        matching_candidate_file_column_indices=[0, 1], metadata_origin='submitted',
        biological_truth_status='Source declarations; biological truth unverified') for position in (2, 3)]
    associations = [dict(source_filename_literal=label, original_source_rows=positions, exact_primary_files=[],
        submitted_factor_rows=[], primary_file_status='unresolved_exact_primary_file',
        filename_stem_primary_candidates=([dict(primary_file_id='one', filename='Sample1.raw'),
                                           dict(primary_file_id='two', filename='Sample1.raw')] if label == 'Sample1' else []),
        filename_stem_submitted_row_candidates=(factors if label == 'Sample1' else []),
        candidate_status='review_required' if label == 'Sample1' else 'no_literal_raw_stem_candidates',
        source_input_table_linkage='unreviewed', biological_independence='not inferred')
        for label, positions in [('Sample1', [0, 1]), ('Sample2', [2])]]
    report = dict(schema_version=1, measurement_manifest_sha256=manifest, input_sha256=info['input_sha256'],
        source_metadata_sha256=info['outputs']['source_metadata.json'], study=info['study'], source_url=info['source_url'],
        source_analysis=info['source_analysis'], features=3, source_filename_field='R.FileName',
        primary_file_metadata_sha256='c' * 64, associations=associations)
    report_path = folder / (manifest + '.json')
    report_path.write_text(json.dumps(report))
    (folder / 'index.json').write_text(json.dumps(dict(schema_version=1, manifests={manifest: report_path.name})))
    (folder / 'validation.json').write_text(json.dumps(dict(schema_version=1,
        status='original_row_filename_associations_reconciled', measurement_packs=1,
        outputs={path.name: digest(path) for path in folder.iterdir()})))
    (root / 'collection.json').write_text(json.dumps(dict(schema_version=1, directory=folder.name,
        validation_sha256=digest(folder / 'validation.json'))))
    monkeypatch.setenv('PARASCAPE_AUDIT', str(tmp_path / 'audit'))
    return pack, root, folder, report_path


def test_repeated_source_rows_keep_independent_original_declarations_and_all_candidates(linked_rows):
    pack, _, _, _ = linked_rows
    original = D.profiles(pack, '001')
    result, record = R.enrich(pack, original)
    pd.testing.assert_frame_equal(result[original.columns], original, check_exact=True)
    assert len(result) == 4 and record['selected_original_source_rows'] == 2
    contexts = [json.loads(value) for value in result.row_design_metadata]
    assert [row['source_row'] for row in contexts] == [0, 1, 0, 1]
    for context in contexts:
        association = context['association']
        assert association['exact_primary_files'] == []
        assert len(association['filename_stem_primary_candidates']) == 2
        assert len(association['filename_stem_submitted_row_candidates']) == 2
        assert association['filename_stem_submitted_row_candidates'][0]['matching_candidate_file_column_indices'] == [0, 1]
        assert association['candidate_status'] == 'review_required'
        assert association['biological_independence'] == 'not inferred'
    assert record['numerical_values_transformed'] is False
    assert record['implementation_sha256'] == R.IMPLEMENTATION_SHA256


@pytest.mark.parametrize('change', ['drop_row', 'swap_row', 'promote_candidate', 'wrong_case', 'drop_duplicate_header', 'changed_factor_row_id', 'infer_independence', 'annotation_context'])
def test_invalid_membership_context_and_candidate_promotion_refused_even_after_rehash(linked_rows, change):
    pack, root, folder, report_path = linked_rows
    report = json.loads(report_path.read_text())
    association = report['associations'][0]
    if change == 'drop_row':
        association['original_source_rows'].pop()
    elif change == 'swap_row':
        association['original_source_rows'][0] = 2
    elif change == 'promote_candidate':
        association['exact_primary_files'] = association['filename_stem_primary_candidates']
        association['primary_file_status'] = 'exact_unique_primary_file'
    elif change == 'wrong_case':
        association['filename_stem_primary_candidates'][0]['filename'] = 'sample1.raw'
    elif change == 'drop_duplicate_header':
        association['filename_stem_submitted_row_candidates'][0]['matching_candidate_file_column_indices'] = [0]
    elif change == 'changed_factor_row_id':
        association['filename_stem_submitted_row_candidates'][0]['source_row_id'] = 'd' * 64
    elif change == 'infer_independence':
        association['biological_independence'] = 'independent samples'
    else:
        report['source_metadata_sha256'] = 'd' * 64
    repin(root, folder, report_path, report)
    with pytest.raises(ValueError):
        R.enrich(pack, D.profiles(pack, '001'))


def test_missing_collection_and_unmatched_query_are_explicit_and_do_not_invent_links(linked_rows, tmp_path):
    pack, _, _, _ = linked_rows
    result, record = R.enrich(pack, D.profiles(pack, '001'), audit_root=tmp_path / 'missing')
    assert 'not installed' in record['association_status']
    assert all('not installed' in json.loads(value)['association']['association_status'] for value in result.row_design_metadata)
    empty, record = R.enrich(pack, D.profiles(pack, 'absent'))
    assert empty.empty and record['selected_original_source_rows'] == 0
    frame = D.profiles(pack, '001')
    frame.loc[0, 'source_row'] = 2
    with pytest.raises(ValueError, match='annotation identity'):
        R.enrich(pack, frame)


def test_selected_measurement_details_and_CSV_keep_candidate_grade_and_duplicate_fields(qtbot, linked_rows, tmp_path):
    pack, _, _, _ = linked_rows
    panel = DepositedValuesPanel(pack.parent)
    qtbot.addWidget(panel)
    panel.features.setPlainText('001')
    panel.show_values()
    assert len(panel.result) == 4
    panel.table.setCurrentCell(0, 0)
    panel.show_measurement_context()
    fields = []
    iterator = QtWidgets.QTreeWidgetItemIterator(panel.measurement_tree)
    while iterator.value() is not None:
        item = iterator.value()
        fields.append((item.text(0), item.text(1)))
        iterator += 1
    assert ('candidate_status', 'review_required') in fields
    association_nodes = panel.measurement_tree.findItems('association',
        QtCore.Qt.MatchFlag.MatchExactly | QtCore.Qt.MatchFlag.MatchRecursive)
    assert len(association_nodes) == 1 and association_nodes[0].isExpanded()
    assert ('1: comment[data file]', 'Sample1.raw') in fields
    assert ('2: comment[data file]', 'Sample1.raw') in fields
    assert ('3: characteristics[biological replicate]', '001') in fields
    assert ('4: characteristics[biological replicate]', 'another declaration') in fields
    output = tmp_path / 'rows.csv'
    assert panel.save(str(output)) == 4
    exported = pd.read_csv(output, dtype=str, keep_default_na=False)
    assert exported.source_token.tolist() == ['-0.000', '1.2500', '0', '-2.500']
    assert 'row_design_metadata' in exported
    record = json.loads(next(tmp_path.glob('rows.csv.treatment-*.json')).read_text())
    assert record['row_design_association_view']['collection']
    assert record['row_design_association_view']['numerical_values_transformed'] is False
    panel.features.setPlainText('CaseID')
    assert panel.measurement_dialog is None and panel.treatment_record is None
