"""Complete exports retain all stored records beyond bounded UI previews."""
import csv
import json
import threading

import pytest

from parascape.sparse_query import prepare_query, QueryCancelled
from test_sparse_values import fixture_collection


@pytest.mark.parametrize('inverse', [False, True])
def test_complete_export_beyond_preview_and_exact_source_metadata(tmp_path, inverse):
    roots, path = fixture_collection(tmp_path, inverse=inverse)
    result = prepare_query(path, ['001', '001', 'absent'], coordinate_root=roots['sparse_coordinates'],
                           axis_root=roots['sparse_axes'], preview_limit=1)
    try:
        assert len(result.preview) == 1 and result.record['exported_records'] == 4
        assert result.record['query_labels'] == ['001', 'absent']
        assert result.record['unmatched_labels'] == ['absent']
        assert len(result.feature_matches) == 2
        assert result.preview[0]['source_cell_label'] == 'bc-1'
        assert result.preview[0]['original_cell_fields'] == ['bc-1']
        destination = tmp_path / 'complete.csv'
        assert result.export(destination) == 4
        with destination.open() as stream:
            rows = list(csv.DictReader(stream))
        assert [r['entry_index'] for r in rows] == ['0', '1', '2', '3']
        assert rows[1]['value'] == str(2**53 + 1)
        assert [r['source_feature_axis_row'] for r in rows] == ['1', '1', '2', '1']
        assert json.loads(rows[2]['original_feature_fields']) == ['001', 'second', 'submitted-genome']
        sidecar = json.loads(destination.with_suffix('.csv.json').read_text())
        assert sidecar['exported_records'] == 4 and sidecar['displayed_records'] == 1
        assert sidecar['descriptor_sha256'] == path.stem
        original = destination.read_bytes()
        with pytest.raises(FileExistsError):
            result.export(destination)
        assert destination.read_bytes() == original
    finally:
        result.close()
    assert not result.csv_path.exists()


def test_empty_query_exports_header_and_explicit_operation_record(tmp_path):
    roots, path = fixture_collection(tmp_path)
    result = prepare_query(path, ['absent'], coordinate_root=roots['sparse_coordinates'], axis_root=roots['sparse_axes'])
    try:
        assert result.preview == [] and result.record['exported_records'] == 0
        target = tmp_path / 'empty.csv'
        result.export(target)
        with target.open() as stream:
            assert list(csv.DictReader(stream)) == []
        assert json.loads(target.with_suffix('.csv.json').read_text())['operations'][0]['positions'] == []
    finally:
        result.close()


def test_cancellation_does_not_publish_partial_result(tmp_path):
    roots, path = fixture_collection(tmp_path)
    event = threading.Event()
    event.set()
    with pytest.raises(QueryCancelled):
        prepare_query(path, ['001'], coordinate_root=roots['sparse_coordinates'], axis_root=roots['sparse_axes'], cancel=event)


def test_changed_spooled_csv_cannot_export(tmp_path):
    roots, path = fixture_collection(tmp_path)
    result = prepare_query(path, ['001'], coordinate_root=roots['sparse_coordinates'], axis_root=roots['sparse_axes'])
    try:
        result.csv_path.write_text('changed')
        with pytest.raises(ValueError, match='changed'):
            result.export(tmp_path / 'changed.csv')
        assert not (tmp_path / 'changed.csv').exists()
    finally:
        result.close()
