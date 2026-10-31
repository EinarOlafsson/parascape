"""Original sparse queries preserve coordinate, token and feature identity."""
import hashlib
import json
import struct

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from parascape import sparse_values as sv


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True) + '\n')
    return sv._sha(path)


def fixture_collection(tmp_path, inverse=False, real=False):
    context = dict(studies='GSE-fixture', source_entities='Sample:GSM-fixture', primary_xml_sha256='a' * 64)
    roots = {k: tmp_path / k for k in ('sparse_coordinates', 'sparse_axes', 'sparse_datasets')}
    for root in roots.values():
        root.mkdir()
    axes = {}
    for role, rows in [('features', [['001', 'first', 'submitted-genome'], ['001', 'second', 'submitted-genome'], ['ABC', '', 'other']]),
                       ('barcodes', [['bc-1'], ['bc-2']])]:
        folder = roots['sparse_axes'] / hashlib.sha256(role.encode()).hexdigest()
        folder.mkdir()
        pq.write_table(pa.table(dict(source_axis_row=pa.array(range(1, len(rows) + 1), type=pa.int64()),
                                     first_field=[row[0] for row in rows],
                                     original_fields=pa.array(rows, type=pa.list_(pa.string())))), folder / 'axis.parquet')
        axis = dict(source_url=role, input_sha256='b' * 64, source_record=context, data_rows=len(rows),
                    original_axis_cells_sha256=hashlib.sha256(role.encode()).hexdigest(),
                    original_headers=['gene_id', 'gene_name', 'genome'] if role == 'features' else [],
                    source_genome_field='genome' if role == 'features' else None,
                    outputs={'axis.parquet': sv._sha(folder / 'axis.parquet')})
        mh = write_json(folder / 'manifest.json', axis)
        validation = dict(state='all_original_axis_cells_reconciled', input_sha256=axis['input_sha256'],
                          manifest_sha256=mh, independent_source_axis_cells_sha256=axis['original_axis_cells_sha256'],
                          all_original_field_cells_and_order_reconciled=True, canonical_gene_or_cell_identity_inferred=False)
        vh = write_json(folder / 'validation.json', validation)
        reference = dict(directory=folder.name, manifest_sha256=mh, validation_sha256=vh,
                         output_name='axis.parquet', output_sha256=axis['outputs']['axis.parquet'],
                         source_url=role, input_sha256=axis['input_sha256'], original_headers=axis['original_headers'],
                         source_genome_field=axis['source_genome_field'])
        coordinate = ('source_column' if inverse else 'source_row') if role == 'features' else ('source_row' if inverse else 'source_column')
        axes[role] = dict(source_coordinate=coordinate, data_rows=len(rows),
                          logical_axis_sha256=axis['original_axis_cells_sha256'], original_axis_artifacts=[reference])
    dims = [2, 3, 4] if inverse else [3, 2, 4]
    # Duplicate coordinates and explicit zero stay separate. No stored entry for feature3.
    rows, columns = ([1, 1, 2, 2], [1, 1, 2, 1]) if inverse else ([1, 1, 2, 1], [1, 1, 2, 2])
    values = [-0.0, 2.5, -3.0, 4.0] if real else [0, 2**53 + 1, -3, 4]
    tokens = ['-0', '2.5', '-3', '4'] if real else ['0', str(2**53 + 1), '-3', '4']
    folder = roots['sparse_coordinates'] / ('c' * 64)
    folder.mkdir()
    pq.write_table(pa.table(dict(entry_index=pa.array(range(4), type=pa.int64()), source_row=pa.array(rows, type=pa.int64()),
                                 source_column=pa.array(columns, type=pa.int64()), row_token=[str(n) for n in rows],
                                 column_token=[str(n) for n in columns], source_token=tokens,
                                 value=pa.array(values, type=pa.float64() if real else pa.int64()))), folder / 'entries.parquet')
    matrix = dict(source_url='matrix', input_sha256='d' * 64, original_source_record=context,
                  declared_dimensions=dims, numeric_field='real' if real else 'integer',
                  outputs={'entries.parquet': sv._sha(folder / 'entries.parquet')})
    mh = write_json(folder / 'manifest.json', matrix)
    receipt = dict(state='all_original_coordinate_records_reconciled', manifest_sha256=mh, input_sha256='d' * 64,
                   original_records=4, matrix_dimensions=dims, all_original_coordinate_tokens_and_order_checked=True,
                   normalization_applied=False, canonical_gene_mapping=False, biological_quality_score=None)
    vh = write_json(folder / 'validation.json', receipt)
    reference = dict(directory=folder.name, manifest_sha256=mh, validation_sha256=vh,
                     output_name='entries.parquet', output_sha256=matrix['outputs']['entries.parquet'])
    descriptor = dict(source_url='matrix', source_record=context, input_sha256='d' * 64, numeric_field=matrix['numeric_field'],
                      declared_dimensions=dims, axis_roles=axes, coordinate_artifact=reference,
                      numeric_payload_verified=True, numeric_transpose_applied=False,
                      biological_or_comparative_admission=False, quantity='original source units unreviewed')
    return roots, save_descriptor(roots, descriptor)


def save_descriptor(roots, descriptor):
    body = json.dumps(descriptor, sort_keys=True) + '\n'
    digest = hashlib.sha256(body.encode()).hexdigest()
    path = roots['sparse_datasets'] / (digest + '.json')
    path.write_text(body)
    write_json(roots['sparse_datasets'] / 'catalogue.json', dict(datasets=[dict(file=path.name, sha256=digest, source_url='matrix')]))
    return path


def open_dataset(roots, path):
    return sv.SparseDataset(path, coordinate_root=roots['sparse_coordinates'], axis_root=roots['sparse_axes'])


@pytest.mark.parametrize('inverse', [False, True])
@pytest.mark.parametrize('real', [False, True])
def test_exact_source_rows_tokens_duplicates_and_numeric_bits(tmp_path, inverse, real):
    roots, path = fixture_collection(tmp_path, inverse, real)
    dataset = open_dataset(roots, path)
    entries, errors = sv.catalogue(roots['sparse_datasets'])
    assert not errors and len(entries) == 1
    matches = [r for batch in dataset.iter_feature_matches('001', batch_size=1) for r in batch.to_pylist()]
    assert [r['source_axis_row'] for r in matches] == [1, 2]
    assert matches[0]['original_fields'] == ['001', 'first', 'submitted-genome']
    assert list(dataset.iter_feature_matches('1')) == []
    observed = [r for batch in dataset.iter_entries([1, 2], batch_size=1) for r in batch.to_pylist()]
    assert [r['entry_index'] for r in observed] == [0, 1, 2, 3]
    assert observed[0]['source_row'] == observed[1]['source_row']
    assert observed[0]['source_column'] == observed[1]['source_column']
    assert observed[0]['source_token'] == ('-0' if real else '0')
    if real:
        assert struct.pack('>d', observed[0]['value']) == struct.pack('>d', -0.0)
    else:
        assert observed[1]['value'] == 2**53 + 1
    assert list(dataset.iter_entries([3])) == []  # No imputed zero.
    assert list(dataset.iter_entries([])) == []
    assert not dataset.info['numeric_transpose_applied']


@pytest.mark.parametrize('positions', [[0], [4], [-1], [1.0], [True]])
def test_invalid_original_positions_refused(tmp_path, positions):
    roots, path = fixture_collection(tmp_path)
    with pytest.raises(ValueError, match='1-based'):
        list(open_dataset(roots, path).iter_entries(positions))


def test_mutated_output_refused_after_open(tmp_path):
    roots, path = fixture_collection(tmp_path)
    dataset = open_dataset(roots, path)
    output = roots['sparse_coordinates'] / ('c' * 64) / 'entries.parquet'
    output.write_bytes(output.read_bytes() + b'changed')
    with pytest.raises(ValueError, match='changed'):
        list(dataset.iter_entries([1]))


@pytest.mark.parametrize('mutation', ['dimensions', 'orientation', 'context', 'axis_digest', 'normalization'])
def test_rehashed_descriptor_cannot_override_receipts(tmp_path, mutation):
    roots, path = fixture_collection(tmp_path)
    descriptor = json.loads(path.read_text())
    if mutation == 'dimensions':
        descriptor['declared_dimensions'][0] += 1
    elif mutation == 'orientation':
        descriptor['axis_roles']['features']['source_coordinate'] = 'source_column'
    elif mutation == 'context':
        descriptor['source_record']['primary_xml_sha256'] = '0' * 64
    elif mutation == 'axis_digest':
        descriptor['axis_roles']['features']['logical_axis_sha256'] = '0' * 64
    else:
        descriptor['numeric_transpose_applied'] = True
    path = save_descriptor(roots, descriptor)
    with pytest.raises(ValueError):
        open_dataset(roots, path)


def test_escaping_artifact_refused(tmp_path):
    roots, path = fixture_collection(tmp_path)
    descriptor = json.loads(path.read_text())
    descriptor['coordinate_artifact']['directory'] = '../escape'
    with pytest.raises(ValueError, match='identity'):
        open_dataset(roots, save_descriptor(roots, descriptor))


def test_each_query_batch_retains_operation_and_source_record(tmp_path):
    roots, path = fixture_collection(tmp_path, inverse=True)
    dataset = open_dataset(roots, path)
    expected = dataset.query_record([1, 2])
    assert expected['descriptor_sha256'] == path.stem
    assert expected['source_sha256'] == 'd' * 64
    assert expected['operations'][0] == dict(operation='select_original_feature_positions',
                                             source_coordinate='source_column', positions=[1, 2])
    assert not any(expected[key] for key in ('normalization_applied', 'imputation_applied',
                                           'coordinate_aggregation_applied', 'numeric_transpose_applied',
                                           'canonical_gene_mapping', 'biological_or_comparative_admission'))
    for batch in dataset.iter_entries([1, 2], batch_size=1):
        assert json.loads(batch.schema.metadata[b'parascape_query_record']) == expected
    assert dataset.query_record([])['operations'][0]['positions'] == []


@pytest.mark.parametrize('inverse', [False, True])
def test_original_cell_and_feature_metadata_by_position(tmp_path, inverse):
    roots, path = fixture_collection(tmp_path, inverse=inverse)
    dataset = open_dataset(roots, path)
    rows = [r for b in dataset.iter_axis_rows('features', [2]) for r in b.to_pylist()]
    assert rows == [dict(source_axis_row=2, first_field='001', original_fields=['001', 'second', 'submitted-genome'])]
    cells = [r for b in dataset.iter_axis_rows('barcodes', [2, 1]) for r in b.to_pylist()]
    assert [r['first_field'] for r in cells] == ['bc-1', 'bc-2']
    with pytest.raises(ValueError):
        list(dataset.iter_axis_rows('barcodes', [3]))
    with pytest.raises(ValueError):
        list(dataset.iter_axis_rows('unknown', [1]))


def test_reader_limits_physical_row_group_scheduling_and_preserves_records(tmp_path, monkeypatch):
    roots, path = fixture_collection(tmp_path)
    descriptor = json.loads(path.read_text())
    reference = descriptor['coordinate_artifact']
    folder = roots['sparse_coordinates']/reference['directory']
    values = pq.read_table(folder/'entries.parquet')
    pq.write_table(values, folder/'entries.parquet', row_group_size=1)
    info = json.loads((folder/'manifest.json').read_text())
    info['outputs']['entries.parquet'] = sv._sha(folder/'entries.parquet')
    reference['manifest_sha256'] = write_json(folder/'manifest.json', info)
    receipt = json.loads((folder/'validation.json').read_text())
    receipt['manifest_sha256'] = reference['manifest_sha256']
    reference['validation_sha256'] = write_json(folder/'validation.json', receipt)
    reference['output_sha256'] = info['outputs']['entries.parquet']
    path = save_descriptor(roots, descriptor)
    original = sv.pq.ParquetFile
    scheduled = []
    class Reader:
        def __init__(self, *args, **kwargs):
            self.reader = original(*args, **kwargs)
        def __getattr__(self, name):
            return getattr(self.reader, name)
        def iter_batches(self, **kwargs):
            assert len(kwargs['row_groups']) == 1 and kwargs['use_threads'] is False
            scheduled.extend(kwargs['row_groups'])
            yield from self.reader.iter_batches(**kwargs)
    monkeypatch.setattr(sv.pq, 'ParquetFile', Reader)
    dataset = open_dataset(roots, path)
    observed = [row for batch in dataset.iter_entries([1, 2], batch_size=3) for row in batch.to_pylist()]
    assert observed == values.to_pylist()
    assert scheduled == [0, 1, 2, 3]
