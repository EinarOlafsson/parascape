"""Browse preserved sparse source records without interpreting absent coordinates.

Feature labels are literal original first fields, not canonical gene identifiers.
Repeated feature labels and repeated coordinates stay distinct. Neither query nor
export sums records, transposes values, normalizes measurements or creates zeros.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path
from typing import Iterable

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from . import paths

_VERIFIED: dict[tuple, bool] = {}


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _checked(path: Path, expected: str) -> None:
    if not re.fullmatch('[a-f0-9]{64}', expected):
        raise ValueError('Invalid sparse artifact digest')
    stat = path.stat()
    key = (str(path.resolve()), stat.st_dev, stat.st_ino, stat.st_size,
           stat.st_mtime_ns, stat.st_ctime_ns, expected)
    if key not in _VERIFIED:
        after = _sha(path)
        current = path.stat()
        if after != expected or any(getattr(current, name) != getattr(stat, name) for name in
                                    ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns')):
            raise ValueError('Sparse artifact changed: ' + path.name)
        # Bound process-local verification memory across different installations.
        if len(_VERIFIED) >= 4096:
            _VERIFIED.clear()
        _VERIFIED[key] = True


def _child(root: Path, name: str) -> Path:
    if not isinstance(name, str) or not re.fullmatch('[a-f0-9]{64}', name):
        raise ValueError('Invalid sparse pack identity')
    child = root / name
    if child.is_symlink() or child.resolve().parent != root.resolve():
        raise ValueError('Sparse artifact lies outside its collection')
    return child


def _batches(path: Path, batch_size: int):
    """Keep one physical row group active, including for billion-record sources."""
    file = pq.ParquetFile(path)
    for group in range(file.num_row_groups):
        yield from file.iter_batches(batch_size=batch_size, row_groups=[group], use_threads=False)


def data_root() -> Path:
    """Locate original sparse dataset descriptors without changing local overrides."""
    return Path(os.environ.get('PARASCAPE_SPARSE_DATASETS',
                              str(Path(paths.data_dir()) / 'sparse_datasets')))


def catalogue(root: str | Path | None = None) -> tuple[list[dict], list[str]]:
    """List bound original-source datasets; missing collections have no entries."""
    root = Path(root) if root is not None else data_root()
    if not (root / 'catalogue.json').exists():
        return [], []
    entries, errors = [], []
    try:
        index = json.loads((root / 'catalogue.json').read_text())
        seen = set()
        for item in index['datasets']:
            try:
                if item['sha256'] in seen:
                    raise ValueError('Repeated sparse dataset descriptor')
                seen.add(item['sha256'])
                descriptor = _child(root, item['sha256']).with_suffix('.json')
                if item['file'] != descriptor.name or descriptor.is_symlink():
                    raise ValueError('Invalid sparse descriptor path')
                _checked(descriptor, item['sha256'])
                info = json.loads(descriptor.read_text())
                if info['source_url'] != item['source_url']:
                    raise ValueError('Sparse catalogue source identity differs')
                entries.append(dict(path=str(descriptor), source_url=info['source_url'],
                                    source_record=info['source_record'],
                                    dimensions=info['declared_dimensions'],
                                    quantity=info['quantity'],
                                    biological_or_comparative_admission=False))
            except (OSError, ValueError, KeyError, TypeError) as error:
                errors.append(str(item) + ': ' + str(error))
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append(str(root) + ': ' + str(error))
    return entries, errors


def coordinate_reconciliation(matrix: dict, receipt: dict) -> dict | None:
    """Require whole-source proof for text coordinates or original binary arrays."""
    source_format = matrix.get('source_format', 'MatrixMarket')
    if source_format in ('HDF5_CSC', 'AnnData_HDF5'):
        storage = matrix['source_storage']
        expected_state = ('all_original_HDF5_CSC_records_reconciled' if source_format == 'HDF5_CSC'
                          else 'all_original_AnnData_layer_records_reconciled')
        if (receipt['state'] != expected_state
                or receipt.get('all_original_binary_values_and_order_checked') is not True
                or receipt.get('derived_coordinate_and_value_tokens_checked') is not True
                or receipt.get('original_text_tokens_present') is not False
                or not re.fullmatch('[a-f0-9]{64}', matrix.get('original_numeric_array_logical_bytes_sha256', ''))
                or receipt.get('independent_original_numeric_array_bytes_sha256') != matrix['original_numeric_array_logical_bytes_sha256']
                or storage.get('format') != source_format
                or storage.get('source_index_base') != 0 or storage.get('output_index_base') != 1
                or storage.get('decimal_tokens_are_derived') is not True
                or storage.get('coordinate_aggregation_applied') is not False):
            raise ValueError('Binary CSC reconstruction lacks a complete source proof')
        if source_format == 'AnnData_HDF5':
            if (storage.get('layout') not in ('csr', 'csc', 'dense')
                    or not isinstance(storage.get('matrix_path'), str) or not storage['matrix_path']
                    or storage.get('source_axis_order') != ['observations', 'variables']
                    or storage.get('original_dtype_roundtrip_verified') is not True
                    or storage.get('layer_aggregation_applied') is not False
                    or receipt.get('all_original_dtype_bits_recoverable') is not True
                    or receipt.get('every_dense_cell_retained') is not (storage['layout'] == 'dense')
                    or (storage['layout'] == 'dense' and matrix['declared_dimensions'][2]
                        != matrix['declared_dimensions'][0] * matrix['declared_dimensions'][1])):
                raise ValueError('AnnData layer lacks its original layout, axes or binary proof')
        return storage
    if (matrix.get('source_format', 'MatrixMarket') != 'MatrixMarket'
            or receipt['state'] != 'all_original_coordinate_records_reconciled'
            or not receipt['all_original_coordinate_tokens_and_order_checked']):
        raise ValueError('Original text coordinate reconciliation differs')
    return None


class SparseDataset:
    """A checked original-coordinate view with explicit source-axis orientation."""

    def __init__(self, descriptor: str | Path, *, coordinate_root: str | Path | None = None,
                 axis_root: str | Path | None = None):
        descriptor = Path(descriptor)
        if descriptor.is_symlink():
            raise ValueError('Sparse descriptor must be a regular collection artifact')
        _checked(descriptor, descriptor.stem)
        self.descriptor_sha256 = descriptor.stem
        self.info = json.loads(descriptor.read_text())
        info = self.info
        if (not info['numeric_payload_verified'] or info['numeric_transpose_applied']
                or info['biological_or_comparative_admission']):
            raise ValueError('Unsupported sparse dataset admission or transformation')
        self.coordinate_root = Path(coordinate_root) if coordinate_root is not None else Path(
            os.environ.get('PARASCAPE_SPARSE_COORDINATES', str(descriptor.parent.parent / 'sparse_coordinates')))
        self.axis_root = Path(axis_root) if axis_root is not None else Path(
            os.environ.get('PARASCAPE_SPARSE_AXES', str(descriptor.parent.parent / 'sparse_axes')))
        dims = info['declared_dimensions']
        if len(dims) != 3 or any(type(n) is not int or n < 0 for n in dims):
            raise ValueError('Invalid sparse dimensions')
        orientation = {role: info['axis_roles'][role]['source_coordinate'] for role in ('features', 'barcodes')}
        if set(orientation.values()) != {'source_row', 'source_column'}:
            raise ValueError('Sparse feature/cell orientations differ')
        self.feature_coordinate = orientation['features']
        self.dimensions = dict(source_row=dims[0], source_column=dims[1])
        self._coordinate_ref = info['coordinate_artifact']
        folder, matrix, receipt = self._artifact(self.coordinate_root, self._coordinate_ref, 'entries.parquet')
        storage = coordinate_reconciliation(matrix, receipt)
        if storage is not None and info.get('source_storage') != storage:
            raise ValueError('Binary CSC descriptor and source storage context differ')
        if (storage is not None and storage.get('format') == 'AnnData_HDF5'
                and orientation != {'features': 'source_column', 'barcodes': 'source_row'}):
            raise ValueError('AnnData must retain observations by variables without transposing')
        if (matrix['input_sha256'] != info['input_sha256'] or matrix['source_url'] != info['source_url']
                or matrix['declared_dimensions'] != dims or receipt['matrix_dimensions'] != dims
                or receipt['original_records'] != dims[2] or matrix['numeric_field'] != info['numeric_field']
                or receipt['normalization_applied'] or receipt['canonical_gene_mapping']
                or receipt['biological_quality_score'] is not None):
            raise ValueError('Sparse numeric validation coverage differs')
        self._context(matrix['original_source_record'])
        self.entries_path = folder / 'entries.parquet'
        schema = pq.ParquetFile(self.entries_path).schema_arrow
        expected = ['entry_index', 'source_row', 'source_column', 'row_token', 'column_token', 'source_token', 'value']
        if (schema.names != expected or any(schema.field(n).type != pa.int64() for n in expected[:3])
                or any(schema.field(n).type != pa.string() for n in expected[3:6])
                or info['numeric_field'] not in ('integer', 'real')
                or schema.field('value').type != (pa.int64() if info['numeric_field'] == 'integer' else pa.float64())
                or pq.ParquetFile(self.entries_path).metadata.num_rows != dims[2]):
            raise ValueError('Sparse numeric storage schema or count differs')
        self._axes = {}
        for role in ('features', 'barcodes'):
            association = info['axis_roles'][role]
            if (association['data_rows'] != self.dimensions[orientation[role]]
                    or not association['original_axis_artifacts']):
                raise ValueError('Sparse axis dimensions or artifacts differ')
            for reference in association['original_axis_artifacts']:
                axis_folder, axis, validation = self._artifact(self.axis_root, reference, 'axis.parquet')
                if (axis['source_url'] != reference['source_url'] or axis['input_sha256'] != reference['input_sha256']
                        or axis['data_rows'] != association['data_rows']
                        or axis['original_axis_cells_sha256'] != association['logical_axis_sha256']
                        or validation['independent_source_axis_cells_sha256'] != association['logical_axis_sha256']
                        or validation['state'] != 'all_original_axis_cells_reconciled'
                        or not validation['all_original_field_cells_and_order_reconciled']
                        or validation['canonical_gene_or_cell_identity_inferred']
                        or axis['original_headers'] != reference['original_headers']
                        or axis['source_genome_field'] != reference['source_genome_field']):
                    raise ValueError('Sparse axis validation or source context differs')
                self._context(axis['source_record'])
                table = pq.ParquetFile(axis_folder / 'axis.parquet')
                if (table.metadata.num_rows != association['data_rows']
                        or table.schema_arrow.names != ['source_axis_row', 'first_field', 'original_fields']
                        or table.schema_arrow.field('source_axis_row').type != pa.int64()
                        or table.schema_arrow.field('first_field').type != pa.string()
                        or table.schema_arrow.field('original_fields').type != pa.list_(pa.string())):
                    raise ValueError('Sparse original axis storage differs')
            # Identical logical aliases remain in provenance; query one table once.
            self._axes[role] = association['original_axis_artifacts'][0]

    def _context(self, record: dict) -> None:
        if any(record[k] != self.info['source_record'][k]
               for k in ('studies', 'source_entities', 'primary_xml_sha256')):
            raise ValueError('Sparse primary study/sample context differs')

    @staticmethod
    def _artifact(root: Path, reference: dict, output: str) -> tuple[Path, dict, dict]:
        folder = _child(root, reference['directory'])
        for name, digest in [('manifest.json', reference['manifest_sha256']),
                             ('validation.json', reference['validation_sha256']),
                             (output, reference['output_sha256'])]:
            if (folder / name).is_symlink():
                raise ValueError('Sparse output must stay inside its artifact')
            _checked(folder / name, digest)
        info = json.loads((folder / 'manifest.json').read_text())
        receipt = json.loads((folder / 'validation.json').read_text())
        if (reference['output_name'] != output or info['outputs'] != {output: reference['output_sha256']}
                or receipt['manifest_sha256'] != reference['manifest_sha256']
                or receipt['input_sha256'] != info['input_sha256']):
            raise ValueError('Sparse artifact receipt differs')
        return folder, info, receipt

    def iter_feature_matches(self, literal_label: str, *, batch_size: int = 65536) -> Iterable[pa.RecordBatch]:
        """Match original first-field labels exactly, retaining every duplicate row."""
        if not isinstance(literal_label, str):
            raise TypeError('Sparse feature queries require a literal source label')
        reference = self._axes['features']
        folder, _, _ = self._artifact(self.axis_root, reference, 'axis.parquet')
        for batch in _batches(folder / 'axis.parquet', batch_size):
            matched = batch.filter(pc.equal(batch.column('first_field'), literal_label))
            if matched.num_rows:
                yield matched

    def iter_axis_rows(self, role: str, source_axis_rows: Iterable[int], *, batch_size: int = 65536) -> Iterable[pa.RecordBatch]:
        """Return unchanged submitted metadata at original feature or cell positions."""
        if role not in ('features', 'barcodes'):
            raise ValueError('Unknown original sparse axis role')
        selected = list(source_axis_rows)
        association = self.info['axis_roles'][role]
        if any(type(n) is not int or not 1 <= n <= association['data_rows'] for n in selected):
            raise ValueError('Axis positions must be original 1-based data rows')
        if not selected:
            return
        reference = self._axes[role]
        folder, _, _ = self._artifact(self.axis_root, reference, 'axis.parquet')
        selected_array = pa.array(selected, type=pa.int64())
        for batch in _batches(folder / 'axis.parquet', batch_size):
            matched = batch.filter(pc.is_in(batch.column('source_axis_row'), value_set=selected_array))
            if matched.num_rows:
                yield matched

    def query_record(self, source_feature_rows: Iterable[int]) -> dict:
        """Describe the exact source-position selection and unchanged value treatment."""
        selected = list(source_feature_rows)
        if any(type(n) is not int or not 1 <= n <= self.dimensions[self.feature_coordinate] for n in selected):
            raise ValueError('Feature positions must be original 1-based axis data rows')
        return dict(schema_version=1, descriptor_sha256=self.descriptor_sha256,
                    source_url=self.info['source_url'], source_sha256=self.info['input_sha256'],
                    source_record=self.info['source_record'],
                    reader_sha256=_sha(Path(__file__)),
                    coordinate_artifact=self._coordinate_ref,
                    axis_roles=self.info['axis_roles'],
                    operations=[dict(operation='select_original_feature_positions',
                                     source_coordinate=self.feature_coordinate, positions=selected),
                                dict(operation='retain_original_entry_order_and_numeric_storage')],
                    normalization_applied=False, imputation_applied=False,
                    coordinate_aggregation_applied=False, numeric_transpose_applied=False,
                    canonical_gene_mapping=False,
                    absent_coordinate_policy='No stored record returned; no zero or missing-value interpretation inferred',
                    biological_or_comparative_admission=False,
                    **({'source_storage': self.info['source_storage']} if 'source_storage' in self.info else {}))

    def iter_entries(self, source_feature_rows: Iterable[int], *, batch_size: int = 65536) -> Iterable[pa.RecordBatch]:
        """Read stored records for original feature positions, preserving entry order.

        Missing coordinates produce no records. Repeated coordinates are never
        summed. Values retain the original exact int64 or real binary64 storage.
        """
        selected = list(source_feature_rows)
        record = self.query_record(selected)
        metadata = {b'parascape_query_record': json.dumps(record, sort_keys=True).encode()}
        if not selected:
            return
        folder, _, _ = self._artifact(self.coordinate_root, self._coordinate_ref, 'entries.parquet')
        selected_array = pa.array(selected, type=pa.int64())
        for batch in _batches(folder / 'entries.parquet', batch_size):
            matched = batch.filter(pc.is_in(batch.column(self.feature_coordinate), value_set=selected_array))
            matched = matched.replace_schema_metadata(metadata)
            if matched.num_rows:
                yield matched
