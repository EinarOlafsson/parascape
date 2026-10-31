"""Inspect preserved workbook source cells without inferring measurement roles."""
from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
import threading
import time

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from . import paths

_VERIFIED = set()


class QueryCancelled(Exception):
    """A source query was cancelled before publishing its result."""


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _check(path, digest, size=None):
    if path.is_symlink() or not re.fullmatch('[a-f0-9]{64}', digest):
        raise ValueError('Invalid workbook artifact identity')
    stat = path.stat()
    signature = stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns
    if size is not None and stat.st_size != size:
        raise ValueError('Workbook artifact size changed: ' + path.name)
    key = str(path.resolve()), signature, digest
    if key not in _VERIFIED:
        if _sha(path) != digest or signature != tuple(getattr(path.stat(), x) for x in
                ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns')):
            raise ValueError('Workbook artifact changed: ' + path.name)
        if len(_VERIFIED) >= 2048:
            _VERIFIED.clear()
        _VERIFIED.add(key)


def data_root():
    """Return the downloaded collection, honoring an explicit local override."""
    return Path(os.environ.get('PARASCAPE_WORKBOOK_VALUES', str(Path(paths.data_dir()) / 'workbook_cells')))


def manifest(folder, *, verify=False):
    """Read original storage facts; admission as comparable measurements stays false."""
    folder = Path(folder)
    if folder.is_symlink() or (folder / 'manifest.json').is_symlink():
        raise ValueError('Workbook collection must use regular artifacts')
    info = json.loads((folder / 'manifest.json').read_text())
    if (info['schema_version'] != 1 or info['operation_version'] != 'source-cells-v1'
            or info['original_bytes_changed'] or info['extraction']['format'] not in ('OOXML', 'BIFF_decoded')):
        raise ValueError('Unsupported workbook source treatment')
    for key in ('normalization', 'aggregation', 'feature_mapping', 'header_inference', 'formula_evaluation'):
        if info['treatment'][key] is not False:
            raise ValueError('Unsupported workbook transformation: ' + key)
    for key in ('biological_quality_score', 'stage_score', 'strain_score'):
        if info['treatment'][key] is not None:
            raise ValueError('Source storage cannot establish biological scores')
    if set(info['outputs']) != {'cells.parquet', 'strings.parquet', 'rows.parquet'}:
        raise ValueError('Workbook storage artifacts differ')
    if verify:
        for name, expected in info['outputs'].items():
            _check(folder / name, expected['sha256'], expected['size'])
        cells = pq.ParquetFile(folder / 'cells.parquet')
        if cells.metadata.num_rows != info['extraction']['cells']:
            raise ValueError('Workbook cell count differs')
        expected_types = {'sheet_index': pa.int32(), 'cell_order': pa.int64(), 'row': pa.int32(),
                          'column': pa.int32(), 'coordinate_inferred': pa.bool_()}
        strings = ('source_address', 'storage_type', 'source_value_token', 'decoded_text',
                   'numeric_binary64_hex', 'formula_xml', 'source_cell_xml', 'source_attributes_json')
        expected_types.update({name: pa.string() for name in strings})
        if set(cells.schema_arrow.names) != set(expected_types) or any(
                cells.schema_arrow.field(name).type != dtype for name, dtype in expected_types.items()):
            raise ValueError('Workbook cell schema differs')
    return info


def validation(folder, info=None):
    """Expose source-bound verification, keeping absent or legacy evidence unknown."""
    folder = Path(folder)
    info = info or manifest(folder)
    path = folder / 'independent_validation.json'
    if not path.exists():
        record = dict(status='independent_verification_pending', independent_whole_source_reconciliation=False,
                      biological_quality_score=None)
    elif path.is_symlink():
        raise ValueError('Workbook validation must be a regular artifact')
    else:
        record = json.loads(path.read_text())
    if path.exists() and (record['measurement_manifest_sha256'] != _sha(folder / 'manifest.json')
            or record['input_sha256'] != info['input_sha256'] or record['measurement_admission']
            or record['biological_quality_score'] is not None):
        raise ValueError('Workbook validation describes a different source or scope')
    if record.get('independent_whole_source_reconciliation'):
        if (record['status'] != 'all_OOXML_source_cells_independently_reconciled'
                or info['extraction']['format'] != 'OOXML' or record['cells'] != info['extraction']['cells']
                or record['rows'] != info['extraction']['rows']
                or record['shared_strings'] != info['extraction']['shared_strings']):
            raise ValueError('Workbook verification coverage differs')
    result = dict(record, decoded_value_checks_agree=True if record.get('independent_whole_source_reconciliation') else None,
                  reported_reader_disagreements=0 if record.get('independent_whole_source_reconciliation') else None)
    if info['extraction']['format'] == 'BIFF_decoded':
        from .workbook_qc import report
        supplementary = report(folder, info)
        if supplementary is not None:
            raw = supplementary['raw_report']
            result.update(status=raw['status'], independent_whole_source_reconciliation=True,
                          decoded_value_checks_agree=raw['independent_Calamine_whole_text_agreement'],
                          reported_reader_disagreements=raw['independent_Calamine_text_disagreements'],
                          supplementary_biff_qc=supplementary, original_validation=record)
    return result


def catalogue(root=None):
    """List all completed source occurrences, including unknown verification grades."""
    root = Path(root) if root is not None else data_root()
    entries, errors = [], []
    if not (root / 'catalogue.json').exists():
        return entries, errors
    try:
        index = json.loads((root / 'catalogue.json').read_text())
        seen = set()
        for item in index['workbooks']:
            if item['status'] != 'complete_source_cell_extraction':
                errors.append(str(item.get('source_url')) + ': extraction requires review')
                continue
            try:
                key = item['directory']
                if not isinstance(key, str) or not re.fullmatch('[a-f0-9]{64}', key) or key in seen:
                    raise ValueError('Invalid or repeated workbook directory')
                seen.add(key)
                folder = root / key
                if folder.is_symlink() or folder.resolve().parent != root.resolve():
                    raise ValueError('Workbook lies outside its collection')
                info = manifest(folder)
                if info['input_sha256'] != item['input_sha256'] or info['source_url'] != item['source_url']:
                    raise ValueError('Workbook catalogue and source identity differ')
                proof = validation(folder, info)
                entries.append(dict(path=str(folder), source_url=info['source_url'],
                                    source_record=info['source_record'], cells=info['extraction']['cells'],
                                    sheets=len(info['extraction']['sheets']), format=info['extraction']['format'],
                                    validation_status=proof['status'],
                                    independent_whole_source_reconciliation=proof['independent_whole_source_reconciliation'],
                                    decoded_value_checks_agree=proof['decoded_value_checks_agree'],
                                    reported_reader_disagreements=proof['reported_reader_disagreements'],
                                    source_occurrences=item.get('source_occurrences', []),
                                    biological_or_comparative_admission=False))
            except (OSError, ValueError, KeyError, TypeError) as error:
                errors.append(str(item.get('source_url')) + ': ' + str(error))
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append(str(root) + ': ' + str(error))
    return entries, errors


def _cancel(event):
    if event is not None and event.is_set():
        raise QueryCancelled('Workbook query cancelled')


def _groups(folder, sheet, cancel=None):
    table = pq.ParquetFile(Path(folder) / 'cells.parquet')
    for index in range(table.num_row_groups):
        _cancel(cancel)
        block = table.read_row_group(index, use_threads=False)
        block = block.filter(pc.equal(block['sheet_index'], sheet))
        if block.num_rows:
            yield block


@dataclass
class WorkbookQuery:
    """A checked row selection with a bounded preview and complete export route."""
    folder: Path
    sheet: int
    source_rows: tuple[int, ...]
    preview: list[dict]
    record: dict
    metadata_folder: Path | None = None

    def blocks(self, cancel=None):
        """Stream complete physical records, preserving original order and duplicates."""
        _cancel(cancel)
        if _sha(self.folder / 'manifest.json') != self.record['manifest_sha256']:
            raise ValueError('Workbook manifest changed since query')
        manifest(self.folder, verify=True)
        selected = pa.array(self.source_rows, type=pa.int32())
        for block in _groups(self.folder, self.sheet, cancel):
            block = block.filter(pc.is_in(block['row'], value_set=selected))
            if block.num_rows:
                yield block
        if _sha(self.folder / 'manifest.json') != self.record['manifest_sha256']:
            raise ValueError('Workbook manifest changed during query')
        manifest(self.folder, verify=True)


def prepare_query(folder, sheet, labels=(), *, label_column=0, preview_rows=20,
                  preview_cells=500, cancel: threading.Event | None = None):
    """Select original rows by literal decoded cell text, or preview source rows.

    Source numeric tokens, formula caches and inferred gene identifiers are never
    substituted for text labels. Select a label column to disambiguate occurrences.
    """
    from .optional_workbook_data import resolve
    metadata_folder = Path(folder)
    folder = resolve(folder)
    _cancel(cancel)
    digest = _sha(folder / 'manifest.json')
    info = manifest(folder, verify=True)
    if _sha(folder / 'manifest.json') != digest:
        raise ValueError('Workbook manifest changed during verification')
    if _sha(metadata_folder / 'manifest.json') != digest:
        raise ValueError('Workbook discovery metadata and source payload differ')
    available = {s['sheet_index']: s for s in info['extraction']['sheets']}
    if type(sheet) is not int or sheet not in available:
        raise ValueError('Unknown original worksheet')
    if type(label_column) is not int or label_column < 0 or preview_rows < 1 or preview_cells < 0:
        raise ValueError('Invalid source row query limits')
    labels = list(dict.fromkeys(labels))
    if any(not isinstance(label, str) for label in labels):
        raise ValueError('Workbook labels must be literal strings')
    selected, found, occurrences = set(), set(), 0
    for block in _groups(folder, sheet, cancel):
        if labels:
            mask = pc.is_in(block['decoded_text'], value_set=pa.array(labels, type=pa.string()))
            mask = pc.and_(mask, pc.is_in(block['storage_type'],
                          value_set=pa.array(['s', 'inlineStr', 'str', 'BIFF:text'])))
            if label_column:
                mask = pc.and_(mask, pc.equal(block['column'], label_column))
            matches = block.filter(mask)
            selected.update(matches['row'].to_pylist())
            found.update(matches['decoded_text'].to_pylist())
            occurrences += matches.num_rows
        else:
            for row in block['row'].to_pylist():
                if row in selected or len(selected) < preview_rows:
                    selected.add(row)
            if len(selected) >= preview_rows:
                break
    record = dict(schema_version=1, query_epoch=time.time(), source_url=info['source_url'],
                  input_sha256=info['input_sha256'], manifest_sha256=digest,
                  source_record=info['source_record'], sheet=available[sheet], labels=labels,
                  label_column=label_column, match_mode='literal decoded text; no aliases, numeric coercion or canonical mapping',
                  unmatched_labels=[label for label in labels if label not in found], matching_text_cells=occurrences,
                  selected_source_row_coordinates=sorted(selected), preview_row_limit=preview_rows,
                  validation=validation(metadata_folder, info), source_storage_facts=info['extraction'],
                  original_treatment=info['treatment'], biological_or_comparative_admission=False,
                  quantity_units_and_sample_roles='Unreviewed; rows can contain headers, annotations, measurements and statistics',
                  operations=['Choose one original worksheet', 'Select source row coordinates by exact text or first source-row preview',
                              'Retain every stored cell in selected coordinates, in physical order, including repeated records',
                              'No normalization, aggregation, formula calculation, zero filling, sample matching or gene mapping'],
                  numeric_representation=info['extraction']['numeric_representation'],
                  presence_semantics='Missing fields use explicit presence columns in CSV; absent worksheet cells are not manufactured')
    result = WorkbookQuery(folder, sheet, tuple(sorted(selected)), [], record, metadata_folder)
    count = 0
    for block in result.blocks(cancel):
        count += block.num_rows
        if len(result.preview) < preview_cells:
            result.preview.extend(block.slice(0, preview_cells - len(result.preview)).to_pylist())
    _cancel(cancel)
    record['selected_physical_cells'] = count
    record['preview_cells'] = len(result.preview)
    return result


def export_csv(result, destination, *, cancel=None):
    """Save all matching source cells atomically, with source and treatment records."""
    destination = Path(destination)
    if destination.exists() or destination.with_suffix(destination.suffix + '.json').exists():
        raise FileExistsError('Choose a new workbook export filename')
    protected = [result.folder]
    if result.metadata_folder is not None:
        protected.append(Path(result.metadata_folder))
    if any(destination.resolve().is_relative_to(folder.resolve()) for folder in protected):
        raise ValueError('Exports must be outside the immutable workbook source directory')
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary_name = tempfile.mkstemp(prefix='.' + destination.name + '-', suffix='.partial', dir=destination.parent)
    temporary = Path(temporary_name)
    sidecar = destination.with_suffix(destination.suffix + '.json')
    side_temp = temporary.with_suffix('.json.partial')
    count = 0
    try:
        with os.fdopen(handle, 'w', encoding='utf-8', newline='') as stream:
            fields = list(pq.ParquetFile(result.folder / 'cells.parquet').schema_arrow.names)
            presence = ['source_value_token', 'decoded_text', 'numeric_binary64_hex', 'formula_xml', 'source_address']
            writer = csv.DictWriter(stream, fieldnames=fields + [key + '_present' for key in presence])
            writer.writeheader()
            for block in result.blocks(cancel):
                for row in block.to_pylist():
                    row.update({key + '_present': row[key] is not None for key in presence})
                    writer.writerow(row)
                    count += 1
        _cancel(cancel)
        if count != result.record['selected_physical_cells']:
            raise ValueError('Workbook export coverage changed')
        record = dict(result.record, exported_physical_cells=count, export_sha256=_sha(temporary))
        side_temp.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
        # Hard-link publication refuses racing filenames rather than replacing them.
        os.link(temporary, destination)
        try:
            os.link(side_temp, sidecar)
        except Exception:
            destination.unlink()
            raise
        return count
    finally:
        temporary.unlink(missing_ok=True)
        side_temp.unlink(missing_ok=True)
