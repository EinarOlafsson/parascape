"""Validate supplementary source-preservation evidence without biological scoring."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import struct


def _digest(value):
    if not isinstance(value, str) or re.fullmatch('[a-f0-9]{64}', value) is None:
        raise ValueError('Invalid workbook QC digest')


def _count(value):
    if type(value) is not int or value < 0:
        raise ValueError('Invalid workbook QC count')
    return value


def _hex(value):
    if not isinstance(value, str) or re.fullmatch('(?:[a-f0-9]{2})*', value) is None:
        raise ValueError('Invalid original formula bytes')
    return bytes.fromhex(value)


def report(folder, info=None, *, root=None):
    """Check core/local QC metadata against its source manifest and coverage.

    Raw audit artifacts need not be downloaded to read a frozen report. Their
    checksums identify backend evidence; this does not claim client reanalysis.
    """
    from . import workbook_values as W
    folder = Path(folder)
    root = Path(root) if root is not None else folder.parent / 'source_qc'
    index_path = root / 'catalogue.json'
    if not index_path.exists():
        return None
    if root.is_symlink() or index_path.is_symlink():
        raise ValueError('Workbook QC metadata must use regular collection files')
    index_bytes = index_path.read_bytes()
    index = json.loads(index_bytes)
    if (index['schema_version'] != 1 or index['measurement_admission'] is not False
            or index['biological_quality_score'] is not None):
        raise ValueError('Unsupported workbook QC scope')
    reference = index['sources'].get(folder.name)
    if reference is None:
        return None
    _digest(folder.name)
    _digest(reference['sha256'])
    if reference['file'] != reference['sha256'] + '.json':
        raise ValueError('Invalid workbook QC metadata filename')
    path = root / reference['file']
    if path.resolve().parent != root.resolve():
        raise ValueError('Workbook QC metadata escapes its collection')
    W._check(path, reference['sha256'])
    wrapped = json.loads(path.read_text())
    info = info or W.manifest(folder)
    manifest_sha = W._sha(folder / 'manifest.json')
    if (info['extraction']['format'] != 'BIFF_decoded' or wrapped['schema_version'] != 1
            or wrapped['source_directory'] != folder.name):
        raise ValueError('Workbook QC format or source identity differs')
    for field, expected in (('input_sha256', info['input_sha256']),
                            ('measurement_manifest_sha256', manifest_sha)):
        if reference[field] != expected or wrapped[field] != expected or wrapped['raw_report'][field] != expected:
            raise ValueError('Workbook QC source binding differs: ' + field)
    for field in ('numerical_values_changed', 'original_receipts_replaced', 'biological_or_comparative_admission'):
        if wrapped[field] is not False:
            raise ValueError('Unsupported workbook QC transformation or admission')
    if wrapped['biological_quality_score'] is not None:
        raise ValueError('Storage checks cannot establish biological quality')
    _digest(wrapped['raw_audit_directory'])
    _digest(wrapped['raw_report_sha256'])
    _digest(wrapped['assembler_sha256'])
    raw = wrapped['raw_report']
    if (raw['schema_version'] != 1 or raw['status'] != 'all_BIFF_physical_cells_and_decoded_values_reconciled'
            or raw['all_physical_cells_and_decoded_values_verified'] is not True
            or raw['source_and_existing_pack_bytes_changed'] is not False
            or raw['formula_bytes_preserved'] is not True or raw['original_formulas_recalculated'] is not False
            or raw['measurement_admission'] is not False or raw['biological_quality_score'] is not None):
        raise ValueError('Unsupported supplementary workbook QC claim')
    _digest(raw['verifier_sha256'])
    count = _count(raw['physical_cells'])
    if _count(raw['decoded_cells_checked']) != count or count != info['extraction']['cells']:
        raise ValueError('Workbook QC cell coverage differs')
    if _count(raw['global_repeated_cell_coordinates']) != 0:
        raise ValueError('Decoded BIFF model does not prove repeated physical cell preservation')
    if set(raw['outputs']) != {'workbook_stream.bin', 'physical_cells.parquet', 'formula_records.parquet'}:
        raise ValueError('Workbook QC artifact bindings differ')
    for artifact in raw['outputs'].values():
        _digest(artifact['sha256'])
        _count(artifact['size'])
    if raw['original_workbook_stream_sha256'] != raw['outputs']['workbook_stream.bin']['sha256']:
        raise ValueError('Original workbook stream digest differs')
    kinds = raw['physical_cell_kinds']
    if set(kinds) - {'text', 'numeric', 'error', 'blank', 'boolean'} or sum(_count(n) for n in kinds.values()) != count:
        raise ValueError('Workbook QC cell-type coverage differs')
    if sum(_count(n) for n in raw['original_error_code_counts'].values()) != kinds.get('error', 0):
        raise ValueError('Original workbook error coverage differs')
    if sum(_count(n) for n in raw['text_confirmation_paths'].values()) != kinds.get('text', 0):
        raise ValueError('Workbook QC text coverage differs')
    differences = _count(raw['independent_Calamine_text_disagreements'])
    if (differences > kinds.get('text', 0)
            or raw['independent_Calamine_whole_text_agreement'] is not (differences == 0)):
        raise ValueError('Workbook QC reader-agreement claim differs')
    sheets = info['extraction']['sheets']
    if len(sheets) != len(raw['sheets']):
        raise ValueError('Workbook QC worksheet coverage differs')
    for source, checked in zip(sheets, raw['sheets']):
        name = source.get('name', source.get('source_attributes', {}).get('name'))
        if (checked['sheet_index'] != source['sheet_index'] or checked['name'] != name
                or _count(checked['cells']) != source['cells']):
            raise ValueError('Workbook QC source worksheet differs')
    if sum(sheet['cells'] for sheet in raw['sheets']) != count:
        raise ValueError('Workbook QC total worksheet count differs')
    total_formulas = _count(raw['formula_records'])
    preview = wrapped['formula_records_preview']
    limit = _count(wrapped['formula_preview_limit'])
    if len(preview) > min(limit, total_formulas) or wrapped['formula_preview_is_complete'] is not (len(preview) == total_formulas):
        raise ValueError('Workbook QC formula-preview coverage differs')
    seen = set()
    for formula in preview:
        payload = _hex(formula['complete_record_payload_hex'])
        cached, tokens = _hex(formula['cached_result_bytes_hex']), _hex(formula['formula_tokens_hex'])
        if (len(payload) < 22 or len(cached) != 8 or payload[6:14] != cached
                or len(tokens) != int.from_bytes(payload[20:22], 'little') or payload[22:22+len(tokens)] != tokens):
            raise ValueError('Original formula payload/cache/token bytes differ')
        if (formula['sheet_index'] not in {s['sheet_index'] for s in sheets}
                or (formula['row'], formula['column']) != tuple(n + 1 for n in struct.unpack_from('<HH', payload))):
            raise ValueError('Original formula coordinates differ')
        position = formula['sheet_index'], formula['row'], formula['column']
        if position in seen:
            raise ValueError('Repeated source formula coordinates')
        seen.add(position)
        offset = _count(formula['record_offset'])
        if offset + 4 + len(payload) > raw['outputs']['workbook_stream.bin']['size']:
            raise ValueError('Original formula offset exceeds its stream')
    return dict(wrapped, metadata_file_sha256=reference['sha256'],
                metadata_index_sha256=hashlib.sha256(index_bytes).hexdigest())


def formula_for(qc, sheet, cell):
    """Return original supplementary formula bytes for one source cell, if previewed."""
    if qc is None or cell is None:
        return None
    for formula in qc['formula_records_preview']:
        if (formula['sheet_index'], formula['row'], formula['column']) == (sheet, cell['row'], cell['column']):
            value = dict(formula, formula_evaluated=False)
            cached = bytes.fromhex(formula['cached_result_bytes_hex'])
            if cached[6:8] == b'\xff\xff':
                value['cached_result_kind'] = {0:'text', 1:'boolean', 2:'error', 3:'empty text'}.get(cached[0], 'unknown')
                if cached[0] in (1, 2):
                    value['cached_result_code'] = cached[2]
            else:
                value['cached_result_kind'] = 'numeric'
            return value
    return None
