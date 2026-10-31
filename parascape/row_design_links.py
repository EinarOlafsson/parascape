"""Attach source-bound filename associations to original quantity rows without inferring sample identity."""
from __future__ import annotations

import hashlib
import json
from numbers import Integral
import os
from pathlib import Path
import re

from . import deposited_values as D
from . import paths

try:
    IMPLEMENTATION_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
except OSError:
    IMPLEMENTATION_SHA256 = None


def enrich(pack, frame, *, audit_root=None):
    """Retain all measurements and separately annotate exact filename links and review-only stem candidates."""
    pack = D.runtime_pack(pack)
    info = D.manifest(pack)
    metadata = json.loads((pack / 'source_metadata.json').read_text())
    if (info['assay'] != 'deposited protein-accession quantities'
            or metadata.get('row_layout', {}).get('layout') != 'original_rowwise_quantities'):
        return frame, None
    if not {'manifest_sha256', 'source_row', 'feature_label', 'source_row_metadata'}.issubset(frame.columns):
        raise ValueError('Row-design views require original measurement and annotation identities')
    manifest_hash = hashlib.sha256((pack / 'manifest.json').read_bytes()).hexdigest()
    if not frame.manifest_sha256.eq(manifest_hash).all():
        raise ValueError('Row-design measurement manifest differs from the original query')
    annotations = metadata['row_annotations']
    if len(annotations) != info['features']:
        raise ValueError('Row-design original annotation coverage differs')
    selected = {}
    for row in frame.itertuples(index=False):
        position = row.source_row
        if not isinstance(position, Integral) or isinstance(position, bool) or not 0 <= position < len(annotations):
            raise ValueError('Row-design source identity is not an original ordinal')
        original = annotations[position]
        if (original.get('source_row') != position or original.get('feature_label') != row.feature_label
                or json.loads(row.source_row_metadata) != original):
            raise ValueError('Row-design original query annotation identity differs')
        selected[int(position)] = original['source_fields']['R.FileName']
    root = Path(audit_root) if audit_root is not None else Path(os.environ.get(
        'PARASCAPE_AUDIT', str(Path(paths.data_dir()) / 'processing_audit')))
    collection = root / 'protein_row_design_links'
    record = dict(operation='original_row_filename_association_view', measurement_manifest_sha256=manifest_hash,
        implementation_sha256=IMPLEMENTATION_SHA256,
        selected_original_source_rows=len(selected), biological_independence='not inferred',
        sample_control_pairing_and_input_table_linkage='unreviewed', numerical_values_transformed=False)
    pointer_path = collection / 'collection.json'
    associations = {}
    if not pointer_path.exists():
        record['association_status'] = 'Source-row association collection is not installed.'
    else:
        pointer = json.loads(pointer_path.read_text())
        directory = pointer['directory']
        if pointer.get('schema_version') != 1 or not re.fullmatch(r'[0-9a-f]{64}', directory):
            raise ValueError('Unsupported row-design collection identity')
        folder = collection / directory
        if folder.resolve().parent != collection.resolve():
            raise ValueError('Row-design collection escapes its installed root')
        receipt_bytes = (folder / 'validation.json').read_bytes()
        if hashlib.sha256(receipt_bytes).hexdigest() != pointer['validation_sha256']:
            raise ValueError('Row-design reconciliation receipt checksum differs')
        receipt = json.loads(receipt_bytes)
        if receipt.get('schema_version') != 1 or receipt['status'] != 'original_row_filename_associations_reconciled':
            raise ValueError('Original row filename associations have not been reconciled')

        def read(name):
            if name != 'index.json' and not re.fullmatch(r'[0-9a-f]{64}\.json', name):
                raise ValueError('Invalid row-design output identity')
            path = folder / name
            if path.resolve().parent != folder.resolve():
                raise ValueError('Row-design output escapes its collection')
            body = path.read_bytes()
            if hashlib.sha256(body).hexdigest() != receipt['outputs'][name]:
                raise ValueError('Row-design output checksum differs')
            return json.loads(body)

        index = read('index.json')
        if index.get('schema_version') != 1 or len(index['manifests']) != receipt['measurement_packs']:
            raise ValueError('Row-design measurement index coverage differs')
        record['collection'] = pointer
        if manifest_hash not in index['manifests']:
            record['association_status'] = 'This exact measurement manifest is not represented in the row-design collection.'
        else:
            name = index['manifests'][manifest_hash]
            report = read(name)
            if (report.get('schema_version') != 1 or report['measurement_manifest_sha256'] != manifest_hash
                    or report['input_sha256'] != info['input_sha256'] or report['source_metadata_sha256'] != info['outputs']['source_metadata.json']
                    or report['study'] != info['study'] or report['source_url'] != info['source_url']
                    or report['source_analysis'] != info['source_analysis'] or report['features'] != info['features']
                    or report['source_filename_field'] != 'R.FileName'):
                raise ValueError('Row-design measurement or original annotation context differs')
            coverage = set()
            for association in report['associations']:
                literal = association['source_filename_literal']
                if literal in associations or not isinstance(literal, str):
                    raise ValueError('Row-design literal filename identities differ')
                exact = association['exact_primary_files']
                stems = association['filename_stem_primary_candidates']
                status = ('source_filename_not_declared' if not literal else 'exact_unique_primary_file' if len(exact) == 1 else
                          'ambiguous_exact_primary_file' if exact else 'unresolved_exact_primary_file')
                if (any(not literal or item['filename'] != literal for item in exact)
                        or any(not literal or item['filename'] != literal + '.raw' for item in stems)
                        or association['primary_file_status'] != status or association['biological_independence'] != 'not inferred'
                        or association['source_input_table_linkage'] != 'unreviewed'):
                    raise ValueError('Row-design exact and candidate filename scope differs')
                candidates = association['filename_stem_submitted_row_candidates']
                expected_status = 'review_required' if stems or candidates else 'no_literal_raw_stem_candidates'
                if association['candidate_status'] != expected_status:
                    raise ValueError('Row-design filename candidates were promoted beyond review scope')
                for key, positions_key, expected in (
                        ('submitted_factor_rows', 'matching_file_column_indices', literal),
                        ('filename_stem_submitted_row_candidates', 'matching_candidate_file_column_indices', literal + '.raw')):
                    for row in association[key]:
                        headers, cells = row['headers'], row['cells']
                        positions = [i for i, (header, value) in enumerate(zip(headers, cells))
                                     if literal and header.strip().lower() == 'comment[data file]' and value == expected]
                        if (len(headers) != len(cells) or any(not isinstance(value, str) for value in [*headers, *cells])
                                or row['source_row_id'] != hashlib.sha256((row['source_sha256'] + ':' + str(row['source_row'])).encode()).hexdigest()
                                or not positions or row[positions_key] != positions):
                            raise ValueError('Row-design original submitted factor occurrence differs')
                for position in association['original_source_rows']:
                    if (not isinstance(position, int) or isinstance(position, bool) or not 0 <= position < len(annotations)
                            or position in coverage or annotations[position]['source_row'] != position
                            or annotations[position]['source_fields']['R.FileName'] != literal):
                        raise ValueError('Row-design original row-to-filename membership differs')
                    coverage.add(position)
                associations[literal] = {key: value for key, value in association.items() if key != 'original_source_rows'}
            if coverage != set(range(info['features'])):
                raise ValueError('Row-design original source-row coverage differs')
            record.update(association_status='Source-bound exact filename links and separately labelled review-only stem candidates.',
                          report_sha256=receipt['outputs'][name], primary_file_metadata_sha256=report['primary_file_metadata_sha256'])
    result = frame.copy()
    provenance = {key: value for key, value in record.items() if key not in ('selected_original_source_rows',)}
    result['row_design_metadata'] = [json.dumps(dict(source_row=int(position), source_filename_field='R.FileName',
        source_filename_literal=selected[int(position)], association=associations.get(selected[int(position)],
            {'association_status': record['association_status'], 'biological_independence': 'not inferred'}),
        association_provenance=provenance), ensure_ascii=False) for position in result.source_row]
    return result, record
