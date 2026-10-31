"""Read source-bound sparse coordinate facts without inferring biological quality."""
from __future__ import annotations

import json
import os
from pathlib import Path
import re

from . import data_components as C, paths


def data_root() -> Path:
    """Locate source-bound coordinate QC in the installed or overridden data."""
    audit = Path(os.environ.get('PARASCAPE_AUDIT', str(Path(paths.data_dir())/'processing_audit')))
    return audit/'sparse_coordinate_qc'


def _json(root, item):
    if not re.fullmatch('[a-f0-9]{64}\\.json', item['file']) or item['file'] != item['sha256']+'.json':
        raise ValueError('Invalid coordinate-QC receipt identity')
    path = root/item['file']
    if path.is_symlink() or path.resolve().parent != root.resolve() or C._hash(path) != item['sha256']:
        raise ValueError('Coordinate-QC receipt checksum/location differs')
    return json.loads(path.read_text())


def report(descriptor, root=None) -> dict | None:
    """Require matching numeric artifacts and internally consistent complete/unknown facts."""
    descriptor = Path(descriptor)
    C._digest(descriptor.stem)
    if C._hash(descriptor) != descriptor.stem:
        raise ValueError('Sparse descriptor changed before coordinate QC')
    info = json.loads(descriptor.read_text())
    root = Path(root) if root is not None else data_root()
    if not (root/'current.json').exists():
        return None
    collection = _json(root, json.loads((root/'current.json').read_text()))
    if collection['schema_version'] != 1:
        raise ValueError('Unsupported coordinate-QC collection')
    matches = [item for item in collection['packs'] if item['source_url'] == info['source_url']]
    if not matches:
        return None
    if len(matches) != 1:
        raise ValueError('Repeated coordinate-QC source identity')
    item = matches[0]
    # Per-source filenames are bindings fingerprints, while sha256 pins the
    # actual receipt bytes independently of that filename.
    if not re.fullmatch('[a-f0-9]{64}\\.json', item['file']):
        raise ValueError('Invalid coordinate-QC pack path')
    file = root/item['file']
    if file.is_symlink() or file.resolve().parent != root.resolve() or C._hash(file) != item['sha256']:
        raise ValueError('Coordinate-QC pack checksum/location differs')
    record = json.loads(file.read_text())
    reference = info['coordinate_artifact']
    expected = dict(directory=reference['directory'], manifest_sha256=reference['manifest_sha256'],
                    validation_sha256=reference['validation_sha256'], output_sha256=reference['output_sha256'])
    if (any(record['bindings'][key] != value for key, value in expected.items())
            or record['source_url'] != info['source_url'] or record['source_input_sha256'] != info['input_sha256']
            or record['original_values_changed'] or record['biological_or_comparative_admission']
            or any(record['source_record'][key] != info['source_record'][key]
                   for key in ('studies', 'source_entities', 'primary_xml_sha256'))):
        raise ValueError('Coordinate-QC source/validation context differs')
    facts = record['result']
    C._digest(record['bindings']['processor_sha256'])
    if collection['processor_sha256'] != record['bindings']['processor_sha256']:
        raise ValueError('Coordinate-QC processor identity differs')
    count = info['declared_dimensions'][2]
    if any(type(facts[key]) is not bool for key in ('source_row_then_column_nondecreasing', 'source_column_then_row_nondecreasing', 'global_duplicate_count_proven')):
        raise ValueError('Coordinate-QC order flags are not explicit booleans')
    ordered = facts['source_row_then_column_nondecreasing'] or facts['source_column_then_row_nondecreasing']
    if (facts['original_records'] != count or facts['global_duplicate_count_proven'] != ordered
            or facts['coordinate_aggregation_applied'] or facts['biological_quality_score'] is not None
            or type(facts['adjacent_repeated_records']) is not int
            or not 0 <= facts['adjacent_repeated_records'] <= max(0, count-1)
            or type(facts['adjacent_duplicate_runs']) is not int
            or not 0 <= facts['adjacent_duplicate_runs'] <= facts['adjacent_repeated_records']):
        raise ValueError('Coordinate-QC count or admission scope differs')
    if ordered:
        duplicate = facts['adjacent_repeated_records']
        if (facts['global_repeated_records'] != duplicate or facts['unique_stored_coordinates'] != count-duplicate
                or facts['global_repeated_coordinate_groups'] != facts['adjacent_duplicate_runs']
                or facts['duplicate_record_fraction'] != (duplicate/count if count else None)):
            raise ValueError('Proven coordinate-QC fractions differ')
    elif any(facts[key] is not None for key in ('global_repeated_records','unique_stored_coordinates',
                                               'global_repeated_coordinate_groups','duplicate_record_fraction')):
        raise ValueError('Unknown global duplication was converted to a known count')
    if item['global_duplicate_count_proven'] != facts['global_duplicate_count_proven'] or item['global_repeated_records'] != facts['global_repeated_records']:
        raise ValueError('Coordinate-QC catalogue summary differs')
    return dict(receipt_sha256=item['sha256'],processor_sha256=record['bindings']['processor_sha256'],
                source_input_sha256=record['source_input_sha256'],output_sha256=reference['output_sha256'],facts=facts)
