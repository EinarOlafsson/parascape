"""Annotate literal source flags and optionally exclude marked original rows without assigning biological quality.

The operation applies only to declared protein-group quantities and exact '+'
tokens in named original fields. Missing declarations and unfamiliar values stay
unreviewed and retained. Rows are keyed by measurement manifest and source ordinal,
so identical labels, repeated samples and different source files never share flags.
This is an explicit view operation; original files and numerical values are untouched.
"""
from __future__ import annotations

import hashlib
import json
import math
from numbers import Integral
from pathlib import Path
import re

import pandas as pd

from . import __version__

FIELDS = ('Reverse', 'Contaminant', 'Potential contaminant', 'Only identified by site')
POLICY = dict(schema_version=1, eligible_assay='deposited protein-group quantities',
              original_field_names=list(FIELDS), marked_literal='+',
              unknown_or_missing_flags='Retain without interpreting as clean or biologically valid',
              quantity_scope='Whole original source rows; no gene, sample or replicate independence inferred',
              numerical_transformations='None; no filtering of source files, normalization or imputation',
              biological_quality_score='Not assigned or calibrated')
POLICY_SHA256 = hashlib.sha256(json.dumps(POLICY, sort_keys=True).encode()).hexdigest()
try:
    IMPLEMENTATION_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
except OSError:
    IMPLEMENTATION_SHA256 = None  # Frozen executables can omit Python source.


def apply(frame: pd.DataFrame, *, exclude_marked: bool = False) -> tuple[pd.DataFrame, dict]:
    """Annotate source flags and return an auditable, optional exclusion of literal source-marked rows."""
    if not isinstance(exclude_marked, bool):
        raise ValueError('Source flag exclusion must be an explicit boolean')
    result = frame.copy()
    report = dict(operation='original_source_flag_view', policy=json.loads(json.dumps(POLICY)), policy_sha256=POLICY_SHA256,
                  implementation_sha256=IMPLEMENTATION_SHA256, application_version=__version__,
                  exclude_source_marked_rows=exclude_marked, input_measurements=len(frame), input_source_rows=0,
                  excluded_measurements=0, excluded_source_rows=0, retained_source_rows=0,
                  source_rows_with_no_flag_fields=0, uninterpreted_source_rows=0, excluded_rows=[])
    if frame.empty:
        for name in ('source_flags', 'source_flag_status', 'source_flag_view_policy'):
            result[name] = pd.Series(dtype=object)
        return result, report
    required = {'manifest_sha256', 'source_row', 'feature_label', 'assay'}
    if not required.issubset(frame.columns):
        raise ValueError('Source flag views require original manifest, row, feature and assay identities')
    records, keys = {}, []
    for row in frame.itertuples(index=False):
        if (not isinstance(row.manifest_sha256, str) or not re.fullmatch(r'[0-9a-f]{64}', row.manifest_sha256)
                or not isinstance(row.source_row, Integral) or isinstance(row.source_row, bool) or row.source_row < 0):
            raise ValueError('Source flag row identities are not original validated ordinals')
        key = row.manifest_sha256, int(row.source_row)
        keys.append(key)
        encoded = getattr(row, 'source_row_metadata', None)
        if isinstance(encoded, float) and math.isnan(encoded):
            encoded = None
        if encoded is not None and not isinstance(encoded, str):
            raise ValueError('Source flag annotations must retain their original JSON representation')
        signature = (row.feature_label, row.assay, encoded)
        if key in records:
            if records[key]['signature'] != signature:
                raise ValueError('Repeated source-row flag declarations or source identities differ')
            continue
        metadata = json.loads(encoded) if encoded else None
        if metadata is not None and (not isinstance(metadata, dict) or metadata.get('source_row') != row.source_row
                                     or metadata.get('feature_label') != row.feature_label):
            raise ValueError('Source flag annotation belongs to a different original row')
        fields = metadata.get('source_fields', {}) if metadata else {}
        if not isinstance(fields, dict):
            raise ValueError('Source flag fields must be original named declarations')
        flags = {name: fields[name] for name in FIELDS if name in fields}
        eligible = row.assay == POLICY['eligible_assay']
        marked = [name for name, value in flags.items() if eligible and value == '+']
        unfamiliar = bool(flags) and (not eligible or any(value not in ('', '+') for value in flags.values()))
        status = ('source_marked' if marked else 'not_interpreted_for_assay' if flags and not eligible else
                  'unreviewed_flag_values' if unfamiliar else 'declared_flags_unmarked' if flags else 'flags_not_declared')
        records[key] = dict(signature=signature, flags=flags, status=status, marked_fields=marked,
                            unfamiliar=unfamiliar, exclude=bool(marked) and exclude_marked)
    report.update(input_source_rows=len(records),
                  excluded_source_rows=sum(record['exclude'] for record in records.values()),
                  source_rows_with_no_flag_fields=sum(not record['flags'] for record in records.values()),
                  uninterpreted_source_rows=sum(record['unfamiliar'] for record in records.values()))
    report['retained_source_rows'] = report['input_source_rows'] - report['excluded_source_rows']
    for (manifest, ordinal), record in records.items():
        if record['exclude']:
            report['excluded_rows'].append(dict(manifest_sha256=manifest, source_row=ordinal,
                feature_label=record['signature'][0], source_marked_fields=record['marked_fields']))
    keep = [not records[key]['exclude'] for key in keys]
    report['excluded_measurements'] = len(keep) - sum(keep)
    result['source_flags'] = [json.dumps(records[key]['flags'], ensure_ascii=False, allow_nan=False) for key in keys]
    result['source_flag_status'] = [records[key]['status'] for key in keys]
    result['source_flag_view_policy'] = json.dumps(dict(policy_sha256=POLICY_SHA256,
        exclude_source_marked_rows=exclude_marked, missing_or_unknown_flags_retained=True,
        biological_quality_score_assigned=False), sort_keys=True)
    return result.loc[keep].reset_index(drop=True), report
