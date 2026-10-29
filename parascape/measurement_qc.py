"""Assay-independent numerical QC without imputing missing measurements or biological labels."""
from __future__ import annotations

import numpy as np
import pandas as pd


def summarize(matrix: pd.DataFrame) -> pd.DataFrame:
    """Describe each sample in a feature-by-sample nonnegative measurement matrix.

    Numeric validity is a format/integrity measure, not calibrated biological quality.
    Missingness and measured zeros remain separate. This function neither transforms
    the input nor infers biological replicates, stage, strain or sample identity.
    """
    if isinstance(matrix.index, pd.MultiIndex) or isinstance(matrix.columns, pd.MultiIndex):
        raise ValueError('QC requires one feature identity and one sample identity per axis')
    if matrix.empty or matrix.index.has_duplicates or matrix.columns.has_duplicates:
        raise ValueError('QC requires nonempty measurements and unique feature/sample identities')
    if matrix.index.isna().any() or matrix.columns.isna().any():
        raise ValueError('QC requires nonmissing feature/sample identities')
    if any(not pd.api.types.is_numeric_dtype(dtype) or pd.api.types.is_bool_dtype(dtype)
           or pd.api.types.is_complex_dtype(dtype)
           for dtype in matrix.dtypes):
        raise ValueError('QC requires real numeric measurements; text, booleans and complex values are not coerced')
    rows = []
    for sample in matrix.columns:
        values = matrix[sample].to_numpy(dtype=float, na_value=np.nan)
        missing = np.isnan(values)
        finite = np.isfinite(values)
        negative = finite & (values < 0)
        valid = finite & ~negative
        observed = int((~missing).sum())
        invalid = int((~finite & ~missing).sum() + negative.sum())
        rows.append(dict(sample=sample, features=len(values), missing=int(missing.sum()),
                         nonfinite_observed=int((~finite & ~missing).sum()),
                         negative=int(negative.sum()), measured_zero=int((valid & (values == 0)).sum()),
                         measured_positive=int((valid & (values > 0)).sum()),
                         observed_coverage_pct=100 * observed / len(values),
                         numeric_validity_observed_pct=100 * int(valid.sum()) / observed if observed else None,
                         status='no_observed_values' if not observed else 'invalid_numeric_values' if invalid else
                                'valid_observed_numeric_values',
                         quality_score=None,
                         quality_score_status='unavailable: no independently calibrated biological quality score',
                         missing_semantics='Unmeasured values are unknown; measured zero is not proof of biological absence',
                         operation='Descriptive numerical QC only; no filtering, imputation, normalization or identity prediction'))
    return pd.DataFrame(rows)
