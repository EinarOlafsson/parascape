"""Numerical QC retains missingness, measured zeros and unsupported quality claims."""
import numpy as np
import pandas as pd
import pytest

from parascape.measurement_qc import summarize


def test_missing_zero_invalid_and_unobserved_are_distinct_without_mutation():
    data = pd.DataFrame({'partial': [0., 2., np.nan, np.inf, -1.],
                         'unobserved': [np.nan] * 5})
    before = data.copy(deep=True)
    rows = summarize(data).set_index('sample')
    a = rows.loc['partial']
    assert (a.missing, a.measured_zero, a.measured_positive, a.nonfinite_observed, a.negative) == (1, 1, 1, 1, 1)
    assert a.observed_coverage_pct == 80 and a.numeric_validity_observed_pct == 50
    assert a.status == 'invalid_numeric_values'
    assert rows.loc['unobserved', 'status'] == 'no_observed_values'
    assert pd.isna(rows.loc['unobserved', 'numeric_validity_observed_pct'])
    assert rows.quality_score.isna().all()
    pd.testing.assert_frame_equal(data, before)


def test_validity_does_not_penalize_missingness_or_treat_zeros_as_missing():
    rows = summarize(pd.DataFrame({'sample': pd.Series([0, pd.NA, 3], dtype='Float64')}))
    assert rows.iloc[0].numeric_validity_observed_pct == 100
    assert rows.iloc[0].measured_zero == 1 and rows.iloc[0].missing == 1
    assert rows.iloc[0].status == 'valid_observed_numeric_values'


@pytest.mark.parametrize('domain', ['signed', 'unspecified_real'])
def test_signed_or_unspecified_quantities_preserve_negative_values(domain):
    data = pd.DataFrame({'sample': [-2., 0., 3., np.nan]})
    result = summarize(data, domain=domain).iloc[0]
    assert result.negative == 1 and not result.negative_values_invalid
    assert result.numeric_validity_observed_pct == 100
    assert result.status == 'valid_observed_numeric_values'
    assert result.quantity_domain == domain and data.iloc[0, 0] == -2


@pytest.mark.parametrize('data', [pd.DataFrame(), pd.DataFrame([[1, 2]], columns=['a', 'a']),
                                    pd.DataFrame({'a': ['1', '2']}), pd.DataFrame({'a': [True]}),
                                    pd.DataFrame({'a': [1 + 2j]}),
                                    pd.DataFrame({'a': [1]}, index=[None])])
def test_ambiguous_identities_and_nonnumeric_measurements_are_rejected(data):
    with pytest.raises(ValueError):
        summarize(data)
