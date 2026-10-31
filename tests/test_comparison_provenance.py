"""Comparison identity, missingness and exported treatment invariants."""
import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from parascape import comparison_provenance as P
from parascape import expression_contrast as X


@pytest.fixture
def counts():
    rng = np.random.default_rng(71)
    return pd.DataFrame(rng.poisson(30, size=(120, 4)),
                        index=['g' + str(i) for i in range(120)], columns=['a1', 'a2', 'b1', 'b2'])


@pytest.mark.parametrize('a,b', [(['a1', 'a1'], ['b1', 'b2']),
                                 (['a1', 'a2'], ['a2', 'b2']),
                                 (['a1', 'absent'], ['b1', 'b2'])])
def test_invalid_run_groups_are_not_silently_repaired(counts, a, b):
    with pytest.raises(ValueError):
        X.contrast(counts, a, b)


@pytest.mark.parametrize('change', ['nan', 'inf', 'negative', 'empty_library', 'duplicate_gene', 'duplicate_run', 'text'])
def test_invalid_model_input_is_refused(counts, change):
    if change in {'nan', 'inf', 'negative'}:
        counts = counts.astype(float)
        counts.iloc[0, 0] = {'nan': np.nan, 'inf': np.inf, 'negative': -1}[change]
    elif change == 'empty_library':
        counts['a1'] = 0
    elif change == 'duplicate_gene':
        counts.index = ['same'] * len(counts)
    elif change == 'duplicate_run':
        counts.columns = ['a1', 'a1', 'b1', 'b2']
    else:
        counts = counts.astype(str)
    with pytest.raises(ValueError):
        X.contrast(counts, ['a1', 'a2'], ['b1', 'b2'])


def test_known_repeated_samples_refused_and_unknown_is_explicit(counts):
    metadata = pd.DataFrame(dict(run_accession=counts.columns, sample_accession=['s1', 's1', 's3', 's4']))
    with pytest.raises(ValueError, match='share a submitted sample'):
        X.contrast(counts, ['a1', 'a2'], ['b1', 'b2'], sample_metadata=metadata)
    metadata.sample_accession = ['s1', '', None, 's4']
    result = X.contrast(counts, ['a1', 'a2'], ['b1', 'b2'], sample_metadata=metadata)
    design = result.processing_record['design']
    assert design['runs_without_sample_identity'] == ['a2', 'b1']
    assert design['biological_independence'].startswith('Not established')
    assert not result.processing_record['control_pairing_verified']
    for broken in [metadata.iloc[:-1], pd.concat([metadata, metadata.iloc[[0]]])]:
        with pytest.raises(ValueError, match='exactly one row'):
            X.contrast(counts, ['a1', 'a2'], ['b1', 'b2'], sample_metadata=broken)
    with pytest.raises(ValueError, match='within one study'):
        X.contrast(counts, ['a1', 'a2'], ['b1', 'b2'], sample_metadata=metadata.assign(study=['s', 's', 't', 't']))


def test_tpm_axes_and_missingness_are_not_silently_repaired(counts):
    tpm = counts / counts.sum() * 1e6
    for broken in [tpm.iloc[::-1], tpm.rename(columns={'a1': 'missing'}), tpm.assign(a1=np.nan), tpm.astype(str)]:
        with pytest.raises(ValueError):
            X.contrast(counts, ['a1', 'a2'], ['b1', 'b2'], normalisation='tpm', tpm=broken)


def test_undefined_normalization_is_not_exported_as_a_valid_comparison(counts):
    counts.iloc[:, :] = 0
    counts.iloc[:30, 0] = 10; counts.iloc[30:60, 1] = 10
    counts.iloc[60:90, 2] = 10; counts.iloc[90:, 3] = 10
    with pytest.raises(ValueError, match='normalization'):
        X.contrast(counts, ['a1', 'a2'], ['b1', 'b2'], normalisation='median_ratio', min_expression=0)


def test_selection_hashes_and_settings_bound_to_export_without_mutating_input(counts, tmp_path):
    before = counts.copy(deep=True)
    result = X.contrast(counts, ['a1', 'a2'], ['b1', 'b2'], normalisation='cpm')
    record = result.processing_record
    assert record['model_inputs']['counts'] == P.matrix_record(counts)
    assert record['normalization'] == 'cpm' and record['tested_genes'] == len(result.table)
    assert record['excluded_genes'] + record['tested_genes'] == len(counts)
    changed = counts.copy(); changed.iloc[0, 0] += 1
    assert P.matrix_record(changed)['sha256'] != record['model_inputs']['counts']['sha256']
    assert P.matrix_record(counts.iloc[::-1])['sha256'] != record['model_inputs']['counts']['sha256']
    assert P.matrix_record(counts.astype(float))['sha256'] != record['model_inputs']['counts']['sha256']
    pd.testing.assert_frame_equal(counts, before)
    path = tmp_path / 'comparison.tsv'
    assert P.export(result.table, path, record) == len(result.table)
    sidecar = json.loads(path.with_suffix('.tsv.json').read_text())
    assert sidecar['export_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert sidecar['design']['group_a'] == ['a1', 'a2']
    with pytest.raises(FileExistsError):
        P.export(result.table, path, record)
    assert not list(tmp_path.glob('*.partial'))
