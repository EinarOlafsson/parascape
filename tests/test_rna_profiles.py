"""Offline profile selection preserves genes, samples, zero/missingness and provenance."""
import json

import numpy as np
import pandas as pd
import pytest

from parascape import rna_profiles as R


@pytest.fixture
def snapshot_root(tmp_path):
    for study in ['study_one', 'study_two']:
        p = tmp_path / 'matrices' / study
        p.mkdir(parents=True)
        tpm = pd.DataFrame({'r1': [0., 20.], 'r2': [10., np.nan]}, index=['TGME49_1', 'TGME49_2'])
        for q in ['salmon', 'kallisto']:
            pd.concat({'count': tpm * 2, 'tpm': tpm, 'length': tpm * 0 + 100}, axis=1).to_parquet(p / f'{q}_ME49_host.parquet')
            pd.DataFrame({'run_accession': ['r1', 'r2'], 'sample_title': ['control', 'treated'],
                          'condition': ['a', 'b']}).to_csv(p / f'{q}_ME49_host_samples.tsv', sep='\t', index=False)
    (tmp_path / 'snapshot.json').write_text(json.dumps({'verified_runs': 4, 'genes': 2}))
    return str(tmp_path)


def test_profiles_preserve_study_run_identity_and_missingness(snapshot_root):
    d = R.profiles('tgme49_1, TGME49_2; absent', root=snapshot_root)
    assert len(d) == 8
    assert d[['gene', 'study', 'run_accession']].duplicated().sum() == 0
    assert (d.tpm == 0).sum() == 2
    assert d.tpm.isna().sum() == 2
    assert 'ABSENT' not in set(d.gene)
    assert R.snapshot(snapshot_root)['verified_runs'] == 4
    assert len(R.profiles('TGME49_1', study='study_one', root=snapshot_root)) == 2


def test_cli_exports_selected_measurements(snapshot_root, tmp_path):
    out = tmp_path / 'profile.csv'
    R.main(['TGME49_1', '--root', snapshot_root, '--output', str(out)])
    assert len(pd.read_csv(out)) == 4
    with pytest.raises(SystemExit):
        R.main(['missing', '--root', snapshot_root])


def test_panel_compares_exports_and_reports_missing(qtbot, snapshot_root, tmp_path):
    from parascape.rna_panel import RNAProfilesPanel
    panel = RNAProfilesPanel(snapshot_root)
    qtbot.addWidget(panel)
    panel.genes.setText('TGME49_1,missing')
    panel.compare()
    assert panel.table.rowCount() == 4
    assert 'MISSING' in panel.summary.text()
    assert panel.save(str(tmp_path / 'out.csv')) == 4
