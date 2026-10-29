"""Offline profile selection preserves genes, samples, zero/missingness and provenance."""
import json
import hashlib
from pathlib import Path

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


def test_qc_preserves_study_run_and_missing_zero_semantics(snapshot_root):
    report = R.technical_qc(root=snapshot_root)
    assert len(report) == 8 and not report[['study', 'run_accession', 'field']].duplicated().any()
    assert set(report.features) == {2} and report.missing.sum() == 4
    assert report.measured_zero.sum() == 4
    assert set(report.numeric_validity_observed_pct) == {100}
    assert report.quality_score.isna().all()
    assert report.source_verification.str.startswith('external matrix').all()
    assert len(R.technical_qc(study='study_one', root=snapshot_root)) == 4
    assert {json.loads(m)['run_accession'] for m in report.sample_metadata} == {'r1', 'r2'}


def test_qc_rejects_changed_declared_outputs_and_misaligned_sample_identities(snapshot_root):
    root = Path(snapshot_root)
    info = R.snapshot(root)
    info['outputs'] = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (root / 'matrices').glob('*/*')}
    (root / 'snapshot.json').write_text(json.dumps(info))
    assert R.technical_qc(root=snapshot_root).source_verification.str.startswith('matches snapshot').all()
    source = root / 'matrices/study_one/salmon_ME49_host_samples.tsv'
    source.write_text('run_accession\nsomething_else\nr2\n')
    with pytest.raises(ValueError, match='checksum'):
        R.technical_qc(root=snapshot_root)
    info['outputs'][str(source.relative_to(root))] = hashlib.sha256(source.read_bytes()).hexdigest()
    (root / 'snapshot.json').write_text(json.dumps(info))
    with pytest.raises(ValueError, match='identities differ'):
        R.technical_qc(root=snapshot_root)


def test_panel_qc_export_uses_full_selection_without_gene_query(qtbot, snapshot_root, tmp_path):
    from parascape.rna_panel import RNAProfilesPanel
    panel = RNAProfilesPanel(snapshot_root)
    qtbot.addWidget(panel)
    path = tmp_path / 'qc.csv'
    assert panel.save_qc(str(path)) == 8
    panel.study.setCurrentIndex(1)
    assert panel.save_qc(str(path)) == 4
    frame = pd.read_csv(path)
    assert frame.study.nunique() == 1 and set(frame.features) == {2}
    assert frame.matrix_sha256.notna().all() and frame.quality_score.isna().all()
