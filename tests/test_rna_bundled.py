"""The released expression snapshot is complete within its explicitly partial scope."""
from pathlib import Path
import hashlib

from parascape import expression_contrast as X
from parascape import rna_profiles as R


def test_bundled_snapshot_manifest_matches_all_matrices():
    root = Path(R.__file__).parent / 'data' / 'rnaseq'
    manifest = R.snapshot(str(root))
    available = X.available(str(root))
    assert set(available.reference) == set(manifest['references'])
    assert 0 < manifest['verified_runs'] <= manifest['planned_runs']
    assert len(manifest['studies']) > 0
    assert sum(s['completed_runs'] for s in manifest['studies']) == manifest['verified_runs']
    for relative, expected in manifest['outputs'].items():
        assert hashlib.sha256((root / relative).read_bytes()).hexdigest() == expected
    for q in X.QUANTIFIERS:
        run_count = 0
        for entry in manifest['studies']:
            c, t, length, samples = X.load(entry['study'], q, 'ME49_host', str(root))
            assert c.columns.equals(t.columns) and c.columns.equals(length.columns)
            assert len(samples) == entry['completed_runs']
            assert not c.index.str.startswith('HOST:').any()
            assert (c >= 0).all().all()
            assert (abs(t.sum() - 1e6) < 2).all()
            run_count += len(samples)
        assert run_count == manifest['verified_runs']
