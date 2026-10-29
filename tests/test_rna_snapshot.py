"""Release snapshots reproduce gene aggregation and reject altered provenance."""
import gzip
import json

import pandas as pd
import pytest

from scripts.build_rna_snapshot import build, sha256
from parascape import expression_contrast as X


@pytest.fixture
def inputs(tmp_path):
    root = tmp_path / 'source'
    ref = root / 'references' / 'ME49_host'
    ref.mkdir(parents=True)
    (ref / 'tx2gene.tsv').write_text('tx1\tTGME49_1\ntx2\tTGME49_1\ntx3\tTGME49_2\nhost\tHOST:1\n')
    (ref / 'reference.json').write_text(json.dumps({'name': 'ME49_host', 'sha256': {}}))
    folder = root / 'rnaseq' / 'study' / 'SRR1'
    folder.mkdir(parents=True)
    row = {'study': 'study', 'run_accession': 'SRR1', 'references': ['ME49_host'],
           'salmon': 'test-salmon', 'kallisto': 'test-kallisto', 'sample_title': 'control rep1'}
    (folder / 'done.json').write_text(json.dumps(row))
    hashes = {}
    for q in X.QUANTIFIERS:
        d = folder / f'{q}_ME49_host'
        d.mkdir()
        p = d / ('quant.sf.gz' if q == 'salmon' else 'abundance.tsv.gz')
        headers = ('Name\tLength\tEffectiveLength\tTPM\tNumReads\n' if q == 'salmon'
                   else 'target_id\tlength\teff_length\test_counts\ttpm\n')
        rows = [('tx1', 100, 100, 100, 2), ('tx2', 300, 300, 300, 3),
                ('tx3', 200, 200, 0, 0), ('host', 500, 500, 600, 5)]
        with gzip.open(p, 'wt') as f:
            f.write(headers)
            for tx, length, effective, tpm, count in rows:
                values = (tx, length, effective, tpm, count) if q == 'salmon' else (tx, length, effective, count, tpm)
                f.write('\t'.join(map(str, values)) + '\n')
        hashes[str(p)] = sha256(p)
    verified = tmp_path / 'verified.json'
    verified.write_text(json.dumps({'study/SRR1': {'sha256': hashes}}))
    plan = tmp_path / 'plan.jsonl'
    plan.write_text(json.dumps(row) + '\n' + json.dumps({**row, 'run_accession': 'SRR2'}) + '\n')
    return root, verified, plan, tmp_path / 'snapshot'


def test_snapshot_aggregates_transcripts_excludes_host_and_keeps_zero(inputs):
    root, verified, plan, out = inputs
    build(root, verified, plan, out)
    manifest = json.loads((out / 'snapshot.json').read_text())
    assert manifest['planned_runs'] == 2 and manifest['verified_runs'] == 1
    assert manifest['studies'][0]['completed_runs'] == 1
    for q in X.QUANTIFIERS:
        c, t, lengths, samples = X.load('study', q, 'ME49_host', str(out))
        assert list(c.index) == ['TGME49_1', 'TGME49_2']
        assert c.loc['TGME49_1', 'SRR1'] == 5
        assert t.loc['TGME49_1', 'SRR1'] == 1e6
        assert t.loc['TGME49_2', 'SRR1'] == 0
        assert lengths.loc['TGME49_1', 'SRR1'] == 250
        assert samples.host_share.iloc[0] == .5
    for relative, digest in manifest['outputs'].items():
        assert sha256(out / relative) == digest
    with pytest.raises(ValueError, match='immutable'):
        build(root, verified, plan, out)


def test_snapshot_rejects_changed_compressed_source(inputs):
    root, verified, plan, out = inputs
    p = root / 'rnaseq/study/SRR1/salmon_ME49_host/quant.sf.gz'
    p.write_bytes(p.read_bytes() + b'changed')
    with pytest.raises(ValueError, match='Changed quantification'):
        build(root, verified, plan, out)
    assert not (out / 'snapshot.json').exists()
