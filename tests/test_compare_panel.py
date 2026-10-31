"""The Compare tab on a synthetic study laid out exactly as the built matrices are."""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parascape import expression_contrast as X  # noqa: E402


@pytest.fixture()
def study_root(tmp_path):
    """A study of 3 wild type + 3 knockout samples, with 300 genes planted 4x up in the knockout."""
    rng = np.random.default_rng(0)
    genes = [f"TGME49_{200000 + i}" for i in range(4000)]
    base = np.exp(rng.normal(3, 1.8, len(genes)))
    up = rng.choice(len(genes), 300, replace=False)
    cols, data = [], []
    for group, mult in (("wt", 1.0), ("ko", None)):
        for i in range(3):
            m = base.copy()
            if mult is None:
                m[up] *= 4
            mu = m / m.sum() * 1e7
            data.append(rng.negative_binomial(12.5, 12.5 / (12.5 + mu)))
            cols.append(f"SRR{group}{i}")
    counts = pd.DataFrame(np.array(data).T, index=genes, columns=cols).astype(float)
    tpm = counts / counts.sum() * 1e6
    d = tmp_path / "matrices" / "GSE_TEST"
    d.mkdir(parents=True)
    pd.concat({"count": counts, "tpm": tpm, "length": counts * 0 + 1500}, axis=1).to_parquet(
        d / "salmon_ME49_host.parquet")
    titles = [f"wild type rep{i + 1}" for i in range(3)] + [f"TRX1 knockout rep{i + 1}" for i in range(3)]
    pd.DataFrame({"run_accession": cols, "sample_title": titles,
                  "condition": ["wild type"] * 3 + ["TRX1 knockout"] * 3,
                  "host_share": [0.6] * 6}).to_csv(d / "salmon_ME49_host_samples.tsv", sep="\t",
                                                   index=False)
    return str(tmp_path), {genes[i] for i in up}


def test_the_matrices_are_listed_and_loaded(study_root):
    root, _up = study_root
    have = X.available(root)
    assert list(have.itertuples(index=False, name=None)) == [("GSE_TEST", "salmon", "ME49_host")]
    counts, tpm, _length, samples = X.load("GSE_TEST", root=root)
    assert counts.shape == (4000, 6) and len(samples) == 6


def test_the_tab_guesses_groups_compares_and_finds_the_planted_genes(study_root, tmp_path):
    from PyQt6 import QtWidgets
    _app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from parascape.compare_panel import ComparePanel
    root, up = study_root
    p = ComparePanel(root=root, titles={"GSE_TEST": "Trx1 knockout"})
    assert p.study.count() == 1 and "Trx1 knockout" in p.study.itemText(0)
    a, b = p.groups()
    assert a == ["SRRwt0", "SRRwt1", "SRRwt2"] and b == ["SRRko0", "SRRko1", "SRRko2"], (a, b)
    res = p.compare()
    assert res is not None
    hit = p.called()
    found = set(hit["gene"])
    assert len(found & up) / len(up) > 0.6, "most planted genes should be called"
    assert len(found - up) / max(len(found), 1) < 0.1, "too many false calls"
    assert (hit.set_index("gene").loc[list(found & up), "log2fc"] > 0).all()
    assert "up" in p.summary.text() and "down" in p.summary.text()
    assert p.save(str(tmp_path / "c.tsv")) == len(res.table)
    import json
    proof = json.loads((tmp_path / 'c.tsv.json').read_text())
    assert proof['study'] == 'GSE_TEST'
    assert proof['design']['group_a'] == a and proof['design']['group_b'] == b
    assert proof['biological_quality_score'] is None
    assert len(proof['design']['runs_without_sample_identity']) == 6
    p.normalisation.setCurrentIndex(p.normalisation.findData('cpm'))
    assert p.result is None and not p.save_btn.isEnabled() and p.result_table.rowCount() == 0
    assert p.save(str(tmp_path / 'stale.csv')) == 0 and not (tmp_path / 'stale.csv').exists()
    assert p.compare() is not None
    p.normalisation.setCurrentIndex(p.normalisation.findData('tpm'))
    assert p.compare() is not None
    assert p.result_table.horizontalHeaderItem(2).text() == 'mean log2(TPM+1)'
    p.set_groups(["SRRwt0"], ["SRRko0", "SRRko1"])
    assert p.compare() is None and "two replicates" in p.summary.text(), "a group of one was compared"
    p.close()
