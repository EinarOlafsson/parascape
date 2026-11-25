"""The Infer benchmark harness (scripts/infer_bench): folds, shuffles, q-values, corrections."""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts", "infer_bench"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, ROOT)

import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e2_confounders as E2  # noqa: E402


def test_a_family_is_never_split_across_folds():
    fam = np.array(["a", "a", "b", "c", "c", "c", "d", "e", "f", "f"], dtype=object)
    folds = H.family_folds(fam, 3, seed=0)
    for f in set(fam):
        assert len(set(folds[fam == f])) == 1, f
    assert set(folds) <= {0, 1, 2}


def test_a_gene_without_an_orthogroup_is_its_own_family():
    nodes = pd.DataFrame({"gene_id": ["g1", "g2", "g3"], "orthogroup": ["OG1", "", None]})
    fam = H.families(nodes)
    assert fam[0] == "OG1" and fam[1] != fam[2]


def test_a_shuffle_keeps_every_class_size_and_the_unlabelled_genes():
    y = np.array(["x", "x", "y", None, "z", "y"], dtype=object)
    s = H.shuffled(y, 1)
    assert s[3] is None
    assert sorted(v for v in s if v is not None) == sorted(v for v in y if v is not None)


def test_benjamini_hochberg_matches_the_textbook_values():
    q = E1.bh(np.array([0.01, 0.04, 0.03, 0.005, np.nan]))
    assert np.allclose(q[:4], [0.02, 0.04, 0.04, 0.02]) and np.isnan(q[4])


def test_a_residual_fit_on_the_training_genes_removes_a_linear_confounder():
    rng = np.random.default_rng(0)
    c = rng.normal(size=400)
    x = 3 * c + rng.normal(scale=0.1, size=400)
    r = E2._residual(x[:200], c[:200], x, c, smooth=False)
    assert abs(np.corrcoef(r[200:], c[200:])[0, 1]) < 0.2, "held-out genes keep the confounder"


def test_a_two_valued_confounder_does_not_need_a_cubic():
    c = np.r_[np.zeros(50), np.ones(50)]
    x = c * 2.0
    r = E2._residual(x, c, x, c, smooth=True)
    assert np.allclose(r, 0, atol=1e-6)


def test_the_final_split_holds_whole_families_and_about_a_fifth_of_them():
    import e7_inference as E7
    fam = np.array([f"f{i // 3}" for i in range(3000)], dtype=object)
    final = E7.final_split(fam)
    for f in pd.unique(fam):
        assert len(set(final[fam == f])) == 1
    assert 0.17 < final.mean() < 0.23


def test_perfect_calibration_has_no_calibration_error():
    import e7_inference as E7
    rng = np.random.default_rng(1)
    conf = rng.uniform(size=20000)
    right = rng.uniform(size=20000) < conf
    assert E7.ece(conf, right) < 0.02
    assert E7.ece(np.full(1000, 0.9), np.zeros(1000, bool)) > 0.8


def test_prepare_keeps_uncovered_genes_empty_and_fills_the_rest():
    X = np.array([[1.0, np.nan], [2.0, 4.0], [3.0, 5.0], [np.nan, np.nan]])
    cover = np.array([True, True, True, False])
    for missing in ("median", "mean", "knn", "mask"):
        E = E3p.prepare(X, cover, missing=missing)
        assert np.isnan(E[3]).all() and not np.isnan(E[:3]).any(), missing


import e3_preparation as E3p  # noqa: E402


def test_one_giant_cluster_recovers_no_more_than_shuffled_labels():
    import e5_clustering as E5
    rng = np.random.default_rng(0)
    y = np.where(rng.uniform(size=600) < 0.1, "rare", "common").astype(object)
    clusters = np.zeros(600, int)                         # one cluster holding every gene
    folds = rng.integers(0, 5, 600)
    real = E5.heldout_recovery(clusters, y, folds)
    null = E5.heldout_recovery(clusters, H.shuffled(y, 1), folds)
    assert abs(real - null) < 0.05


def test_a_cluster_per_class_recovers_the_held_out_genes():
    import e5_clustering as E5
    y = np.array(["a"] * 100 + ["b"] * 100, dtype=object)
    clusters = np.r_[np.zeros(100, int), np.ones(100, int)]
    folds = np.tile(np.arange(5), 40)
    assert E5.heldout_recovery(clusters, y, folds) == 1.0
