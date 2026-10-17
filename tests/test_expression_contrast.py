"""The expression contrast: a moderated test that holds its error rates, and normalisations that do what
their papers say -- checked on simulated counts where the truth is known."""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parascape import expression_contrast as X  # noqa: E402


def _simulate(n_genes=6000, n=3, effect=None, seed=0, depth=(8e6, 12e6), dispersion=0.08):
    """Negative-binomial counts: gene means log-normal, per-gene fold changes in `effect` (B over A)."""
    rng = np.random.default_rng(seed)
    base = np.exp(rng.normal(3.0, 1.8, n_genes))
    fold = np.ones(n_genes) if effect is None else effect
    cols, data = [], []
    for g, mult in (("a", np.ones(n_genes)), ("b", fold)):
        for i in range(n):
            mu = base * mult
            mu = mu / mu.sum() * rng.uniform(*depth)
            r = 1 / dispersion
            data.append(rng.negative_binomial(r, r / (r + mu)))
            cols.append(f"{g}{i}")
    return pd.DataFrame(np.array(data).T, index=[f"g{i}" for i in range(n_genes)], columns=cols)


def test_null_p_values_are_uniform_and_false_positives_at_their_nominal_rate():
    counts = _simulate(seed=1)
    for norm in ("tmm", "median_ratio", "upper_quartile", "cpm"):
        c = X.contrast(counts, ["a0", "a1", "a2"], ["b0", "b1", "b2"], normalisation=norm)
        p = c.table["p"].to_numpy()
        assert 0.03 < (p < 0.05).mean() < 0.07, (norm, (p < 0.05).mean())
        assert stats.kstest(p, "uniform").pvalue > 1e-3, norm
        assert (c.table["fdr"] < 0.05).sum() <= 5, "the null produced discoveries"


def test_planted_changes_are_found_and_the_fdr_holds():
    rng = np.random.default_rng(2)
    n = 6000
    effect = np.ones(n)
    changed = rng.choice(n, 600, replace=False)
    effect[changed] = np.where(rng.random(600) < 0.5, 4.0, 0.25)
    counts = _simulate(n_genes=n, effect=effect, seed=3)
    c = X.contrast(counts, ["a0", "a1", "a2"], ["b0", "b1", "b2"], normalisation="tmm")
    calls = set(c.table.loc[c.table["fdr"] < 0.05, "gene"])
    truth = {f"g{i}" for i in changed}
    tested = set(c.table["gene"])
    false = len(calls - truth) / max(len(calls), 1)
    power = len(calls & truth) / len(truth & tested)
    assert false < 0.08, f"observed false discovery proportion {false:.3f}"
    assert power > 0.5, f"power {power:.2f}"
    up = c.table.set_index("gene").loc[[f"g{i}" for i in changed if effect[i] > 1 and f"g{i}" in tested]]
    assert (up["log2fc"] > 0).mean() > 0.95, "fold changes have the wrong sign"


def test_tmm_and_median_ratio_correct_a_composition_bias():
    """Five genes take a third of one group's reads: library-size scaling alone shifts everything else."""
    counts = _simulate(seed=4)
    biased = counts.copy()
    top = biased.index[:5]
    for c in ("b0", "b1", "b2"):
        biased.loc[top, c] = biased[c].sum() // 10
    rest = biased.index[5:]
    for norm, tol in (("tmm", 0.15), ("median_ratio", 0.15)):
        c = X.contrast(biased, ["a0", "a1", "a2"], ["b0", "b1", "b2"], normalisation=norm)
        lfc = c.table.set_index("gene").reindex(rest)["log2fc"].dropna()
        assert abs(lfc.median()) < tol, (norm, lfc.median())
    naive = X.contrast(biased, ["a0", "a1", "a2"], ["b0", "b1", "b2"], normalisation="cpm")
    lfc = naive.table.set_index("gene").reindex(rest)["log2fc"].dropna()
    assert lfc.median() < -0.3, "the bias was not there to correct"


def test_median_ratio_is_the_deseq2_definition():
    counts = pd.DataFrame({"s1": [10, 20, 30, 0], "s2": [20, 40, 60, 5], "s3": [5, 10, 15, 1]})
    f = X.norm_factors(counts, "median_ratio")
    lib = counts.sum()
    size = f * lib
    keep = (counts > 0).all(axis=1)
    geo = np.exp(np.log(counts[keep]).mean(axis=1))
    expected = np.exp(np.log(counts[keep].div(geo, axis=0)).median())
    ratio = size / expected
    assert np.allclose(ratio / ratio.iloc[0], 1.0), "size factors differ from DESeq2's up to a constant"


def test_control_relative_puts_the_controls_at_zero():
    counts = _simulate(n_genes=500, seed=5)
    e = X.normalised(counts, method="tmm")
    rel = X.control_relative(e, ["a0", "a1", "a2"])
    assert np.allclose(rel[["a0", "a1", "a2"]].mean(axis=1), 0.0)
    with pytest.raises(ValueError):
        X.control_relative(e, ["nope"])


def test_moderation_borrows_strength_but_not_from_nothing():
    rng = np.random.default_rng(6)
    s2 = rng.chisquare(4, 5000) / 4 * 0.2
    post, d0, s0 = X.squeeze_variances(s2, 4.0)
    assert d0 > 20, "homogeneous variances should give a strong prior"
    assert np.std(post) < np.std(s2), "moderation did not shrink"
    s2h = np.exp(rng.normal(np.log(0.2), 1.5, 5000))
    _post, d0h, _ = X.squeeze_variances(s2h, 4.0)
    assert d0h < d0, "very heterogeneous variances should give a weaker prior"


def test_a_group_of_one_is_refused():
    counts = _simulate(n_genes=200, seed=7)
    with pytest.raises(ValueError):
        X.contrast(counts, ["a0"], ["b0", "b1"])
