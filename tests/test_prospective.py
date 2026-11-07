"""Frozen prospective claims are graded unchanged, or not at all (instructions 63, 64.38)."""
from __future__ import annotations

import hashlib
import json

import pandas as pd
import pytest

from parascape import prospective as P


def _freeze(tmp_path, n=300):
    claims = pd.DataFrame({"gene_id": [f"g{i}" for i in range(n)], "question": "q",
                           "probability_needed": [1 - i / n for i in range(n)]})
    claims["rank"] = range(1, n + 1)
    path = tmp_path / "claims.tsv"
    claims.to_csv(path, sep="\t", index=False)
    (tmp_path / "manifest.json").write_text(json.dumps({
        "question": "q", "frozen_at": "2026-11-06T00:00:00",
        "claims_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}))
    return claims


def test_a_ranking_that_puts_the_truth_on_top_grades_well(tmp_path):
    _freeze(tmp_path)
    truth = pd.Series({f"g{i}": i < 60 for i in range(0, 250)})   # g250..g299 unmeasured
    g = P.grade(tmp_path, truth)
    assert g["tested"] == 250 and g["unmeasured"] == 50
    assert g["top50_precision"] == 1.0 and g["top100_correct"] == 60
    assert g["auc"] == pytest.approx(1.0)


def test_unmeasured_genes_are_counted_never_scored_as_wrong(tmp_path):
    _freeze(tmp_path)
    truth = pd.Series({"g0": True, "g1": False})
    g = P.grade(tmp_path, truth)
    assert g["top50_tested"] == 2 and g["top50_correct"] == 1 and g["unmeasured"] == 298


def test_a_set_changed_after_the_freeze_is_refused(tmp_path):
    claims = _freeze(tmp_path)
    claims.loc[0, "probability_needed"] = 0.0
    claims.to_csv(tmp_path / "claims.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="changed after the freeze"):
        P.grade(tmp_path, pd.Series({"g0": True}))


def test_the_shipped_freeze_still_matches_its_hash():
    import os
    folder = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results", "prospective_claims_2026_11_06")
    if not os.path.isdir(folder):
        pytest.skip("results are not shipped with an installed copy")
    claims, manifest = P.load(folder)
    assert len(claims) == manifest["n_claims"]
