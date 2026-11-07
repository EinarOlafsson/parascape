"""The inference atlas manifest and query cache (instructions 64.26-64.29)."""
from __future__ import annotations

import pandas as pd

from parascape import atlas as AT, strategies as S


def _ledger():
    return pd.DataFrame({"strategy": ["feature_knn"] * 3, "target": ["compartment"] * 3, "mode": ["together"] * 3,
                         "gene_id": ["g1", "g2", "g3"]})


def test_every_method_has_a_role_and_a_specific_status_for_known_and_unknown_genes():
    claims = pd.DataFrame({"organism": ["Tg"], "generator": ["random_forest"], "target": ["compartment"]})
    m = AT.manifest("Tg", ledger=_ledger(), claims=claims, calibration=lambda key: None)
    assert len(m) == len(S.catalog()) and m["role"].notna().all()
    knn = m.set_index("strategy").loc["feature_knn"]
    assert knn.heldout_genes == 3 and "held-out rows" in knn.known_genes and "on demand" in knn.unknown_genes
    rf = m.set_index("strategy").loc["random_forest"]
    assert rf.deployment_claims == 1 and "deployment claims" in rf.unknown_genes
    assert AT.coverage(m)["every_method_has_a_status"]


def test_a_query_is_run_once_and_then_served_from_the_cache(tmp_path):
    calls = []

    def run():
        calls.append(1)
        return pd.DataFrame({"gene_id": ["g1"], "score": [0.5]})
    first, how1 = AT.cached_or_run("Tg", "feature_knn", run, target="compartment", genes=["g2", "g1"],
                                   directory=tmp_path)
    second, how2 = AT.cached_or_run("Tg", "feature_knn", run, target="compartment", genes=["g1", "g2"],
                                    directory=tmp_path)
    assert (how1, how2) == ("ran", "cached") and len(calls) == 1 and second.equals(first)
    _f, how3 = AT.cached_or_run("Tg", "feature_knn", run, target="compartment", genes=["g1"], directory=tmp_path)
    assert how3 == "ran"                     # a different list is a different query
