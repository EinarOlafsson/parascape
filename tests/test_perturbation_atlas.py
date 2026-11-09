"""The perturbation atlas (instruction 71 product 2): query and recovery on a planted store."""
from __future__ import annotations

import numpy as np
import pandas as pd

from parascape import perturbation_atlas as P


def _store(n=2000, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n)]
    effect = rng.normal(size=n)
    cols = {"a1": effect + rng.normal(scale=0.5, size=n),       # study A, the perturbation
            "b1": effect + rng.normal(scale=0.5, size=n),       # study B, the same perturbation
            "a2": rng.normal(size=n), "b2": rng.normal(size=n), "c1": rng.normal(size=n)}
    lfc = pd.DataFrame(cols, index=genes)
    reg = pd.DataFrame({"study": ["A", "B", "A", "B", "C"], "contrast_label": list(cols)},
                       index=list(cols))
    return {"lfc": lfc, "registry": reg, "release": "planted"}


def test_the_same_perturbation_in_another_study_comes_first():
    store = _store()
    hits = P.query(store["lfc"]["a1"], store, exclude_study="A")
    assert hits.iloc[0]["signature"] == "b1" and hits.iloc[0]["rank"] == 1
    assert "a2" not in set(hits["signature"])          # its own study is held out


def test_recovery_reports_the_partner_rank_against_chance():
    out = P.recovery(_store(), {"planted": ["a1", "b1"]})
    assert set(out["query"]) == {"a1", "b1"}
    assert out["top1"].all() and (out["best_partner_rank"] < out["chance_rank"]).all()


def test_too_few_shared_genes_are_not_compared():
    store = _store(n=500)
    assert P.query(store["lfc"]["a1"], store).empty


def test_host_replicate_groups_need_two_studies_and_a_named_perturbation():
    import pandas as pd
    from parascape import perturbation_atlas as P
    reg = pd.DataFrame({"study": ["A", "B", "A", "C", "D"], "species": ["human"] * 5, "axis": ["genotype"] * 5,
                        "treated": ["deltaBFD1", "delta-BFD1", "unstated", "unstated", "MYR1"]},
                       index=["a1", "b1", "a2", "c2", "d1"])
    groups = P.replicate_groups({"registry": reg})
    assert list(groups) == ["human|genotype|deltabfd1"]          # MYR1 is one study; "unstated" names nothing
    assert all(len(set(reg.loc[v, "study"])) >= 2 for v in groups.values())


def test_the_symbol_map_is_one_to_one_within_each_species():
    """A symbol that names several genes in a species is dropped, not guessed (instruction 71)."""
    from parascape import perturbation_atlas as PA

    mapping = PA.host_symbol_map()
    assert set(mapping.species) == {"human", "mouse"}
    for species, group in mapping.groupby("species"):
        assert group.symbol.is_unique, f"{species}: a symbol names more than one gene"
        assert group.gene_id.is_unique


def test_crossing_the_species_puts_human_and_mouse_on_one_index():
    """`cross_species` re-indexes both species on shared symbols and says what the crossing cost."""
    from parascape import perturbation_atlas as PA

    store = PA.host_store()
    crossed = PA.cross_species(store)
    assert list(crossed["lfc"].columns) == list(store["lfc"].columns)
    assert crossed["lfc"].index.is_unique
    assert not any(str(i).startswith("HOST:") for i in crossed["lfc"].index)
    assert crossed["dropped"]["rows_without_symbol"] > 0
    registry = crossed["registry"]
    human = [c for c in crossed["lfc"].columns if registry.loc[c, "species"] == "human"]
    mouse = [c for c in crossed["lfc"].columns if registry.loc[c, "species"] == "mouse"]
    shared = crossed["lfc"][human].notna().any(axis=1) & crossed["lfc"][mouse].notna().any(axis=1)
    assert int(shared.sum()) > 1000, "the species should now meet on thousands of symbols"
