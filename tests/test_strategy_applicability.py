"""Strategy applicability per space (instruction 53, WP8)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from parascape import strategy_applicability as SA


@pytest.fixture
def bare():
    """A space with genes, one categorical and one numeric column, and nothing else."""
    nodes = pd.DataFrame({"gene_id": ["a", "b", "c"], "product": ["x", "y", "z"],
                          "compartment": ["nucleus", "cytosol", "nucleus"],
                          "length": [100.0, 200.0, 300.0]})
    return SA.evidence_of("Xx", nodes, (), partner=False)


def test_identity_columns_do_not_count_as_evidence(bare):
    """`gene_id` and `product` must not satisfy a strategy asking for a categorical column."""
    assert bare.categorical == 1 and bare.numeric == 1
    assert bare.satisfies("a categorical column")[0]
    assert not bare.satisfies("two categorical columns")[0]
    assert not bare.satisfies("two numeric columns")[0]


def test_a_space_without_a_graph_cannot_run_network_strategies(bare):
    """Every edge requirement fails, and says which layer is missing rather than failing silently."""
    for need in ("an edge layer", "at least one edge layer", "two or more edge layers",
                 "an interaction layer", "a co-mention layer", "a structural-similarity layer"):
        can, why = bare.satisfies(need)
        assert not can and why


def test_the_layer_names_of_the_shipped_graphs_are_recognised():
    """`xlms` and `ip_ms` are interaction evidence: a check that missed them would call Plasmodium bare."""
    nodes = pd.DataFrame({"gene_id": ["a", "b"], "compartment": ["nucleus", "cytosol"]})
    evidence = SA.evidence_of("Xx", nodes, ("xlms__a", "ip_ms__w", "struct__b", "comention__a"), partner=False)
    assert evidence.satisfies("an interaction layer")[0]
    assert evidence.satisfies("a crosslink or IP-MS layer")[0]
    assert evidence.satisfies("a structural-similarity layer")[0]
    assert evidence.satisfies("a co-mention layer")[0]


def test_both_shipped_parasite_spaces_can_run_every_strategy():
    """The strategies were written against these two spaces, so a block here is a bug in the check."""
    import parascape.strategy_catalog  # noqa: F401  (registers them)
    import parascape.strategy_graph  # noqa: F401
    import parascape.strategy_learning  # noqa: F401
    from parascape import paths, strategies

    assert len(strategies.REGISTRY) > 30
    for code, table, graph in (("Tg", "nodes.parquet", "graph.npz"),
                               ("Pf", "pf_nodes.parquet", "pf_graph.npz")):
        nodes = pd.read_parquet(paths.cache_file(table))
        layers = tuple(np.load(paths.cache_file(graph), allow_pickle=True).files)
        frame = SA.applicability(SA.evidence_of(code, nodes, layers))
        assert frame.applicable.all(), f"{code}: {frame[~frame.applicable].missing.tolist()}"


def test_an_unknown_requirement_is_allowed_rather_than_refused(bare):
    """A requirement this module has never seen must not silently block a strategy."""
    can, why = bare.satisfies("a sequencing depth of at least 30x")
    assert can and not why
