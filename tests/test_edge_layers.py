"""THE EDGE LAYERS AS AN OBJECT (instruction 101, C1).

Einar's correction: the network is not co-expression drawn as a graph. Co-expression is one layer
among thirteen, each worth a different amount, and these tests hold the two rules that came out of
measuring that in figure 7.
"""
import numpy as np
import pytest

from parascape import edge_layers as EL


def test_a_layer_that_recovers_nothing_carries_no_weight():
    """cofitness has 88,997 edges and does worse than guessing, so it must not vote."""
    assert EL.WEIGHT["cofitness"] == 0.0
    assert not EL.Layer("cofitness", 88997, EL.WEIGHT["cofitness"]).carries_information
    assert EL.Layer("xlms", 2842, EL.WEIGHT["xlms"]).carries_information


def test_size_is_not_worth():
    """The ordering by weight must not be the ordering by number of edges, or the weighting is
    doing nothing that taking the union would not."""
    by_weight = sorted(EL.WEIGHT, key=lambda k: -EL.WEIGHT[k])
    assert by_weight[0] == "comention", "the smallest useful layer is the most informative"
    assert by_weight[-1] == "cofitness", "the largest layer is the least informative"


def test_the_circular_layer_is_refused_for_the_label_it_is_built_from():
    assert "compartment" in EL.CIRCULAR_FOR
    assert "lopit_map" in EL.CIRCULAR_FOR["compartment"]


def test_a_weighted_vote_is_not_carried_by_the_layer_with_the_most_edges():
    """One worthless layer shouting ten times must lose to one good layer speaking once."""
    labels = np.array(["a", "wrong", "right", "wrong", "wrong"])
    known = np.array([False, True, True, True, True])
    adjacency = {"cofitness": [np.array([1, 3, 4]), *[np.array([])] * 4],
                 "xlms": [np.array([2]), *[np.array([])] * 4]}
    assert EL.weighted_vote(adjacency, labels, known, 0) == "right"


def test_a_gene_with_no_weighted_neighbour_gets_no_label():
    labels = np.array(["a", "b"])
    known = np.array([False, True])
    adjacency = {"cofitness": [np.array([1]), np.array([])]}
    assert EL.weighted_vote(adjacency, labels, known, 0) == "", "a zero-weight layer never votes"


def test_the_weights_say_where_they_came_from():
    assert EL.MEASURED_IN.startswith("network_layers_")
