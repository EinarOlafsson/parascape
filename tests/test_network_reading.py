"""Reading the network (instruction 87): links of the chosen layers, link steps, hit rates."""
import numpy as np
import pandas as pd

from parascape import network_reading as NR


def _links():
    # 0-1-2 in co-fitness (a chain), 3-4 by crosslink, 5 alone
    return pd.DataFrame({"a": [0, 1, 3, 0], "b": [1, 2, 4, 2],
                         "source": ["cofitness", "cofitness", "xlms", "comention"],
                         "strength": [0.9, 0.8, 0.99, 0.2]})


def test_only_the_chosen_layers_above_the_strength_count():
    links = NR.chosen(_links(), ["cofitness"], 0.85)
    assert links[["a", "b"]].values.tolist() == [[0, 1]]


def test_link_steps_from_a_gene():
    adj = NR.adjacency(NR.chosen(_links(), ["cofitness", "xlms"]), 6)
    d = NR.hops(adj, [0])
    assert d[:3].tolist() == [0, 1, 2]
    assert np.isinf(d[3]) and np.isinf(d[5])


def test_a_genes_links_strongest_first():
    mine = NR.links_of(NR.chosen(_links(), ["cofitness", "comention"]), 0)
    assert mine["partner"].tolist() == [1, 2]
    assert NR.linkage(NR.chosen(_links(), ["cofitness"]), 1).per_layer == [("cofitness", 2)]


def test_agreement_reads_direct_links_against_chance():
    links = pd.DataFrame({"a": [0, 2, 4], "b": [1, 3, 5], "source": ["x"] * 3,
                          "strength": [1.0] * 3})
    adj = NR.adjacency(links, 6)
    result = NR.agreement(adj, ["a", "a", "b", "b", "a", "b"])
    assert result.n == 6 and abs(result.accuracy - 4 / 6) < 1e-9 and result.chance == 0.5


def test_unlinked_genes_ring_the_outside_of_a_link_layout():
    rng = np.random.default_rng(0)
    a = rng.integers(0, 40, 200)
    b = rng.integers(0, 40, 200)
    links = pd.DataFrame({"a": a, "b": b, "source": "x", "strength": 1.0})
    links = links[links.a != links.b]
    xy = NR.spectral_layout(NR.adjacency(links, 50))
    assert xy.shape == (50, 2)
    outside = np.linalg.norm(xy[40:], axis=1)
    assert np.allclose(outside, 1.6)


def test_map_layout_is_two_dimensional_and_centred():
    xyz = np.random.default_rng(1).normal(size=(100, 3))
    xy = NR.map_layout(xyz)
    assert xy.shape == (100, 2) and abs(np.median(xy[:, 0])) < 0.5
