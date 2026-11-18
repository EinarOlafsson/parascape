"""The edge layers as an object, with what each one is measured to be worth.

Instruction 101, item C1, and Einar's correction of 2026-11-18: the network is not co-expression
drawn as a graph. Nodes carry traits and edges carry derived relationships of many KINDS, and
co-expression is one layer among thirteen. That makes the edge layer a first-class thing rather
than a view of the map, and this module is where it lives.

WHY A WEIGHT, AND WHERE IT COMES FROM. Figure 7 measured what each layer is worth by hiding a
compartment and recovering it from a gene's neighbours in that layer alone. The result is not a
ranking for its own sake: it is what makes the union usable. Taking every layer equally lets the
biggest one decide, and the biggest, `cofitness` with 88,997 edges, is WORSE than guessing the
commonest compartment, so unioning it with the 2,842 crosslinks buries them. Weighting by the
measured worth gave 54.8 per cent against 49.1 for the same genes unweighted, fitted on half the
labelled genes and scored on the other half.

The weights are RECORDED HERE WITH THEIR PROVENANCE rather than recomputed at import: they came
from one run against one release, and a number whose origin is not written down is a number
nobody can check.
"""
from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np

#: The run these weights came from, under `results/`.
MEASURED_IN = "network_layers_2026_11_18"

#: What each layer recovers above chance, measured by hiding a compartment and reading the
#: neighbours in that layer alone. Zero means it did no better than guessing the commonest answer.
WEIGHT = {
    "comention": 0.4195, "ip_ms": 0.3545, "xlms": 0.3482, "unwritten_interaction": 0.3200,
    "comention_ft": 0.3064, "orthogroup": 0.3002, "struct": 0.2103, "domain": 0.2091,
    "coexpression": 0.0963, "structural_hole": 0.0633, "cotranslation": 0.0394,
    "cofitness": 0.0000,
}

#: Built FROM the compartment labels, so it must never be used to recover a compartment.
CIRCULAR_FOR = {"compartment": ("compartment", "compartment_best", "lopit_map", "lopit_unified")}


@dataclass(frozen=True)
class Layer:
    """One kind of relationship between genes."""

    name: str
    edges: int
    weight: float

    @property
    def carries_information(self) -> bool:
        """Whether it recovered anything above guessing the commonest answer."""
        return self.weight > 0


def load(cache: str) -> dict:
    """Every edge layer of the parasite graph, as name -> (gene index a, gene index b, weight)."""
    store = np.load(os.path.join(cache, "graph.npz"), allow_pickle=True)
    names = sorted({k.rsplit("__", 1)[0] for k in store.files if "__" in k})
    return {n: (store[f"{n}__a"], store[f"{n}__b"], store[f"{n}__w"]) for n in names}


def layers(cache: str) -> list:
    """The layers with their size and their measured worth, worth first."""
    out = [Layer(name, int(len(a)), WEIGHT.get(name, 0.0))
           for name, (a, _b, _w) in load(cache).items()]
    return sorted(out, key=lambda layer: (-layer.weight, -layer.edges))


def usable_for(target: str, cache: str) -> list:
    """The layers that may be used to recover `target`, excluding any built from it.

    A layer that joins genes BECAUSE they share the label being recovered returns the label from
    itself. Figure 7 reported a perfect score from exactly that before the layer was excluded, and
    it is the one mistake here that looks like the best result.
    """
    banned = {name for name, built_from in CIRCULAR_FOR.items() if target in built_from}
    return [layer for layer in layers(cache) if layer.name not in banned]


def weighted_vote(adjacency: dict, labels, known, gene: int) -> str:
    """The label a gene's neighbours give it, each layer weighted by what it is worth.

    :param adjacency: layer name -> a list of neighbour arrays, one per gene.
    :param labels: one label per gene; only those where `known` is true are read.
    :returns: the winning label, or "" when no weighted layer reaches a labelled neighbour.
    """
    tally: dict = {}
    for name, neighbours in adjacency.items():
        weight = WEIGHT.get(name, 0.0)
        if weight <= 0:
            continue
        for other in neighbours[gene]:
            if other != gene and known[other]:
                tally[labels[other]] = tally.get(labels[other], 0.0) + weight
    if not tally:
        return ""
    return max(tally.items(), key=lambda kv: kv[1])[0]
