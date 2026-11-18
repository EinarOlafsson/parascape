"""What the whole network says that no single kind of evidence does.

Figure 7 of the ParaScape paper (instruction 101). Einar's correction of 2026-11-18 defines the
object: this is not co-expression drawn as a graph. Nodes carry TRAITS and edges carry DERIVED
RELATIONSHIPS of many kinds -- co-expression, shared fitness behaviour, crosslinks, pull-downs,
being translated together, sharing a domain, being named in the same paper -- and co-expression is
one layer among thirteen.

THE TEST. Take a trait a gene carries, hide it, and try to recover it from the gene's neighbours in
one layer at a time, then from all layers at once. A layer that recovers it is carrying information
about that trait; the question the figure answers is whether the union carries more than its best
single layer, which is the only honest way to say that combining evidence was worth doing.

THE CIRCULAR LAYER IS EXCLUDED BY NAME. The `compartment` layer joins genes BECAUSE they share a
compartment, so using it to predict compartment would recover the label from itself and report 100
per cent. It is dropped, and the drop is stated rather than assumed, because this is the single
mistake that would make the figure meaningless while looking like its best result.

Writes `results/network_layers_<date>/`.
"""
from __future__ import annotations

import json
import os
from collections import Counter
from datetime import date

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: The trait recovered, and the edge layer that is built from it and must not be used.
TARGET, CIRCULAR = "compartment_best", "compartment"

#: Labels that are not a compartment.
NOT_A_LABEL = ("unassigned", "nan", "", "none")

#: A layer that reaches fewer genes than this is reported as too sparse rather than scored.
MIN_REACHED = 50


def layers(path: str) -> dict:
    """Every edge layer as (gene index a, gene index b, weight)."""
    store = np.load(path, allow_pickle=True)
    names = sorted({k.rsplit("__", 1)[0] for k in store.files if "__" in k})
    return {n: (store[f"{n}__a"], store[f"{n}__b"], store[f"{n}__w"]) for n in names}


def neighbours(a, b, n_genes: int) -> list:
    """Undirected adjacency as a list of neighbour arrays."""
    out = [[] for _ in range(n_genes)]
    for i, j in zip(a, b):
        out[i].append(j)
        out[j].append(i)
    return [np.array(v, dtype=int) for v in out]


def recover_on(adjacency, labels: np.ndarray, known: np.ndarray, score_on: np.ndarray) -> dict:
    """:func:`recover`, restricted to the genes in `score_on`."""
    right = answered = 0
    for i in np.where(score_on)[0]:
        votes = [labels[j] for j in adjacency[i] if known[j] and j != i]
        if not votes:
            continue
        answered += 1
        right += int(Counter(votes).most_common(1)[0][0] == labels[i])
    return {"answered": answered, "right": right,
            "accuracy": right / answered if answered else float("nan")}


def recover(adjacency, labels: np.ndarray, known: np.ndarray) -> dict:
    """Predict each known label from its neighbours' labels, never from its own.

    The gene being scored is hidden: only the labels of its neighbours are read, which is what
    makes this a held-out measurement rather than a description of the graph.
    """
    right = answered = 0
    for i in np.where(known)[0]:
        votes = [labels[j] for j in adjacency[i] if known[j] and j != i]
        if not votes:
            continue
        answered += 1
        right += int(Counter(votes).most_common(1)[0][0] == labels[i])
    return {"answered": answered, "right": right,
            "accuracy": right / answered if answered else float("nan")}


def weighted_recover(adjacencies: dict, weights: dict, labels: np.ndarray,
                     known: np.ndarray, score_on: np.ndarray) -> dict:
    """Combine the layers by voting, each layer weighted by what it is worth.

    The naive union lets the biggest layer decide: cofitness carries 88,997 edges and is WORSE
    than guessing the commonest compartment, so unioning it with the crosslinks buries them. A
    vote weighted by each layer's own measured accuracy keeps the sparse, reliable evidence
    audible.
    """
    right = answered = 0
    for i in np.where(score_on)[0]:
        tally = {}
        for name, adjacency in adjacencies.items():
            weight = weights.get(name, 0.0)
            if weight <= 0:
                continue
            for j in adjacency[i]:
                if known[j] and j != i:
                    tally[labels[j]] = tally.get(labels[j], 0.0) + weight
        if not tally:
            continue
        answered += 1
        right += int(max(tally.items(), key=lambda kv: kv[1])[0] == labels[i])
    return {"answered": answered, "right": right,
            "accuracy": right / answered if answered else float("nan")}


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    nodes = pd.read_parquet(os.path.join(cache, "nodes.parquet"))
    labels = nodes[TARGET].astype(str).str.strip().str.lower().to_numpy()
    known = ~np.isin(labels, NOT_A_LABEL)
    n_genes = len(nodes)
    print(f"{n_genes:,} genes, {int(known.sum()):,} with a {TARGET}")

    every = layers(os.path.join(cache, "graph.npz"))
    print(f"{len(every)} edge layers; dropping {CIRCULAR!r}, which is built from the label itself")
    usable = {k: v for k, v in every.items() if k != CIRCULAR}

    baseline = float(pd.Series(labels[known]).value_counts(normalize=True).max())
    rows = []
    for name, (a, b, _w) in usable.items():
        adjacency = neighbours(a, b, n_genes)
        score = recover(adjacency, labels, known)
        reached = sum(1 for i in np.where(known)[0] if len(adjacency[i]))
        rows.append({"layer": name, "edges": int(len(a)), "reached": reached, **score,
                     "lift": score["accuracy"] - baseline if score["answered"] else np.nan,
                     "enough": reached >= MIN_REACHED})
    # The union: every layer at once, which is the graph a person actually looks at.
    all_a = np.concatenate([v[0] for v in usable.values()])
    all_b = np.concatenate([v[1] for v in usable.values()])
    union = neighbours(all_a, all_b, n_genes)
    score = recover(union, labels, known)
    rows.append({"layer": "ALL LAYERS AT ONCE", "edges": int(len(all_a)),
                 "reached": sum(1 for i in np.where(known)[0] if len(union[i])), **score,
                 "lift": score["accuracy"] - baseline, "enough": True})

    # WEIGHTS FITTED ON ONE HALF OF THE GENES, SCORED ON THE OTHER. Weighting each layer by its
    # own accuracy and then scoring on the same genes would be fitting and reporting on one set.
    rng = np.random.default_rng(0)
    with_label = np.where(known)[0]
    first = np.zeros(n_genes, dtype=bool)
    first[rng.choice(with_label, size=len(with_label) // 2, replace=False)] = True
    second = known & ~first
    adjacencies = {name: neighbours(v[0], v[1], n_genes) for name, v in usable.items()}
    weights = {}
    for name, adjacency in adjacencies.items():
        fitted = recover_on(adjacency, labels, known, first)
        weights[name] = max(0.0, fitted["accuracy"] - baseline) if fitted["answered"] else 0.0
    weighted = weighted_recover(adjacencies, weights, labels, known, second)
    plain = weighted_recover(adjacencies, {k: 1.0 for k in adjacencies}, labels, known, second)
    rows.append({"layer": "ALL LAYERS, WEIGHTED", "edges": int(len(all_a)),
                 "reached": weighted["answered"], **weighted,
                 "lift": weighted["accuracy"] - baseline, "enough": True})

    table = pd.DataFrame(rows).sort_values("accuracy", ascending=False)
    print(f"\ncommonest-compartment baseline: {baseline:.1%}\n")
    print(f"  {'layer':24s} {'edges':>8s} {'answered':>9s} {'right':>7s} {'lift':>7s}")
    for row in table.itertuples():
        if not row.enough:
            print(f"  {row.layer:24s} {row.edges:8,}  too sparse: reached {row.reached} genes")
            continue
        print(f"  {row.layer:24s} {row.edges:8,} {row.answered:9,} {row.accuracy:7.1%} "
              f"{row.lift:+7.1%}")

    scored = table[table["enough"] & (table["layer"] != "ALL LAYERS AT ONCE")]
    best = scored.iloc[0]
    union_row = table[table["layer"] == "ALL LAYERS AT ONCE"].iloc[0]
    print(f"\nBest single layer: {best['layer']} at {best['accuracy']:.1%} on "
          f"{best['answered']:,} genes.")
    print(f"All layers at once: {union_row['accuracy']:.1%} on {union_row['answered']:,} genes "
          f"-- {union_row['answered'] - best['answered']:+,} genes reached and "
          f"{union_row['accuracy'] - best['accuracy']:+.1%} accuracy.")
    print(f"\nWeighting each layer by what it is worth, fitted on half the labelled genes and "
          f"scored on\nthe other half: {weighted['accuracy']:.1%} against {plain['accuracy']:.1%} "
          f"for the same genes unweighted, on {weighted['answered']:,} genes.")
    print("Reach and accuracy are reported together: a layer that answers for a tenth of the "
          "genes\nhas not beaten one that answers for all of them by being right more often on "
          "its tenth.")

    out = os.path.join(ROOT, "results", f"network_layers_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "layers.tsv"), sep="\t", index=False)
    json.dump({"target": TARGET, "excluded_layer": CIRCULAR, "genes": n_genes,
               "genes_with_label": int(known.sum()), "baseline": baseline,
               "layers": table.to_dict("records"),
               "layer_weights": {k: float(v) for k, v in weights.items()},
               "weighted_on_held_out": weighted, "unweighted_same_genes": plain},
              open(os.path.join(out, "summary.json"), "w"), indent=2, default=float)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
