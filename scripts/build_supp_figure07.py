"""Supplementary figure 7a: the design choices behind figure 7, each tested with data.

  A  why the `compartment` layer is excluded: with it, recovery is near-perfect and circular
  B  the half-split used to fit layer weights: five splits, vote vs best single layer
  C  the layer weights across those five splits
  D  weighted vs unweighted vote, same held-out genes, five splits
  E  accuracy by layers reaching a gene, five splits (the trust curve is not one split)
  F  per compartment: held-out accuracy against how common the compartment is
"""
from __future__ import annotations

import os
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402
import build_network_layers as NL                       # noqa: E402
import build_figure07_v2 as F7                          # noqa: E402
import paper_panels as PP                               # noqa: E402

F = "supplementary_figure_7a"
CACHE = F7.CACHE
SEEDS = range(5)


def run(adj, labels, known, baseline, seed):
    n = len(labels)
    rng = np.random.default_rng(seed)
    wl = np.where(known)[0]
    first = np.zeros(n, bool); first[rng.choice(wl, len(wl) // 2, replace=False)] = True
    second = known & ~first
    w = {}
    for k, a in adj.items():
        f = NL.recover_on(a, labels, known, first)
        w[k] = max(0.0, f["accuracy"] - baseline) if f["answered"] else 0.0
    live = [k for k in adj if w[k] > 0]
    rows = []
    for i in np.where(second)[0]:
        tally, plain, reach = {}, {}, []
        for k in live:
            v = [labels[j] for j in adj[k][i] if known[j] and j != i]
            if not v:
                continue
            reach.append(k)
            for lab in v:
                tally[lab] = tally.get(lab, 0) + w[k]; plain[lab] = plain.get(lab, 0) + 1
        if not tally:
            continue
        best = max(reach, key=lambda k: w[k])
        bv = Counter(labels[j] for j in adj[best][i] if known[j] and j != i).most_common(1)[0][0]
        rows.append({"truth": labels[i], "vote": max(tally, key=tally.get),
                     "plain": max(plain, key=plain.get), "best": bv, "n_layers": len(reach)})
    d = pd.DataFrame(rows)
    return d, w


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    nodes = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    labels = nodes[NL.TARGET].astype(str).str.strip().str.lower().to_numpy()
    known = ~np.isin(labels, NL.NOT_A_LABEL)
    every = NL.layers(os.path.join(CACHE, "graph.npz"))
    n = len(nodes)
    baseline = float(pd.Series(labels[known]).value_counts(normalize=True).max())
    srcs = [os.path.join(CACHE, "graph.npz"), (os.path.join(CACHE, "nodes.parquet"), [NL.TARGET])]

    # A: the circular layer
    circ = NL.recover(NL.neighbours(*every[NL.CIRCULAR][:2], n), labels, known)
    usable = {k: v for k, v in every.items() if k != NL.CIRCULAR}
    adj = {k: NL.neighbours(v[0], v[1], n) for k, v in usable.items()}
    d0, _ = run(adj, labels, known, baseline, 0)
    a = pd.DataFrame([{"graph": "with compartment layer (circular)",
                       "correct_pct": circ["accuracy"] * 100, "answered": circ["answered"]},
                      {"graph": "without (figure 7)", "correct_pct": (d0.vote == d0.truth).mean()
                       * 100, "answered": len(d0)}])
    FS.use()
    fa, ax = FS.figure(2, vertical=True)
    FS.bars(ax, a.correct_pct.to_numpy(), ["with", "without"], highlight=["without"])
    ax.axhline(baseline * 100, color=FS.INK, lw=0.6, ls="--", label="commonest compartment")
    ax.set_ylabel("Compartment correct (%)"); ax.set_ylim(0, 105); ax.set_xlim(-0.6, 1.4)
    ax.set_xlabel("Compartment layer")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "a", fig=fa, data=a, stats=a, sources=srcs, what="the circular layer, with and without")

    runs = {s: run(adj, labels, known, baseline, s) for s in SEEDS}
    rb = pd.DataFrame([{"seed": s, "vote": (d.vote == d.truth).mean() * 100,
                        "best_layer": (d.best == d.truth).mean() * 100,
                        "unweighted": (d.plain == d.truth).mean() * 100, "genes": len(d)}
                       for s, (d, _) in runs.items()])

    def seed_rows(cols, labels_, xlabel, legend):
        FS.use()
        nr = len(rb); rh = 0.2
        w, h = 0.9 + 1.6 + 1.7, 0.15 + nr * rh + 0.6
        fig = plt.figure(figsize=(w, h)); ax = fig.add_axes([0.9 / w, 0.6 / h, 1.6 / w, nr * rh / h])
        y = np.arange(nr)
        ax.hlines(y, rb[cols[0]], rb[cols[1]], color=FS.CONTROL, lw=0.6)
        ax.plot(rb[cols[0]], y, "o", color=FS.CONTROL, ms=3.5, alpha=FS.ALPHA, label=legend[0])
        ax.plot(rb[cols[1]], y, "o", color=FS.ACCENT, ms=3.5, alpha=FS.ALPHA, label=legend[1])
        ax.set_yticks(y); ax.set_yticklabels([f"split {s}" for s in rb.seed])
        ax.set_ylim(nr - 0.5, -0.5); ax.tick_params(axis="y", length=0); ax.set_xlabel(xlabel)
        FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
        return fig

    PP.panel(F, "b", fig=seed_rows(["best_layer", "vote"], None, "Compartment correct (%)",
                                   ["best single layer (ref.)", "weighted vote"]),
             data=rb, stats={"vote minus best layer, mean (points)": float((rb.vote -
                             rb.best_layer).mean()), "splits": len(rb)}, sources=srcs,
             what="five half-splits: vote vs best single layer")
    W = pd.DataFrame({s: w for s, (_, w) in runs.items()})
    FS.use()
    nr = len(W); rh = 0.2
    w_, h_ = 1.45 + 1.6 + 0.3, 0.15 + nr * rh + 0.6
    fc = plt.figure(figsize=(w_, h_)); ax = fc.add_axes([1.45 / w_, 0.6 / h_, 1.6 / w_, nr * rh / h_])
    order = list(W.mean(axis=1).sort_values(ascending=False).index)
    for yy, k in enumerate(order):
        ax.plot(W.loc[k] * 100, [yy] * len(W.columns), "o", color=FS.ACCENT, ms=3, alpha=FS.ALPHA)
    ax.set_yticks(range(nr)); ax.set_yticklabels([F7.LAYER_NAMES.get(k, k) for k in order])
    ax.set_ylim(nr - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Layer weight (points over guessing)")
    PP.panel(F, "c", fig=fc, data=(W * 100).rename_axis("layer"),
             stats={"weight sd across splits, max (points)": float((W * 100).std(axis=1).max())},
             sources=srcs, what="layer weights across five splits")
    PP.panel(F, "d", fig=seed_rows(["unweighted", "vote"], None, "Compartment correct (%)",
                                   ["unweighted vote (ref.)", "weighted vote"]),
             data=rb, stats={"weighted minus unweighted, mean (points)": float((rb.vote -
                             rb.unweighted).mean())}, sources=srcs,
             what="weighted vs unweighted vote")
    curves = []
    FS.use()
    fe, ax = FS.square_figure()
    for s, (d, _) in runs.items():
        k = d.n_layers.clip(upper=6)
        c = (d.vote == d.truth).groupby(k).mean() * 100
        curves += [{"split": s, "layers": int(i), "correct_pct": v} for i, v in c.items()]
        ax.plot(c.index, c.values, "-", color=FS.ACCENT, lw=0.8, alpha=FS.ALPHA)
    ax.axhline(baseline * 100, color=FS.INK, lw=0.6, ls="--", label="commonest compartment")
    ax.plot([], [], "-", color=FS.ACCENT, label="one half-split")
    ax.set_xticks(range(1, 7)); ax.set_xticklabels(["1", "2", "3", "4", "5", "6+"])
    ax.set_xlabel("Layers reaching the gene"); ax.set_ylabel("Compartment correct (%)")
    ax.set_ylim(0, 100)
    FS.legend(ax, loc="lower right")
    PP.panel(F, "e", fig=fe, data=pd.DataFrame(curves), stats={"splits": len(runs)},
             sources=srcs, what="accuracy by layers reaching a gene, five splits")
    freq = pd.Series(labels[known]).value_counts(normalize=True) * 100
    acc = (d0.vote == d0.truth).groupby(d0.truth).mean() * 100
    fr = pd.DataFrame({"share_pct": freq, "correct_pct": acc}).dropna()
    FS.use()
    ff, ax = FS.square_figure()
    ax.plot(fr.share_pct, fr.correct_pct, "o", color=FS.ACCENT, ms=4, alpha=FS.ALPHA, lw=0,
            label="compartment")
    ax.plot([0, 100], [0, 100], color=FS.BACKGROUND, lw=0.6, label="accuracy = share")
    ax.set_xlim(0, 60); ax.set_ylim(0, 100)
    ax.set_xlabel("Labelled genes in the compartment (%)"); ax.set_ylabel("Compartment correct (%)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "f", fig=ff, data=fr.rename_axis("compartment"), stats={"compartments": len(fr)},
             sources=srcs, what="per compartment: accuracy against how common it is")
    print(a.to_string(index=False)); print(rb.round(1).to_string(index=False))
    PP.compose(F, "supplementary", "Every design choice behind figure 7, varied with the rest "
               "held fixed.")


if __name__ == "__main__":
    main()
