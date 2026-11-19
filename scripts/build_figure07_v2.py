"""Figure 7, one claim: in a graph whose genes carry traits and whose edges carry many kinds of
derived relationship, a placement is as trustworthy as the number of independent kinds of
evidence behind it.

Result (2026-11-19, held-out half): the weighted vote is barely better than the best single layer
reaching each gene (54.8 % against 53.7 % correct), so combining does not beat the best layer
gene by gene. What combining buys is a measure of trust: 39.8 % correct when one layer reaches
the gene, 79.8 % when six or more do (commonest-compartment baseline 32.7 %).

Einar, 2026-11-19: placeholders unhelpful, graph choices poor. All panels real. The test (as in
build_network_layers.py): hide a gene's compartment, recover it from its neighbours. The
`compartment` layer is excluded (built from the label). Layer weights are fitted on one half of
the labelled genes and everything below is scored on the OTHER half.

  A  per layer: accuracy (x) and genes answered (dot area), with the weighted vote (the trade-off)
  B  which genes each layer reaches: overlap between layers (map)                   (why combine)
  C  per compartment: best single layer reaching the gene vs weighted vote          (where it helps)
  D  accuracy by how many layers reach the gene                                     (more = better)
  E  accuracy by the winning share of the weighted vote                             (confidence)
  F  genes placed against accuracy, most layers first                             (what you get)
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

OUT = os.path.join(ROOT, "results", "figure07_v2")
PNL = os.path.join(OUT, "panels")
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
LAYER_NAMES = {"comention": "co-mentioned (curated)", "comention_ft": "co-mentioned (full text)",
               "ip_ms": "pull-down (IP-MS)", "xlms": "crosslink (XL-MS)",
               "unwritten_interaction": "predicted interaction", "orthogroup": "orthogroup",
               "struct": "similar structure", "domain": "shared domain",
               "structural_hole": "structural hole", "coexpression": "co-expression",
               "cotranslation": "co-translation", "cofitness": "co-fitness"}
COMP_NAMES = {"ribosome_rnp": "ribosome", "apical_secretory": "apical secretory",
              "pellicle_pm_imc": "pellicle", "golgi_secretory_vesicular": "Golgi", "er": "ER"}


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def vote(adj, i, labels, known):
    v = [labels[j] for j in adj[i] if known[j] and j != i]
    return Counter(v).most_common(1)[0][0] if v else None


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    nodes = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    labels = nodes[NL.TARGET].astype(str).str.strip().str.lower().to_numpy()
    known = ~np.isin(labels, NL.NOT_A_LABEL)
    n = len(nodes)
    every = NL.layers(os.path.join(CACHE, "graph.npz"))
    usable = {k: v for k, v in every.items() if k != NL.CIRCULAR}
    adj = {k: NL.neighbours(v[0], v[1], n) for k, v in usable.items()}
    baseline = float(pd.Series(labels[known]).value_counts(normalize=True).max())

    rng = np.random.default_rng(0)                       # the same split as the builder
    with_label = np.where(known)[0]
    first = np.zeros(n, dtype=bool)
    first[rng.choice(with_label, size=len(with_label) // 2, replace=False)] = True
    second = known & ~first
    weights = {}
    for k, a in adj.items():
        f = NL.recover_on(a, labels, known, first)
        weights[k] = max(0.0, f["accuracy"] - baseline) if f["answered"] else 0.0
    live = [k for k in adj if weights[k] > 0]
    print(f"{len(adj)} layers, {len(live)} with positive weight; baseline {baseline:.3f}")

    rows = []
    for i in np.where(second)[0]:
        tally, reach = {}, []
        per = {}
        for k in live:
            p = vote(adj[k], i, labels, known)
            if p is None:
                continue
            reach.append(k); per[k] = p
            for j in adj[k][i]:
                if known[j] and j != i:
                    tally[labels[j]] = tally.get(labels[j], 0.0) + weights[k]
        if not tally:
            continue
        win, wv = max(tally.items(), key=lambda kv: kv[1])
        best = max(reach, key=lambda k: weights[k])
        rows.append({"gene": i, "truth": labels[i], "pred": win,
                     "share": wv / sum(tally.values()), "n_layers": len(reach),
                     "best_layer": best, "best_pred": per[best], "layers": reach})
    d = pd.DataFrame(rows)
    d["right"] = d["pred"] == d["truth"]
    d["best_right"] = d["best_pred"] == d["truth"]
    print(f"held-out genes answered {len(d)}; weighted {d.right.mean():.3f}, best single layer "
          f"reaching each gene {d.best_right.mean():.3f}")

    # A: per-layer answered vs accuracy on held-out half
    per_layer = []
    for k in adj:
        r = NL.recover_on(adj[k], labels, known, second)
        if r["answered"] >= 25:
            per_layer.append((k, r["answered"], r["accuracy"]))
    pl = pd.DataFrame(per_layer, columns=["layer", "answered", "acc"])
    # Named rows, not a labelled scatter: twelve labels piled up where the accurate sparse
    # layers sit. Accuracy along x; how many genes the layer answers is the dot area.
    order_a = sorted(adj, key=lambda k: -weights[k])
    pl = pl.set_index("layer").reindex([k for k in order_a if k in set(pl["layer"])])
    pl.loc["weighted vote"] = [len(d), d.right.mean()]
    FS.use()
    rh, nr = 0.2, len(pl)
    w, h = 1.55 + 1.6 + 1.7, 0.15 + nr * rh + 0.6
    figa = plt.figure(figsize=(w, h))
    ax = figa.add_axes([1.55 / w, 0.6 / h, 1.6 / w, nr * rh / h])
    y = np.arange(nr)
    size = lambda a: 4 + 60 * np.sqrt(np.asarray(a, float) / len(d))
    cols = [FS.ACCENT if k == "weighted vote" else FS.CONTROL for k in pl.index]
    ax.scatter(pl["acc"] * 100, y, s=size(pl["answered"]), c=cols, alpha=FS.ALPHA, lw=0)
    ax.axvline(baseline * 100, color=FS.INK, lw=0.6, ls="--")
    ax.set_yticks(y); ax.set_yticklabels([LAYER_NAMES.get(k, k) for k in pl.index])
    ax.set_ylim(nr - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    ax.set_xlim(20, 85); ax.set_xlabel("Compartment correct (%)")
    for a_ in (100, 1000):
        ax.scatter([], [], s=size(a_), c=FS.CONTROL, alpha=FS.ALPHA, lw=0,
                   label=f"{a_:,} genes answered")
    ax.plot([], [], color=FS.INK, lw=0.6, ls="--", label="commonest compartment")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    # B: overlap of genes reached (Jaccard), layers ordered by weight
    order = sorted(adj, key=lambda k: -weights[k])
    reached = {k: set(np.where([len(x) > 0 for x in adj[k]])[0]) for k in order}
    J = np.array([[len(reached[a] & reached[b]) / max(1, len(reached[a] | reached[b]))
                   for b in order] for a in order]) * 100
    FS.use()
    nl = len(order); cell = 0.16
    w, h = 1.6 + nl * cell + 0.9, 0.15 + nl * cell + 1.6
    figb = plt.figure(figsize=(w, h))
    ax = figb.add_axes([1.6 / w, 1.6 / h, nl * cell / w, nl * cell / h])
    im = ax.imshow(J, cmap="Blues", vmin=0, vmax=100, interpolation="nearest")
    names = [LAYER_NAMES.get(k, k) for k in order]
    ax.set_xticks(range(nl)); ax.set_xticklabels(names, rotation=90)
    ax.set_yticks(range(nl)); ax.set_yticklabels(names)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = figb.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
    cb.set_label("Genes reached by both (% of either)"); cb.outline.set_linewidth(0.6)

    # C: per compartment, best single layer vs weighted vote
    comp = (d.groupby("truth")[["best_right", "right"]].mean() * 100)
    comp["n"] = d.groupby("truth").size()
    comp = comp[comp["n"] >= 20].sort_values("right", ascending=False)
    FS.use()
    rh, nr = 0.2, len(comp)
    w, h = 1.25 + 1.6 + 1.6, 0.15 + nr * rh + 0.6
    figc = plt.figure(figsize=(w, h))
    ax = figc.add_axes([1.25 / w, 0.6 / h, 1.6 / w, nr * rh / h])
    y = np.arange(nr)
    ax.hlines(y, comp["best_right"], comp["right"], color=FS.CONTROL, lw=0.6)
    ax.plot(comp["best_right"], y, "o", color=FS.CONTROL, ms=3.5, alpha=FS.ALPHA,
            label="best single layer (ref.)")
    ax.plot(comp["right"], y, "o", color=FS.ACCENT, ms=3.5, alpha=FS.ALPHA,
            label="weighted vote")
    ax.set_yticks(y); ax.set_yticklabels([COMP_NAMES.get(c, c) for c in comp.index])
    ax.set_ylim(nr - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 100); ax.set_xlabel("Compartment correct (%)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    # D: accuracy by number of layers reaching the gene
    d["k"] = d["n_layers"].clip(upper=6)
    byk = d.groupby("k")["right"].agg(["mean", "size"])
    FS.use()
    figd, ax = FS.square_figure()
    ax.axhline(baseline * 100, color=FS.INK, lw=0.6, ls="--", label="commonest compartment")
    ax.plot(byk.index, byk["mean"] * 100, "o-", color=FS.ACCENT, ms=4, lw=1.2,
            label="weighted vote")
    ax.set_xticks(list(byk.index))
    ax.set_xticklabels([str(k) if k < 6 else "6+" for k in byk.index])
    ax.set_xlabel("Layers reaching the gene"); ax.set_ylabel("Compartment correct (%)")
    ax.set_ylim(0, 100)
    FS.legend(ax, loc="lower right")

    # E: accuracy by winning share
    bins = np.linspace(0, 1, 11)
    d["sb"] = pd.cut(d["share"], bins, include_lowest=True)
    bys = d.groupby("sb", observed=True)["right"].agg(["mean", "size"])
    mids = [b.mid for b in bys.index]
    FS.use()
    fige, ax = FS.square_figure()
    ax.plot([0, 100], [0, 100], color=FS.BACKGROUND, lw=0.6, label="share = accuracy")
    ax.plot(np.array(mids) * 100, bys["mean"] * 100, "o-", color=FS.ACCENT, ms=4, lw=1.2,
            label="weighted vote")
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.set_xlabel("Winning share of the vote (%)"); ax.set_ylabel("Compartment correct (%)")
    FS.legend(ax, loc="upper left")

    # F: genes placed vs accuracy, most confident first, union vs best single layer
    # Ranked by how many kinds of evidence reach the gene, then by vote share. Ranking by share
    # alone put single-neighbour votes (share 1.0, 61 % right) first.
    s = d.sort_values(["n_layers", "share"], ascending=False)
    cum = s["right"].cumsum().to_numpy() / np.arange(1, len(s) + 1) * 100
    # single layers: rank the best layer's genes by that layer's weight then own vote agreement
    t = d.assign(wt=d["best_layer"].map(weights)).sort_values("wt", ascending=False)
    cum_t = t["best_right"].cumsum().to_numpy() / np.arange(1, len(t) + 1) * 100
    FS.use()
    figf, ax = FS.square_figure()
    ax.plot(np.arange(1, len(t) + 1), cum_t, color=FS.CONTROL, lw=1.2,
            label="best single layer (ref.)")
    ax.plot(np.arange(1, len(s) + 1), cum, color=FS.ACCENT, lw=1.2, label="weighted vote")
    ax.axhline(baseline * 100, color=FS.INK, lw=0.6, ls="--", label="commonest compartment")
    ax.set_xlabel("Held-out genes placed"); ax.set_ylabel("Compartment correct (%)")
    ax.set_ylim(0, 100)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    at80 = int((cum >= 80).sum()); at80_t = int((cum_t >= 80).sum())
    print(f"genes placed at >= 80% accuracy: vote {at80}, best single layer {at80_t}")
    print(byk.to_string()); print(bys.to_string())
    print(comp.round(1).to_string())

    panels = [("a", save(figa, "a_tradeoff")), ("b", save(figb, "b_overlap")),
              ("c", save(figc, "c_compartments")), ("d", save(figd, "d_layers")),
              ("e", save(fige, "e_share")), ("f", save(figf, "f_placed"))]
    import assemble_paper_figures as A_
    out = A_.compose("07", [{"png": p, "letter": l} for l, p in panels], OUT, rows=[3, 3])
    print("figure 7 ->", out)


if __name__ == "__main__":
    main()
