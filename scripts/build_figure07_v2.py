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
    # A: per layer, held-out accuracy (bars); the weighted vote last, in the accent; guessing
    # the commonest compartment as a dotted line. Genes answered are in data.csv and the legend.
    order_a = sorted(adj, key=lambda k: -weights[k])
    pl = pl.set_index("layer").reindex([k for k in order_a if k in set(pl["layer"])])
    pl.loc["weighted vote"] = [len(d), d.right.mean()]
    SHORT_L = {"comention": "co-mention", "comention_ft": "co-mention FT", "ip_ms": "IP-MS",
               "xlms": "XL-MS", "unwritten_interaction": "pred. interact.",
               "orthogroup": "orthogroup", "struct": "structure", "domain": "domain",
               "structural_hole": "struct. hole", "coexpression": "co-expression",
               "cotranslation": "co-translation", "cofitness": "co-fitness",
               "weighted vote": "all combined"}
    FS.use()
    figa, ax = FS.ref_bars(pl["acc"].to_numpy() * 100, [SHORT_L.get(k, k) for k in pl.index],
                           highlight=["all combined"], ylabel="Correctly localised (%)")
    ax.axhline(baseline * 100, color=FS.INK, lw=0.6, ls=(0, (1.5, 1.5)),
               label="majority-class baseline")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02))
    ax.set_ylim(0, 100)

    # B: overlap of genes reached (Jaccard), layers ordered by weight
    order = sorted(adj, key=lambda k: -weights[k])
    reached = {k: set(np.where([len(x) > 0 for x in adj[k]])[0]) for k in order}
    J = np.array([[len(reached[a] & reached[b]) / max(1, len(reached[a] | reached[b]))
                   for b in order] for a in order]) * 100
    FS.use()
    nl = len(order); cell = 0.13
    w, h = 1.0 + nl * cell + 0.9, 0.15 + nl * cell + 1.0
    figb = plt.figure(figsize=(w, h))
    ax = figb.add_axes([1.0 / w, 1.0 / h, nl * cell / w, nl * cell / h])
    im = ax.imshow(J, cmap="Blues", vmin=0, vmax=100, interpolation="nearest")
    names = [SHORT_L.get(k, k) for k in order]
    ax.set_xticks(range(nl)); ax.set_xticklabels(names, rotation=90)
    ax.set_yticks(range(nl)); ax.set_yticklabels(names)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = figb.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
    cb.set_label("Genes in common (%)"); cb.outline.set_linewidth(0.6)

    # C: per compartment, weighted vote (bars) against the best single layer (tick)
    comp = (d.groupby("truth")[["best_right", "right"]].mean() * 100)
    comp["n"] = d.groupby("truth").size()
    comp = comp[comp["n"] >= 20].sort_values("right", ascending=False)
    CN = {"ribosome_rnp": "ribosome", "apical_secretory": "apical", "pellicle_pm_imc": "pellicle",
          "golgi_secretory_vesicular": "Golgi", "er": "ER", "mitochondrion": "mito."}
    FS.use()
    names_c = [CN.get(c, c) for c in comp.index]
    figc, ax = FS.ref_bars(comp["right"].to_numpy(), names_c, refs=comp["best_right"].to_numpy(),
                           highlight=names_c, ref_label="best single evidence type",
                           ylabel="Correctly localised (%)")
    ax.set_ylim(0, 100)

    # D: accuracy by number of layers reaching the gene
    d["k"] = d["n_layers"].clip(upper=6)
    byk = d.groupby("k")["right"].agg(["mean", "size"])
    FS.use()
    figd, ax = FS.square_figure()
    ax.axhline(baseline * 100, color=FS.INK, lw=0.6, ls="--", label="majority-class baseline")
    ax.plot(byk.index, byk["mean"] * 100, "o-", color=FS.ACCENT, ms=4, lw=1.2,
            label="all evidence combined")
    ax.set_xticks(list(byk.index))
    ax.set_xticklabels([str(k) if k < 6 else "6+" for k in byk.index])
    ax.set_xlabel("Evidence types reaching the gene"); ax.set_ylabel("Correctly localised (%)")
    ax.set_ylim(0, 100)
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2)

    # E: accuracy by winning share
    bins = np.linspace(0, 1, 11)
    d["sb"] = pd.cut(d["share"], bins, include_lowest=True)
    bys = d.groupby("sb", observed=True)["right"].agg(["mean", "size"])
    mids = [b.mid for b in bys.index]
    FS.use()
    fige, ax = FS.square_figure()
    ax.plot([0, 100], [0, 100], color=FS.BACKGROUND, lw=0.6, label="perfect calibration")
    ax.plot(np.array(mids) * 100, bys["mean"] * 100, "o-", color=FS.ACCENT, ms=4, lw=1.2,
            label="all evidence combined")
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.set_xlabel("Vote share of the predicted niche (%)"); ax.set_ylabel("Correctly localised (%)")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2)

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
            label="best single evidence type")
    ax.plot(np.arange(1, len(s) + 1), cum, color=FS.ACCENT, lw=1.2, label="all evidence combined")
    ax.axhline(baseline * 100, color=FS.INK, lw=0.6, ls="--", label="majority-class baseline")
    ax.set_xlabel("Predictions, most confident first"); ax.set_ylabel("Correctly localised (%)")
    ax.set_ylim(0, 100)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    at80 = int((cum >= 80).sum()); at80_t = int((cum_t >= 80).sum())
    print(f"genes placed at >= 80% accuracy: vote {at80}, best single layer {at80_t}")
    print(byk.to_string()); print(bys.to_string())
    print(comp.round(1).to_string())

    import paper_panels as PP
    from scipy.stats import spearmanr as _sp, binomtest
    GRAPH = os.path.join(CACHE, "graph.npz")
    NODES = (os.path.join(CACHE, "nodes.parquet"), [NL.TARGET])
    NLR = [os.path.join(ROOT, "results", "network_layers_2026_11_18", f)
           for f in ("layers.tsv", "summary.json")]
    pred = d.drop(columns=["layers"]).assign(layers=d["layers"].map(";".join))
    common = {"held-out genes answered": len(d), "baseline (commonest compartment)": baseline,
              "weights fitted on": "one half of labelled genes (seed 0)",
              "scored on": "the other half"}
    # ---- Figure 6: each evidence type, tested on held-out genes (storyboard 2026-11-20) ----
    F = "supplementary_figure_6b__final"     # the new Fig. 6 is the predictive matrix
    # A: one gene and its neighbours, the picture of guilt by association
    gid_index = {g: k for k, g in enumerate(nodes["gene_id"].astype(str))}
    cand = d[(d.n_layers >= 4) & d.right].copy()
    tot = {i_: sum(len(adj[k][i_]) for k in live) for i_ in cand.gene}
    cand["deg"] = cand.gene.map(tot)
    cand = cand[(cand.deg >= 8) & (cand.deg <= 22)].sort_values("n_layers", ascending=False)
    ego = int(cand.gene.iloc[0])
    lay_cols = {}
    pal = [FS.ACCENT, FS.TEAL, FS.NAVY, FS.MAGENTA, FS.AMBER, FS.ACCENT_LIGHT]
    nb = []
    for k in live:
        for j_ in adj[k][ego]:
            if j_ != ego:
                nb.append((k, int(j_)))
    layers_here = sorted({k for k, _ in nb}, key=lambda k: -weights[k])
    for q, k in enumerate(layers_here):
        lay_cols[k] = pal[q % len(pal)]
    uniq = sorted({j_ for _, j_ in nb})
    ang = {j_: 2 * np.pi * q / len(uniq) for q, j_ in enumerate(uniq)}
    comp_names = sorted({labels[j_] for j_ in uniq if known[j_]} | {labels[ego]})
    house = [FS.ACCENT, FS.TEAL, FS.MAGENTA, FS.AMBER, FS.NAVY, FS.ACCENT_LIGHT, FS.DARK_GREY]
    ccol = {c: house[q % len(house)] for q, c in enumerate(comp_names)}
    FS.use()
    fig0, ax = FS.square_figure()
    for k, j_ in nb:
        off = (layers_here.index(k) - len(layers_here) / 2) * 0.02
        ax.plot([0, np.cos(ang[j_]) + off], [0, np.sin(ang[j_]) + off], color=lay_cols[k],
                lw=0.9, alpha=0.8, zorder=1)
    for j_ in uniq:
        ax.scatter(np.cos(ang[j_]), np.sin(ang[j_]), s=40, zorder=3, lw=0.4, edgecolor=FS.INK,
                   color=ccol[labels[j_]] if known[j_] else "white")
    ax.scatter(0, 0, s=130, zorder=4, lw=1.0, edgecolor=FS.INK, color=ccol[labels[ego]],
               linestyle=(0, (1.5, 1.0)))
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    for k in layers_here:
        ax.plot([], [], color=lay_cols[k], lw=1.2, label=SHORT_L.get(k, k))
    for c in comp_names:
        ax.scatter([], [], s=25, color=ccol[c], lw=0.4, edgecolor=FS.INK,
                   label=COMP_NAMES.get(c, c))
    ax.scatter([], [], s=25, color="white", lw=0.4, edgecolor=FS.INK, label="unlabelled")
    ax.scatter([], [], s=40, color="white", lw=1.0, edgecolor=FS.INK, linestyle=(0, (1.5, 1.0)),
               label="held-out gene (true niche)")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=FS.NOTE_PT,
              handlelength=1.0, labelspacing=0.2)
    egoname = str(nodes["gene_id"].iloc[ego])
    PP.panel(F, "a", fig=fig0, data=pd.DataFrame(
        [{"gene": egoname, "neighbour": str(nodes["gene_id"].iloc[j_]), "evidence_type": k,
          "neighbour_niche": labels[j_] if known[j_] else ""} for k, j_ in nb]),
             stats={"gene": egoname, "true niche": labels[ego], "evidence types": len(layers_here),
                    "neighbours": len(uniq)}, sources=[GRAPH, NODES],
             what="one held-out gene, its neighbours and the evidence linking them")
    PP.panel(F, "b", fig=figa, data=pl.rename_axis("layer").reset_index(),
             stats={**common, **{f"weight, {k}": v for k, v in weights.items()}},
             sources=NLR + [GRAPH, NODES], what="held-out accuracy per evidence type")
    PP.panel(F, "c", fig=figb, data=pd.DataFrame(J, index=order, columns=order)
             .rename_axis("layer"), stats={"layers": len(order)}, sources=[GRAPH],
             what="genes reached by two evidence types")
    b = binomtest(int(d.right.sum()), len(d), d.best_right.mean())
    PP.panel(F, "d", fig=figc, data=comp.rename_axis("compartment"),
             stats={**common, "weighted vote correct": d.right.mean(),
                    "best single layer correct": d.best_right.mean(),
                    "binomial test, vote vs best-layer rate, p": b.pvalue},
             sources=[GRAPH, NODES], what="per niche: combined vote vs best single evidence type")
    out = PP.compose(F, "supplementary", "Each evidence type in the gene graph, tested on held-out genes.")

    # ---- Figure 9: inference with an uncertainty estimate ----
    F = "figure_9__final"
    rk = _sp(d["n_layers"], d["right"])
    PP.panel(F, "a", fig=figd, data=byk.rename(columns={"mean": "correct", "size": "genes"})
             .rename_axis("layers_reaching_gene"),
             stats={**common, "Spearman rho, layers vs correct": rk.statistic, "p": rk.pvalue},
             sources=[GRAPH, NODES], what="accuracy by the number of agreeing evidence types")
    PP.panel(F, "b", fig=fige, data=bys.rename(columns={"mean": "correct", "size": "genes"})
             .assign(bin=[str(b_) for b_ in bys.index]).reset_index(drop=True),
             stats=common, sources=[GRAPH, NODES], what="calibration: vote share vs accuracy")
    PP.panel(F, "c", fig=figf, data=pred, stats={**common,
             "genes placed at >= 80 % accuracy, vote (most layers first)": at80,
             "genes placed at >= 80 % accuracy, best single layer": at80_t},
             sources=[GRAPH, NODES], what="accuracy against the number of predictions made")
    # D: what inference adds, per target that beats the baseline: measured / predicted / unknown
    FMb = os.path.join(ROOT, "results", "filling_the_map_2026_11_18", "best_per_target.tsv")
    fmb = pd.read_csv(FMb, sep="\t")
    pf_ = pd.read_parquet(os.path.join(CACHE, "pf_nodes.parquet"))
    tgt = {("Tg", "lopit_unified"): "Tg localisation",
           ("Pf", "lopit_pf_location"): "Pf localisation", ("Tg", "cellcycle_phase"): "Tg phase"}
    rows = []
    for (org, t_), nm in tgt.items():
        r_ = fmb[(fmb.organism == org) & (fmb.target == t_)].iloc[0]
        tab = nodes if org == "Tg" else pf_
        col = tab[t_].astype(str).str.lower()
        miss = tab[t_].isna() | col.isin(["nan", "none", "", "unassigned", "unknown", "<na>"])
        n_all = len(tab); n_known = int((~miss).sum()); n_miss = int(miss.sum())
        n_pred = n_miss * r_.answer_rate
        rows.append({"target": nm, "measured": n_known / n_all * 100,
                     "predicted": n_pred / n_all * 100,
                     "unknown": (n_miss - n_pred) / n_all * 100,
                     "accuracy_when_predicted": r_.right_when_answered})
    st_ = pd.DataFrame(rows)
    FS.use()
    figg, ax = FS.plot_figure(len(st_) * FS.SLOT_IN * 25.4 * 1.2, FS.PLOT_H_MM)
    x_ = np.arange(len(st_)); bot = np.zeros(len(st_))
    for col_, c_ in (("measured", FS.DARK_GREY), ("predicted", FS.ACCENT),
                     ("unknown", "#E6E6E6")):
        ax.bar(x_, st_[col_], bottom=bot, width=0.6, color=c_, lw=0, label=col_)
        bot += st_[col_].to_numpy()
    ax.set_xticks(x_); ax.set_xticklabels(st_.target, rotation=45, ha="right",
                                          rotation_mode="anchor")
    ax.tick_params(axis="x", length=0); ax.set_ylim(0, 100); ax.set_ylabel("Genes (%)")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=3)
    PP.panel(F, "d", fig=figg, data=st_, stats=st_, sources=[FMb], what="measured, predicted and "
             "still unknown, per target")
    out = PP.compose(F, "main", "ParaScape infers missing labels with an uncertainty estimate: "
                     "confidence grows with agreeing evidence and is calibrated.")
    print("figure 7 ->", out)


if __name__ == "__main__":
    main()
