"""Figure 6 (Einar, 2026-11-20): the predictive value of each information class on every label and
on every other class, tested on held-out genes.

  A  information class x label: gain over a shuffled-label null (held-out kNN, 2 seeds x 2 folds);
     cells not above the replicate noise floor are blank
  B  controls: a planted label per class (must be recovered by its own class) and a random label
     (must be at chance everywhere)
  C  information class x information class: shared nearest neighbours (how much two classes agree
     about which genes are alike), e.g. secondary structure vs 3Di
  D  label x label association over shared genes (does a sparse label track a well-covered one?)
  E  binding partners: how well each class separates crosslinked/pulled-down pairs from random pairs
  F  combining classes: held-out score as the best classes are added, per label
"""
from __future__ import annotations

import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402
import build_information_matrix as IM                   # noqa: E402
import paper_panels as PP                               # noqa: E402

F = "figure_6__final"
VN = {"transcript": "mRNA", "translation": "translation", "chromatin": "chromatin",
      "protein_abundance": "protein", "fitness": "CRISPR screens", "place": "localisation",
      "interaction": "interactions", "sequence": "sequence", "structure_computed": "structure",
      "conservation": "conservation", "esm": "ESM-2", "secondary_structure": "sec. structure",
      "3di": "3Di", "domain": "domains", "host_response": "host response",
      "phenotype": "phenotype"}
LN = {"effect on host": "effect on host", "enzymatic function (EC class)": "EC class",
      "hyperLOPIT niche": "hyperLOPIT niche", "phenotype score, fibroblasts": "fibroblast PS",
      "phenotype score, IFN-gamma macrophages": "IFN-gamma PS",
      "phenotype score, in vivo (brain)": "in vivo PS", "cell-cycle phase": "cell-cycle phase",
      "nutrient dependence (glucose)": "glucose dependence", "imaging phenotype": "imaging phenotype"}


def heat(M, rows, cols, cbl, cmap="Blues", vmin=0, vmax=None, cell=3.6, mask=None):
    FS.use()
    fig, ax = FS.plot_figure(len(cols) * cell, len(rows) * cell, extra_bottom_mm=10)
    A = np.ma.masked_invalid(M if mask is None else np.where(mask, M, np.nan))
    cm = matplotlib.colormaps[cmap].copy(); cm.set_bad("#EEEEEE")
    im = ax.imshow(A, cmap=cm, vmin=vmin, vmax=vmax, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows); ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(cbl); cb.outline.set_linewidth(0.6)
    return fig, ax


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    d = sorted(glob.glob(os.path.join(ROOT, "results", "information_matrix_*")))[-1]
    cells = pd.read_csv(os.path.join(d, "cells.tsv"), sep="\t")
    real = cells[~cells.label.str.startswith(("PLANTED", "RANDOM"))]
    views = [v for v in VN if v in set(cells.view)]
    labs = [l for l in LN if l in set(real.label)]
    G = real.pivot_table(index="label", columns="view", values="gain").reindex(index=labs,
                                                                              columns=views)
    C = real.pivot_table(index="label", columns="view", values="called").reindex(index=labs,
                                                                                columns=views)
    figa, ax = heat(G.to_numpy(float), [LN[l] for l in labs], [VN[v] for v in views],
                    "Gain over shuffled labels", vmin=0, vmax=0.6,
                    mask=C.fillna(False).to_numpy(bool))
    PP.panel(F, "a", fig=figa, data=real, stats=real[["label", "view", "genes", "score",
             "null_mean", "null_sd", "replicate_diff", "gain", "z", "called"]],
             sources=[os.path.join(d, "cells.tsv")],
             what="each information class predicting each label on held-out genes")

    ctl = cells[cells.label.str.startswith(("PLANTED", "RANDOM"))].copy()
    ctl["row"] = ctl.label.str.replace("PLANTED: ", "", regex=False)
    rows_ = [v for v in views if v in set(ctl.row)] + ["RANDOM (must be chance)"]
    Pm = ctl.pivot_table(index="row", columns="view", values="gain").reindex(index=rows_,
                                                                             columns=views)
    figb, ax = heat(Pm.to_numpy(float), [VN.get(r, "random label") for r in rows_],
                    [VN[v] for v in views], "Gain over shuffled labels", vmin=0, vmax=0.5)
    ax.set_ylabel("Planted label")
    PP.panel(F, "b", fig=figb, data=ctl, stats={"diagonal mean gain": float(np.nanmean(
             [Pm.loc[v, v] for v in views if v in Pm.index])), "random label max |gain|":
             float(Pm.loc["RANDOM (must be chance)"].abs().max())},
             sources=[os.path.join(d, "cells.tsv")], what="planted positives and a random label")

    # C: shared nearest neighbours between classes
    nodes = pd.read_parquet(os.path.join(IM.CACHE, "nodes.parquet"))
    vw = IM.build_views(nodes); ex = IM.extra_views(nodes)
    emb = {}
    for name, cols in vw.items():
        X = nodes[cols].to_numpy(float); cov = (~np.isnan(X)).mean(axis=1) >= 0.5
        if cov.sum() >= 100:
            emb[name] = IM.embed(X, cov)
    for name, X in ex.items():
        cov = ~np.isnan(X).any(axis=1)
        if cov.sum() >= 100:
            emb[name] = IM.embed(X, cov)
    vv = [v for v in views if v in emb]
    k = 15
    nbr = {}
    for v in vv:
        ok = ~np.isnan(emb[v]).any(axis=1)
        nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine").fit(emb[v][ok])
        _, ind = nn.kneighbors(emb[v][ok])
        gi = np.flatnonzero(ok)
        nbr[v] = {g: set(gi[r[1:]]) for g, r in zip(gi, ind)}
    S = np.full((len(vv), len(vv)), np.nan)
    for i, a in enumerate(vv):
        for j, b in enumerate(vv):
            common = sorted(set(nbr[a]) & set(nbr[b]))
            if len(common) < 100:
                continue
            cs = set(common)
            vals = [len((nbr[a][g] & cs) & (nbr[b][g] & cs)) / k for g in common]
            exp = k / len(common)
            S[i, j] = np.mean(vals) / exp            # fold over random neighbours
    Sl = np.log2(S)
    figc, ax = heat(Sl, [VN[v] for v in vv], [VN[v] for v in vv],
                    "Shared neighbours (log$_2$ fold over random)", cmap="Blues", vmin=0,
                    vmax=float(np.nanpercentile(Sl[~np.eye(len(vv), dtype=bool)], 98)))
    PP.panel(F, "c", fig=figc, data=pd.DataFrame(S, index=vv, columns=vv), stats={
             "k": k, "measure": "mean fraction of 15 nearest neighbours shared, / random"},
             sources=[os.path.join(IM.CACHE, "nodes.parquet"), IM.STRUCT],
             what="how much two information classes agree on which genes are alike")

    # D: label x label association
    a = pd.read_csv(os.path.join(d, "label_association.tsv"), sep="\t")
    D = np.full((len(labs), len(labs)), np.nan); N = np.zeros_like(D)
    for r in a.itertuples():
        if r.label_a in labs and r.label_b in labs:
            i, j = labs.index(r.label_a), labs.index(r.label_b)
            D[i, j] = D[j, i] = abs(r.value); N[i, j] = N[j, i] = r.shared_genes
    figd, ax = heat(D, [LN[l] for l in labs], [LN[l] for l in labs],
                    "Association (Cramer V, $\\epsilon^2$ or |$\\rho$|)", vmin=0, vmax=0.6)
    PP.panel(F, "d", fig=figd, data=a, stats=a, sources=[os.path.join(d, "label_association.tsv")],
             what="label x label association over shared genes")

    # E: binding partners per class
    pr = pd.read_csv(os.path.join(d, "partners.tsv"), sep="\t").sort_values("auc", ascending=False)
    FS.use()
    fige, ax = FS.ref_bars(pr.auc.to_numpy(), [VN.get(v, v) for v in pr.view],
                           highlight=[VN.get(v, v) for v in pr.view[:3]],
                           ylabel="AUC, partners vs random pairs")
    ax.axhline(0.5, color=FS.INK, lw=0.6, ls=(0, (1.5, 1.5)), label="chance")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02)); ax.set_ylim(0.4, 1.0)
    PP.panel(F, "e", fig=fige, data=pr, stats=pr, sources=[os.path.join(d, "partners.tsv")],
             what="binding partners separated from random pairs, per class")

    # F: combining classes
    cb_ = pd.read_csv(os.path.join(d, "combinations.tsv"), sep="\t")
    FS.use()
    figf, ax = FS.square_figure()
    for q, l in enumerate(labs):
        s_ = cb_[cb_.label == l]
        if s_.empty:
            continue
        ax.plot(s_.step, s_.score, "-o", ms=3, lw=1.0, color=FS.SERIES[q % len(FS.SERIES)]
                if q < 6 else FS.CONTROL, label=LN[l])
    ax.set_xlabel("Information classes combined"); ax.set_ylabel("Held-out score (AUC or $\\rho$)")
    ax.set_xticks([1, 2, 3, 4]); ax.set_ylim(0, 1)
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2)
    PP.panel(F, "f", fig=figf, data=cb_, stats=cb_, sources=[os.path.join(d, "combinations.tsv")],
             what="held-out score as the best classes are combined")
    print(PP.compose(F, "main", "Which information predicts which: every information class "
                     "tested against every label and every other class on held-out genes."))


if __name__ == "__main__":
    main()
