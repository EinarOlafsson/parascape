"""Map of the E8 "all views" clustering with the host-facing genes on it (Einar, 2026-11-25).

Same embedding and Louvain clustering as `show_host_clusters.py` (all views, cosine kNN 15, seed 0, resolution
1.0; the most host-facing cluster re-clustered at resolution 4.0). UMAP (cosine, 15 neighbours, min_dist 0.3,
seed 0) of every clustered gene. Panels, all the same size:

A. every gene coloured by its cluster (13 clusters); the most host-facing cluster in the accent colour;
B. effect on host (measured, 65 genes) on the map;
C. host-facing (papers): host-facing (72) and not host-facing (679) of the genes papers studied;
D. a UMAP of the most host-facing cluster only, its four most host-facing sub-clusters coloured,
   measured host-facing genes outlined in black.

    python scripts/infer_bench/show_host_map.py OUT_DIR
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e5_clustering as E5  # noqa: E402
import e8_selected_clustering as E8  # noqa: E402
from parascape import figstyle as F  # noqa: E402

RES, SUB_RES = 1.0, 4.0
MM = 1 / 25.4


def umap2(E, seed=0):
    import umap
    return umap.UMAP(n_neighbors=15, min_dist=0.3, metric="cosine", random_state=seed).fit_transform(E)


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if len(argv) != 1:
        raise SystemExit(__doc__)
    out = argv[0]
    os.makedirs(out, exist_ok=True)
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    eff = np.asarray(L["effect on host"][0], dtype=object) == "host-facing"
    pap = np.asarray(L["host-facing (papers)"][0], dtype=object)
    E, rows = E8.embed(nodes, views, extra, {v: c for v, c in views.items() if len(c) >= 2}, list(extra))
    c = E5.louvain(E, RES, 0)
    e_r, p_r = eff[rows], pap[rows]
    shares = pd.Series(e_r).groupby(c).mean().sort_values(ascending=False)
    top = int(shares.index[0])
    order = {int(k): i for i, k in enumerate(shares.index)}          # cluster number = rank by host-facing share
    U = umap2(E)
    inside = np.flatnonzero(c == top)
    Es = E[inside]
    s = E5.louvain(Es, SUB_RES, 0)
    sub_share = pd.Series(e_r[inside]).groupby(s).mean().sort_values(ascending=False)
    best4 = [int(k) for k in sub_share.index[:4]]
    Us = umap2(Es)
    pd.DataFrame({"gene_id": nodes["gene_id"].to_numpy()[rows], "product": nodes["product"].to_numpy()[rows],
                  "umap_1": U[:, 0], "umap_2": U[:, 1], "cluster": [order[int(k)] + 1 for k in c],
                  "effect_on_host": e_r, "host_facing_papers": np.where(pd.notna(p_r), p_r, "")}
                 ).to_csv(os.path.join(out, "map_all_genes.csv"), index=False)
    sub_rank = {k: i + 1 for i, k in enumerate(sub_share.index)}
    pd.DataFrame({"gene_id": nodes["gene_id"].to_numpy()[rows][inside],
                  "product": nodes["product"].to_numpy()[rows][inside], "umap_1": Us[:, 0], "umap_2": Us[:, 1],
                  "subcluster": [sub_rank[int(k)] for k in s], "effect_on_host": e_r[inside],
                  "host_facing_papers": np.where(pd.notna(p_r[inside]), p_r[inside], "")}
                 ).to_csv(os.path.join(out, "map_cluster1.csv"), index=False)
    # ---------------------------------------------------------------- figure: four equal square maps
    F.use()
    side, gap, top_m, left = 52.0, 10.0, 11.0, 6.0
    W = left + 4 * side + 3 * gap + 4
    Hh = top_m + side + 16
    fig = plt.figure(figsize=(W * MM, Hh * MM))
    axes = [fig.add_axes([(left + i * (side + gap)) / W, 1 - (top_m + side) / Hh, side / W, side / Hh])
            for i in range(4)]
    dot = 1.2
    greys = ["#B4B4B4", "#D2D2D2"]
    ax = axes[0]
    for k in np.unique(c):
        m = c == k
        col = F.ACCENT if k == top else greys[order[int(k)] % 2]
        ax.scatter(U[m, 0], U[m, 1], s=dot, c=col, alpha=F.ALPHA, lw=0, rasterized=True)
    ax = axes[1]
    ax.scatter(U[:, 0], U[:, 1], s=dot, c=F.BACKGROUND, alpha=F.ALPHA, lw=0, rasterized=True)
    ax.scatter(U[e_r, 0], U[e_r, 1], s=6, c=F.ACCENT, alpha=F.ALPHA, lw=0)
    ax = axes[2]
    ax.scatter(U[:, 0], U[:, 1], s=dot, c=F.BACKGROUND, alpha=F.ALPHA, lw=0, rasterized=True)
    no, yes = p_r == "not host-facing", p_r == "host-facing"
    ax.scatter(U[no, 0], U[no, 1], s=4, c=F.NAVY, alpha=F.ALPHA, lw=0)
    ax.scatter(U[yes, 0], U[yes, 1], s=6, c=F.ACCENT, alpha=F.ALPHA, lw=0)
    ax = axes[3]
    cols4 = (F.ACCENT, F.TEAL, F.MAGENTA, F.AMBER)
    ax.scatter(Us[:, 0], Us[:, 1], s=3, c=F.BACKGROUND, alpha=F.ALPHA, lw=0)
    for k, col in zip(best4, cols4):
        m = s == k
        ax.scatter(Us[m, 0], Us[m, 1], s=7, c=col, alpha=F.ALPHA, lw=0)
    e_in = e_r[inside]
    ax.scatter(Us[e_in, 0], Us[e_in, 1], s=7, facecolors="none", edgecolors=F.INK, lw=0.4)
    for i, ax in enumerate(axes):
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP 1", labelpad=2); ax.set_ylabel("UMAP 2", labelpad=2)
        ax.set_aspect("auto")
        fig.text((left + i * (side + gap) - 4) / W, 1 - (top_m - 5) / Hh, "ABCD"[i], fontsize=F.LETTER_PT,
                 weight="bold", va="top")
    from matplotlib.lines import Line2D
    mk = lambda col, lab, **kw: Line2D([], [], ls="", marker="o", ms=3, color=col, alpha=F.ALPHA, label=lab, **kw)
    hs = [mk(F.ACCENT, "cluster 1 (A) / host-facing (B, C) / sub-cluster 1 (D)"), mk(F.NAVY, "not host-facing, papers (C)"),
          mk(F.TEAL, "sub-cluster 2 (D)"), mk(F.MAGENTA, "sub-cluster 3 (D)"), mk(F.AMBER, "sub-cluster 4 (D)"),
          Line2D([], [], ls="", marker="o", ms=3, mfc="none", mec=F.INK, mew=0.4, label="effect on host (D)")]
    fig.legend(handles=hs, loc="lower center", ncol=3, frameon=False, fontsize=F.FONT_PT)
    fig.text(0.99, 0.995, "DEVELOPMENT FIGURE - made by a script (show_host_map.py)", ha="right", va="top",
             color=F.DARK_GREY, fontsize=F.NOTE_PT)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"fig_host_map.{ext}"), dpi=300)
    print("top cluster", top, "size", len(inside), "best sub-clusters", best4,
          [int((s == k).sum()) for k in best4], [int(e_in[s == k].sum()) for k in best4])
    return 0


if __name__ == "__main__":
    sys.exit(main())
