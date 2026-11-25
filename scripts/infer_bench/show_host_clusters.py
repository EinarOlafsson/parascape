"""Which all-views clusters hold the host-facing genes? (the E8 reference clustering, opened up)

Rebuilds the E8 "all views" embedding (every view, nothing left out: no feature is built from the host labels),
clusters it with Louvain (cosine kNN 15, seed 0) at resolution 0.5 and 1.0 exactly as E8 does, and for every
cluster counts:
  - effect on host: measured host-facing genes (65 in the genome; every other gene is 'other');
  - host-facing (papers): host-facing among the genes some paper studied (72 of 751).
Enrichment is a one-sided hypergeometric test against the genome (effect on host) or against all studied genes
(papers), Benjamini-Hochberg over clusters. The untouched 20 % of gene families (E7 final split) are counted
separately, as a held-out check.

Writes clusters.csv (one row per resolution x cluster), genes.csv (every gene in an enriched cluster, with its
labels), fig_host_clusters.png/.svg and a legend .md.

    python scripts/infer_bench/show_host_clusters.py OUT_DIR
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e5_clustering as E5  # noqa: E402
import e7_inference as E7  # noqa: E402
import e8_selected_clustering as E8  # noqa: E402
from parascape import figstyle as F  # noqa: E402

IM = H.IM
RES = (0.5, 1.0)
SUB_RES = 4.0
MM = 1 / 25.4


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if len(argv) != 1:
        raise SystemExit(__doc__)
    out = argv[0]
    os.makedirs(out, exist_ok=True)
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    eff = np.asarray(L["effect on host"][0], dtype=object)
    pap = np.asarray(L["host-facing (papers)"][0], dtype=object)
    is_eff = eff == "host-facing"
    studied = pd.notna(pap) & (pap.astype(str) != "None")
    is_pap = pap == "host-facing"
    fin = E7.final_split(H.families(nodes))
    allf = {v: c for v, c in views.items() if len(c) >= 2}
    E, rows = E8.embed(nodes, views, extra, allf, list(extra))
    n = len(nodes)
    crow, grow, labs = [], [], {}
    for res in RES:
        c = np.full(n, -1)
        c[rows] = E5.louvain(E, res, 0)
        on = c >= 0
        N, K = int(on.sum()), int((is_eff & on).sum())
        Ns, Ks = int((studied & on).sum()), int((is_pap & on).sum())
        for k in np.unique(c[on]):
            m = c == k
            x, s = int((is_eff & m).sum()), int(m.sum())
            xs, ss = int((is_pap & m).sum()), int((studied & m).sum())
            crow.append({"resolution": res, "cluster": int(k), "genes": s,
                         "effect_on_host": x, "effect_share": x / s, "effect_genome_share": K / N,
                         "effect_fold": (x / s) / (K / N), "effect_p": hypergeom.sf(x - 1, N, K, s),
                         "effect_untouched": int((is_eff & m & fin).sum()),
                         "papers_studied": ss, "papers_host_facing": xs,
                         "papers_share": xs / ss if ss else np.nan, "papers_studied_share": Ks / Ns,
                         "papers_p": hypergeom.sf(xs - 1, Ns, Ks, ss) if ss else np.nan,
                         "papers_untouched": int((is_pap & m & fin).sum())})
        labs[res] = c
    C = pd.DataFrame(crow)
    for res in RES:
        m = C["resolution"] == res
        for col in ("effect", "papers"):
            C.loc[m, f"{col}_q"] = E1.bh(C.loc[m, f"{col}_p"].fillna(1).to_numpy())
        for _, r in C[m & ((C["effect_q"] < 0.05) | (C["papers_q"] < 0.05))].iterrows():
            for i in np.flatnonzero(labs[res] == r["cluster"]):
                grow.append({"resolution": res, "cluster": int(r["cluster"]), "gene_id": nodes["gene_id"].iat[i],
                             "product": nodes["product"].iat[i], "effect_on_host": bool(is_eff[i]),
                             "papers": pap[i] if studied[i] else "", "untouched": bool(fin[i])})
    C.to_csv(os.path.join(out, "clusters.csv"), index=False)
    G = pd.DataFrame(grow)
    G.to_csv(os.path.join(out, "genes.csv"), index=False)
    # ---------------------------------------------------------------- inside the most enriched cluster (res 1.0)
    top = C[C["resolution"] == 1.0].sort_values("effect_share", ascending=False)["cluster"].iat[0]
    members = np.flatnonzero(labs[1.0] == top)
    Esub = E[np.searchsorted(rows, members)]
    sub = E5.louvain(Esub, SUB_RES, 0)
    Kp, Ksp = int(is_eff[members].sum()), int((is_pap[members]).sum())
    Nsp = int(studied[members].sum())
    srow, brow = [], []
    for j in np.unique(sub):
        m = members[sub == j]
        x, xs = int(is_eff[m].sum()), int(is_pap[m].sum())
        ss = int(studied[m].sum())
        prod = nodes["product"].iloc[m].fillna("").str.lower()
        srow.append({"subcluster": int(j), "genes": len(m), "effect_on_host": x, "effect_share": x / len(m),
                     "effect_parent_share": Kp / len(members),
                     "effect_p": hypergeom.sf(x - 1, len(members), Kp, len(m)),
                     "papers_studied": ss, "papers_host_facing": xs, "papers_share": xs / ss if ss else np.nan,
                     "papers_parent_share": Ksp / Nsp,
                     "papers_p": hypergeom.sf(xs - 1, Nsp, Ksp, ss) if ss else np.nan,
                     "hypothetical": int(prod.str.contains("hypothetical").sum()),
                     "untouched_effect": int((is_eff[m] & fin[m]).sum())})
        for i in m:
            brow.append({"subcluster": int(j), "gene_id": nodes["gene_id"].iat[i], "product": nodes["product"].iat[i],
                         "effect_on_host": bool(is_eff[i]), "papers": pap[i] if studied[i] else "",
                         "untouched": bool(fin[i])})
    S = pd.DataFrame(srow)
    for col in ("effect", "papers"):
        S[f"{col}_q"] = E1.bh(S[f"{col}_p"].fillna(1).to_numpy())
    S.to_csv(os.path.join(out, "subclusters.csv"), index=False)
    pd.DataFrame(brow).to_csv(os.path.join(out, "subcluster_genes.csv"), index=False)
    # ---------------------------------------------------------------- figure
    F.use()
    cell = F.SLOT_IN * 25.4
    ph, labw, gapx, rowgap = 30.0, 16.0, 20.0, 20.0
    C1 = C[C["resolution"] == 1.0].sort_values("effect_share", ascending=False)
    S1 = S.sort_values("effect_share", ascending=False)
    wA, wC = cell * len(C1), cell * len(S1)
    W = max(labw + 2 * wA + gapx, labw + wC) + 4
    Hh = 6 + 3 * ph + 3 * rowgap + 6
    fig = plt.figure(figsize=(W * MM, Hh * MM))
    ymax = 100 * float(np.nanmax(np.r_[C1[["effect_share", "papers_share"]].to_numpy().ravel(),
                                       S1[["effect_share", "papers_share"]].to_numpy().ravel()])) * 1.08
    specs = [("A", C1, "effect", labw, 0, wA, "effect_genome_share", "Cluster"),
             ("B", C1, "papers", labw + wA + gapx, 0, wA, "papers_studied_share", "Cluster"),
             ("C", S1, "effect", labw, 1, wC, "effect_parent_share", "Sub-cluster of cluster 1"),
             ("D", S1, "papers", labw, 2, wC, "papers_parent_share", "Sub-cluster of cluster 1")]
    for let, T, col, x0, r_, w, refcol, xl in specs:
        ytop = 6 + r_ * (ph + rowgap)
        ax = fig.add_axes([x0 / W, 1 - (ytop + ph) / Hh, w / W, ph / Hh])
        share = 100 * T[f"{col}_share"].to_numpy()
        sig = T[f"{col}_q"].to_numpy() < 0.05
        xs = np.arange(len(T))
        ax.bar(xs, np.nan_to_num(share), width=F.BAR_FRACTION,
               color=[F.ACCENT if q else F.CONTROL for q in sig], alpha=F.ALPHA)
        ax.axhline(100 * T[refcol].iat[0], color=F.INK, lw=0.6, ls=":", alpha=0.6)
        ax.set_xlim(-0.6, len(T) - 0.4)
        ax.set_ylim(0, ymax)
        ax.set_xticks(xs)
        ax.set_xticklabels([str(k + 1) for k in range(len(T))])
        ax.set_xlabel(xl, labelpad=2)
        ax.set_ylabel("Effect on host (%)" if col == "effect" else "Host-facing, papers (%)", labelpad=2)
        fig.text((x0 - 14) / W, 1 - (ytop - 2) / Hh, let, fontsize=F.LETTER_PT, weight="bold", va="top")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=F.ACCENT, alpha=F.ALPHA, label="enriched (q < 0.05)"),
                        Patch(color=F.CONTROL, alpha=F.ALPHA, label="not enriched"),
                        Line2D([], [], color=F.INK, lw=0.6, ls=":", label="reference rate")],
               loc="lower center", ncol=3, frameon=False, fontsize=F.FONT_PT)
    fig.text(0.99, 0.995, "DEVELOPMENT FIGURE - made by a script (show_host_clusters.py)", ha="right", va="top",
             color=F.DARK_GREY, fontsize=F.NOTE_PT)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"fig_host_clusters.{ext}"), dpi=300)
    print(S.sort_values("effect_share", ascending=False).round(4).head(8).to_string())
    print(C.round(4).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
