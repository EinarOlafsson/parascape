"""Figure of E10: the spaCR OPS clustering ledger repeated on ParaScape gene data.

DEVELOPMENT FIGURE made by a script from the E10 output folders (one per parallel group, each with its own
baseline). Columns: k = 100 and k = 300 clusters (matched across arms). Rows:
A. held-out CORUM complexes kept together (>= 50 % of members in one cluster), excess over a size-matched null;
B. held-out STRING pairs called (same cluster), recall excess over the null;
C. planted modules recovered (positive control, 20 x 15 genes at the calibrated strength).
Bars = mean of 3 seeds, points = seeds. Bar colour (A, B): blue = better than the baseline (paired bootstrap 95 %
interval above 0 in the mean over seeds), magenta = worse, light blue = no detectable difference; baseline and
shuffled-gene control grey (control hatched); dotted line = baseline mean.

    python scripts/infer_bench/fig_e10.py E10_DIR OUT_DIR
"""
from __future__ import annotations

import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parascape import figstyle as F  # noqa: E402
import e10_ops_ledger as E10  # noqa: E402

SHORT = {"baseline": "baseline", "robust z": "robust z", "rank-normal": "rank-normal", "KNN impute": "KNN fill",
         "no per-view PCA": "no view PCA", "equal view weight": "equal views", "whitened (TVN analogue)": "whitened",
         "measured views only": "measured only", "fused graph (SNF-like)": "fused graph", "k-means": "k-means",
         "Ward": "Ward", "average linkage, cosine": "average link", "spectral, cosine kNN 15": "spectral",
         "consensus, seeds": "cons. seeds", "consensus, feature bootstrap": "cons. features",
         "whitened + consensus, seeds": "white.+cons.", "no per-view PCA + whitened + consensus, seeds": "noPCA+wh.+cons.",
         "shuffled genes (negative control)": "shuffled"}
MM = 1 / 25.4


def load(src):
    L, B = [], []
    for d in sorted(glob.glob(os.path.join(src, "*"))):
        if os.path.exists(os.path.join(d, "ledger.parquet")):
            L.append(pd.read_parquet(os.path.join(d, "ledger.parquet")).assign(group=os.path.basename(d)))
            B.append(pd.read_parquet(os.path.join(d, "bootstrap.parquet")).assign(group=os.path.basename(d)))
    L = pd.concat(L)
    B = pd.concat(B) if B else pd.DataFrame()
    # the baseline is in every group (deterministic); keep one copy
    base = L[L["arm"] == "baseline"]
    L = pd.concat([base[base["group"] == base["group"].iloc[0]], L[L["arm"] != "baseline"]])
    return L, B


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if len(argv) != 2:
        raise SystemExit(__doc__)
    src, out = argv
    os.makedirs(out, exist_ok=True)
    F.use()
    L, B = load(src)
    order = [a for a, _ in E10.ARMS if a in set(L["arm"])]
    T = L[L["split"] == "test"].copy()
    ks = sorted(T["k"].unique())
    verdict = {}
    if len(B):
        m = B.groupby(["arm", "k"])[["complexes_diff", "complexes_lo", "complexes_hi", "string_recall_diff",
                                     "string_lo", "string_hi"]].mean()
        for (a, k), r in m.iterrows():
            verdict[(a, k, "complexes")] = ("better" if r["complexes_lo"] > 0 else
                                            ("worse" if r["complexes_hi"] < 0 else "same"))
    for met, colname in (("string", "string_f1_excess"), ("planted", "planted")):
        base_ = T[T["arm"] == "baseline"].set_index(["k", "seed"])[colname]
        for (a, k), g in T.groupby(["arm", "k"]):
            d = np.array([v - base_.get((k, sd), np.nan) for sd, v in zip(g["seed"], g[colname])])
            d = d[np.isfinite(d)]
            if len(d):
                verdict[(a, k, met)] = "better" if (d > 0).all() else ("worse" if (d < 0).all() else "same")
    sizes = pd.read_csv(os.path.join(src, "sizes.csv")) if os.path.exists(os.path.join(src, "sizes.csv")) else None
    if sizes is not None:
        for r in sizes.itertuples():
            T.loc[(T["arm"] == r.arm) & (T["k"] == r.k), "largest_share"] = r.largest_share
    rows = (("complexes_excess", "Complexes kept (excess)", "complexes"),
            ("string_f1_excess", "STRING F1 (excess)", "string"),
            ("planted", "Planted modules (%)", "planted"),
            ("largest_share", "Largest cluster (%)", None))
    rows = tuple(r for r in rows if r[0] in T)
    cell = F.SLOT_IN * 25.4
    pw, ph, labw, gapx, gapy = cell * len(order), 30.0, 22.0, 14.0, 6.0
    W = labw + len(ks) * (pw + gapx) + 2
    Hh = 8 + len(rows) * (ph + gapy) + 36
    fig = plt.figure(figsize=(W * MM, Hh * MM))
    colour = {"better": F.ACCENT, "worse": F.MAGENTA, "same": F.ACCENT_LIGHT}
    for ri, (col, ylab, met) in enumerate(rows):
        pct = col in ("planted", "largest_share")
        vals_all = 100 * T[col] if pct else T[col]
        lo, hi = float(np.nanmin(vals_all)), float(np.nanmax(vals_all))
        pad = 0.08 * (hi - lo if hi > lo else 1)
        for ci, k in enumerate(ks):
            x0, yt = labw + ci * (pw + gapx), 8 + ri * (ph + gapy)
            ax = fig.add_axes([x0 / W, 1 - (yt + ph) / Hh, pw / W, ph / Hh])
            for i, a in enumerate(order):
                v = T[(T["arm"] == a) & (T["k"] == k)][col].to_numpy(float)
                v = 100 * v if pct else v
                if not np.isfinite(v).any():
                    ax.text(i, lo + pad, "n/a", ha="center", va="bottom", color=F.CONTROL, fontsize=F.FONT_PT,
                            rotation=90)
                    continue
                ctrl = a == "baseline" or a.startswith("shuffled")
                c = F.CONTROL if ctrl else (colour[verdict.get((a, k, met), "same")] if met else F.NAVY)
                ax.bar(i, np.nanmean(v), width=F.BAR_FRACTION, color=c, alpha=F.ALPHA,
                       hatch="////" if a.startswith("shuffled") else None, edgecolor=F.INK if a.startswith("shuffled") else None,
                       lw=0)
                if col != "largest_share":          # sizes were measured on seed 0 only: no seed points
                    ax.scatter(np.full(len(v), i) + np.linspace(-0.12, 0.12, len(v)), v, s=4, c=F.INK,
                               alpha=F.ALPHA, lw=0, zorder=3)
            b = T[(T["arm"] == "baseline") & (T["k"] == k)][col].mean()
            ax.axhline(100 * b if pct else b, color=F.INK, lw=0.6, ls=":", alpha=0.6)
            ax.axhline(0, color=F.INK, lw=0.6)
            ax.set_ylim(min(0, lo) - pad, hi + pad)
            ax.set_xlim(-0.6, len(order) - 0.4)
            ax.set_xticks(range(len(order)))
            ax.set_xticklabels([SHORT.get(a, a) for a in order] if ri == len(rows) - 1 else [], rotation=90)
            if ci == 0:
                ax.set_ylabel(ylab, labelpad=2)
            if ri == 0:
                ax.set_title(f"k = {k} clusters", fontsize=F.FONT_PT, pad=3)
            if ci == 0:
                fig.text((x0 - 20) / W, 1 - (yt - 3) / Hh, "ABCD"[ri], fontsize=F.LETTER_PT, weight="bold", va="top")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    hs = [Patch(color=F.ACCENT, alpha=F.ALPHA, label="better than baseline"),
          Patch(color=F.MAGENTA, alpha=F.ALPHA, label="worse than baseline"),
          Patch(color=F.ACCENT_LIGHT, alpha=F.ALPHA, label="no detectable difference"),
          Patch(color=F.CONTROL, alpha=F.ALPHA, label="baseline (ref.)"),
          Patch(facecolor=F.CONTROL, alpha=F.ALPHA, hatch="////", edgecolor=F.INK, lw=0, label="shuffled genes"),
          Patch(color=F.NAVY, alpha=F.ALPHA, label="cluster size (D)"),
          Line2D([], [], color=F.INK, lw=0.6, ls=":", label="baseline mean")]
    fig.legend(handles=hs, loc="lower center", ncol=3, frameon=False, fontsize=F.FONT_PT)
    fig.text(0.995, 0.998, "DEVELOPMENT FIGURE - made by a script (fig_e10.py)", ha="right", va="top",
             color=F.DARK_GREY, fontsize=F.NOTE_PT)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"fig_e10.{ext}"), dpi=300)
    agg = (T.groupby(["arm", "k"])[["complexes_excess", "string_recall_excess", "string_f1_excess", "planted"]]
           .agg(["mean", "std"]).round(4))
    agg.to_csv(os.path.join(out, "e10_table.csv"))
    if len(B):
        B.groupby(["arm", "k"]).mean(numeric_only=True).round(4).to_csv(os.path.join(out, "e10_bootstrap_mean.csv"))
    print(os.path.join(out, "fig_e10.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
