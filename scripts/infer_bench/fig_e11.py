"""Figure of E11: per-label optimised clustering against the references and the shuffled-label search.

DEVELOPMENT FIGURE made by a script from the E11 output folder (one sub-folder per label, `<name>` real search and
`<name>s` the same search on shuffled development labels). Panels:
A. untouched-gene class recovery F1 per label: E8 reference, E10 default, the per-label optimised choice, and
   the choice of the search run on shuffled labels (negative control);
B. stability of the chosen clustering (mean best Jaccard over 5 column redraws) per label;
Tables (not drawn): e11_summary.csv - what the optimiser chose per label; e11_positive_controls.csv - per-class
F1 on untouched genes for the ribosome and proteasome classes.

    python scripts/infer_bench/fig_e11.py E11_DIR OUT_DIR
"""
from __future__ import annotations

import glob
import json
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
import fig_e9 as FE9  # noqa: E402

MM = 1 / 25.4
CONTROL = ("40S", "60S", "ribosom", "proteasome", "PA700", "PA28")


def load(src):
    rows, per = [], []
    for d in sorted(glob.glob(os.path.join(src, "*"))):
        f = os.path.join(d, "final.json")
        if not os.path.exists(f):
            continue
        j = json.load(open(f))
        name = os.path.basename(d)
        rows.append({"run": name, "label": "planted modules (8 SD)" if j["planted"] else j["label"],
                     "shuffled": j["shuffled"], "test_f1": j["test_f1"], "cv_f1": j["cv_f1"],
                     "stability": j["stability"], "clusters": j["test_clusters"], "largest": j["test_largest_share"],
                     **{f"ref_{k}": v["test_f1"] for k, v in j["references"].items()},
                     **{f"cfg_{k}": v for k, v in j["chosen_config"].items()}})
        pc = os.path.join(d, "per_class.csv")
        if os.path.exists(pc) and not j["shuffled"]:
            p = pd.read_csv(pc)
            p["label"] = rows[-1]["label"]
            per.append(p)
    return pd.DataFrame(rows), (pd.concat(per) if per else pd.DataFrame())


def describe(r) -> str:
    m = r["cfg_method"]
    clus = {"louvain": f"Louvain {r['cfg_resolution']:.1f}" + (" cons." if r["cfg_consensus"] > 1 else ""),
            "kmeans": f"k-means {r['cfg_k']}", "ward": f"Ward {r['cfg_k']}",
            "hdbscan": f"HDBSCAN {r['cfg_min_size']}"}[m]
    sup = "" if r["cfg_supervised"] == "none" else f", {r['cfg_supervised'].upper()}"
    return f"{r['cfg_extract']}{sup}; {clus}"


def main(argv=None) -> int:
    src, out = (argv if argv is not None else sys.argv[1:])
    os.makedirs(out, exist_ok=True)
    F.use()
    T, P = load(src)
    real = T[~T["shuffled"]].set_index("label")
    shuf = T[T["shuffled"]].set_index("label")
    labels = list(real.index)
    short = [FE9.FE8_SHORT.get(x, x) for x in labels]
    refcols = [c for c in T.columns if c.startswith("ref_")]
    series = [(refcols[0], "E8 reference", F.CONTROL), (refcols[1], "E10 default", F.ACCENT_LIGHT),
              ("chosen", "optimised (this label)", F.ACCENT), ("shuffled", "optimised on shuffled labels", F.MAGENTA)]
    cell = F.SLOT_IN * 25.4
    n = len(labels)
    gw = cell * len(series)                                  # one slot per series inside each label group
    pw = gw * n + cell * (n - 1) * 0.5
    W, labw = max(pw + 30, 170), 22
    Hh = 8 + 40 + 40 + 26 + 8 * n + 30
    fig = plt.figure(figsize=(W * MM, Hh * MM))
    ax = fig.add_axes([labw / W, 1 - 48 / Hh, pw / W, 40 / Hh])
    x0s = np.arange(n) * (len(series) + 0.5)
    for j, (col, name, colr) in enumerate(series):
        if col == "chosen":
            v = real["test_f1"].to_numpy(float)
        elif col == "shuffled":
            v = np.array([shuf["test_f1"].get(l, np.nan) for l in labels], float)
        else:
            v = real[col].to_numpy(float)
        ax.bar(x0s + j, np.nan_to_num(v), width=F.BAR_FRACTION, color=colr, alpha=F.ALPHA, label=name,
               hatch="////" if col == "shuffled" else None, edgecolor=F.INK if col == "shuffled" else None, lw=0)
    ax.set_xticks(x0s + (len(series) - 1) / 2)
    ax.set_xticklabels([])
    ax.set_ylabel("Class recovery F1", labelpad=2)
    ax.set_xlim(-0.6, x0s[-1] + len(series) - 0.4)
    fig.text(2 / W, 1 - 6 / Hh, "A", fontsize=F.LETTER_PT, weight="bold", va="top")
    ax2 = fig.add_axes([labw / W, 1 - 92 / Hh, pw / W, 34 / Hh])
    ax2.bar(x0s + 2, real["stability"].to_numpy(float), width=F.BAR_FRACTION, color=F.ACCENT, alpha=F.ALPHA)
    ax2.bar(x0s + 3, [shuf["stability"].get(l, np.nan) for l in labels], width=F.BAR_FRACTION, color=F.MAGENTA,
            alpha=F.ALPHA, hatch="////", edgecolor=F.INK, lw=0)
    ax2.axhline(0.75, color=F.INK, lw=0.6, ls=":", alpha=0.6)
    ax2.set_ylim(0, 1.05)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(x0s + (len(series) - 1) / 2)
    ax2.set_xticklabels(short, rotation=90)
    ax2.set_ylabel("Stability (Jaccard)", labelpad=2)
    fig.text(2 / W, 1 - 50 / Hh, "B", fontsize=F.LETTER_PT, weight="bold", va="top")
    fig.legend(*ax.get_legend_handles_labels(), loc="upper center", ncol=4, frameon=False, fontsize=F.FONT_PT,
               bbox_to_anchor=(0.5, 0.0))
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"fig_e11.{ext}"), dpi=300, bbox_inches="tight")
    tab = real.reset_index()[["label", "test_f1", "stability", "clusters", "largest"] + refcols].copy()
    tab["shuffled_search_f1"] = [shuf["test_f1"].get(l, np.nan) for l in tab["label"]]
    tab["chosen"] = [describe(r) for _, r in real.reset_index().iterrows()]
    tab.round(3).to_csv(os.path.join(out, "e11_summary.csv"), index=False)
    if len(P):
        ctrl = P[P["class"].astype(str).str.contains("|".join(CONTROL))]
        ctrl.round(3).to_csv(os.path.join(out, "e11_positive_controls.csv"), index=False)
    print(tab.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
