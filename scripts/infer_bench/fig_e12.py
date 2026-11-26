"""Figure of E12: pair recovery (AUROC against degree-matched decoys) per pair set and embedding.

    python scripts/infer_bench/fig_e12.py E12_DIR OUT_DIR
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from parascape import figstyle as F  # noqa: E402

MM = 1 / 25.4
SHORT = {"crosslink (XL-MS, measured)": "crosslink", "STRING, measured in Toxoplasma": "STRING,\nmeasured",
         "STRING, copied by orthology": "STRING,\northology", "CORUM, copied by orthology": "CORUM,\northology"}


def main(argv=None) -> int:
    src, out = (argv if argv is not None else sys.argv[1:])
    os.makedirs(out, exist_ok=True)
    F.use()
    t = pd.read_csv(os.path.join(src, "summary.csv"))
    sets = list(dict.fromkeys(t["pairs"]))
    emb = [("E8 reference", F.ACCENT, None), ("E10 default", F.TEAL, None),
           ("E10 default, genes shuffled (control)", F.CONTROL, "////")]
    cell = F.SLOT_IN * 25.4
    pw = cell * (len(sets) * (len(emb) + 0.5))
    W, H_ = pw + 30, 74
    fig = plt.figure(figsize=(W * MM, H_ * MM))
    ax = fig.add_axes([20 / W, 26 / H_, pw / W, 38 / H_])
    x0 = np.arange(len(sets)) * (len(emb) + 0.5)
    for j, (e, c, h) in enumerate(emb):
        v = [float(t[(t.pairs == s) & (t.embedding == e)]["auroc"].iloc[0]) for s in sets]
        ax.bar(x0 + j, v, width=F.BAR_FRACTION, color=c, alpha=F.ALPHA, hatch=h,
               edgecolor=F.INK if h else None, lw=0, label=e.replace(", genes shuffled (control)", ", shuffled genes"))
    ax.axhline(0.5, color=F.INK, lw=0.6, ls=":", alpha=0.6)
    ax.set_ylim(0.4, 0.9)
    ax.set_xticks(x0 + 1)
    ax.set_xticklabels([SHORT[s] for s in sets])
    ax.set_ylabel("AUROC (pairs vs decoys)", labelpad=2)
    fig.legend(*ax.get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False, fontsize=F.FONT_PT)
    fig.text(0.995, 0.995, "DEVELOPMENT FIGURE - made by a script (fig_e12.py)", ha="right", va="top",
             color=F.DARK_GREY, fontsize=F.NOTE_PT)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"fig_e12.{ext}"), dpi=300)
    open(os.path.join(out, "fig_e12_legend.md"), "w").write(
        "**Figure E12. Interacting genes have similar profiles.** PRELIMINARY (development figure). AUROC of the "
        "cosine similarity of the two genes' profiles for known pairs against degree-matched decoys (5 per pair: "
        "the same first gene with a partner of the same number of known partners and no recorded interaction). "
        "Pair sets: StarPath XL-MS crosslinks (2,810), STRING v12.0 pairs measured in Toxoplasma (2,555), STRING "
        "pairs copied by orthology (score ≥ 700; 28,623), CORUM human complexes copied by orthology (17,917). "
        "Profiles: all views except the interaction counts (which track degree), as the E8 reference (per view "
        "PCA 30) or the E10 default (no per view PCA, whitening); hatched grey: E10 default with gene rows "
        "shuffled (negative control). The y axis starts at 0.4; the dotted line is chance (0.5).\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
