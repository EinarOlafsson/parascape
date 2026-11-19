"""Figure 8, one claim: secreted effectors vary between strains far more than housekeeping genes,
and the most variable of them are the most disordered -- the part AlphaFold cannot model.

Einar, 2026-11-19: point unclear, placeholders unhelpful. All panels real:
results/strain_variation_by_class_2026_11_19 (nonsynonymous / synonymous SNPs across strains,
each class against 2,000 coding-length-matched random sets), AlphaFold DB mean pLDDT per gene
(`mean_plddt`), and the AF3 screen bait foldability (208 host-facing genes with an AFDB entry).

  A  per class: nonsynonymous / synonymous ratio vs length-matched random genes  (the claim)
  B  per class: the per-gene ratios (median, IQR)                                (not a few genes)
  C  per class: AlphaFold DB mean pLDDT (median, IQR)                            (secreted = disordered)
  D  per gene: strain variation vs pLDDT, dense granules marked                  (variable = disordered)
  E  the 208 host-facing baits: residues with pLDDT >= 70                        (what AF3 can see)
  F  the most variable dense-granule proteins                                    (who they are)
"""
from __future__ import annotations

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402

OUT = os.path.join(ROOT, "results", "figure08_v2")
PNL = os.path.join(OUT, "panels")
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
SV = os.path.join(ROOT, "results", "strain_variation_by_class_2026_11_19")
BAITS = ("/media/carruthers/mnt3/claude/toxoplasma_projects/af3_folding/reports/"
         "bait_folded_fraction.json")
SECRETED = ("dense granules (hyperLOPIT)", "carries a signal peptide",
            "apical secretory (rhoptries, micronemes)")
SHORT = {"dense granules (hyperLOPIT)": "dense granules",
         "carries a signal peptide": "signal peptide",
         "apical secretory (rhoptries, micronemes)": "rhoptries, micronemes",
         "surface and pellicle": "surface, pellicle",
         "required in vitro (lowest fitness tenth)": "essential (lowest tenth)"}
MYR1 = "TGME49_254470"
CUT = 0.30


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def rows_fig(nr, left=1.45, value=1.6, right=1.6, rh=0.2):
    w, h = left + value + right, 0.15 + nr * rh + 0.6
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([left / w, 0.6 / h, value / w, nr * rh / h])
    ax.set_ylim(nr - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    return fig, ax


def iqr_rows(ax, groups, order, accent):
    for y, k in enumerate(order):
        v = np.asarray(groups[k], float); v = v[~np.isnan(v)]
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        col = FS.ACCENT if k in accent else FS.CONTROL
        ax.hlines(y, q1, q3, color=col, lw=1.2, alpha=FS.ALPHA)
        ax.plot(med, y, "o", color=col, ms=3.5)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([SHORT.get(k, k) for k in order])


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    bc = pd.read_csv(os.path.join(SV, "by_class.tsv"), sep="\t").sort_values("ratio",
                                                                            ascending=False)
    pg = pd.read_csv(os.path.join(SV, "per_gene.tsv"), sep="\t")
    n = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    n = n.set_index(n["gene_id"].astype(str))
    order = list(bc["class"])

    # A
    FS.use()
    figa, ax = rows_fig(len(order))
    y = np.arange(len(order))
    ax.hlines(y, bc["length_matched_baseline"], bc["ratio"], color=FS.CONTROL, lw=0.6)
    ax.plot(bc["length_matched_baseline"], y, "o", color=FS.CONTROL, ms=3.5, alpha=FS.ALPHA,
            label="length-matched random (ref.)")
    cols = [FS.ACCENT if c in SECRETED else FS.DARK_GREY for c in order]
    ax.scatter(bc["ratio"], y, s=14, c=cols, alpha=FS.ALPHA, lw=0, zorder=3)
    ax.scatter([], [], s=14, c=FS.ACCENT, alpha=FS.ALPHA, lw=0, label="secreted class")
    ax.scatter([], [], s=14, c=FS.DARK_GREY, alpha=FS.ALPHA, lw=0, label="other class")
    ax.set_yticks(y); ax.set_yticklabels([SHORT.get(k, k) for k in order])
    ax.set_xlabel("Nonsynonymous / synonymous SNPs")
    FS.pad_limits(ax, list(bc["ratio"]) + list(bc["length_matched_baseline"]), axis="x")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    # B
    groups = {c: pg.loc[pg["class"] == c, "ratio"].to_numpy() for c in order}
    FS.use()
    figb, ax = rows_fig(len(order), right=0.3)
    iqr_rows(ax, groups, order, SECRETED)
    allmed = float(np.nanmedian(pg.drop_duplicates("gene_id")["ratio"]))
    ax.axvline(allmed, color=FS.INK, lw=0.6, ls="--")
    ax.set_xlabel("Nonsynonymous / synonymous, per gene")

    # C
    pl = n["mean_plddt"]
    groups_c = {c: pl.reindex(pg.loc[pg["class"] == c, "gene_id"].astype(str)).to_numpy()
                for c in order}
    FS.use()
    figc, ax = rows_fig(len(order), right=0.3)
    iqr_rows(ax, groups_c, order, SECRETED)
    ax.axvline(70, color=FS.INK, lw=0.6, ls="--")
    ax.set_xlabel("AlphaFold DB mean pLDDT")
    for k in order:
        v = groups_c[k]; v = v[~np.isnan(v)]
        print(f"  pLDDT {SHORT.get(k, k):28s} median {np.median(v):.1f}  n={len(v)}")

    # D
    g = pg.drop_duplicates("gene_id").set_index("gene_id")
    dd = pd.DataFrame({"ratio": g["ratio"], "plddt": pl.reindex(g.index)}).dropna()
    dd = dd[dd["ratio"] > 0]
    dg = set(pg.loc[pg["class"] == "dense granules (hyperLOPIT)", "gene_id"])
    isdg = dd.index.isin(dg)
    rho = spearmanr(dd["ratio"], dd["plddt"]).statistic
    FS.use()
    figd, ax = FS.square_figure()
    ax.plot(dd.loc[~isdg, "plddt"], dd.loc[~isdg, "ratio"], "o", color=FS.CONTROL, ms=1.2,
            alpha=0.3, lw=0, rasterized=True, label="other genes (ref.)")
    ax.plot(dd.loc[isdg, "plddt"], dd.loc[isdg, "ratio"], "o", color=FS.ACCENT, ms=2.2,
            alpha=0.8, lw=0, label="dense granules")
    ax.plot([], [], " ", label=f"Spearman $\\rho$ = {rho:.2f}")
    ax.set_yscale("log")
    ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels(["0.1", "1", "10"]); ax.minorticks_off()
    ax.set_xlabel("AlphaFold DB mean pLDDT")
    ax.set_ylabel("Nonsynonymous / synonymous")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    print(f"  ratio vs pLDDT: Spearman {rho:.3f} over {len(dd)} genes")

    # E
    raw = json.load(open(BAITS))
    vals = np.array([v for v in raw.values() if isinstance(v, (int, float))], float)
    FS.use()
    fige, ax = FS.square_figure()
    ax.hist(vals, bins=np.linspace(0, 1, 21), color=FS.CONTROL, alpha=FS.ALPHA, lw=0,
            label="host-facing proteins")
    ax.axvline(CUT, color=FS.INK, lw=0.6, ls="--", label="co-folding admission cut")
    ax.axvline(float(raw[MYR1]), color=FS.ACCENT, lw=1.4, label="MYR1")
    ax.set_xlabel("Residues with pLDDT $\\geq$ 70 (fraction)"); ax.set_ylabel("Proteins")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    print(f"  baits {len(vals)}, below cut {(vals < CUT).mean() * 100:.1f} %")

    # F
    dgr = g.loc[g.index.isin(dg)].copy()
    dgr["product"] = n["product"].reindex(dgr.index).astype(str)
    top = dgr.sort_values("ratio", ascending=False).head(15)

    #: Published names where the product line is a phrase rather than a name.
    KNOWN = {"TGME49_240060": "IST", "TGME49_239010": "TEEGR"}

    def name(gid, prod):
        """The protein name with its gene number: 'GRA37 (236890)'."""
        import re
        short = gid.replace("TGME49_", "")
        if gid in KNOWN:
            return f"{KNOWN[gid]} ({short})"
        tok = re.findall(r"\b([A-Z]{2,}[0-9]*[A-Z0-9]*)\b", prod)
        if tok:
            return f"{tok[-1]} ({short})"
        return f"hypothetical ({short})" if "hypothetical" in prod.lower() else f"{short}"
    FS.use()
    figf, ax = rows_fig(len(top), left=2.0, right=0.3)
    yy = np.arange(len(top))
    ax.plot(top["ratio"], yy, "o", color=FS.ACCENT, ms=3.5, alpha=FS.ALPHA)
    ax.axvline(allmed, color=FS.INK, lw=0.6, ls="--")
    ax.set_yticks(yy); ax.set_yticklabels([name(a, b) for a, b in zip(top.index, top["product"])])
    ax.set_xscale("log"); ax.set_xlabel("Nonsynonymous / synonymous")
    ax.set_xticks([1, 3, 10]); ax.set_xticklabels(["1", "3", "10"]); ax.minorticks_off()

    panels = [("a", save(figa, "a_classes")), ("b", save(figb, "b_per_gene")),
              ("c", save(figc, "c_plddt")), ("d", save(figd, "d_variation_vs_order")),
              ("e", save(fige, "e_baits")), ("f", save(figf, "f_top_dense_granules"))]
    import assemble_paper_figures as A_
    out = A_.compose("08", [{"png": p, "letter": l} for l, p in panels], OUT, rows=[3, 3])
    print("figure 8 ->", out)


if __name__ == "__main__":
    main()
