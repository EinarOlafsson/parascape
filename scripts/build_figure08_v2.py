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
SHORT = {"dense granules (hyperLOPIT)": "dense gran.",
         "carries a signal peptide": "sig. peptide",
         "apical secretory (rhoptries, micronemes)": "rhop./mic.",
         "surface and pellicle": "pellicle", "mitochondrion": "mito.",
         "apicoplast": "apicoplast", "nucleus": "nucleus", "cytosol": "cytosol",
         "required in vitro (lowest fitness tenth)": "essential"}
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

    # A: bars, a tick for the length-matched random genes of each class
    FS.use()
    sn = [SHORT.get(k, k) for k in order]
    figa, ax = FS.ref_bars(bc["ratio"].to_numpy(), sn, refs=bc["length_matched_baseline"].to_numpy(),
                           highlight=[SHORT.get(c, c) for c in SECRETED],
                           ref_label="length-matched random genes",
                           ylabel="pN/pS")

    # B: per-gene distribution per class
    groups = {c: pg.loc[pg["class"] == c, "ratio"].to_numpy() for c in order}
    allmed = float(np.nanmedian(pg.drop_duplicates("gene_id")["ratio"]))
    FS.use()
    figb, ax = FS.dist_bars([np.clip(groups[c], None, 8) for c in order], sn,
                            highlight=[SHORT.get(c, c) for c in SECRETED],
                            ylabel="pN/pS", ref_line=allmed,
                            ref_label="genome median")

    # C: per-gene pLDDT per class
    pl = n["mean_plddt"]
    groups_c = {c: pl.reindex(pg.loc[pg["class"] == c, "gene_id"].astype(str)).to_numpy()
                for c in order}
    FS.use()
    figc, ax = FS.dist_bars([groups_c[c] for c in order], sn,
                            highlight=[SHORT.get(c, c) for c in SECRETED],
                            ylabel="AlphaFold DB mean pLDDT", ref_line=70,
                            ref_label="pLDDT 70")
    ax.set_ylim(20, 100)

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
    ax.set_ylabel("pN/pS")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    print(f"  ratio vs pLDDT: Spearman {rho:.3f} over {len(dd)} genes")

    # E
    raw = json.load(open(BAITS))
    vals = np.array([v for v in raw.values() if isinstance(v, (int, float))], float)
    FS.use()
    fige, ax = FS.square_figure()
    xs_ = np.sort(vals)
    ax.step(xs_, np.arange(1, len(xs_) + 1) / len(xs_) * 100, where="post", color=FS.ACCENT,
            lw=1.2, label="host-facing proteins")
    ax.axvline(CUT, color=FS.INK, lw=0.6, ls="--", label="co-folding admission cut")
    ax.scatter([float(raw[MYR1])], [float((vals <= raw[MYR1]).mean() * 100)], s=24,
               color=FS.MAGENTA, zorder=4, label="MYR1")
    ax.set_xlabel("Residues with pLDDT > 70 (fraction)"); ax.set_ylabel("Proteins (cumulative %)")
    ax.set_ylim(0, 100)
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2)
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
            return KNOWN[gid]
        tok = re.findall(r"\b([A-Z]{2,}[0-9]*[A-Z0-9]*)\b", prod)
        if tok:
            return tok[-1]
        return short          # gene number for unnamed proteins; full ids in data.csv
    FS.use()
    figf, ax = FS.ref_bars(top["ratio"].to_numpy(), [name(a, b) for a, b in
                                                     zip(top.index, top["product"])],
                           highlight=[name(a, b) for a, b in zip(top.index, top["product"])],
                           ylabel="pN/pS", log=True)
    ax.axhline(allmed, color=FS.INK, lw=0.6, ls=(0, (1.5, 1.5)), label="genome median")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02))
    ax.set_yticks([1, 3, 10]); ax.set_yticklabels(["1", "3", "10"]); ax.minorticks_off()

    import paper_panels as PP
    F = "figure_11__final"
    from scipy.stats import kruskal, mannwhitneyu
    SVF = [os.path.join(SV, f) for f in ("by_class.tsv", "per_gene.tsv", "summary.json")]
    NOD = (os.path.join(CACHE, "nodes.parquet"), ["mean_plddt", "product", "snp_nonsynonymous",
                                                    "snp_synonymous"])
    kw = kruskal(*[v[~np.isnan(v)] for v in groups.values()])
    r_ = spearmanr(dd["ratio"], dd["plddt"])
    PP.panel(F, "a", fig=figb, data=pg, stats=pd.DataFrame(
        [{"class": c, "genes": int(np.sum(~np.isnan(v))), "median": float(np.nanmedian(v))}
         for c, v in groups.items()] + [{"class": "Kruskal-Wallis across classes",
                                         "H": kw.statistic, "p": kw.pvalue}]),
             sources=SVF, what="pN/pS per gene, per protein class")
    sec_v = np.concatenate([v[~np.isnan(v)] for c, v in groups_c.items() if c in SECRETED])
    oth_v = np.concatenate([v[~np.isnan(v)] for c, v in groups_c.items() if c not in SECRETED])
    u = mannwhitneyu(sec_v, oth_v, alternative="less")
    PP.panel(F, "b", fig=figd, data=dd.assign(dense_granule=isdg).rename_axis("gene_id"),
             stats={"Spearman rho": r_.statistic, "p": r_.pvalue, "genes": len(dd)},
             sources=SVF + [NOD], what="per gene: strain variation vs pLDDT")
    PP.panel(F, "c", fig=fige, data=pd.Series(raw, name="folded_fraction").rename_axis(
             "gene_id").reset_index(), stats={"proteins": len(vals), "admission cut": CUT,
             "below cut (%)": float((vals < CUT).mean() * 100), "MYR1": float(raw[MYR1])},
             sources=[BAITS], what="host-facing baits: residues with pLDDT >= 70")
    out = PP.compose(F, "main", "Secreted effectors vary between strains far more than housekeeping genes, and the most variable are the most disordered.", rows=[3, 3])
    print("figure 8 ->", out)


if __name__ == "__main__":
    main()
