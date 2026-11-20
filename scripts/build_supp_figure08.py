"""Supplementary figure 8a: the design choices behind figure 8, each tested with data.

Figure 8 compares the nonsynonymous / synonymous SNP ratio (N/S) of each class with random genes
MATCHED ON CODING LENGTH, and relates N/S to AlphaFold DB pLDDT. Here:

  A  why match on length: N/S against coding length, per gene
  B  the class ratios against the genome-wide ratio (unmatched) and the length-matched one
  C  the dense-granule ratio against its 2,000 length-matched random sets
  D  pLDDT against coding length: is the N/S-pLDDT link length in disguise?
  E  N/S against pLDDT within length quintiles (the link holds at every length)
  F  a second measure of variation: SNPs per kb, per class (the ranking does not hinge on N/S)
"""
from __future__ import annotations

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
import paper_panels as PP                               # noqa: E402

F = "supplementary_figure_8a"
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
SV = os.path.join(ROOT, "results", "strain_variation_by_class_2026_11_19")
SHORT = {"dense granules (hyperLOPIT)": "dense granules", "carries a signal peptide":
         "signal peptide", "apical secretory (rhoptries, micronemes)": "rhoptries, micronemes",
         "surface and pellicle": "surface, pellicle",
         "required in vitro (lowest fitness tenth)": "essential (lowest tenth)"}
SECRETED = ("dense granules (hyperLOPIT)", "carries a signal peptide",
            "apical secretory (rhoptries, micronemes)")


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    pg = pd.read_csv(os.path.join(SV, "per_gene.tsv"), sep="\t")
    bc = pd.read_csv(os.path.join(SV, "by_class.tsv"), sep="\t").sort_values("ratio",
                                                                            ascending=False)
    n = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    n = n.set_index(n["gene_id"].astype(str))
    g = pg.drop_duplicates("gene_id").set_index("gene_id")
    g["plddt"] = n["mean_plddt"].reindex(g.index)
    g["snp_per_kb"] = (n["snp_total_all_strains"] / (n["length"] * 3 / 1000)).reindex(g.index)
    srcs = [os.path.join(SV, f) for f in ("per_gene.tsv", "by_class.tsv", "summary.json")] + \
        [(os.path.join(CACHE, "nodes.parquet"), ["mean_plddt", "snp_total_all_strains", "length"])]
    gg = g[(g.ratio > 0) & g.length.notna()]

    # A
    r = spearmanr(gg.length, gg.ratio)
    FS.use()
    fa, ax = FS.square_figure()
    ax.plot(gg.length, gg.ratio, "o", color=FS.ACCENT, ms=1.2, alpha=0.3, lw=0, rasterized=True,
            label="gene")
    ax.plot([], [], " ", label=f"Spearman $\\rho$ = {r.statistic:.2f}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks([100, 1000, 10000]); ax.set_xticklabels(["100", "1,000", "10,000"])
    ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels(["0.1", "1", "10"]); ax.minorticks_off()
    ax.set_xlabel("Coding length (aa)"); ax.set_ylabel("Nonsynonymous / synonymous")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "a", fig=fa, data=gg[["length", "ratio"]],
             stats={"Spearman rho": r.statistic, "p": r.pvalue, "genes": len(gg)},
             sources=srcs, what="N/S against coding length")

    # B
    genome = float(np.nanmedian(g.ratio))
    FS.use()
    nr = len(bc); rh = 0.2
    w, h = 1.45 + 1.6 + 1.7, 0.15 + nr * rh + 0.6
    fb = plt.figure(figsize=(w, h)); ax = fb.add_axes([1.45 / w, 0.6 / h, 1.6 / w, nr * rh / h])
    y = np.arange(nr)
    ax.axvline(genome, color=FS.INK, lw=0.6, ls="--", label="genome median (unmatched)")
    ax.plot(bc.length_matched_baseline, y, "o", color=FS.CONTROL, ms=3.5, alpha=FS.ALPHA,
            label="length-matched random")
    ax.plot(bc.ratio, y, "o", color=FS.ACCENT, ms=3.5, alpha=FS.ALPHA, label="class")
    ax.set_yticks(y); ax.set_yticklabels([SHORT.get(c, c) for c in bc["class"]])
    ax.set_ylim(nr - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Nonsynonymous / synonymous")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "b", fig=fb, data=bc.assign(genome_median=genome), stats={
        "genome median N/S": genome, "median length-matched baseline":
        float(bc.length_matched_baseline.median())}, sources=srcs,
             what="class ratios against unmatched and length-matched baselines")

    # C: dense granules against length-matched random sets (redrawn here, 2,000 sets)
    dg = g[g.index.isin(pg.loc[pg["class"] == "dense granules (hyperLOPIT)", "gene_id"])]
    pool = g.dropna(subset=["length", "ratio"]).sort_values("length")
    pos = {gid: i for i, gid in enumerate(pool.index)}
    half = int(len(pool) * 0.05)
    rng = np.random.default_rng(0)
    ratios = pool.ratio.to_numpy()
    obs = float(np.nanmedian(dg.ratio))
    centres = np.array([pos[x] for x in dg.index if x in pos])
    null = []
    for _ in range(2000):
        pick = np.clip(centres + rng.integers(-half, half + 1, len(centres)), 0, len(pool) - 1)
        null.append(float(np.nanmedian(ratios[pick])))
    null = np.array(null)
    FS.use()
    fc, ax = FS.square_figure()
    ax.hist(null, bins=40, color=FS.CONTROL, alpha=FS.ALPHA, lw=0, label="length-matched sets")
    ax.axvline(obs, color=FS.ACCENT, lw=1.4, label="dense granules")
    ax.set_xlabel("Median nonsynonymous / synonymous"); ax.set_ylabel("Random sets")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "c", fig=fc, data=pd.DataFrame({"null_median_ratio": null}),
             stats={"observed median": obs, "null mean": null.mean(), "null sd": null.std(ddof=1),
                    "permutation p": (1 + int((null >= obs).sum())) / (len(null) + 1),
                    "genes": len(centres), "sets": len(null),
                    "window": "+/- 5 % of genes ranked by length"},
             sources=srcs, what="dense granules against length-matched random sets")

    # D: pLDDT vs length
    gl = g.dropna(subset=["plddt", "length"])
    r2 = spearmanr(gl.length, gl.plddt)
    FS.use()
    fd, ax = FS.square_figure()
    ax.plot(gl.length, gl.plddt, "o", color=FS.ACCENT, ms=1.2, alpha=0.3, lw=0, rasterized=True,
            label="gene")
    ax.plot([], [], " ", label=f"Spearman $\\rho$ = {r2.statistic:.2f}")
    ax.set_xscale("log"); ax.set_xticks([100, 1000, 10000])
    ax.set_xticklabels(["100", "1,000", "10,000"]); ax.minorticks_off()
    ax.set_xlabel("Coding length (aa)"); ax.set_ylabel("AlphaFold DB mean pLDDT")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "d", fig=fd, data=gl[["length", "plddt"]], stats={"Spearman rho": r2.statistic,
             "p": r2.pvalue, "genes": len(gl)}, sources=srcs, what="pLDDT against length")

    # E: N/S vs pLDDT within length quintiles
    ge = g[(g.ratio > 0)].dropna(subset=["plddt", "length"]).copy()
    ge["q"] = pd.qcut(ge.length, 5, labels=False)
    rows = []
    for q in range(5):
        s = ge[ge.q == q]; rr = spearmanr(s.ratio, s.plddt)
        rows.append({"length_quintile": q + 1, "min_length": int(s.length.min()),
                     "max_length": int(s.length.max()), "rho": rr.statistic, "p": rr.pvalue,
                     "genes": len(s)})
    e = pd.DataFrame(rows)
    FS.use()
    fe, ax = FS.figure(5, vertical=True)
    FS.bars(ax, e.rho.to_numpy(), [f"{a}-{b}" for a, b in zip(e.min_length, e.max_length)],
            highlight=[f"{a}-{b}" for a, b in zip(e.min_length, e.max_length)])
    ax.axhline(0, color=FS.INK, lw=0.6)
    ax.set_ylabel("Spearman $\\rho$, N/S vs pLDDT"); ax.set_xlabel("Coding length (aa), fifths")
    ax.set_xlim(-0.6, 4.4)
    PP.panel(F, "e", fig=fe, data=ge[["length", "q", "ratio", "plddt"]], stats=e, sources=srcs,
             what="N/S against pLDDT within length fifths")

    # F: SNPs per kb per class
    classes = list(bc["class"])
    groups = {c: g.loc[g.index.isin(pg.loc[pg["class"] == c, "gene_id"]), "snp_per_kb"].dropna()
              for c in classes}
    FS.use()
    nr = len(classes)
    w, h = 1.45 + 1.6 + 0.3, 0.15 + nr * rh + 0.6
    ff = plt.figure(figsize=(w, h)); ax = ff.add_axes([1.45 / w, 0.6 / h, 1.6 / w, nr * rh / h])
    for yy, c in enumerate(classes):
        v = groups[c].to_numpy(); q1, med, q3 = np.percentile(v, [25, 50, 75])
        col = FS.ACCENT if c in SECRETED else FS.CONTROL
        ax.hlines(yy, q1, q3, color=col, lw=1.2, alpha=FS.ALPHA); ax.plot(med, yy, "o", color=col,
                                                                          ms=3.5)
    ax.set_yticks(range(nr)); ax.set_yticklabels([SHORT.get(c, c) for c in classes])
    ax.set_ylim(nr - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("SNPs per kb of coding sequence")
    med = pd.Series({c: groups[c].median() for c in classes})
    rk = spearmanr(bc.set_index("class").ratio.reindex(classes), med.reindex(classes))
    PP.panel(F, "f", fig=ff, data=pd.concat([pd.DataFrame({"class": c, "snp_per_kb": v})
                                             for c, v in groups.items()]),
             stats={"Spearman rho, class N/S vs class SNPs per kb": rk.statistic, "p": rk.pvalue,
                    **{f"median SNPs/kb, {SHORT.get(c, c)}": float(m) for c, m in med.items()}},
             sources=srcs, what="SNPs per kb per class")
    print(e.round(3).to_string(index=False))
    PP.compose(F, "supplementary", "Every design choice behind figure 8, varied or checked: "
               "length matching is needed, and the variation-disorder link is not length.")


if __name__ == "__main__":
    main()
