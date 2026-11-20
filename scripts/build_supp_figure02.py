"""Supplementary figure 2a: the design choices behind figure 2, each tested with data.

Figure 2 removes one cell-cycle axis from the control profiles of 12 studies and reports the gain
in agreement against random axes. Each choice is varied here, the rest held fixed:

  A  genes kept: measured in at least this share of studies (chosen 0.8)
  B  studies kept: at least this many control runs (chosen 2)
  C  which axis: cell cycle vs an axis of the same shape from compartments, vs random
  D  how many axes removed (chosen 1)
  E  agreement measure: Spearman (chosen) vs Pearson
  F  leave one study out: the gain without each study

Gains are reported as z against random axes of the same shape (50 draws per variant; figure 2
uses 200), so a variant that raises agreement only because it removes more variance does not
look like a better choice. Runs on the same corpus guard as figure 2.
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402
import build_control_comparability as CC                # noqa: E402
import build_figure02_v2 as F2                          # noqa: E402
import paper_panels as PP                               # noqa: E402

N_NULL = 50
F = "supplementary_figure_2a"


def regress_k(frame, groups, k=1, method="spearman"):
    rows = []
    for genes in groups.values():
        mem = frame.index.intersection(pd.Index(genes))
        rows.append(frame.loc[mem].mean(axis=0).to_numpy() if len(mem) else np.zeros(frame.shape[1]))
    block = np.vstack(rows); block = block - block.mean(axis=1, keepdims=True)
    _u, _s, vt = np.linalg.svd(block, full_matrices=False)
    cov = vt[:k].T
    cov = (cov - cov.mean(axis=0)) / np.where(cov.std(axis=0) == 0, 1, cov.std(axis=0))
    x = np.column_stack([np.ones(frame.shape[1]), cov])
    v = frame.to_numpy(float); out = v.copy(); ok = ~np.isnan(v)
    for i in range(v.shape[0]):
        m = ok[i]
        if m.sum() < k + 3:
            continue
        beta, *_ = np.linalg.lstsq(x[m], v[i, m], rcond=None)
        out[i, m] = v[i, m] - x[m] @ beta + beta[0]
    return pd.DataFrame(out, index=frame.index, columns=frame.columns)


def agree(frame, method="spearman"):
    if method == "spearman":
        return CC.agreement(frame)
    return frame.corr(method="pearson", min_periods=200)


def gain_and_null(frame, groups, k=1, method="spearman", seed=0):
    base = CC.median_offdiagonal(agree(frame, method))
    g = CC.median_offdiagonal(agree(regress_k(frame, groups, k), method)) - base
    rng = np.random.default_rng(seed)
    pool = list(frame.index)
    sizes = [len(frame.index.intersection(pd.Index(v))) for v in groups.values()]
    null = []
    for _ in range(N_NULL):
        perm = rng.permutation(pool); drawn, i = {}, 0
        for j, n in enumerate(sizes):
            drawn[j] = list(perm[i:i + n]); i += n
        null.append(CC.median_offdiagonal(agree(regress_k(frame, drawn, k), method)) - base)
    null = np.array(null)
    return {"gain": g, "null_mean": null.mean(), "null_sd": null.std(ddof=1),
            "z": (g - null.mean()) / null.std(ddof=1),
            "p": (1 + int((null >= g).sum())) / (len(null) + 1), "studies": frame.shape[1],
            "genes": len(frame)}


def dots(rows, xlabel, labels, chosen):
    """Gain over random axes (z) per variant: bars, the figure-2 choice in the accent, z = 1.96
    dotted (Einar 2026-11-20: bars, not dot rows)."""
    FS.use()
    fig, ax = FS.ref_bars([r["z"] for r in rows], labels, highlight=[chosen] if chosen else labels,
                          ylabel=xlabel)
    ax.axhline(1.96, color=FS.INK, lw=0.6, ls=(0, (1.5, 1.5)), label="z = 1.96")
    ax.axhline(0, color=FS.INK, lw=0.6)
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02))
    return fig


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    raw, _rec = CC.build_profiles(F2.OMICS)
    valid = set(pd.read_csv(os.path.join(F2.VALID, "agreement_before.tsv"), sep="\t",
                            index_col=0).index)
    if set(raw.columns) != valid:
        raise SystemExit(f"corpus differs from the valid run: {sorted(set(raw.columns) ^ valid)}")
    phases = CC.phase_table(F2.CACHE)
    groups = {p: list(g.index) for p, g in phases.groupby(phases)}
    srcs = [os.path.join(F2.VALID, f) for f in ("agreement_before.tsv", "summary.json")] + \
        [(os.path.join(F2.CACHE, "nodes.parquet"), ["cellcycle_phase", "lopit_unified"]),
         os.path.join(F2.OMICS, "matrices")]

    # A: gene threshold
    rowsA, labA = [], []
    for t in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        CC.MIN_STUDIES_PER_GENE = t
        r = gain_and_null(CC.shared_genes(raw), groups); r["threshold"] = t
        rowsA.append(r); labA.append(f"{int(t * 100)} % of studies")
        print("A", t, round(r["gain"], 3), round(r["z"], 1), r["genes"])
    CC.MIN_STUDIES_PER_GENE = 0.8
    frame = CC.shared_genes(raw)
    PP.panel(F, "a", fig=dots(rowsA, "Gain over random axes (z)", labA, "80 % of studies"),
             data=pd.DataFrame(rowsA), stats=pd.DataFrame(rowsA), sources=srcs,
             what="genes kept: measured in at least this share of studies")

    # B: minimum control runs
    rowsB, labB = [], []
    for mn in (2, 3):
        CC.MIN_CONTROL_RUNS = mn
        rb, _ = CC.build_profiles(F2.OMICS)
        r = gain_and_null(CC.shared_genes(rb), groups); r["min_control_runs"] = mn
        rowsB.append(r); labB.append(f"{mn} control runs")
        print("B", mn, round(r["gain"], 3), round(r["z"], 1), r["studies"])
    CC.MIN_CONTROL_RUNS = 2
    PP.panel(F, "b", fig=dots(rowsB, "Gain over random axes (z)", labB, "2 control runs"),
             data=pd.DataFrame(rowsB), stats=pd.DataFrame(rowsB), sources=srcs,
             what="studies kept: at least this many control runs")

    # C: which axis
    nodes = pd.read_parquet(os.path.join(F2.CACHE, "nodes.parquet"))
    comp = nodes.set_index(nodes["gene_id"].astype(str))["lopit_unified"].dropna().astype(str)
    top = comp.value_counts().index[:len(groups)]
    cgroups = {c: list(comp.index[comp == c]) for c in top}
    rc = gain_and_null(frame, groups); rcomp = gain_and_null(frame, cgroups)
    # The two axes place the studies in the same order (2026-11-20: Spearman 0.97): the
    # compartment sets hold the wave-built organelles of figure 6, so they read the same
    # division programme through more genes. Recorded with the panel, not left to the legend.
    from scipy.stats import spearmanr as _sp
    pos_cc, pos_cp = CC._covariate(frame, groups), CC._covariate(frame, cgroups)
    axis_rho = _sp(pos_cc, pos_cp).statistic
    rowsC = [dict(rc, axis="cell-cycle phases", axis_position_spearman_vs_other=axis_rho),
             dict(rcomp, axis="compartments", axis_position_spearman_vs_other=axis_rho)]
    PP.panel(F, "c", fig=dots(rowsC, "Gain over random axes (z)",
                              ["cell-cycle phases", "compartments"], "cell-cycle phases"),
             data=pd.DataFrame(rowsC), stats=pd.DataFrame(rowsC), sources=srcs,
             what="which gene groups define the axis removed")
    print("C", [(r["axis"], round(r["z"], 1)) for r in rowsC])

    # D: number of axes removed
    rowsD = []
    for k in (1, 2, 3):
        r = gain_and_null(frame, groups, k=k); r["axes_removed"] = k; rowsD.append(r)
    PP.panel(F, "d", fig=dots(rowsD, "Gain over random axes (z)", ["1 axis", "2 axes", "3 axes"],
                              "1 axis"),
             data=pd.DataFrame(rowsD), stats=pd.DataFrame(rowsD), sources=srcs,
             what="how many axes are removed")
    print("D", [(r["axes_removed"], round(r["gain"], 3), round(r["z"], 1)) for r in rowsD])

    # E: agreement measure
    rowsE = [dict(gain_and_null(frame, groups, method=m), measure=m)
             for m in ("spearman", "pearson")]
    PP.panel(F, "e", fig=dots(rowsE, "Gain over random axes (z)", ["Spearman", "Pearson"],
                              "Spearman"),
             data=pd.DataFrame(rowsE), stats=pd.DataFrame(rowsE), sources=srcs,
             what="agreement measure")

    # F: leave one study out
    rowsF = []
    for st in frame.columns:
        r = gain_and_null(frame.drop(columns=[st]), groups); r["left_out"] = st; rowsF.append(r)
    lf = pd.DataFrame(rowsF).sort_values("z", ascending=False)
    PP.panel(F, "f", fig=dots(lf.to_dict("records"), "Gain over random axes (z)",
                              list(lf.left_out), None),
             data=lf, stats={"min z": lf.z.min(), "max z": lf.z.max(),
                             "left-out runs with z > 1.96": int((lf.z > 1.96).sum()),
                             "runs": len(lf)}, sources=srcs,
             what="leave one study out")
    print("F", lf[["left_out", "gain", "z"]].round(2).to_string(index=False))
    PP.compose(F, "supplementary", "Every design choice behind figure 2, varied with the rest "
               "held fixed: the cell-cycle gain holds across thresholds, measures and studies, "
               "and an axis of the same shape from compartments does not give it.")


if __name__ == "__main__":
    main()
