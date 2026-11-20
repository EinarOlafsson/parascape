"""Figure 2, one claim: published Toxoplasma studies disagree mainly because their controls sit at
different points of the cell cycle; removing that axis brings them into agreement, and a random
axis of the same shape does not.

Einar, 2026-11-19: the point was unclear, the placeholders (cond. 1, group 1) said nothing, and
much could be real already. Every panel here is real:

  A  control-to-control agreement between studies, before                     (the problem)
  B  the same after removing the cell-cycle axis, same scale and order        (the fix)
  C  agreement per study, before -> after (named rows)                        (who moved)
  D  the gain against 200 random axes of the same shape                       (it is not any axis)
  E  where each study's controls sit in the cell cycle: phase-gene levels     (the cause)
  F  per pair of studies: cell-cycle distance vs agreement before             (cause -> problem)

THE CORPUS GUARD. A re-run while the quantification queue rebuilt matrices used 8 studies instead
of 12 and moved the headline gain from +0.055 to +0.031 (results/INVALID_incomplete_corpus_*).
This script refuses to draw unless the studies match the valid 2026-11-18 run exactly.
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
import build_control_comparability as CC                # noqa: E402

OUT = os.path.join(ROOT, "results", "figure02_v2")
PNL = os.path.join(OUT, "panels")
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
OMICS = os.environ.get("PARASCAPE_OMICS", "/nas_mnt/carruthers/Einar/parascape_omics")
VALID = os.path.join(ROOT, "results", "control_comparability_2026_11_18")
PHASE_ORDER = ["G1a", "G1b", "S", "M", "C"]


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def heat(scores, order, title_cb, vmin, vmax):
    s = scores.loc[order, order]
    FS.use()
    n = len(order)
    cell = 0.16
    w, h = 1.15 + n * cell + 0.85, 0.2 + n * cell + 0.95
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([1.15 / w, 0.95 / h, n * cell / w, n * cell / h])
    im = ax.imshow(s.to_numpy(), cmap="Blues", vmin=vmin, vmax=vmax, aspect="equal",
                   interpolation="nearest")
    ax.set_xticks(range(n)); ax.set_xticklabels(order, rotation=90)
    ax.set_yticks(range(n)); ax.set_yticklabels(order)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.07, pad=0.04)
    cb.set_label(title_cb); cb.outline.set_linewidth(0.6)
    return fig


def per_study(scores):
    return scores.mask(np.eye(len(scores), dtype=bool)).median(axis=1)


def panel_c(before, after):
    """Median agreement of each study after removing the axis (bars), before it (tick)."""
    b, a = per_study(before), per_study(after)
    order = list((a - b).sort_values(ascending=False).index)
    FS.use()
    fig, ax = FS.ref_bars(a[order].to_numpy(), order, refs=b[order].to_numpy(),
                          highlight=order, ref_label="before removing the cell-cycle axis",
                          ylabel="Median agreement (Spearman $\\rho$)")
    FS.pad_limits(ax, list(b) + list(a), frac=0.3)
    return fig


def panel_d(null_gain, gain):
    """The gain from the cell-cycle axis against random axes of the same shape (mean, 95 %)."""
    FS.use()
    lo, hi = np.percentile(null_gain, [2.5, 97.5])
    m = float(null_gain.mean())
    fig, ax = FS.ref_bars([m, gain], ["random axes", "cell cycle"], highlight=["cell cycle"],
                          ylabel="Gain in median agreement ($\\Delta\\rho$)")
    ax.errorbar([0], [m], yerr=[[m - lo], [hi - m]], fmt="none", ecolor=FS.INK, lw=0.6)
    ax.set_ylim(0, gain * 1.15)
    return fig


def phase_levels(frame, groups):
    rows = {}
    for p in PHASE_ORDER:
        g = frame.index.intersection(pd.Index(groups.get(p, [])))
        if len(g):
            rows[p] = frame.loc[g].mean(axis=0)
    m = pd.DataFrame(rows)                 # study x phase
    return m.sub(m.mean(axis=0), axis=1)   # each phase centred across studies


def panel_e(levels, position):
    order = list(position.sort_values().index)
    m = levels.loc[order]
    FS.use()
    nr, nc = m.shape
    cell, cw = 0.16, 0.32                # phase names need a full slot each, or "G1a" runs into "G1b"
    w, h = 1.15 + nc * cw + 0.95, 0.2 + nr * cell + 0.55
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([1.15 / w, 0.55 / h, nc * cw / w, nr * cell / h])
    v = float(np.nanmax(np.abs(m.to_numpy())))
    im = ax.imshow(m.to_numpy(), cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto",
                   interpolation="nearest")
    ax.set_xticks(range(nc)); ax.set_xticklabels(m.columns)
    ax.set_yticks(range(nr)); ax.set_yticklabels(order)
    ax.set_xlabel("Cell-cycle phase genes")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.09, pad=0.05)
    cb.set_label("Relative level (log$_2$)"); cb.outline.set_linewidth(0.6)
    return fig


def panel_f(before, position):
    studies = list(before.index)
    xs, ys = [], []
    for i, a in enumerate(studies):
        for b in studies[i + 1:]:
            xs.append(abs(position[a] - position[b])); ys.append(before.loc[a, b])
    xs, ys = np.array(xs), np.array(ys)
    rho = spearmanr(xs, ys).statistic
    FS.use()
    fig, ax = FS.square_figure()
    ax.plot(xs, ys, "o", color=FS.ACCENT, ms=3.5, alpha=FS.ALPHA, lw=0, label="pair of studies")
    k = np.polyfit(xs, ys, 1)
    gx = np.linspace(xs.min(), xs.max(), 50)
    ax.plot(gx, np.polyval(k, gx), color=FS.INK, lw=0.8, label="least-squares line")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.set_xlabel("Cell-cycle distance (SD)"); ax.set_ylabel("Spearman $\\rho$, before correction")
    FS.pad_limits(ax, ys, frac=0.08)
    return fig, rho, len(xs)


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    raw, record = CC.build_profiles(OMICS)
    valid = set(pd.read_csv(os.path.join(VALID, "agreement_before.tsv"), sep="\t",
                            index_col=0).index)
    have = set(raw.columns)
    if have != valid:
        raise SystemExit(f"corpus differs from the valid 2026-11-18 run: missing "
                         f"{sorted(valid - have)}, extra {sorted(have - valid)}. The quant queue "
                         f"is probably rebuilding matrices; wait and re-run.")
    frame = CC.shared_genes(raw)
    phases = CC.phase_table(CACHE)
    groups = {p: list(g.index) for p, g in phases.groupby(phases)}
    before = CC.agreement(frame)
    after = CC.agreement(CC.correct_for(frame, groups))
    base, post = CC.median_offdiagonal(before), CC.median_offdiagonal(after)
    gain = post - base
    rng = np.random.default_rng(0)
    pool = list(frame.index)
    sizes = [len(frame.index.intersection(pd.Index(g))) for g in groups.values()]
    null = []
    for _ in range(CC.N_NULL):
        perm = rng.permutation(pool)
        drawn, i = {}, 0
        for k, n in enumerate(sizes):
            drawn[k] = list(perm[i:i + n]); i += n
        null.append(CC.median_offdiagonal(CC.agreement(CC.correct_for(frame, drawn))) - base)
    null = np.array(null)
    p = (1.0 + int((null >= gain).sum())) / (len(null) + 1.0)
    os.makedirs(OUT, exist_ok=True)
    pd.DataFrame({"null_gain": null}).to_csv(os.path.join(OUT, "null_gains.tsv"), sep="\t",
                                            index=False)
    print(f"{frame.shape[1]} studies, {len(frame)} genes; median agreement {base:.3f} -> "
          f"{post:.3f}, gain {gain:+.3f}; null {null.mean():+.3f} (sd {null.std(ddof=1):.3f}); "
          f"p = {p:.4f}")

    position = pd.Series(CC._covariate(frame, groups), index=frame.columns)
    levels = phase_levels(frame, groups)
    order = list(position.sort_values().index)          # studies by cell-cycle position
    lo = float(np.nanmin(before.where(~np.eye(len(before), dtype=bool)).to_numpy()))
    figf, rho, npairs = panel_f(before, position)
    print(f"cell-cycle distance vs agreement before: Spearman {rho:.3f} over {npairs} pairs")
    import paper_panels as PP
    from scipy.stats import wilcoxon, spearmanr
    F = "figure_2__final"
    summary = os.path.join(VALID, "summary.json")
    tables = [os.path.join(VALID, f) for f in ("agreement_before.tsv",
                                                "agreement_after_cellcycle.tsv", "studies.tsv")]
    nodes_src = (os.path.join(CACHE, "nodes.parquet"), ["cellcycle_phase"])
    matrices = os.path.join(OMICS, "matrices")
    pb, pa = per_study(before), per_study(after)
    w = wilcoxon(pa - pb)
    common = {"studies": frame.shape[1], "genes": len(frame),
              "median agreement before": base, "median agreement after": post}
    # Storyboard 2026-11-20: A agreement before (lower) and after (upper) in one matrix;
    # B what differs: each study by its G1 and S/M subtranscriptome (Behnke et al. 2010);
    # C that difference explains disagreement; D each study improves (slope graph);
    # E it is the cell cycle, not any axis (strip of random axes).
    n_ = len(order)
    split = np.full((n_, n_), np.nan)
    B_, A_ = before.loc[order, order].to_numpy(), after.loc[order, order].to_numpy()
    for r_ in range(n_):
        for c_ in range(n_):
            split[r_, c_] = B_[r_, c_] if r_ > c_ else (A_[r_, c_] if r_ < c_ else np.nan)
    FS.use()
    figa, ax = FS.plot_figure(n_ * 3.4, n_ * 3.4, extra_bottom_mm=12)
    cmap = matplotlib.colormaps["Blues"].copy(); cmap.set_bad("white")
    im = ax.imshow(np.ma.masked_invalid(split), cmap=cmap, vmin=lo, vmax=1,
                   interpolation="nearest")
    ax.plot([-0.5, n_ - 0.5], [-0.5, n_ - 0.5], color=FS.INK, lw=0.6)
    ax.set_xticks(range(n_)); ax.set_xticklabels(order, rotation=90)
    ax.set_yticks(range(n_)); ax.set_yticklabels(order); ax.tick_params(length=0)
    cb = figa.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("Spearman $\\rho$"); cb.outline.set_linewidth(0.6)
    ax.plot([], [], " ", label="below diagonal: before correction")
    ax.plot([], [], " ", label="above diagonal: after correction")
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=FS.NOTE_PT,
              handlelength=0, borderaxespad=0.1)
    PP.panel(F, "a", fig=figa, data=pd.DataFrame(split, index=order, columns=order),
             stats=common, sources=tables + [nodes_src, matrices],
             what="agreement of control profiles: before correction (below the diagonal), after (above)")

    g1 = levels[[c for c in ("G1a", "G1b") if c in levels]].mean(axis=1)
    sm = levels[[c for c in ("S", "M", "C") if c in levels]].mean(axis=1)
    FS.use()
    figb, ax = FS.square_figure()
    ax.scatter(g1, sm, s=18, color=FS.ACCENT, alpha=FS.ALPHA, lw=0, zorder=3, label="study")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02))
    ax.axhline(0, color=FS.BACKGROUND, lw=0.6); ax.axvline(0, color=FS.BACKGROUND, lw=0.6)
    ax.set_xlabel("G1 subtranscriptome (log$_2$, centred)")
    ax.set_ylabel("S/M subtranscriptome (log$_2$, centred)")
    rg = spearmanr(g1, sm)
    PP.panel(F, "b", fig=figb, data=pd.DataFrame({"G1": g1, "S_M": sm, "axis_position": position}),
             stats={"Spearman rho, G1 vs S/M": rg.statistic, "p": rg.pvalue,
                    "subtranscriptomes": "Behnke et al. 2010: G1 (G1a, G1b), S/M (S, M, C)"},
             sources=[nodes_src, matrices], what="cell-cycle state of each study's controls")

    xs = []
    sts = list(before.index)
    for i_, a_ in enumerate(sts):
        for b_ in sts[i_ + 1:]:
            xs.append({"study_a": a_, "study_b": b_,
                       "cell_cycle_distance_sd": abs(position[a_] - position[b_]),
                       "agreement_before": before.loc[a_, b_]})
    fd = pd.DataFrame(xs)
    rr = spearmanr(fd.cell_cycle_distance_sd, fd.agreement_before)
    PP.panel(F, "c", fig=figf, data=fd, stats={"Spearman rho": rr.statistic, "p": rr.pvalue,
                                               "pairs": len(fd)},
             sources=tables + [matrices], what="cell-cycle distance vs agreement, per pair")

    FS.use()
    figd, ax = FS.plot_figure(26, FS.PLOT_H_MM)
    for st in pb.index:
        up = pa[st] > pb[st]
        ax.plot([0, 1], [pb[st], pa[st]], "-o", ms=3, lw=0.8,
                color=FS.ACCENT if up else FS.CONTROL, alpha=0.8)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["before", "after"]); ax.set_xlim(-0.3, 1.3)
    ax.set_ylabel("Median Spearman $\\rho$"); ax.tick_params(axis="x", length=0)
    PP.panel(F, "d", fig=figd,
             data=pd.DataFrame({"before": pb, "after": pa, "gain": pa - pb}).rename_axis("study"),
             stats={"Wilcoxon signed-rank, after vs before": w.statistic, "p": w.pvalue,
                    "studies improved": int((pa > pb).sum()), "studies": len(pa)},
             sources=tables, what="each study, agreement before and after correction")

    FS.use()
    fige, ax = FS.plot_figure(22, FS.PLOT_H_MM)
    rng_ = np.random.default_rng(0)
    ax.scatter(rng_.uniform(-0.18, 0.18, len(null)), null, s=4, color=FS.CONTROL,
               alpha=FS.ALPHA, lw=0, label="random gene sets")
    ax.scatter([0], [gain], s=30, color=FS.ACCENT, zorder=3, label="cell-cycle genes")
    ax.set_xlim(-0.5, 0.5); ax.set_xticks([])
    ax.set_ylabel("Gain in median Spearman $\\rho$"); ax.set_ylim(0, gain * 1.15)
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02))
    PP.panel(F, "e", fig=fige, data=pd.DataFrame({"null_gain": null}),
             stats={"observed gain": gain, "null mean": null.mean(),
                    "null sd": null.std(ddof=1), "null draws": len(null),
                    "permutation p": p, "z vs null": (gain - null.mean()) / null.std(ddof=1)},
             sources=[os.path.join(OUT, "null_gains.tsv"), summary, matrices],
             what="the cell-cycle axis against random axes of the same shape")
    out = PP.compose(F, "main", "Control samples of published RNA-seq studies differ in "
                     "cell-cycle state, and correcting for it brings them into agreement.")
    print("figure 2 ->", out)


if __name__ == "__main__":
    main()
