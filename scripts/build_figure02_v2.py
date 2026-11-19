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
    b, a = per_study(before), per_study(after)
    order = list((a - b).sort_values(ascending=False).index)
    FS.use()
    fig, ax = FS.figure(len(order), vertical=False, value_mm=40.0)
    y = np.arange(len(order))
    ax.hlines(y, b[order], a[order], color=FS.CONTROL, lw=0.6)
    ax.plot(b[order], y, "o", color=FS.CONTROL, ms=4, alpha=FS.ALPHA, label="before (ref.)")
    ax.plot(a[order], y, "o", color=FS.ACCENT, ms=4, alpha=FS.ALPHA,
            label="cell cycle removed")
    ax.set_yticks(y); ax.set_yticklabels(order); ax.set_ylim(len(order) - 0.5, -0.5)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Median agreement (Spearman $\\rho$)")
    FS.pad_limits(ax, list(b) + list(a), axis="x", frac=0.12)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return fig


def panel_d(null_gain, gain):
    FS.use()
    fig, ax = FS.square_figure()
    ax.hist(null_gain, bins=25, color=FS.CONTROL, alpha=FS.ALPHA, lw=0,
            label="random axes (ref.)")
    ax.axvline(gain, color=FS.ACCENT, lw=1.4, label="cell-cycle axis")
    ax.set_xlabel("Gain in median agreement ($\\Delta\\rho$)"); ax.set_ylabel("Random axes")
    lo = min(null_gain.min(), 0)
    ax.set_xlim(lo - 0.005, max(gain, null_gain.max()) * 1.15)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
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
    cell, cw = 0.16, FS.SLOT_IN          # phase names need a full slot each, or "G1a" runs into "G1b"
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
    ax.set_xlabel("Cell-cycle distance (SD)"); ax.set_ylabel("Agreement before (Spearman $\\rho$)")
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
    panels = [("a", save(heat(before, order, "Agreement (Spearman $\\rho$)", lo, 1), "a_before")),
              ("b", save(heat(after, order, "Agreement (Spearman $\\rho$)", lo, 1), "b_after")),
              ("c", save(panel_c(before, after), "c_per_study")),
              ("d", save(panel_d(null, gain), "d_null")),
              ("e", save(panel_e(levels, position), "e_phase_levels")),
              ("f", save(figf, "f_distance"))]
    import assemble_paper_figures as A
    out = A.compose("02", [{"png": q, "letter": l} for l, q in panels], OUT, rows=[3, 3])
    print("figure 2 ->", out)


if __name__ == "__main__":
    main()
