"""The eight panels whose results are measured but whose figures were not drawn.

Figures 2, 3, 4, 7, 8, 9, 10 and 11 of the ParaScape paper (instruction 101). Every panel reads
the result table its own script wrote and draws nothing it did not measure. Geometry, palette and
opacity come from :mod:`parascape.figstyle`, so a panel with five categories here is the same
width as a panel with five categories anywhere else in the paper.

The form of each panel is chosen by what the measurement is, per the `publication-figures` skill:
a paired comparison gets a parity diagonal, a distribution of more than twenty values gets a box
rather than a bar with whiskers, and a count across named categories gets bars of uniform width.
"""
from __future__ import annotations

import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import figstyle as FS      # noqa: E402

OUT = os.path.join(ROOT, "results", "paper_panels_2026_11_18")


def newest(pattern: str) -> str:
    found = sorted(glob.glob(os.path.join(ROOT, "results", pattern)))
    if not found:
        raise SystemExit(f"missing: {pattern}")
    return found[-1]


def save(fig, name: str) -> None:
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=600 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.02)
    print("  wrote", name)


def bars(values, labels, xlabel, highlight=(), reference=None, letter="a", height=52.0,
         fmt="{:.2f}", ref_label=""):
    """Horizontal bars of uniform width, control grey, the highlighted ones in series colour."""
    FS.use()
    fig, ax = FS.figure(len(values), height_mm=height)
    y = list(range(len(values)))
    colours = [FS.SERIES[0] if l in highlight else FS.CONTROL for l in labels]
    ax.barh(y, values, height=FS.bar_width(), linewidth=0, color=colours, alpha=FS.ALPHA)
    if reference is not None:
        ax.axvline(reference, color=FS.INK, linewidth=0.6, linestyle=(0, (3, 2)), zorder=4)
        ax.annotate(ref_label or "matched baseline", (reference, len(values) - 0.35),
                    xytext=(3, 0), textcoords="offset points", fontsize=5.0, color=FS.INK,
                    ha="left", va="center")
    for i, (v, l) in enumerate(zip(values, labels)):
        ax.annotate(fmt.format(v), (v, i), xytext=(3, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=5.6,
                    color=FS.SERIES[0] if l in highlight else FS.INK)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.8)
    ax.set_xlabel(xlabel, fontsize=6.0, labelpad=2)
    ax.set_xlim(0, max(values) * 1.20)
    ax.set_ylim(-0.6, len(values) - 0.4)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", labelsize=5.5, pad=1.5)
    ax.spines["left"].set_visible(False)
    FS.letter(fig, letter)
    return fig


def figure2():
    """Each study before and after the cell-cycle correction: a paired comparison, so a parity line."""
    c = pd.read_csv(newest("control_comparability_*/clusters.tsv"), sep="\t", index_col=0)
    c = c.dropna(subset=["median_agreement_before", "median_agreement_after"])
    FS.use()
    fig, ax = FS.figure(6, height_mm=54.0)
    lo, hi = 0.72, 0.95
    ax.plot([lo, hi], [lo, hi], "-", color="#c9c9c9", linewidth=0.7, zorder=1)
    better = c["median_agreement_after"] > c["median_agreement_before"]
    ax.plot(c.loc[~better, "median_agreement_before"], c.loc[~better, "median_agreement_after"],
            "o", color=FS.CONTROL, markersize=3.6, markeredgewidth=0, alpha=FS.ALPHA, zorder=3)
    ax.plot(c.loc[better, "median_agreement_before"], c.loc[better, "median_agreement_after"],
            "o", color=FS.SERIES[0], markersize=3.6, markeredgewidth=0, alpha=FS.ALPHA, zorder=3)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("agreement with the other studies, before", fontsize=6.0, labelpad=2)
    ax.set_ylabel("after removing the cell-cycle axis", fontsize=6.0, labelpad=2)
    ax.tick_params(labelsize=5.5, pad=1.5)
    ax.annotate(f"{int(better.sum())} of {len(c)} studies agree better\nwith the rest once cell "
                f"cycle is removed", (hi - 0.005, lo + 0.008), fontsize=5.2, color="#7d7d7d",
                ha="right", va="bottom", linespacing=1.4)
    FS.letter(fig, "a")
    return fig


def figure3():
    m = pd.read_csv(newest("coexpression_modules_*/modules.tsv"), sep="\t")
    named = m[m["named"]].sort_values("share").tail(10)
    labels = [f"{r.compartment} (module {int(r.module)})" for r in named.itertuples()]
    fig = bars(list(named["share"]), labels,
               "share of the module assigned to that compartment",
               highlight={labels[-1]}, letter="a", height=56.0, fmt="{:.0%}")
    ax = fig.axes[0]
    ax.plot(named["matched_baseline"], range(len(named)), "|", color=FS.INK, markersize=4,
            markeredgewidth=0.8, zorder=5)
    ax.annotate("| matched baseline", (0.98, 0.04), xycoords="axes fraction", fontsize=5.0,
                color=FS.INK, ha="right")
    return fig


def figure4():
    """Distributions of pair similarity: more than twenty values, so a box and never a bar."""
    p = pd.read_csv(newest("what_structures_contrasts_*/pairs.tsv"), sep="\t")
    p["abs_rho"] = p["rho"].abs()
    groups = [("same\nstudy", p[p["same_study"]]["abs_rho"]),
              ("across studies,\nsame kind", p[(~p["same_study"]) & p["same_axis"]]["abs_rho"]),
              ("across studies,\nother kind",
               p[(~p["same_study"]) & (~p["same_axis"])]["abs_rho"])]
    FS.use()
    fig, ax = FS.figure(len(groups), height_mm=52.0, vertical=True)
    box = ax.boxplot([g.to_numpy() for _n, g in groups], widths=FS.bar_width(),
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color=FS.INK, linewidth=0.9),
                     whiskerprops=dict(color=FS.INK, linewidth=0.5),
                     capprops=dict(color=FS.INK, linewidth=0.5),
                     boxprops=dict(linewidth=0))
    for patch, name in zip(box["boxes"], [n for n, _g in groups]):
        # match on the first word, not the whole label: renaming "same study" to
        # "same\nstudy" for the axis silently turned the highlighted box grey.
        patch.set_facecolor(FS.SERIES[0] if name.split()[0] == "same" else FS.CONTROL)
        patch.set_alpha(FS.ALPHA)
    for i, (_n, g) in enumerate(groups, start=1):
        ax.annotate(f"n = {len(g):,}\n{g.median():.3f}", (i, 0.97),
                    xycoords=("data", "axes fraction"), fontsize=5.0, ha="center", va="top",
                    color=FS.INK, linespacing=1.35)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([n for n, _g in groups], fontsize=5.6)
    ax.set_ylabel("similarity of two contrasts (|Spearman|)", fontsize=6.0, labelpad=2)
    ax.set_ylim(0, 0.72)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", labelsize=5.5, pad=1.5)
    FS.letter(fig, "a")
    return fig


def figure7():
    """Reach against accuracy per edge layer: both matter, so both are axes."""
    L = pd.read_csv(newest("network_layers_*/layers.tsv"), sep="\t")
    L = L[L["enough"].astype(bool)].copy()
    FS.use()
    fig, ax = FS.figure(6, height_mm=54.0)
    union = L[L["layer"].str.contains("ALL LAYERS")]
    single = L[~L["layer"].str.contains("ALL LAYERS")]
    ax.plot(single["answered"], single["accuracy"], "o", color=FS.CONTROL, markersize=3.6,
            markeredgewidth=0, alpha=FS.ALPHA, zorder=3)
    ax.plot(union["answered"], union["accuracy"], "o", color=FS.SERIES[0], markersize=4.6,
            markeredgewidth=0, alpha=FS.ALPHA, zorder=4)
    FS.place_labels(ax, fig,
                    [(r.answered, r.accuracy) for r in L.itertuples()],
                    [r.layer.replace(" AT ONCE", "").lower() for r in L.itertuples()],
                    [FS.INK if "ALL" in r.layer else FS.CONTROL for r in L.itertuples()])
    base = float(L["accuracy"].iloc[0] - L["lift"].iloc[0])
    ax.axhline(base, color=FS.INK, linewidth=0.6, linestyle=(0, (3, 2)), zorder=2)
    ax.annotate("guessing the commonest compartment", (0.98, base), xycoords=("axes fraction", "data"),
                xytext=(0, 3), textcoords="offset points", fontsize=5.0, color=FS.INK, ha="right")
    ax.set_xscale("log")
    # Explicit ticks: the automatic minor labels on a log axis wrote 2x10^2, 3x10^2 and 4x10^2
    # through each other along the bottom.
    from matplotlib.ticker import FixedLocator, NullFormatter
    ax.xaxis.set_major_locator(FixedLocator([100, 300, 1000, 3000]))
    ax.set_xticklabels(["100", "300", "1,000", "3,000"], fontsize=5.5)
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(90, 7000)
    ax.set_xlabel("genes the layer can answer for", fontsize=6.0, labelpad=2)
    ax.set_ylabel("right when it answers", fontsize=6.0, labelpad=2)
    ax.set_ylim(0.25, 0.80)
    ax.tick_params(labelsize=5.5, pad=1.5)
    FS.letter(fig, "a")
    return fig


def figure8():
    c = pd.read_csv(newest("strain_variation_by_class_*/by_class.tsv"), sep="\t")
    c = c.sort_values("ratio")
    labels = [s.split(" (")[0] for s in c["class"]]
    top = labels[-2:]
    fig = bars(list(c["ratio"]), labels,
               "nonsynonymous over synonymous SNPs between strains",
               highlight=set(top), reference=float(c["length_matched_baseline"].median()),
               letter="a", height=54.0, ref_label="length-matched baseline")
    return fig


def figure9():
    L = pd.read_csv(newest("species_transfer_*/layers.tsv"), sep="\t")
    L = L[L["measurable"].astype(bool)].sort_values("transfer")
    short = {"how the message is written (codon bias)": "how the message is written",
             "how long the protein is": "how long the protein is (control)"}
    labels = [short.get(s, s) for s in L["layer"]]
    # The CONTROL is grey, by the house rule: protein length must transfer between true
    # orthologues, so it is the reference and not the result. The finding is that where the
    # protein is transfers and how the message is written does not.
    fig = bars(list(L["transfer"]), labels,
               "how much of the layer transfers between species",
               highlight={"where the protein is"}, letter="a", height=48.0, fmt="{:.3f}")
    ax = fig.axes[0]
    ax.plot(L["baseline"], range(len(L)), "|", color=FS.INK, markersize=4, markeredgewidth=0.8,
            zorder=5)
    ax.annotate("| chance", (0.98, 0.06), xycoords="axes fraction", fontsize=5.0, color=FS.INK,
                ha="right")
    return fig


def figure10():
    e = pd.read_csv(newest("host_response_shared_*/per_effector.tsv"), sep="\t")
    e = e.sort_values("host_genes_moved").tail(12)
    FS.use()
    fig, ax = FS.figure(len(e), height_mm=56.0)
    y = list(range(len(e)))
    own, shared = e["its_own_only"].to_numpy(), (e["host_genes_moved"] - e["its_own_only"]).to_numpy()
    ax.barh(y, own, height=FS.bar_width(), linewidth=0, color=FS.SERIES[0], alpha=FS.ALPHA,
            label="moved by this effector only")
    ax.barh(y, shared, left=own, height=FS.bar_width(), linewidth=0, color=FS.CONTROL,
            alpha=FS.ALPHA, label="also moved by others")
    ax.set_yticks(y)
    # 222100 is the accession fragment the source uses for the effector the paper discovered.
    names = {"222100": "SOS1"}
    ax.set_yticklabels([f"{names.get(r.effector, r.effector)}{' *' if r.translocon else ''}"
                        for r in e.itertuples()], fontsize=5.8)
    ax.set_xlabel("host genes moved", fontsize=6.0, labelpad=2)
    ax.set_ylim(-0.6, len(e) - 0.4)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", labelsize=5.5, pad=1.5)
    ax.spines["left"].set_visible(False)
    leg = ax.legend(fontsize=5.2, frameon=False, loc="lower right", handlelength=1.0,
                    borderaxespad=0.2)
    for text in leg.get_texts():
        text.set_color(FS.INK)
    ax.annotate("* translocon", (0.98, 0.30), xycoords="axes fraction", fontsize=5.0,
                color=FS.INK, ha="right")
    FS.letter(fig, "a")
    return fig


def figure11():
    b = pd.read_csv(newest("filling_the_map_*/best_per_target.tsv"), sep="\t")
    b = b.sort_values("lift_over_baseline")
    FS.use()
    fig, ax = FS.figure(len(b), height_mm=56.0)
    y = list(range(len(b)))
    wins = b["lift_over_baseline"] > 0
    ax.barh(y, b["lift_over_baseline"], height=FS.bar_width(), linewidth=0,
            color=[FS.SERIES[0] if w else FS.CONTROL for w in wins], alpha=FS.ALPHA)
    ax.axvline(0, color=FS.INK, linewidth=0.6)
    for i, row in enumerate(b.itertuples()):
        # ALWAYS to the right of zero. Writing the negative ones leftwards put them straight
        # through the target names on the y axis.
        ax.annotate(f"{row.right_when_answered:.0%} vs {row.commonest_answer_baseline:.0%}, "
                    f"answers {row.answer_rate:.0%}",
                    (max(row.lift_over_baseline, 0.0), i), xytext=(4, 0),
                    textcoords="offset points", va="center", ha="left",
                    fontsize=4.9, color=FS.INK if row.lift_over_baseline > 0 else FS.CONTROL)
    ax.set_yticks(y)
    ax.set_yticklabels([t[:26] for t in b["target"]], fontsize=5.6)
    ax.set_xlabel("accuracy above guessing the commonest answer", fontsize=6.0, labelpad=2)
    ax.set_xlim(-0.12, 0.80)
    ax.set_ylim(-0.6, len(b) - 0.4)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", labelsize=5.5, pad=1.5)
    ax.spines["left"].set_visible(False)
    FS.letter(fig, "a")
    return fig


def main() -> None:
    for name, fn in (("figure02_comparability", figure2), ("figure03_modules", figure3),
                     ("figure04_what_structures", figure4), ("figure07_network_layers", figure7),
                     ("figure08_strain_variation", figure8), ("figure09_species_transfer", figure9),
                     ("figure10_host_shared", figure10), ("figure11_filling_the_map", figure11)):
        save(fn(), name)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
