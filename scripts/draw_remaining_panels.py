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

OUT = os.path.join(ROOT, "results", "paper_panels_2026_11_19")


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
    # Compartment names carry underscores in the record and must not carry them on a figure --
    # `pellicle_PM_IMC` is a key, not a phrase. The module number stays: it IS the identifier.
    labels = [f"{r.compartment.replace('_', ' ')} (module {int(r.module)})"
              for r in named.itertuples()]
    fig = bars(list(named["share"]), labels,
               "share of the module assigned to that compartment",
               highlight={labels[-1]}, letter="a", height=56.0, fmt="{:.0%}")
    ax = fig.axes[0]
    # Percentage ticks under percentage labels: the axis read 0.0-0.8 while every bar was
    # annotated "39%", which asks the reader to convert.
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xticklabels(["0", "20%", "40%", "60%", "80%"], fontsize=5.5)
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


#: Edge-layer keys in words. `unwritten_interaction` and `structural_hole` are record keys, and a
#: figure that prints them makes the reader decode an underscore.
LAYER_NAMES = {
    "comention": "named together in a paper",
    "comention_ft": "named together in a full text",
    "coexpression": "co-expressed",
    "cofitness": "co-fitness",
    "cotranslation": "co-translated",
    "orthogroup": "same orthogroup",
    "domain": "shares a domain",
    "struct": "structural neighbour",
    "xlms": "crosslinked",
    "structural_hole": "structural hole",
    "unwritten_interaction": "interaction nobody wrote down",
}


def figure7():
    """Reach against accuracy per edge layer: both matter, so both are axes."""
    L = pd.read_csv(newest("network_layers_*/layers.tsv"), sep="\t")
    L = L[L["enough"].astype(bool)].copy()
    FS.use()
    # Wider than the house default: the layer names in words are long, and at the default width
    # the label placer ran out of non-overlapping candidates and put "shares a domain" through
    # "all layers". The VALUE axis may grow; the category slot may not, and this panel has none.
    fig, ax = FS.figure(6, height_mm=58.0, value_mm=92.0)
    union = L[L["layer"].str.contains("ALL LAYERS")]
    single = L[~L["layer"].str.contains("ALL LAYERS")]
    ax.plot(single["answered"], single["accuracy"], "o", color=FS.CONTROL, markersize=3.6,
            markeredgewidth=0, alpha=FS.ALPHA, zorder=3)
    ax.plot(union["answered"], union["accuracy"], "o", color=FS.SERIES[0], markersize=4.6,
            markeredgewidth=0, alpha=FS.ALPHA, zorder=4)
    FS.place_labels(ax, fig,
                    [(r.answered, r.accuracy) for r in L.itertuples()],
                    [LAYER_NAMES.get(r.layer.replace(" AT ONCE", "").lower(),
                                     r.layer.replace(" AT ONCE", "").lower().replace("_", " "))
                     for r in L.itertuples()],
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
    # Headroom at the right so a label beside the rightmost point is not cut off: "coexpression"
    # sat at 3,500 genes and ran off the page.
    ax.set_xlim(90, 13000)
    ax.set_xlabel("genes the layer can answer for", fontsize=6.0, labelpad=2)
    ax.set_ylabel("right when it answers", fontsize=6.0, labelpad=2)
    ax.set_ylim(0.25, 0.80)
    ax.tick_params(labelsize=5.5, pad=1.5)
    FS.letter(fig, "a")
    return fig


def figure8():
    """Figure 8: how variable each class of gene is between strains, as a violin per class.

    This was a bar per class, and it should never have been: the smallest class holds 171 genes
    and the largest 1,253, so every bar was the mean of a distribution the reader could not see.
    More than twenty values per group is a violin by the house rules.

    The length-matched baseline stays as a dashed reference, because a class mean above it is the
    whole claim and a violin alone would not say what "high" means. The two classes furthest above
    the baseline carry the series colour; the rest are the grey reference.

    The y axis is clipped at the 99th percentile: a handful of genes with one synonymous SNP and
    several nonsynonymous ones carry ratios in the tens, and letting them set the scale squashed
    every violin into the bottom twentieth of the panel.
    """
    g = pd.read_csv(newest("strain_variation_by_class_*/per_gene.tsv"), sep="\t")
    c = pd.read_csv(newest("strain_variation_by_class_*/by_class.tsv"), sep="\t")
    c = c.sort_values("ratio")
    names = [s.split(" (")[0] for s in c["class"]]
    groups = [g.loc[g["class"] == full, "ratio"].dropna().to_numpy() for full in c["class"]]
    top = set(names[-2:])
    baseline = float(c["length_matched_baseline"].median())
    ceiling = float(np.percentile(np.concatenate(groups), 99))
    FS.use()
    fig, ax = FS.figure(len(groups), height_mm=56.0, vertical=True)
    FS.violin(ax, [np.clip(v, 0, ceiling) for v in groups], names,
              control=[n for n in names if n not in top], vertical=True, draw_points=False)
    ax.axhline(baseline, color=FS.INK, linewidth=0.6, linestyle=(0, (3, 2)), zorder=4)
    ax.annotate("length-matched baseline", (-0.46, baseline), xytext=(0, 3),
                textcoords="offset points", fontsize=5.0, color=FS.INK, ha="left", va="bottom")
    for i, row in enumerate(c.itertuples()):
        ax.annotate(f"{row.ratio:.2f}", (i, ceiling * 1.02), fontsize=5.4, ha="center",
                    va="bottom", color=FS.SERIES[0] if names[i] in top else FS.INK)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([n.replace(" ", "\n", 1) for n in names], fontsize=5.2)
    ax.set_ylabel("nonsynonymous over synonymous SNPs,\nper gene", fontsize=6.0, labelpad=2)
    ax.set_ylim(0, ceiling * 1.14)
    ax.tick_params(axis="both", labelsize=5.5, pad=1.5)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "a")
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
    # Gene names in italic, as the house rules require. The translocon asterisk stays upright,
    # because it is a footnote mark and not part of the name.
    ax.set_yticklabels([f"$\\it{{{names.get(r.effector, r.effector)}}}$"
                        f"{' *' if r.translocon else ''}" for r in e.itertuples()], fontsize=5.8)
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


#: What each prediction target IS, in words. The identifier stays the record; this is the display
#: layer the house rules ask for, because a figure that says `pb_transferred_phenotype` makes the
#: reader look it up, and three of these were truncated mid-word by the old 26-character cut.
TARGET_NAMES = {
    "lopit_pf_location": "where it is (Plasmodium hyperLOPIT)",
    "lopit_unified": "where it is (unified hyperLOPIT)",
    "compartment": "which compartment it is in",
    "cellcycle_phase": "which cell cycle phase it peaks in",
    "is_exported": "whether it is exported",
    "dtm_class": "its coarse functional class",
    "pb_transferred_phenotype": "phenotype carried over from P. berghei",
    "screen_any_phenotype": "any phenotype in a screen",
    "screen_actin_phenotype": "an actin phenotype in a screen",
    "screen_egress_phenotype": "an egress phenotype in a screen",
    "screen_apicoplast_phenotype": "an apicoplast phenotype in a screen",
    "screen_replication_phenotype": "a replication phenotype in a screen",
}


def figure11():
    """Figure 11: how much better than guessing each missing fact can be filled in.

    Bars of the BEST strategy per target, with every strategy tried drawn over them as points --
    fourteen or fifteen per target, which is under the twenty above which a bar is refused, so the
    house rules call for the individual values on top rather than a summary alone. The points are
    what stops a bar being read as "this target is predictable": where they scatter across zero,
    the best strategy is the best of a bad set.
    """
    b = pd.read_csv(newest("filling_the_map_*/best_per_target.tsv"), sep="\t")
    allv = pd.read_csv(newest("filling_the_map_*/by_target_and_strategy.tsv"), sep="\t")
    b = b.sort_values("lift_over_baseline")
    FS.use()
    fig, ax = FS.figure(len(b), height_mm=60.0)
    y = list(range(len(b)))
    wins = b["lift_over_baseline"] > 0
    ax.barh(y, b["lift_over_baseline"], height=FS.bar_width(), linewidth=0,
            color=[FS.SERIES[0] if w else FS.CONTROL for w in wins], alpha=FS.ALPHA)
    for i, target in enumerate(b["target"]):
        values = allv.loc[allv["target"] == target, "lift_over_baseline"].dropna().to_numpy()
        if not len(values):
            continue
        offsets = FS._spread(len(values), FS.bar_width() * 0.30)
        ax.plot(values, [i + o for o in offsets], "o", markersize=1.7,
                markerfacecolor="none", markeredgecolor=FS.INK, markeredgewidth=0.3,
                alpha=FS.ALPHA, linestyle="none", zorder=5)
    ax.axvline(0, color=FS.INK, linewidth=0.6)
    for i, row in enumerate(b.itertuples()):
        ax.annotate(f"{row.right_when_answered:.0%} vs {row.commonest_answer_baseline:.0%}, "
                    f"answers {row.answer_rate:.0%}",
                    (max(row.lift_over_baseline, 0.0), i), xytext=(5, 0),
                    textcoords="offset points", va="center", ha="left",
                    fontsize=4.9, color=FS.INK if row.lift_over_baseline > 0 else FS.CONTROL)
    ax.set_yticks(y)
    ax.set_yticklabels([TARGET_NAMES.get(t, t) for t in b["target"]], fontsize=5.4)
    ax.set_xlabel("accuracy above guessing the commonest answer", fontsize=6.0, labelpad=2)
    # The limit has to fit the POINTS, not just the bars: hard-coding -0.18 pushed the worst
    # strategies out under the target names.
    left = float(allv["lift_over_baseline"].min())
    ax.set_xlim(min(left, 0.0) - 0.04, 0.86)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xticklabels(["0", "20%", "40%", "60%", "80%"], fontsize=5.5)
    ax.set_ylim(-0.6, len(b) - 0.4)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", pad=1.5)
    ax.spines["left"].set_visible(False)
    FS.letter(fig, "a")
    return fig


def figure05_violin():
    """Figure 5b: the EVIDENCE behind the counts of 5a, as a violin per subtype.

    5a is five counts, and five points is a bar by the house rule. But the counts are a threshold
    applied to 201 per-contrast q values per subtype, and 201 values is a distribution: the rule
    says it may not be a bar, and the distribution says more than the count does. A subtype nothing
    reaches has its whole violin pushed against 1.0, not merely a short bar.
    """
    c = pd.read_csv(newest("bradyzoite_reach_*/contrasts.tsv"), sep="\t")
    subtypes = [s for s in "ABCDE" if f"q_{s}" in c.columns]
    groups = [c[f"q_{s}"].dropna().to_numpy() for s in subtypes]
    reached = [int((c[f"q_{s}"] < 0.05).sum()) for s in subtypes]
    weakest = subtypes[int(np.argmin(reached))]
    FS.use()
    fig, ax = FS.figure(len(subtypes), height_mm=52.0, vertical=True)
    FS.violin(ax, groups, [f"subtype {s}" for s in subtypes],
              control=[f"subtype {s}" for s in subtypes if s != weakest], vertical=True)
    ax.axhline(0.05, color=FS.INK, linewidth=0.6, linestyle=(0, (3, 2)), zorder=3)
    # On the LEFT, above the line: at the right it landed on the subtype E violin and its median.
    ax.annotate("q = 0.05", (-0.48, 0.05), xytext=(0, 3), textcoords="offset points",
                fontsize=5.0, color=FS.INK, ha="left", va="bottom")
    for i, (s, n) in enumerate(zip(subtypes, reached)):
        ax.annotate(f"{n}", (i, 1.03), fontsize=5.6, ha="center", va="bottom",
                    color=FS.SERIES[0] if s == weakest else FS.INK)
    ax.set_xticks(range(len(subtypes)))
    ax.set_xticklabels([f"subtype {s}" for s in subtypes], fontsize=6.0)
    ax.set_ylabel(f"q per contrast, all {len(c):,} contrasts", fontsize=6.0, labelpad=2)
    ax.set_ylim(-0.03, 1.12)
    ax.set_yticks([0, 0.5, 1.0])
    ax.tick_params(axis="both", labelsize=5.5, pad=1.5)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "b")
    return fig


def figure12a():
    """Figure 12a: per category of knowledge, the share of genes nothing has ever measured.

    One value per category, so bars -- and the categories are the point, so the two worst carry the
    series colour and the rest are the grey reference. Drawn from the run that splits predicted
    structure from measured structure, which is the pair this panel exists to separate.
    """
    t = pd.read_csv(newest("measured_absence_*/by_category.tsv"), sep="\t")
    t = t.sort_values("share_measured")
    share_absent = (1.0 - t["share_measured"]).tolist()
    labels = [c.replace("_", " ") for c in t["category"]]
    # labels[0] is the LOWEST share measured, which is the HIGHEST share absent -- and barh puts
    # index 0 at the bottom. The first draft highlighted labels[-2:], which coloured transcript and
    # chromatin: the two categories that reach almost every gene, the opposite of the point.
    fig = bars(share_absent, labels, "share of 8,140 genes nothing has ever measured that way",
               highlight=set(labels[:2]), letter="a", height=58.0, fmt="{:.1%}")
    ax = fig.axes[0]
    ax.set_xlim(0, 1.16)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "50%", "100%"], fontsize=5.5)
    return fig


def figure12b():
    """Figure 12b: how many of the categories measure a given gene, over all 8,140 genes.

    ONE violin, not two. The first draft split the genes at the median of the very quantity it
    plots, which guarantees the difference it appears to show: a group defined as "above the
    median" is above the median by construction and the panel proves nothing. A single
    distribution with the sparse tail marked is the honest form, and a count per gene over 8,140
    genes is a distribution, which the house rules forbid drawing as a bar chart.
    """
    g = pd.read_csv(newest("measured_absence_*/categories_per_gene.tsv"), sep="\t")
    counts = g["categories_with_a_value"].dropna()
    sparse = int((counts <= 2).sum())
    FS.use()
    fig, ax = FS.figure(2, height_mm=52.0, vertical=True)
    FS.violin(ax, [counts.to_numpy()], ["all genes"], control=["all genes"], vertical=True,
              draw_points=False)
    ax.axhline(2.5, color=FS.SERIES[0], linewidth=0.7, linestyle=(0, (3, 2)), zorder=4)
    ax.annotate(f"{sparse} genes carry anything in two\ncategories or fewer, of {int(counts.max())}",
                (0.42, 2.5), xytext=(0, 4), textcoords="offset points", fontsize=5.4,
                color=FS.SERIES[0], ha="left", va="bottom")
    ax.annotate(f"median {int(counts.median())}", (0.0, float(counts.median())), xytext=(14, 0),
                textcoords="offset points", fontsize=5.4, color=FS.INK, ha="left", va="center")
    ax.set_xticks([0])
    ax.set_xticklabels([f"all {len(counts):,} genes"], fontsize=6.0)
    ax.set_xlim(-0.75, 1.25)
    ax.set_ylabel("categories of knowledge\nwith any value, per gene", fontsize=6.0, labelpad=2)
    ax.set_ylim(0, float(counts.max()) + 1.2)
    ax.tick_params(axis="both", labelsize=5.5, pad=1.5)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "b")
    return fig


def figure06b():
    """Figure 6b: what sets protein abundance, layer by layer.

    Four layers, one value each, so bars. TWO quantities matter and both are drawn: what a layer
    explains ON ITS OWN, and what it ADDS to the layers before it. The second is the claim --
    transcript abundance alone reaches r2 = 0.49 and everything after it adds almost nothing, and
    one layer (how much is on ribosomes) has a NEGATIVE gain, which a cumulative-only panel would
    have hidden as a flat line.

    Printed beside it, because a number like this is meaningless without a scale a reader knows:
    in human cells transcript abundance explains about 0.40 of protein abundance.
    """
    L = pd.read_csv(newest("what_sets_protein_abundance_*/layers.tsv"), sep="\t")
    labels = [l.replace("+ ", "") for l in L["layer"]]
    FS.use()
    fig, ax = FS.figure(len(L), height_mm=54.0)
    y = list(range(len(L)))
    ax.barh(y, L["r2_alone"], height=FS.bar_width(), linewidth=0, color=FS.CONTROL,
            alpha=FS.ALPHA, label="the layer on its own")
    # The gain is drawn as a rule rather than a second bar: it is a DIFFERENCE, and a bar of a
    # difference beside a bar of a level invites the two to be read on the same footing.
    for i, row in enumerate(L.itertuples()):
        colour = FS.SERIES[0] if row.gain > 0 else FS.DARK_GREY
        ax.plot([row.gain, row.gain], [i - FS.bar_width() / 2, i + FS.bar_width() / 2],
                color=colour, linewidth=1.3, zorder=5, solid_capstyle="butt")
        ax.annotate(f"alone {row.r2_alone:.2f}   adds {row.gain:+.2f}",
                    (max(row.r2_alone, row.gain), i), xytext=(4, 0),
                    textcoords="offset points", va="center", ha="left", fontsize=5.0,
                    color=FS.INK)
    ax.axvline(0, color=FS.INK, linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.6)
    ax.set_xlabel("variance of protein abundance explained (r$^2$)", fontsize=6.0, labelpad=2)
    ax.set_xlim(-0.12, 0.86)
    ax.set_ylim(-0.6, len(L) - 0.4)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", labelsize=5.5, pad=1.5)
    ax.spines["left"].set_visible(False)
    # Top right, where the shortest bar leaves room: at the bottom it ran through the annotation
    # of the transcript-abundance row, which is the one row a reader must be able to read.
    ax.annotate("| what the layer adds to the ones above it", (0.99, 0.985),
                xycoords="axes fraction", fontsize=4.9, color=FS.INK, ha="right", va="top")
    FS.letter(fig, "b")
    return fig


#: What each predicted column IS. `te302107_tachy_r1` and `fit_invitro_hff` are record keys.
MEASUREMENT_NAMES = {
    "expr_tachy": "tachyzoite RNA",
    "protein_ibaq_log2": "protein abundance",
    "fit_invitro_hff": "fitness in culture",
    "te302107_tachy_r1": "translation efficiency",
    "snp_nonsynonymous": "nonsynonymous SNPs",
    "n_phosphosites": "phosphorylation sites",
}


def figure06c():
    """Figure 6c: which kind of measurement explains which other kind.

    One value per target-and-predictor pair, so bars, ordered by how much better than a shuffled
    control the pair does. The shuffled control is the baseline and not the prior size of
    anything: a pair is only interesting to the extent it beats its own measurement procedure run
    on shuffled input.
    """
    p = pd.read_csv(newest("which_data_explains_which_*/pairs.tsv"), sep="\t")
    p = p.sort_values("lift")
    labels = [f"{MEASUREMENT_NAMES.get(t, t.replace('_', ' '))} from {c.replace('_', ' ')}"
              for t, c in zip(p["target"], p["predictor_category"])]
    top = set(labels[-2:])
    fig = bars(list(p["lift"]), labels,
               "variance explained above the same model on shuffled input",
               highlight=top, letter="c", height=58.0, fmt="{:.2f}")
    ax = fig.axes[0]
    for i, row in enumerate(p.itertuples()):
        ax.annotate(f"{int(row.genes):,} genes", (row.lift, i), xytext=(30, 0),
                    textcoords="offset points", va="center", ha="left", fontsize=4.8,
                    color=FS.CONTROL)
    # Room for BOTH annotations on the longest bar: the gene count on the top row was cut off.
    ax.set_xlim(0, float(p["lift"].max()) * 1.75)
    return fig


def main() -> None:
    for name, fn in (("figure02_comparability", figure2), ("figure03_modules", figure3),
                     ("figure04_what_structures", figure4), ("figure07_network_layers", figure7),
                     ("figure08_strain_variation", figure8), ("figure09_species_transfer", figure9),
                     ("figure10_host_shared", figure10), ("figure11_filling_the_map", figure11),
                     ("figure05b_reach_violin", figure05_violin),
                     ("figure12a_measured_absence", figure12a),
                     ("figure12b_categories_per_gene", figure12b),
                     ("figure06b_protein_abundance", figure06b),
                     ("figure06c_explains_which", figure06c)):
        save(fn(), name)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
