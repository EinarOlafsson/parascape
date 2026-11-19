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


def newest_or_none(pattern: str):
    """The newest match, or None -- for a result that may not have been produced yet."""
    found = sorted(glob.glob(os.path.join(ROOT, "results", pattern)))
    return found[-1] if found else None


def save(fig, name: str) -> None:
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=600 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.02)
    print("  wrote", name)


def bars(values, labels, ylabel, highlight=(), reference=None, letter="a", height=None,
         fmt="{:.2f}", ref_label="", show_values=True):
    """VERTICAL bars in the house geometry, with the value written above each one.

    Vertical by default now. Nine panels of this paper were horizontal, which is many more than
    "rare occasion": a long category name fits beside a horizontal bar without any thought, and
    that is a reason to shorten the name rather than to turn the graph on its side.

    `height` is the length of the value axis in millimetres; left alone it is the house default,
    so two panels of the same category count match exactly.
    """
    FS.use()
    fig, ax = FS.figure(len(values), height_mm=height or FS.VALUE_MM, vertical=True)
    FS.bars(ax, list(values), list(labels), highlight=set(highlight))
    if reference is not None:
        ax.axhline(reference, color=FS.INK, linewidth=0.6, linestyle=(0, (3, 2)), zorder=4)
        ax.annotate(ref_label or "baseline", (len(values) - 0.5, reference), xytext=(0, 2),
                    textcoords="offset points", fontsize=FS.NOTE_PT, color=FS.INK,
                    ha="right", va="bottom")
    if show_values:
        for i, (v, l) in enumerate(zip(values, labels)):
            ax.annotate(fmt.format(v), (i, v), xytext=(0, 2), textcoords="offset points",
                        ha="center", va="bottom", fontsize=FS.NOTE_PT,
                        color=FS.SERIES[0] if l in highlight else FS.INK)
    ax.set_ylabel(ylabel)
    ax.set_ylim(min(0.0, float(min(values)) * 1.25), float(max(values)) * 1.18)
    ax.set_xlim(-0.6, len(values) - 0.4)
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
    ax.set_xlabel("agreement with the other studies, before")
    ax.set_ylabel("after removing the cell-cycle axis")
    ax.tick_params(labelsize=5.5, pad=1.5)
    ax.annotate(f"{int(better.sum())} of {len(c)} studies agree better\nwith the rest once cell "
                f"cycle is removed", (hi - 0.005, lo + 0.008), fontsize=FS.NOTE_PT, color="#7d7d7d",
                ha="right", va="bottom", linespacing=1.4)
    FS.letter(fig, "a")
    return fig


def figure3():
    m = pd.read_csv(newest("coexpression_modules_*/modules.tsv"), sep="\t")
    named = m[m["named"]].sort_values("share").tail(10)
    # Short labels, because they are now written under a vertical bar. The module number is the
    # identifier and stays; the compartment name loses its underscores.
    labels = [f"{r.compartment.replace('_', ' ')} {int(r.module)}" for r in named.itertuples()]
    fig = bars(list(named["share"]), labels, "share of the module",
               highlight={labels[-1]}, letter="a", fmt="{:.0%}", show_values=False)
    ax = fig.axes[0]
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_yticklabels(["0", "20%", "40%", "60%", "80%"])
    ax.plot(range(len(named)), named["matched_baseline"], "_", color=FS.INK, markersize=5,
            markeredgewidth=0.9, zorder=5)
    ax.annotate("_ matched baseline", (0.99, 0.98), xycoords="axes fraction",
                fontsize=FS.NOTE_PT, color=FS.INK, ha="right", va="top")
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
                    xycoords=("data", "axes fraction"), fontsize=FS.NOTE_PT, ha="center", va="top",
                    color=FS.INK, linespacing=1.35)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([n for n, _g in groups], fontsize=FS.NOTE_PT)
    ax.set_ylabel("pair similarity (|Spearman|)")
    ax.set_ylim(0, 0.72)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", pad=1.5)
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
                xytext=(0, 3), textcoords="offset points", fontsize=FS.NOTE_PT, color=FS.INK, ha="right")
    ax.set_xscale("log")
    # Explicit ticks: the automatic minor labels on a log axis wrote 2x10^2, 3x10^2 and 4x10^2
    # through each other along the bottom.
    from matplotlib.ticker import FixedLocator, NullFormatter
    ax.xaxis.set_major_locator(FixedLocator([100, 300, 1000, 3000]))
    ax.set_xticklabels(["100", "300", "1,000", "3,000"], fontsize=FS.NOTE_PT)
    ax.xaxis.set_minor_formatter(NullFormatter())
    # Headroom at the right so a label beside the rightmost point is not cut off: "coexpression"
    # sat at 3,500 genes and ran off the page.
    ax.set_xlim(90, 13000)
    ax.set_xlabel("genes the layer can answer for")
    ax.set_ylabel("right when it answers")
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
                textcoords="offset points", fontsize=FS.NOTE_PT, color=FS.INK, ha="left", va="bottom")
    for i, row in enumerate(c.itertuples()):
        ax.annotate(f"{row.ratio:.2f}", (i, ceiling * 1.02), fontsize=FS.NOTE_PT, ha="center",
                    va="bottom", color=FS.SERIES[0] if names[i] in top else FS.INK)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_ylabel("nonsynonymous / synonymous")
    ax.set_ylim(0, ceiling * 1.14)
    ax.tick_params(axis="both", pad=1.5)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "a")
    return fig


def figure9():
    L = pd.read_csv(newest("species_transfer_*/layers.tsv"), sep="\t")
    L = L[L["measurable"].astype(bool)].sort_values("transfer")
    short = {"how the message is written (codon bias)": "how it is written",
             "how long the protein is": "protein length (control)",
             "whether losing it costs fitness": "fitness cost",
             "where the protein is": "where it is"}
    labels = [short.get(s, s) for s in L["layer"]]
    # Protein length is the CONTROL and therefore grey: it must transfer between true orthologues.
    fig = bars(list(L["transfer"]), labels, "share that transfers",
               highlight={"where it is"}, letter="a", fmt="{:.2f}")
    ax = fig.axes[0]
    ax.plot(range(len(L)), L["baseline"], "_", color=FS.INK, markersize=5, markeredgewidth=0.9,
            zorder=5)
    ax.annotate("_ chance", (0.99, 0.98), xycoords="axes fraction", fontsize=FS.NOTE_PT,
                color=FS.INK, ha="right", va="top")
    return fig


def figure10():
    """Vertical stacked bars: a composition that adds to the whole, per effector."""
    e = pd.read_csv(newest("host_response_shared_*/per_effector.tsv"), sep="\t")
    e = e.sort_values("host_genes_moved").tail(12)
    names = {"222100": "SOS1"}          # the accession fragment the source uses
    labels = [f"{names.get(r.effector, r.effector)}{'*' if r.translocon else ''}"
              for r in e.itertuples()]
    own = e["its_own_only"].to_numpy()
    shared = (e["host_genes_moved"] - e["its_own_only"]).to_numpy()
    FS.use()
    fig, ax = FS.figure(len(e), vertical=True)
    x = list(range(len(e)))
    ax.bar(x, own, width=FS.bar_width(), linewidth=0, color=FS.SERIES[0], alpha=FS.ALPHA,
           label="this effector only")
    ax.bar(x, shared, bottom=own, width=FS.bar_width(), linewidth=0, color=FS.CONTROL,
           alpha=FS.ALPHA, label="also moved by others")
    ax.set_xticks(x)
    # Gene names italic, the translocon mark upright.
    ax.set_xticklabels([f"$\\it{{{l.rstrip('*')}}}$" + ("*" if l.endswith("*") else "")
                        for l in labels], rotation=30, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("host genes moved")
    ax.set_xlim(-0.6, len(e) - 0.4)
    leg = ax.legend(frameon=False, loc="upper left", handlelength=0.9, borderaxespad=0.2,
                    labelspacing=0.3)
    for text in leg.get_texts():
        text.set_color(FS.INK)
    ax.annotate("* translocon", (0.99, 0.98), xycoords="axes fraction", fontsize=FS.NOTE_PT,
                color=FS.INK, ha="right", va="top")
    FS.letter(fig, "a")
    return fig


#: Shorter still, for a tilted category axis under a vertical bar.
SHORT_TARGETS = {
    "lopit_pf_location": "location (Pf)",
    "lopit_unified": "location (unified)",
    "compartment": "compartment",
    "cellcycle_phase": "cell cycle phase",
    "is_exported": "exported",
    "dtm_class": "functional class",
    "pb_transferred_phenotype": "phenotype from P. berghei",
    "screen_any_phenotype": "any phenotype",
    "screen_actin_phenotype": "actin phenotype",
    "screen_egress_phenotype": "egress phenotype",
    "screen_apicoplast_phenotype": "apicoplast phenotype",
    "screen_replication_phenotype": "replication phenotype",
}

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
    """Vertical bars of the best strategy per target, with every strategy tried drawn over them.

    The long per-bar annotation is gone: it said "77% vs 24%, answers 3%" on every bar, which is
    three numbers twelve times over. Those belong in the legend. What the panel has to show is
    that the best strategy is often the best of a set that is mostly below zero, and that is the
    points.
    """
    b = pd.read_csv(newest("filling_the_map_*/best_per_target.tsv"), sep="\t")
    allv = pd.read_csv(newest("filling_the_map_*/by_target_and_strategy.tsv"), sep="\t")
    b = b.sort_values("lift_over_baseline", ascending=False)
    labels = [SHORT_TARGETS.get(t, TARGET_NAMES.get(t, t)) for t in b["target"]]
    FS.use()
    fig, ax = FS.figure(len(b), vertical=True)
    x = list(range(len(b)))
    wins = (b["lift_over_baseline"] > 0).to_numpy()
    ax.bar(x, b["lift_over_baseline"], width=FS.bar_width(), linewidth=0,
           color=[FS.SERIES[0] if w else FS.CONTROL for w in wins], alpha=FS.ALPHA)
    for i, target in enumerate(b["target"]):
        values = allv.loc[allv["target"] == target, "lift_over_baseline"].dropna().to_numpy()
        if not len(values):
            continue
        offsets = FS._spread(len(values), FS.bar_width() * 0.30)
        ax.plot([i + o for o in offsets], values, "o", markersize=1.8, markerfacecolor="none",
                markeredgecolor=FS.INK, markeredgewidth=0.3, alpha=FS.ALPHA, linestyle="none",
                zorder=5)
    ax.axhline(0, color=FS.INK, linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("above guessing")
    ax.set_xlim(-0.6, len(b) - 0.4)
    low = float(allv["lift_over_baseline"].min())
    ax.set_ylim(min(low, 0.0) * 1.1, float(b["lift_over_baseline"].max()) * 1.2)
    ax.set_yticks([0, 0.2, 0.4, 0.6])
    ax.set_yticklabels(["0", "20%", "40%", "60%"])
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
                fontsize=FS.NOTE_PT, color=FS.INK, ha="left", va="bottom")
    for i, (s, n) in enumerate(zip(subtypes, reached)):
        ax.annotate(f"{n}", (i, 1.03), fontsize=FS.NOTE_PT, ha="center", va="bottom",
                    color=FS.SERIES[0] if s == weakest else FS.INK)
    ax.set_xticks(range(len(subtypes)))
    ax.set_xticklabels([f"subtype {s}" for s in subtypes])
    ax.set_ylabel("q per contrast")
    ax.set_ylim(-0.03, 1.12)
    ax.set_yticks([0, 0.5, 1.0])
    ax.tick_params(axis="both", pad=1.5)
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
    fig = bars(share_absent, labels, "never measured",
               highlight=set(labels[:2]), letter="a", fmt="{:.0%}", show_values=False)
    ax = fig.axes[0]
    ax.set_ylim(0, 1.08)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(["0", "50%", "100%"])
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
                (0.42, 2.5), xytext=(0, 4), textcoords="offset points", fontsize=FS.NOTE_PT,
                color=FS.SERIES[0], ha="left", va="bottom")
    ax.annotate(f"median {int(counts.median())}", (0.0, float(counts.median())), xytext=(14, 0),
                textcoords="offset points", fontsize=FS.NOTE_PT, color=FS.INK, ha="left", va="center")
    ax.set_xticks([0])
    ax.set_xticklabels([f"all {len(counts):,} genes"])
    ax.set_xlim(-0.75, 1.25)
    ax.set_ylabel("categories per gene")
    ax.set_ylim(0, float(counts.max()) + 1.2)
    ax.tick_params(axis="both", pad=1.5)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "b")
    return fig


def figure06b():
    """What sets protein abundance, layer by layer: what each explains alone, and what it adds.

    Transcript abundance alone reaches r2 = 0.49 and everything after it adds almost nothing --
    and one layer ADDS -0.05, which a cumulative-only panel would have drawn as a flat line.
    """
    L = pd.read_csv(newest("what_sets_protein_abundance_*/layers.tsv"), sep="\t")
    short = {"transcript abundance": "transcript",
             "+ how fast it is made and destroyed": "turnover",
             "+ how much is on ribosomes": "ribosomes",
             "+ how the message is written": "message"}
    labels = [short.get(l, l.replace("+ ", "")) for l in L["layer"]]
    FS.use()
    fig, ax = FS.figure(len(L), vertical=True)
    x = list(range(len(L)))
    ax.bar(x, L["r2_alone"], width=FS.bar_width(), linewidth=0, color=FS.CONTROL,
           alpha=FS.ALPHA, label="alone")
    # The gain is a DIFFERENCE and gets a rule, not a second bar beside a level.
    for i, row in enumerate(L.itertuples()):
        ax.plot([i - FS.bar_width() / 2, i + FS.bar_width() / 2], [row.gain, row.gain],
                color=FS.SERIES[0] if row.gain > 0 else FS.DARK_GREY, linewidth=1.3, zorder=5,
                solid_capstyle="butt", label="adds" if i == 0 else None)
    ax.axhline(0, color=FS.INK, linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("protein r$^2$")
    ax.set_xlim(-0.6, len(L) - 0.4)
    leg = ax.legend(frameon=False, loc="upper right", handlelength=0.9, borderaxespad=0.2,
                    labelspacing=0.3)
    for text in leg.get_texts():
        text.set_color(FS.INK)
    FS.letter(fig, "b")
    return fig


#: What each predicted column IS. `te302107_tachy_r1` and `fit_invitro_hff` are record keys.
MEASUREMENT_NAMES = {
    "expr_tachy": "tachyzoite RNA",
    "protein_ibaq_log2": "protein abundance",
    "fit_invitro_hff": "fitness in culture",
    "te302107_tachy_r1": "translation efficiency",
    "snp_nonsynonymous": "nonsynonymous SNPs",
    "n_phosphosites": "phosphosites",
}


def figure06c():
    """Which kind of measurement explains which other kind, ordered by lift over a shuffled run."""
    p = pd.read_csv(newest("which_data_explains_which_*/pairs.tsv"), sep="\t")
    p = p.sort_values("lift", ascending=False)
    labels = [f"{MEASUREMENT_NAMES.get(t, t.replace('_', ' '))}\nfrom {c.replace('_', ' ')}"
              for t, c in zip(p["target"], p["predictor_category"])]
    fig = bars(list(p["lift"]), labels, "above shuffled",
               highlight=set(labels[:2]), letter="c", fmt="{:.2f}", show_values=False)
    return fig


def figure10b():
    """Figure 10b: when an effector goes, does the PARASITE change in the same cell?

    This matters because every host effect attributed to an effector assumes the parasite is
    otherwise unchanged. If knocking out an effector also moves the parasite, part of the host
    response is downstream of a sick parasite rather than of the missing effector.

    One point per knockout, 85 of them, so a scatter: the parasite evidence against the host
    evidence, each as -log10 of the corrected q. Grey for everything and colour for what passes,
    which is the house form for one value per gene against another. The axes are evidence, not
    effect size, and the panel says so -- a target can reach q here on a few hundred cells.

    ONLY ONE SAMPLE IS DRAWN. The other 10x run in this deposit compares against a named neutral
    control (UPRT) and this one against every other assigned cell, so their q values do not mean
    the same thing and pooling them would be the mistake figure 4 is about. The other run is
    reported in the legend instead.
    """
    import glob as _glob
    directory = sorted(_glob.glob(os.path.join(ROOT, "results", "dual_perturb_single_cell_*")))[-1]
    wide = sorted(_glob.glob(os.path.join(directory, "*_targets.tsv")),
                  key=lambda f: -len(open(f).readlines()))[0]
    t = pd.read_csv(wide, sep="\t")
    t = t.dropna(subset=["parasite_q", "host_q"])
    # The AXES ARE THE UNCORRECTED p, and the CORRECTION decides the colour. Drawn on q, 83 of the
    # 85 points sat at exactly q = 1 and the panel was a pile at the origin with two outliers --
    # true, and showing nothing about the other 83. On p the 85 spread out, and a reader can see
    # that they lie where chance would put them, which is the actual finding.
    smallest = 1e-8     # so a p of exactly 0 is still drawable
    x = -np.log10(t["parasite_p"].clip(lower=smallest))
    y = -np.log10(t["host_p"].clip(lower=smallest))
    hit = (t["parasite_q"] < 0.05) | (t["host_q"] < 0.05)
    FS.use()
    fig, ax = FS.figure(6, height_mm=58.0, value_mm=86.0)
    ax.plot(x[~hit], y[~hit], "o", color=FS.CONTROL, markersize=3.2, markeredgewidth=0,
            alpha=FS.ALPHA, zorder=3)
    ax.plot(x[hit], y[hit], "o", color=FS.SERIES[0], markersize=4.2, markeredgewidth=0,
            alpha=FS.ALPHA, zorder=4)
    cut = -np.log10(0.05)
    for value, axis in ((cut, ax.axvline), (cut, ax.axhline)):
        axis(value, color=FS.INK, linewidth=0.5, linestyle=(0, (3, 2)), zorder=2)
    if hit.any():
        FS.place_labels(ax, fig, list(zip(x[hit], y[hit])),
                        [f"{g} ({int(c)} cells)" for g, c in
                         zip(t.loc[hit, "target"], t.loc[hit, "cells"])],
                        [FS.INK] * int(hit.sum()))
    ax.set_xlabel("parasite changed (-log$_{10}$ p)")
    ax.set_ylabel("host changed (-log$_{10}$ p)")
    # TOP right. Top left ran through the host hit's label and bottom right through the parasite
    # hit's: the label placer measures the points, not an annotation added afterwards, so a note
    # has to go where neither hit is -- and both hits are at the left and the bottom right.
    ax.annotate(f"{len(t)} knockouts tested in one 10x run\n"
                f"coloured: q < 0.05 after correction across all {len(t)}\n"
                f"dashed: uncorrected p = 0.05",
                (0.99, 0.97), xycoords="axes fraction", fontsize=FS.NOTE_PT, color=FS.INK,
                ha="right", va="top", linespacing=1.5)
    ax.tick_params(labelsize=5.5, pad=1.5)
    FS.letter(fig, "b")
    return fig


def figure08b():
    """Figure 8b: AF3 interface confidence for host partners against matched controls.

    THE DATA MAY BE SIMULATED. The 24 co-folds are prepared (`prepare_af3_host_cofolds.py`) but
    the GPU step has not been run. Einar, 2026-11-19: "for graphs that have not data yet, put in
    place holder data so i can evaluate the choice of graph." So when the AF3 scores are absent
    this draws the SAME graph on seeded synthetic values and says so on its face, which lets the
    form be judged now and cannot be mistaken for a result.

    THE FORM. Two arms of twelve, which is under the twenty above which a bar is refused, so the
    individual co-folds are drawn as points over a bar of the median. Twelve points is few enough
    that the reader should see every one of them -- with an n this small, the spread is the
    finding or there is no finding.
    """
    pairs = pd.read_csv(newest("af3_host_cofolds_*/pairs.tsv"), sep="\t")
    scored = newest_or_none("af3_host_cofolds_*/interface_scores.tsv")
    simulated = scored is None
    if simulated:
        rng = np.random.default_rng(8)
        # Deliberately a MODEST separation, not a clean one: a placeholder that shows a dramatic
        # effect invites the figure to be designed around a result nobody has.
        centres = {"only in the bait": 0.42, "in both, not enriched": 0.30}
        pairs = pairs.assign(iptm=[float(np.clip(rng.normal(centres[a], 0.13), 0.05, 0.95))
                                   for a in pairs["arm"]])
    else:
        pairs = pairs.merge(pd.read_csv(scored, sep="\t"), on="name", how="inner")
    arms = ["only in the bait", "in both, not enriched"]
    groups = [pairs.loc[pairs["arm"] == a, "iptm"].to_numpy() for a in arms]
    FS.use()
    fig, ax = FS.figure(len(arms), vertical=True)
    medians = [float(np.median(g)) if len(g) else np.nan for g in groups]
    ax.bar(range(len(arms)), medians, width=FS.bar_width(), linewidth=0,
           color=[FS.SERIES[0], FS.CONTROL], alpha=FS.ALPHA)
    for i, values in enumerate(groups):
        offsets = FS._spread(len(values), FS.bar_width() * 0.32)
        ax.plot([i + o for o in offsets], values, "o", markersize=2.6, markerfacecolor="none",
                markeredgecolor=FS.INK, markeredgewidth=0.4, alpha=FS.ALPHA, linestyle="none",
                zorder=5)
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels(["only in bait", "in both, level"], rotation=30, ha="right",
                       rotation_mode="anchor")
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("interface confidence (ipTM)")
    # Top tick at 0.8, not 1.0: a tick label in the top-left corner sits under the panel letter.
    ax.set_ylim(0, 0.95)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xlim(-0.6, len(arms) - 0.4)
    if simulated:
        ax.annotate("SIMULATED", (0.5, 0.99), xycoords="axes fraction",
                    fontsize=FS.NOTE_PT, color=FS.SERIES[0], ha="center", va="top",
                    fontweight="bold")
    FS.letter(fig, "b")
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
                     ("figure06c_explains_which", figure06c),
                     ("figure10b_single_cell", figure10b),
                     ("figure08b_cofold", figure08b)):
        save(fn(), name)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
