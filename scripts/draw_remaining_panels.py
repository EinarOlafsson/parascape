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
import matplotlib.ticker
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
         ref_label="", percent=False):
    """Vertical bars in the house geometry. NOTHING is written on the plotting area.

    Einar, 2026-11-19: "there should only be text on x and y axes or in legends". So there is no
    value above a bar, no baseline label beside its line, no n and no note in a corner. A
    reference line is identified by a LEGEND ENTRY; every count belongs in the figure legend.
    """
    FS.use()
    fig, ax = FS.figure(len(values), height_mm=height or FS.VALUE_MM, vertical=True)
    FS.bars(ax, list(values), list(labels), highlight=set(highlight))
    if reference is not None:
        ax.axhline(reference, color=FS.INK, linewidth=0.6, linestyle=(0, (3, 2)), zorder=4,
                   label=ref_label or "baseline")
        FS.legend(ax, loc="upper right")
    ax.set_ylabel(ylabel)
    span = list(values) + ([reference] if reference is not None else []) + [0.0]
    FS.pad_limits(ax, span)
    ax.set_xlim(-0.6, len(values) - 0.4)
    if percent:
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _p: f"{v * 100:.0f}"))
    FS.letter(fig, letter)
    return fig


def figure2():
    """Each study before and after the cell-cycle correction: a paired comparison, so a parity line."""
    c = pd.read_csv(newest("control_comparability_*/clusters.tsv"), sep="\t", index_col=0)
    c = c.dropna(subset=["median_agreement_before", "median_agreement_after"])
    FS.use()
    fig, ax = FS.square_figure()
    lo, hi = 0.72, 0.95
    ax.plot([lo, hi], [lo, hi], "-", color="#c9c9c9", linewidth=0.7, zorder=1)
    better = c["median_agreement_after"] > c["median_agreement_before"]
    ax.plot(c.loc[~better, "median_agreement_before"], c.loc[~better, "median_agreement_after"],
            "o", color=FS.CONTROL, markersize=3.6, markeredgewidth=0, alpha=FS.ALPHA, zorder=3)
    ax.plot(c.loc[better, "median_agreement_before"], c.loc[better, "median_agreement_after"],
            "o", color=FS.SERIES[0], markersize=3.6, markeredgewidth=0, alpha=FS.ALPHA, zorder=3)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("Spearman $\\rho$, all axes")
    ax.set_ylabel("Spearman $\\rho$, cell cycle removed")
    FS.letter(fig, "")          # figure 2 is a single panel, so it carries no letter

    return fig


def figure3():
    m = pd.read_csv(newest("coexpression_modules_*/modules.tsv"), sep="\t")
    named = m[m["named"]].sort_values("share").tail(10)
    # Short labels, because they are now written under a vertical bar. The module number is the
    # identifier and stays; the compartment name loses its underscores.
    labels = [f"{r.compartment.replace('_', ' ')} {int(r.module)}" for r in named.itertuples()]
    fig = bars(list(named["share"]), labels, "Genes in compartment (%)",
               highlight={labels[-1]}, letter="", percent=True)
    ax = fig.axes[0]
    ax.plot(range(len(named)), named["matched_baseline"], "_", color=FS.INK, markersize=5,
            markeredgewidth=0.9, zorder=5, linestyle="none", label="Size-matched random set")
    FS.legend(ax, loc="upper left")
    return fig


def figure4():
    """Distributions of pair similarity: more than twenty values, so a box and never a bar."""
    p = pd.read_csv(newest("what_structures_contrasts_*/pairs.tsv"), sep="\t")
    p["abs_rho"] = p["rho"].abs()
    # SHORT labels, tilted. Spelled out in the figure legend, per the skill: a 3-condition panel
    # is 38 mm wide and three two-line phrases written flat ran straight through each other.
    groups = [("Same study", p[p["same_study"]]["abs_rho"]),
              ("Same kind", p[(~p["same_study"]) & p["same_axis"]]["abs_rho"]),
              ("Other kind", p[(~p["same_study"]) & (~p["same_axis"])]["abs_rho"])]
    FS.use()
    fig, ax = FS.figure(len(groups), height_mm=52.0, vertical=True)
    box = ax.boxplot([g.to_numpy() for _n, g in groups], widths=FS.bar_width(),
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color=FS.INK, linewidth=0.9),
                     whiskerprops=dict(color=FS.INK, linewidth=0.5),
                     capprops=dict(color=FS.INK, linewidth=0.5),
                     boxprops=dict(linewidth=0))
    # The accent is on the FIRST group by position, not by matching its label text. Matching on
    # the label has now silently turned this box grey TWICE, once for a line break and once for a
    # capital letter; position cannot drift when a label is reworded.
    for i, patch in enumerate(box["boxes"]):
        patch.set_facecolor(FS.ACCENT if i == 0 else FS.CONTROL)
        patch.set_alpha(FS.ALPHA)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([n for n, _g in groups], rotation=30, ha="right",
                       rotation_mode="anchor")
    ax.set_ylabel("|Spearman $\\rho$|")
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
    fig, ax = FS.square_figure()
    union = L[L["layer"].str.contains("ALL LAYERS")]
    single = L[~L["layer"].str.contains("ALL LAYERS")]
    ax.plot(single["answered"], single["accuracy"], "o", color=FS.CONTROL, markersize=3.6,
            markeredgewidth=0, alpha=FS.ALPHA, zorder=3)
    ax.plot(union["answered"], union["accuracy"], "o", color=FS.SERIES[0], markersize=4.6,
            markeredgewidth=0, alpha=FS.ALPHA, zorder=4)
    # Label ONLY what the reader should look at: the two union points, the most accurate layer
    # and the broadest one. Thirteen labels at 8 pt in a 46 mm square ran through each other, and
    # the skill is explicit that everything else is the grey cloud and belongs in the legend.
    show = set(L.loc[L["layer"].str.contains("ALL LAYERS"), "layer"])
    show.add(L.loc[L["accuracy"].idxmax(), "layer"])
    show.add(L.loc[L["answered"].idxmax(), "layer"])
    marked = L[L["layer"].isin(show)]
    FS.place_labels(ax, fig,
                    [(r.answered, r.accuracy) for r in marked.itertuples()],
                    [LAYER_NAMES.get(r.layer.replace(" AT ONCE", "").lower(),
                                     r.layer.replace(" AT ONCE", "").lower().replace("_", " "))
                     for r in marked.itertuples()],
                    [FS.INK if "ALL" in r.layer else FS.CONTROL for r in marked.itertuples()])
    base = float(L["accuracy"].iloc[0] - L["lift"].iloc[0])
    ax.axhline(base, color=FS.INK, linewidth=0.6, linestyle=(0, (3, 2)), zorder=2)
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
    ax.set_xlabel("Genes answered")
    ax.set_ylabel("Accuracy (%)")
    # The title says per cent, so the ticks must read per cent and not 0.3-0.8.
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _p: f"{v * 100:.0f}"))
    ax.set_ylim(0.25, 0.80)
    ax.tick_params(pad=1.5)
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
    ax.set_ylabel("Nonsynonymous / synonymous SNPs")
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
    # CHANCE-CORRECTED, so the axis is ONE quantity. `transfer` is Spearman for the three
    # continuous layers and categorical agreement for "where it is", whose chance baseline is
    # 0.374 rather than 0 -- so the raw values are not the same measurement and no objective axis
    # title was possible. (observed - chance) / (1 - chance) is Cohen's kappa: 0 is chance for
    # every layer, 1 is perfect, and the four bars become comparable.
    kappa = (L["transfer"] - L["baseline"]) / (1.0 - L["baseline"])
    fig = bars(list(kappa), labels, "Chance-corrected agreement ($\\kappa$)",
               highlight={"where it is"}, letter="")
    ax = fig.axes[0]
    ax.axhline(0, color=FS.INK, linewidth=0.6)      # chance, which kappa puts at zero
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
    ax.set_ylabel("Host genes")
    ax.set_xlim(-0.6, len(e) - 0.4)
    leg = ax.legend(frameon=False, loc="upper left", handlelength=0.9, borderaxespad=0.2,
                    labelspacing=0.3)
    for text in leg.get_texts():
        text.set_color(FS.INK)
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
    ax.set_ylabel("Accuracy above baseline (%)")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _p: f"{v * 100:.0f}"))
    ax.set_xlim(-0.6, len(b) - 0.4)
    # Cut at the 5th percentile of the strategy points, and SAY SO in the figure legend. A single
    # strategy at -0.9 was setting the bottom of the axis, so two thirds of the panel was empty
    # white space below the data -- which the house rules forbid.
    low = float(np.percentile(allv["lift_over_baseline"].dropna(), 5))
    ax.set_ylim(min(low, -0.05) * 1.15, float(b["lift_over_baseline"].max()) * 1.18)
    ax.set_yticks([-0.2, 0, 0.2, 0.4, 0.6])
    # The unit is in the axis title, so the ticks are plain numbers.
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _p: f"{v * 100:.0f}"))
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
    ax.set_xticks(range(len(subtypes)))
    ax.set_xticklabels(list(subtypes))
    ax.set_ylabel("$q$")
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
    fig = bars(share_absent, labels, "Genes never measured (%)",
               highlight=set(labels[:2]), letter="a", percent=True)
    ax = fig.axes[0]
    ax.set_ylim(0, 1.04)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
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
    ax.set_xticks([0])
    ax.set_xticklabels(["All genes"])
    ax.set_xlim(-0.75, 1.25)
    ax.set_ylabel("Categories per gene")
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
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
    ax.set_ylabel("Protein abundance explained ($r^2$)")
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
    # ONE line, tilted. Fourteen two-line labels in a 14-condition panel overlapped each other
    # 64 times by the checker and were unreadable; the pairing is spelled out in the legend.
    labels = [f"{MEASUREMENT_NAMES.get(t, t.replace('_', ' '))} / {c.replace('_', ' ')}"
              for t, c in zip(p["target"], p["predictor_category"])]
    fig = bars(list(p["lift"]), labels, "above shuffled",
               highlight=set(labels[:2]), letter="c")
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
    fig, ax = FS.square_figure()
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
    ax.set_xlabel("$-$log$_{10}$ $p$, parasite")
    ax.set_ylabel("$-$log$_{10}$ $p$, host")
    ax.tick_params(pad=1.5)
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
    ax.set_ylabel("ipTM")
    # Top tick at 0.8, not 1.0: a tick label in the top-left corner sits under the panel letter.
    ax.set_ylim(0, 0.95)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xlim(-0.6, len(arms) - 0.4)
    if simulated:
        # In the LEGEND, not on the plotting area -- but it must appear somewhere on the figure,
        # because a panel of invented numbers that looks like a result is the worst thing here.
        ax.plot([], [], " ", label="Simulated: co-folds prepared, not run")
        FS.legend(ax, loc="upper center")
    FS.letter(fig, "b")
    return fig


def figure08b_foldability():
    """Figure 8b, redesigned: which effectors AlphaFold can see at all.

    The co-folding panel this replaces could not work. MYR1 has a folded fraction of 0.0 -- zero
    of 810 residues above AFDB pLDDT 70 -- and the AF3 screen (session repo-1f) established over
    ~400 PDB complexes that AF3 maps the folded, domain-mediated fraction (median DockQ 0.92) and
    is AT CHANCE for disordered, motif-mediated binding (iLIS AUC 0.515). Co-folding a bait with
    no ordered residues returns noise, which is what ipTM 0.19 was.

    So the honest panel is the limit itself: the distribution of foldability across the 208 baits,
    the cut below which co-folding is not attempted, and where MYR1 falls. 208 values is far more
    than twenty, so a violin; the cut and the bait are legend entries, not text on the plot.
    """
    import json as _json
    path = ("/media/carruthers/mnt3/claude/toxoplasma_projects/af3_folding/reports/"
            "bait_folded_fraction.json")
    raw = _json.load(open(path))
    values = np.array([v for v in raw.values() if isinstance(v, (int, float))], dtype=float)
    myr1 = float(raw["TGME49_254470"])
    cut = 0.30
    FS.use()
    fig, ax = FS.figure(2, vertical=True)
    # 208 host-facing genes WITH an AFDB entry. Eight of the 216 have none and are absent from
    # the file rather than scored as zero, so the denominator is 208 -- not 216, and not the
    # proteome. Without that, 21 per cent below the cut reads as a proteome-wide disorder rate.
    FS.violin(ax, [values], ["Host-facing"], control=["Host-facing"], vertical=True,
              draw_points=False)
    # "Screen admission cut", NOT "disorder threshold": 0.30 is the AF3 screen's operational
    # threshold for admitting a bait to co-folding, not an established definition of disorder.
    ax.axhline(cut, color=FS.INK, linewidth=0.6, linestyle=(0, (2, 2)), zorder=4,
               label="Screen admission cut")
    ax.plot([0], [myr1], "o", markersize=4.5, markerfacecolor=FS.ACCENT,
            markeredgecolor=FS.ACCENT, alpha=FS.ALPHA, zorder=6, linestyle="none",
            label="MYR1")
    # Centre right: the upper right sat over the top of the violin, and a legend never goes
    # over data.
    FS.legend(ax, loc="center right")
    ax.set_xticks([0])
    ax.set_xticklabels(["Host-facing"])
    # Name the QUANTITY: this is not an AF3 output, it is the fraction of residues whose
    # AlphaFold DB pLDDT is at or above 70.
    ax.set_ylabel("Residues with pLDDT $\\geq$ 70 (fraction)")
    ax.set_ylim(-0.06, 1.06)
    ax.set_xlim(-0.75, 1.25)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "b")
    return fig


def figure3b():
    """3b: every module against its own size-matched baseline, not just the named ten.

    3a shows ten named modules. This shows ALL 24 against the chance share for a set of that
    size, so a reader can see how many clear it and how many do not -- which 3a cannot say
    because it only draws the winners. A paired comparison, so a parity diagonal.
    """
    m = pd.read_csv(newest("coexpression_modules_*/modules.tsv"), sep="\t")
    FS.use()
    fig, ax = FS.square_figure()
    hi = float(max(m["share"].max(), m["matched_baseline"].max())) * 1.1
    ax.plot([0, hi], [0, hi], "-", color=FS.BACKGROUND, linewidth=0.7, zorder=1)
    beats = m["share"] > m["matched_baseline"]
    ax.plot(m.loc[~beats, "matched_baseline"], m.loc[~beats, "share"], "o", color=FS.CONTROL,
            markersize=3.4, markeredgewidth=0, alpha=FS.ALPHA, zorder=3)
    ax.plot(m.loc[beats, "matched_baseline"], m.loc[beats, "share"], "o", color=FS.ACCENT,
            markersize=3.4, markeredgewidth=0, alpha=FS.ALPHA, zorder=4)
    ax.set_xlabel("Size-matched baseline (fraction)")
    ax.set_ylabel("Module in one compartment (fraction)")
    ax.set_xlim(0, hi); ax.set_ylim(0, hi)
    FS.letter(fig, "b")
    return fig


def figure3c():
    """3c: how big the modules are. 24 values, so a violin, with the genes they hold."""
    m = pd.read_csv(newest("coexpression_modules_*/modules.tsv"), sep="\t")
    FS.use()
    fig, ax = FS.figure(2, vertical=True)
    FS.violin(ax, [m["genes"].to_numpy()], ["Modules"], control=["Modules"], vertical=True)
    ax.set_xticks([0]); ax.set_xticklabels(["Modules"])
    ax.set_ylabel("Genes per module")
    ax.set_xlim(-0.75, 1.25)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "c")
    return fig


def figure3d():
    """3d: how far each module sits above the DETECTION LIMIT of the matched null.

    A share above baseline is not yet a finding: the question is whether the gap exceeds what the
    null's own spread admits. This is that gap per module, which 3a and 3b both leave implicit.
    """
    m = pd.read_csv(newest("coexpression_modules_*/modules.tsv"), sep="\t")
    m = m.sort_values("over_detection_limit", ascending=False)
    # An UNNAMED module has no compartment (NaN), and it still belongs on this panel: whether a
    # module clears the detection limit is a fact about the module, not about whether we named it.
    labels = [f"{str(c).replace('_', ' ')} {int(n)}" if isinstance(c, str) else f"module {int(n)}"
              for c, n in zip(m["compartment"], m["module"])]
    fig = bars(list(m["over_detection_limit"]), labels,
               "Share above the detection limit", highlight=set(labels[:2]), letter="d",
               reference=1.0, ref_label="Detection limit")
    return fig


def main() -> None:
    for name, fn in (("figure02_comparability", figure2), ("figure03_modules", figure3), ("figure03b_vs_baseline", figure3b),
                     ("figure03c_sizes", figure3c), ("figure03d_detection", figure3d),
                     ("figure04_what_structures", figure4), ("figure07_network_layers", figure7),
                     ("figure08_strain_variation", figure8), ("figure09_species_transfer", figure9),
                     ("figure10_host_shared", figure10), ("figure11_filling_the_map", figure11),
                     ("figure05b_reach_violin", figure05_violin),
                     ("figure12a_measured_absence", figure12a),
                     ("figure12b_categories_per_gene", figure12b),
                     ("figure06b_protein_abundance", figure06b),
                     ("figure06c_explains_which", figure06c),
                     ("figure10b_single_cell", figure10b),
                     ("figure08b_cofold", figure08b),
                     ("figure08b_foldability", figure08b_foldability)):
        save(fn(), name)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
