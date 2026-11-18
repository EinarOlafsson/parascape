"""A paper panel: which evidence places a protein in its compartment.

For the ParaScape paper, figure 5 -- which measurement explains which other measurement. This panel
is the localisation instance of it: the 13 coarse hyperLOPIT niches recovered from RNA co-expression
alone, against the same recovery from all allowed evidence together.

WHY A SLOPE GRAPH. The quantity is a PAIRED difference -- the same niche scored two ways -- and the
thing worth seeing is which niches cross over, where RNA alone matches or beats everything else. A
bar chart of two series hides the pairing, and a heatmap hides it twice over. Two columns of points
joined per niche shows the crossing directly, and the four that cross are the claim.

Reads `results/coexpression_components_localization_2026_11_06/scores.tsv`, written by the notebook
of instruction 71 product 1. Nothing is recomputed here: a figure script that re-derives its numbers
is a second copy of the analysis.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager

#: Where the scored table is, relative to the repository root.
SOURCE = ("results/coexpression_components_localization_2026_11_06/scores.tsv")

#: Millimetres, because a paper panel is specified in the units the page is in.
WIDTH_MM, HEIGHT_MM = 105.0, 58.0

#: The two columns compared, in the order they are drawn.
LEFT, RIGHT = "coexp_cv_auc", "all_evidence_cv_auc"

#: Ink. The crossing niches are the claim, so they are the only ones that carry colour.
INK = "#1a1a1a"
MUTED = "#9a9a9a"
WINS = "#1f6f8b"


def _font() -> None:
    """Open Sans everywhere, which is the typeface of the rest of the figures."""
    # Only Regular, Italic and Bold of Open Sans are installed on this machine, so asking for
    # Light or SemiBold silently falls back and prints a warning for every text object. Regular
    # for the body and Bold for the panel letter is the nearest honest pair.
    if "Open Sans" in {f.name for f in font_manager.fontManager.ttflist}:
        plt.rcParams["font.family"] = "Open Sans"
        plt.rcParams["font.weight"] = "normal"
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                         "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                         "ytick.major.width": 0.5, "font.size": 6.2})


def load(root: str) -> pd.DataFrame:
    """The scored niches, best-recovered first, with the paired difference."""
    frame = pd.read_csv(os.path.join(root, SOURCE), sep="\t")
    for column in (LEFT, RIGHT):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=[LEFT, RIGHT]).copy()
    frame["gain"] = frame[RIGHT] - frame[LEFT]
    return frame.sort_values(LEFT, ascending=False).reset_index(drop=True)


def _spread(values, low: float, high: float):
    """Label positions in value order, evenly over the whole height.

    Repelling labels locally left them bunched at the top -- most compartments score above 0.9 --
    with the leader lines crossing each other and the lower half of the panel empty. Spreading them
    evenly in VALUE ORDER cannot cross: if label i is above label j exactly when point i is above
    point j, no two leaders can intersect. The same argument as the cluster labels of figure 1.
    """
    order = sorted(range(len(values)), key=lambda i: values[i])
    step = (high - low) / max(1, len(values) - 1)
    out = [0.0] * len(values)
    for rank, i in enumerate(order):
        out[i] = low + rank * step
    return out


def draw(frame: pd.DataFrame, letter: str = "a"):
    """The slope graph. Returns the figure.

    One line per niche from RNA-alone to all-evidence. A line that falls is a niche where RNA alone
    does better, and those are drawn in ink while the rest recede, so the eye lands on the crossing
    without a legend saying which to look at.

    Names sit at the RIGHT end only. Labelling both ends put the niche names where the y-axis label
    already was and the two ran through each other; one end also halves the number of collisions to
    repel. The claim belongs in the legend, not printed under the panel, which is what made the
    first version taller than it was wide with the text wider than the figure.
    """
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    ax = fig.add_axes([0.155, 0.13, 0.40, 0.84])
    x = (0.0, 1.0)
    rights = _spread(list(frame[RIGHT]), 0.655, 1.005)
    for row, right_y in zip(frame.itertuples(), rights):
        rna, allev = getattr(row, LEFT), getattr(row, RIGHT)
        rna_wins = rna >= allev
        colour = WINS if rna_wins else MUTED
        ax.plot(x, (rna, allev), "-", color=colour, linewidth=1.1 if rna_wins else 0.7,
                solid_capstyle="round", zorder=3 if rna_wins else 2)
        ax.plot(x[0], rna, "o", color=colour, markersize=2.8, markeredgewidth=0, zorder=4)
        ax.plot(x[1], allev, "o", color=colour, markersize=2.8, markeredgewidth=0, zorder=4)
        # a hairline from the point to its label where the repeller had to move it
        ax.plot([1.0, 1.085], [allev, right_y], "-", color=colour, linewidth=0.35,
                zorder=1, clip_on=False, alpha=0.75)
        ax.annotate(f"{rna:.2f}  {row.niche}  {allev:.2f}", (1.10, right_y),
                    xycoords=("axes fraction", "data"), ha="left", va="center",
                    fontsize=5.5, color=INK if rna_wins else MUTED,
                    annotation_clip=False)

    ax.set_xlim(-0.03, 1.03)
    ax.set_xticks(x)
    ax.set_xticklabels(["RNA\nco-expression", "all evidence\ntogether"], fontsize=5.8)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_ylabel("compartment recovered (cross-validated AUC)", fontsize=5.8, labelpad=2)
    ax.set_ylim(0.645, 1.015)
    ax.set_yticks([0.7, 0.8, 0.9, 1.0])
    ax.tick_params(axis="y", labelsize=5.5, pad=1.5)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(INK)
    fig.text(0.012, 0.98, letter, fontsize=7.5, fontweight="bold", color=INK,
             ha="left", va="top")
    return fig


def legend_text(frame: pd.DataFrame) -> str:
    """The claim, for the figure legend rather than printed on the panel."""
    won = frame[frame[LEFT] >= frame[RIGHT]]
    return (f"Each line is one of the {len(frame)} coarse hyperLOPIT compartments, recovered from "
            f"RNA co-expression alone (left) and from all allowed evidence together (right), as "
            f"cross-validated AUC. {len(won)} of {len(frame)} are placed as well or better by "
            f"transcription alone ({', '.join(won.niche)}); each is built in a wave during "
            f"division, so its genes share a schedule. The compartments transcription cannot place "
            f"are the ones made continuously: plasma membrane {frame.set_index('niche').loc['plasma membrane', LEFT]:.2f} "
            f"and ER {frame.set_index('niche').loc['ER', LEFT]:.2f}.")


def main() -> None:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    _font()
    frame = load(root)
    fig = draw(frame)
    out = os.path.join(root, "results", "coexpression_components_localization_2026_11_06")
    for extension in ("pdf", "png"):
        path = os.path.join(out, f"panel_which_evidence_places_a_protein.{extension}")
        fig.savefig(path, dpi=600 if extension == "png" else None,
                    bbox_inches="tight", pad_inches=0.02)
        print("wrote", path)
    print()
    print(legend_text(frame))
    print(f"median RNA alone {frame[LEFT].median():.3f}, all evidence {frame[RIGHT].median():.3f}")


if __name__ == "__main__":
    main()
