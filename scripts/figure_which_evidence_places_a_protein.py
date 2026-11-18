"""A paper panel: which evidence places a protein in its compartment.

For the ParaScape paper, figure 5 -- which measurement explains which other measurement. This panel
is the localisation instance of it: the 13 coarse hyperLOPIT niches recovered from RNA co-expression
alone, against the same recovery from all allowed evidence together.

WHY A PARITY SCATTER. The quantity is a PAIRED comparison and the claim is which compartments do as
well on transcription alone, so a diagonal of parity reads the answer directly. Three forms were
drawn before this one was chosen; `draw` records why the other two lose.

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
from matplotlib import text as mtext

#: Where the scored table is, relative to the repository root.
SOURCE = ("results/coexpression_components_localization_2026_11_06/scores.tsv")

#: Millimetres, because a paper panel is specified in the units the page is in.
WIDTH_MM, HEIGHT_MM = 88.0, 62.0

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


#: Where a label may sit relative to its point, in points, tried in this order.
CANDIDATES = ((4.5, 0.0, "left", "center"), (-4.5, 0.0, "right", "center"),
              (4.5, 5.0, "left", "bottom"), (-4.5, 5.0, "right", "bottom"),
              (4.5, -5.5, "left", "top"), (-4.5, -5.5, "right", "top"),
              (0.0, 6.5, "center", "bottom"), (0.0, -7.0, "center", "top"))


def _box(text, renderer):
    """The bounding box of a label, measured on the TEXT and nothing else.

    `Annotation.get_window_extent` includes the leader line, which measures a box several times
    too large and makes every candidate position look occupied. This is the trap recorded against
    figure 1 of the bradyzoite paper, and it is why this goes through `matplotlib.text.Text`.
    """
    return mtext.Text.get_window_extent(text, renderer=renderer)


def _place(ax, fig, points, labels, colours):
    """Label every point, each in the first candidate position that hits nothing already placed.

    Guessed offsets collided -- apicoplast over cytosol, micronemes over apical -- because an
    offset that works for one arrangement of points works for no other. Measuring is cheap at
    thirteen labels and is the only version that stays correct when the numbers change.
    """
    renderer = fig.canvas.get_renderer()
    taken = []
    for (x, y) in points:
        taken.append(ax.transData.transform((x, y)))
    placed_boxes = []
    for (x, y), name, colour in zip(points, labels, colours):
        chosen = None
        for dx, dy, ha, va in CANDIDATES:
            text = ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                               fontsize=5.3, color=colour, ha=ha, va=va, annotation_clip=False)
            fig.canvas.draw()
            box = _box(text, renderer).expanded(1.06, 1.16)
            # Inside the page first: the left-hand candidate put "plasma membrane" half off the
            # figure, which no overlap test would ever have caught.
            width, height = fig.get_size_inches() * fig.dpi
            clash = box.x0 < 2 or box.y0 < 2 or box.x1 > width - 2 or box.y1 > height - 2
            clash = clash or any(box.overlaps(other) for other in placed_boxes)
            if not clash:
                clash = any(box.x0 - 1 <= px <= box.x1 + 1 and box.y0 - 1 <= py <= box.y1 + 1
                            for px, py in taken)
            if not clash:
                chosen, placed_boxes = text, placed_boxes + [box]
                break
            text.remove()
        if chosen is None:                       # nowhere clear: put it right and accept it
            placed_boxes.append(_box(
                ax.annotate(name, (x, y), xytext=(4.5, 0.0), textcoords="offset points",
                            fontsize=5.3, color=colour, ha="left", va="center"), renderer))


def draw(frame: pd.DataFrame, letter: str = "a"):
    """A parity scatter: each compartment placed from RNA alone against from all evidence.

    WHY THIS FORM, after three others were drawn. The quantity is a PAIRED comparison and the claim
    is which compartments do as well on transcription alone, so a diagonal of parity reads the
    answer directly -- on or below the line means RNA alone is enough, and the distance from the
    line is the gain from everything else. A heatmap of the same numbers has to print the values
    into the cells, because colour cannot separate 0.94 from 0.93, and once they are printed the
    colour is decoration on a table in which the crossings are invisible. A slope graph keeps the
    pairing but is a form the reader has to be taught. A ranked difference states the claim most
    directly and throws away the absolute performance, so it cannot distinguish a compartment that
    gains nothing because transcription is enough from one that gains nothing because both fail.
    """
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    ax = fig.add_axes([0.145, 0.165, 0.825, 0.80])
    lo, hi = 0.645, 1.03
    ax.plot([lo, hi], [lo, hi], "-", color="#c9c9c9", linewidth=0.7, zorder=1)
    wins = (frame[LEFT] >= frame[RIGHT]).to_numpy()
    for (x, y), won in zip(zip(frame[LEFT], frame[RIGHT]), wins):
        ax.plot(x, y, "o", color=WINS if won else MUTED, markersize=3.6,
                markeredgewidth=0, zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_ylim(0.775, 1.025)
    ax.set_xlabel("placed from RNA co-expression alone (AUC)", fontsize=6.0, labelpad=2)
    ax.set_ylabel("placed from all evidence (AUC)", fontsize=6.0, labelpad=2)
    ax.tick_params(labelsize=5.5, pad=1.5)
    ax.set_xticks([0.7, 0.8, 0.9, 1.0])
    ax.set_yticks([0.8, 0.9, 1.0])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    # Horizontal, in the empty corner. Rotated along the diagonal it ran through the data, which
    # is the one place on this panel a note must never be.
    ax.annotate("on or below the line:\ntranscription alone is enough", (0.995, 0.787),
                fontsize=5.1, color="#7d7d7d", ha="right", va="bottom", linespacing=1.35)
    _place(ax, fig, list(zip(frame[LEFT], frame[RIGHT])), list(frame.niche),
           [INK if w else MUTED for w in wins])
    fig.text(0.013, 0.985, letter, fontsize=7.5, fontweight="bold", color=INK,
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
