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
import sys

import matplotlib
matplotlib.use("Agg")
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from parascape import figstyle as FS   # noqa: E402

#: Where the scored table is, relative to the repository root.
SOURCE = ("results/coexpression_components_localization_2026_11_06/scores.tsv")

#: Size comes from `figstyle.figure` now, so this panel and figure 2 -- the other parity scatter --
#: are the same size by construction rather than by two numbers that happened to agree.

#: The two columns compared, in the order they are drawn.
LEFT, RIGHT = "coexp_cv_auc", "all_evidence_cv_auc"

#: Ink and series come from the house style (`parascape.figstyle`), so this panel matches every
#: other panel of the paper. The compartments all-evidence wins are the REFERENCE here and are
#: therefore grey; the crossings carry the first series colour.
INK, MUTED, WINS = FS.INK, FS.CONTROL, FS.SERIES[0]


def _font() -> None:
    """The house style: palette, thin axes, no gridlines, Open Sans where it is installed."""
    FS.use()


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
    # The house figure, not a hand-rolled one. This panel had its own figsize and its own axes
    # fractions, so it did not match figure 2 -- the OTHER parity scatter of the same paper -- and
    # rule 4 says two panels of the same shape are the same size. `figure()` with the same
    # arguments figure 2 uses makes them identical.
    # Wider than the default, for the same reason figure 7 is: the compartment names are long and
    # "plasma membrane" sits at the far left of the data, so at the default width the placer had
    # nowhere to put it but off the page. A value axis may grow; a category slot may not.
    fig, ax = FS.figure(6, height_mm=58.0, value_mm=92.0)
    # `lo` is the LEFT LIMIT and `hi` the right; the parity line is drawn across the same span.
    # The limit was 0.645 with a point at 0.67, so "plasma membrane" had nowhere to go but off the
    # page and through its own marker.
    lo, hi = 0.62, 1.03
    ax.plot([lo, hi], [lo, hi], "-", color="#c9c9c9", linewidth=0.7, zorder=1)
    wins = (frame[LEFT] >= frame[RIGHT]).to_numpy()
    for (x, y), won in zip(zip(frame[LEFT], frame[RIGHT]), wins):
        # 60 per cent, like every other drawn element in the paper. These markers were opaque,
        # which is why this panel read brighter than figure 2 beside it.
        ax.plot(x, y, "o", color=WINS if won else MUTED, markersize=3.6,
                markeredgewidth=0, alpha=FS.ALPHA, zorder=3)
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
    # The house label placer and the house panel letter, so a label here is measured and rejected
    # the same way as on every other panel. The two local helpers this panel carried were a
    # second copy of `figstyle.place_labels`, and a second copy is a second set of bugs.
    FS.place_labels(ax, fig, list(zip(frame[LEFT], frame[RIGHT])), list(frame.niche),
                    [INK if w else MUTED for w in wins])
    FS.letter(fig, letter)
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
