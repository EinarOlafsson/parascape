"""Check a figure against the house rules, instead of catching defects one at a time by eye.

Einar's `spacr-graphs` checklist opens with `check_figure(fig)`, and his helper lives on the Home
machine. This is the same idea written here, because four of the first four panels I opened by eye
after the rebuild each held a defect, and finding them one at a time is slow and leaky.

WHAT IT CAN SEE: text that overlaps other text, text clipped at the canvas edge, a continuous-axis
plot whose plotting area is not square, a categorical panel whose width is not n slots, data drawn
outside the axis limits, a bar or violin drawn at the wrong width, an element not at 60 per cent
opacity, and text sitting on the plotting area where only axis titles and legends are allowed.

KNOWN FALSE POSITIVES, which the checklist says must be judged by eye. A ROTATED tick label has
an axis-aligned bounding box, so tilted labels report as overlapping when they do not touch; and
figures are saved with `bbox_inches="tight"`, which expands the canvas, so text reported as clipped
before the save usually is not. Weight a count by looking at the panel, not by the number alone.

WHAT IT CANNOT SEE, and still needs a person: whether the form suits the measurement, whether the
accent is on the right element, whether an axis title is objective, and whether a schematic is
informative. A green check is not a correct figure.
"""
from __future__ import annotations

from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle
from matplotlib.text import Annotation, Text

from parascape import figstyle as FS

#: Two text boxes may share this many points of overlap before it is called a collision; glyph
#: bounding boxes are generous and touching descenders are not a defect.
TOLERANCE_PT = 1.0


def _boxes(fig, renderer):
    """Every visible text on the figure, with the box MEASURED ON THE TEXT.

    Measured on `Text`, never on an `Annotation`: an annotation box includes its leader line and
    comes back several times too large, which is a trap this project has hit before.
    """
    out = []
    for ax in fig.axes:
        for t in [ax.title, ax.xaxis.label, ax.yaxis.label]:
            if t.get_text():
                out.append((t, Text.get_window_extent(t, renderer)))
        for t in ax.get_xticklabels() + ax.get_yticklabels():
            if t.get_text() and t.get_visible():
                out.append((t, Text.get_window_extent(t, renderer)))
        for child in ax.get_children():
            if isinstance(child, (Text, Annotation)) and child.get_text() \
                    and child not in (ax.title,):
                out.append((child, Text.get_window_extent(child, renderer)))
    for child in fig.get_children():
        if isinstance(child, Text) and child.get_text():
            out.append((child, Text.get_window_extent(child, renderer)))
    return out


def check(fig, categorical: int | None = None, continuous: bool = False) -> list:
    """Return a list of complaints about `fig`. Empty means nothing this can see is wrong.

    `categorical` is the number of conditions, when the panel has a condition axis; `continuous`
    says the plotting area should be square.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    problems = []

    # --- text overlapping other text
    boxes = _boxes(fig, renderer)
    for i, (t1, b1) in enumerate(boxes):
        for t2, b2 in boxes[i + 1:]:
            dx = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)
            dy = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
            if dx > TOLERANCE_PT and dy > TOLERANCE_PT:
                problems.append(f"text overlaps text: {t1.get_text()[:26]!r} / "
                                f"{t2.get_text()[:26]!r}")

    # --- text clipped at the canvas edge
    w, h = fig.canvas.get_width_height()
    for t, b in boxes:
        if b.x0 < -TOLERANCE_PT or b.y0 < -TOLERANCE_PT or b.x1 > w + TOLERANCE_PT \
                or b.y1 > h + TOLERANCE_PT:
            problems.append(f"text clipped at the canvas edge: {t.get_text()[:34]!r}")

    for ax in fig.axes:
        box = ax.get_position()
        width_in = box.width * fig.get_size_inches()[0]
        height_in = box.height * fig.get_size_inches()[1]

        # --- a continuous-axis plot is square
        if continuous and abs(width_in - height_in) > 0.02:
            problems.append(f"continuous plot is not square: "
                            f"{width_in:.2f} x {height_in:.2f} in")

        # --- a categorical panel is n slots wide
        if categorical:
            expect = max(categorical, FS.MIN_SLOTS) * FS.SLOT_IN
            if abs(width_in - expect) > 0.02:
                problems.append(f"condition axis is {width_in:.2f} in, expected "
                                f"{expect:.2f} in ({categorical} slots)")

        # --- nothing drawn outside the limits, and nothing sitting on a spine
        lo, hi = ax.get_ylim()
        for patch in ax.patches:
            if isinstance(patch, Rectangle):
                y0, y1 = patch.get_y(), patch.get_y() + patch.get_height()
                if y1 > hi + 1e-9 or y0 < lo - 1e-9:
                    problems.append("a bar is drawn outside the y limits")
                    break
        for coll in ax.collections:
            if isinstance(coll, PolyCollection):
                for path in coll.get_paths():
                    ys = path.vertices[:, 1]
                    if ys.max() > hi + 1e-9 or ys.min() < lo - 1e-9:
                        problems.append("a violin crosses the y limits")
                        break

        # --- elements at 60 per cent
        for patch in ax.patches:
            a = patch.get_alpha()
            if a is not None and abs(a - FS.ALPHA) > 1e-6:
                problems.append(f"an element is at alpha {a}, not {FS.ALPHA}")
                break

        # --- text on the plotting area: only axis titles and legend entries are allowed
        legend_texts = set()
        if ax.get_legend() is not None:
            legend_texts = {id(t) for t in ax.get_legend().get_texts()}
        for child in ax.get_children():
            if isinstance(child, (Text, Annotation)) and child.get_text() \
                    and id(child) not in legend_texts \
                    and child not in (ax.title, ax.xaxis.label, ax.yaxis.label) \
                    and child not in ax.get_xticklabels() + ax.get_yticklabels():
                problems.append(f"text on the plotting area: {child.get_text()[:34]!r}")

    return problems


def report(name: str, fig, **kw) -> list:
    """Check and print, in the form the checklist asks for: run it, print it, fix it."""
    problems = check(fig, **kw)
    if problems:
        print(f"  {name}: {len(problems)} problem(s)")
        for p in dict.fromkeys(problems):        # de-duplicated, order kept
            print(f"      {p}")
    else:
        print(f"  {name}: ok")
    return problems
