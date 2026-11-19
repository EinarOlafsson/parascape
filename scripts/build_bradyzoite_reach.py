"""Which bradyzoite subtype of a brain cyst does anything in the literature reach?

Figure 5 of the ParaScape paper (instruction 101). Five subtypes of in vivo bradyzoite (Ulu et al.,
Nat Commun 2026) are marked from the published block of average expression per subtype, and each is
tested against every within-study contrast ParaScape stores -- drugs, knockouts, pH switches, time
points -- by :mod:`parascape.cluster_axis`, which is the engine behind the Reach module.

THE QUESTION THIS ANSWERS. A compound that moves only the genes of one subtype will only clear that
subtype, and a tissue cyst is not one transcriptional state. So: is there a subtype that nothing
published reaches?

The three refusals of the engine apply throughout: the baseline is a random gene set matched to each
subtype in size and expression and never the number of genes in it; a subtype the contrast barely
measured reads unknown rather than unaffected; and every tested subtype carries the smallest effect
its own comparison could have resolved.

Writes `results/bradyzoite_reach_<date>/`.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import date

import matplotlib
matplotlib.use("Agg")
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import cluster_axis as CA        # noqa: E402
from parascape import figstyle as FS           # noqa: E402

#: The published block: one average expression per gene per subtype.
BLOCK = "deposit_brady_subtypes.tsv"
SUBTYPES = tuple("ABCDE")

#: Matched random sets behind each subtype, per contrast.
N_NULL = 1000



def sweep(cache: str, n_null: int = N_NULL) -> tuple:
    """Every stored contrast against all five subtypes. Returns the per-contrast table and markers."""
    block = pd.read_csv(os.path.join(cache, BLOCK), sep="\t")
    columns = {g: f"bzsub_{g}_expr" for g in SUBTYPES}
    markers = CA.markers_from_blocks(block, columns)
    expression = block[list(columns.values())].mean(axis=1)
    expression.index = block["gene_id"].astype(str)
    expression = expression.dropna()

    catalogue = CA.axis_catalogue()
    rows = []
    for i, column in enumerate(catalogue.index, 1):
        try:
            result = CA.test_clusters_against_axis(
                markers, CA.axis_values(column), expression, n_null=n_null, axis_name=str(column))
        except (ValueError, KeyError) as error:
            rows.append({"contrast": column, "study": catalogue.loc[column, "study"],
                         "reached": "", "n_reached": 0, "error": str(error)})
            continue
        rows.append({"contrast": column,
                     "study": catalogue.loc[column, "study"],
                     "label": str(catalogue.loc[column, "contrast_label"]),
                     "reached": ",".join(result.reached),
                     "n_reached": len(result.reached),
                     "n_unknown": len(result.unknown),
                     **{f"q_{s}": float(result.table.loc[s, "q"])
                        if s in result.table.index else float("nan") for s in SUBTYPES}})
        if i % 25 == 0:
            print(f"  {i}/{len(catalogue)} contrasts", flush=True)
    return pd.DataFrame(rows), markers


def counts(table: pd.DataFrame) -> pd.Series:
    """How many contrasts reach each subtype."""
    return pd.Series({s: int(table["reached"].fillna("").str.split(",").apply(
        lambda r, s=s: s in r).sum()) for s in SUBTYPES})


def draw(reached: pd.Series, total: int, markers: dict, letter: str = "a"):
    """How many of the stored contrasts reach each subtype, in the house style.

    A bar per subtype, because the claim is a comparison of five counts with a meaningful zero and
    a reader has to see that one bar is a twentieth of another. Five counts is five points, which
    is well under the twenty above which a bar is refused in favour of a violin.

    Geometry, palette and opacity come from :mod:`parascape.figstyle`, so this panel is the same
    width per category as every other panel of the paper and its bars are the same width and
    spacing. The subtype nothing reaches carries the first series colour; the rest are grey,
    because they are the reference against which it is remarkable.
    """
    FS.use()
    order = reached.sort_values(ascending=False)
    lowest = order.index[-1]
    fig, ax = FS.figure(len(order), vertical=True)
    names = [f"subtype {s}" for s in order.index]
    FS.bars(ax, order.to_numpy(), names, highlight={f"subtype {lowest}"}, rotate=0)
    for i, value in enumerate(order.to_numpy()):
        ax.annotate(f"{value}", (i, value), xytext=(0, 2), textcoords="offset points",
                    ha="center", va="bottom", fontsize=FS.NOTE_PT,
                    color=FS.SERIES[0] if order.index[i] == lowest else FS.INK)
    ax.set_ylabel(f"contrasts reaching it (of {total})")
    ax.set_ylim(0, float(order.max()) * 1.16)
    ax.set_xlim(-0.6, len(order) - 0.4)
    FS.letter(fig, letter)
    return fig, lowest


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    print("sweeping every stored contrast against the five subtypes...")
    table, markers = sweep(cache)
    reached = counts(table)
    total = len(table)
    any_reached = int((table["n_reached"] > 0).sum())
    print(f"\n{any_reached} of {total} contrasts reach at least one subtype")
    print(reached.to_string())

    out = os.path.join(ROOT, "results", f"bradyzoite_reach_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "contrasts.tsv"), sep="\t", index=False)
    fig, lowest = draw(reached, total, markers)
    for extension in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"panel_reach_by_subtype.{extension}"),
                    dpi=600 if extension == "png" else None, bbox_inches="tight", pad_inches=0.02)
    by_study = table[table["n_reached"] > 0]["study"].value_counts()
    json.dump({"contrasts": total, "reaching_any": any_reached,
               "by_subtype": {k: int(v) for k, v in reached.items()},
               "markers_per_subtype": {k: len(v) for k, v in markers.items()},
               "least_reached": lowest, "n_null": N_NULL,
               "top_studies": {str(k): int(v) for k, v in by_study.head(6).items()}},
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print(f"\nleast reached: subtype {lowest} at {reached[lowest]} of {total}")
    sizes = {len(v) for v in markers.values()}
    print(f"legend: each subtype is marked by "
          f"{'the same ' + str(sizes.pop()) if len(sizes) == 1 else 'its own number of'} genes, "
          f"most specific first; a contrast reaches a subtype when it moves those markers more "
          f"than expression-matched random genes, corrected across the five.")
    print("wrote", out)


if __name__ == "__main__":
    main()
