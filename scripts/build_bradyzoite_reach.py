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
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import cluster_axis as CA        # noqa: E402

#: The published block: one average expression per gene per subtype.
BLOCK = "deposit_brady_subtypes.tsv"
SUBTYPES = tuple("ABCDE")

#: Matched random sets behind each subtype, per contrast.
N_NULL = 1000

WIDTH_MM, HEIGHT_MM = 88.0, 52.0
INK, MUTED, MARK = "#1a1a1a", "#b4b4b4", "#b3452c"


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
    """How many of the stored contrasts reach each subtype.

    A bar per subtype, because the claim is a comparison of five counts and nothing else: the
    quantity has a meaningful zero, the categories are few and unordered, and a reader needs to
    see that one bar is a twentieth of another. The subtype nothing reaches carries the only
    colour, so the figure says what it is about without a legend.
    """
    if "Open Sans" in {f.name for f in font_manager.fontManager.ttflist}:
        plt.rcParams["font.family"] = "Open Sans"
    plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 6.2,
                         "axes.linewidth": 0.5, "xtick.major.width": 0.5})
    order = reached.sort_values(ascending=True)
    lowest = order.index[0]
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    ax = fig.add_axes([0.125, 0.21, 0.825, 0.70])
    y = range(len(order))
    ax.barh(list(y), order.to_numpy(),
            color=[MARK if s == lowest else MUTED for s in order.index],
            height=0.62, linewidth=0)
    for i, (name, value) in enumerate(order.items()):
        ax.annotate(f"{value}", (value, i), xytext=(3, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=5.8,
                    color=MARK if name == lowest else INK)
    ax.set_yticks(list(y))
    # Just the name. Every subtype keeps the same number of markers, so printing it five times
    # says nothing and costs width; it goes in the legend instead.
    ax.set_yticklabels([f"subtype {s}" for s in order.index], fontsize=6.0)
    ax.set_xlabel(f"published contrasts that reach the subtype (of {total})", fontsize=6.0,
                  labelpad=2)
    ax.set_xlim(0, float(order.max()) * 1.14)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", labelsize=5.5, pad=1.5)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    fig.text(0.013, 0.98, letter, fontsize=7.5, fontweight="bold", color=INK, ha="left", va="top")
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
