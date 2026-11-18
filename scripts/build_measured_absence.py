"""What nothing in the literature measures, and which experiment would close the most of it.

Figure 12 of the ParaScape paper (instruction 101), and the close Einar asked for: *"explaining
which datasets would be most valuable to the field."* So this is prescriptive rather than a caveat
section.

THE ARGUMENT. One dataset that does not mention a gene has not shown the gene is absent -- it has
shown that experiment did not reach it. The union of 177 datasets is different: where no assay of a
KIND has ever carried a value for a gene, "we do not know" becomes a measured statement with a
denominator, and the size of that gap per kind of knowledge is what tells a laboratory where an
experiment is worth the most.

WHAT IS COUNTED. Each gene, for each category of knowledge (where the protein is, whether losing it
costs fitness, how much RNA there is, and so on), either carries a value in at least one column of
that category or does not. The result is per category: how many genes nothing has ever measured.
The categories are read from `parascape.information_categories`, which is the same table the
information map uses, so this cannot drift from what the program says it holds.

WHAT IT DOES NOT CLAIM. A gene with no measurement is not a gene with no function, and the figure
never says so. It says nobody has looked, which is the one thing an empty cell is allowed to mean.

Writes `results/measured_absence_<date>/`.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import date

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import information_categories as IC    # noqa: E402

#: Categories that say nothing about the biology of a gene, so a gap in them is not a gap worth
#: filling. `identity` is names, `attention` is how much somebody has already published.
NOT_KNOWLEDGE = ("identity", "attention")


def is_a_measurement(series: pd.Series) -> bool:
    """Whether a column records something that was measured, as opposed to computed for every gene.

    THE RULE THAT MAKES THIS FIGURE MEAN ANYTHING. A column with no missing value at all is not a
    measurement of 8,140 genes; it is a count, a flag or a derivation in which absence was written
    as zero or as False. `n_bioid_studies` is 0 for a gene nobody pulled down, `has_phospho` is
    False for a gene nobody assayed, and `compartment_best` carries a best guess for everything.
    Counting those as measurements said that every gene in the genome has had its localisation,
    its interactions and its modifications measured, which is the exact error -- absence written
    as zero -- that this figure exists to expose.
    """
    return bool(series.isna().any())


def category_columns(frame: pd.DataFrame) -> dict:
    """Category -> the columns that record a measurement, not a derivation or a count of zero."""
    derived = set(getattr(IC, "DERIVED_FROM", {}))
    out = {}
    for column in frame.columns:
        if column in derived or not is_a_measurement(frame[column]):
            continue
        for category in IC.categories_of(column):
            if category in NOT_KNOWLEDGE:
                continue
            out.setdefault(category, []).append(column)
    return {k: v for k, v in sorted(out.items())}


def coverage(frame: pd.DataFrame, columns: dict) -> pd.DataFrame:
    """Per category: how many genes carry at least one value, and how many nothing has measured."""
    total = len(frame)
    rows = []
    for category, cols in columns.items():
        present = frame[cols].notna().any(axis=1)
        rows.append({"category": category,
                     "description": IC.CATEGORIES.get(category, ""),
                     "columns": len(cols),
                     "genes_measured": int(present.sum()),
                     "genes_never_measured": int(total - present.sum()),
                     "share_measured": float(present.mean())})
    return pd.DataFrame(rows).sort_values("genes_never_measured", ascending=False)


def darkest_genes(frame: pd.DataFrame, columns: dict) -> pd.Series:
    """How many categories of knowledge each gene carries anything in."""
    counts = pd.Series(0, index=frame.index, dtype=int)
    for cols in columns.values():
        counts += frame[cols].notna().any(axis=1).astype(int)
    return counts


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    nodes = pd.read_parquet(os.path.join(cache, "nodes.parquet"))
    nodes = nodes.set_index(nodes["gene_id"].astype(str))
    columns = category_columns(nodes)
    dropped = len(nodes.columns) - sum(len(v) for v in columns.values())
    print(f"{len(nodes):,} genes, {len(columns)} categories of knowledge, "
          f"{sum(len(v) for v in columns.values()):,} measurement columns "
          f"({dropped:,} columns set aside as derivations, flags or counts of zero)")

    table = coverage(nodes, columns)
    print("\nwhat nothing has ever measured, per category:")
    for row in table.itertuples():
        print(f"  {row.category:18s} {row.genes_never_measured:5,} genes never measured "
              f"({1 - row.share_measured:5.1%})   {row.columns:3d} columns")

    counts = darkest_genes(nodes, columns)
    distribution = counts.value_counts().sort_index()
    dark = counts[counts <= 2]
    print(f"\n{len(dark):,} genes carry anything in two categories or fewer, of {len(columns)}")
    print(f"{int((counts == 0).sum()):,} carry nothing at all")

    # The prescriptive line: the category where one new assay would give the most genes their
    # first measurement of that kind. Ranked by genes, not by how interesting the assay sounds.
    best = table.iloc[0]
    print(f"\nThe kind of experiment that would give the most genes a first measurement is "
          f"{best['category']}:\n{best['genes_never_measured']:,} genes "
          f"({1 - best['share_measured']:.1%}) have never been measured that way.")

    out = os.path.join(ROOT, "results", f"measured_absence_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "by_category.tsv"), sep="\t", index=False)
    counts.rename("categories_with_a_value").to_csv(
        os.path.join(out, "categories_per_gene.tsv"), sep="\t")
    dark.rename("categories_with_a_value").to_csv(os.path.join(out, "darkest_genes.tsv"), sep="\t")
    json.dump({"genes": int(len(nodes)), "categories": len(columns),
               "columns": int(sum(len(v) for v in columns.values())),
               "genes_with_nothing": int((counts == 0).sum()),
               "genes_with_two_or_fewer": int(len(dark)),
               "distribution": {int(k): int(v) for k, v in distribution.items()},
               "by_category": table.to_dict("records"),
               "largest_gap": best["category"]},
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
