"""The co-expression modules of Toxoplasma, and what each one is.

Figure 3 of the ParaScape paper (instruction 101). Genes that rise and fall together across 57
independent studies form modules; the question the figure answers is what those modules ARE, which
is only answerable because the compartment, the fitness cost and the function of their members are
in the same table.

HOW A MODULE IS NAMED. For each module and each compartment, how many of its members are assigned
to that compartment, against what the SAME procedure gives on random gene sets matched to the
module in size and in how well measured its members are. Never against the prior size of the
compartment: a compartment holding a third of the genome is a third of any module by shape alone,
and dividing by that prior is the error this project refuses by standing instruction.

A module whose best enrichment does not clear its own matched baseline is reported as unnamed. That
is a result -- genes that move together without sharing a compartment -- and not a failure to try
harder.

Writes `results/coexpression_modules_<date>/`.
"""
from __future__ import annotations

import json
import os
from datetime import date

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

COMPONENTS = os.path.join(ROOT, "results", "coexpression_components_localization_2026_11_06",
                          "coexpression_components.parquet")

#: How many modules the genes are cut into. Stated, not tuned to make the answer better.
N_MODULES = 24

#: A module smaller than this is not named.
MIN_MODULE = 20

#: Matched random sets behind each module-compartment pair.
N_NULL = 1000

#: The window a matched gene is drawn from, as a fraction of all genes.
MATCH_WINDOW = 0.10


def modules(components: pd.DataFrame, n: int = N_MODULES) -> pd.Series:
    """Cut the genes into modules on their co-expression coordinates, by Ward linkage."""
    values = components.to_numpy(dtype=float)
    values = (values - values.mean(axis=0)) / (values.std(axis=0) + 1e-12)
    link = linkage(pdist(values, metric="euclidean"), method="ward")
    return pd.Series(fcluster(link, t=n, criterion="maxclust"), index=components.index,
                     name="module")


def name_modules(assignment: pd.Series, labels: pd.Series, measured: pd.Series,
                 seed: int = 0) -> pd.DataFrame:
    """For every module, the compartment it holds more of than matched random genes do.

    `measured` is what the null is matched on -- how many categories of knowledge each gene
    carries -- because a module of well-studied genes would otherwise look enriched for every
    compartment that is easy to assign.
    """
    rng = np.random.default_rng(seed)
    order = measured.sort_values().index.to_numpy()
    position = {g: i for i, g in enumerate(order)}
    labels_by_position = labels.reindex(order).to_numpy()
    half = max(1, int(round(len(order) * MATCH_WINDOW / 2)))

    rows = []
    for module, members in assignment.groupby(assignment):
        genes = [g for g in members.index if g in position]
        if len(genes) < MIN_MODULE:
            continue
        centres = np.array([position[g] for g in genes], dtype=int)
        picks = np.empty((N_NULL, len(centres)), dtype=int)
        for i in range(N_NULL):
            for j, centre in enumerate(centres):
                lo, hi = max(0, centre - half), min(len(order), centre + half + 1)
                picks[i, j] = rng.integers(lo, hi)
        drawn = labels_by_position[picks]
        observed = labels.reindex(genes)
        best = None
        for compartment in observed.dropna().unique():
            if compartment in ("unassigned", "nan"):
                continue
            share = float((observed == compartment).mean())
            null = (drawn == compartment).mean(axis=1)
            null_mean, null_sd = float(null.mean()), float(null.std(ddof=1))
            if null_sd <= 0:
                continue
            p = (1.0 + int((null >= share).sum())) / (len(null) + 1.0)
            over = (share - null_mean) / (1.96 * null_sd)
            if best is None or over > best["over_detection_limit"]:
                best = {"compartment": compartment, "share": share,
                        "matched_baseline": null_mean,
                        "over_detection_limit": over, "p": p}
        row = {"module": int(module), "genes": len(genes)}
        row.update(best or {"compartment": "", "share": np.nan, "matched_baseline": np.nan,
                            "over_detection_limit": np.nan, "p": np.nan})
        rows.append(row)
    table = pd.DataFrame(rows)
    table["q"] = (table["p"] * len(table)).clip(upper=1.0)
    table["named"] = (table["q"] < 0.05) & (table["over_detection_limit"] >= 1.0)
    return table.sort_values("over_detection_limit", ascending=False)


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    components = pd.read_parquet(COMPONENTS)
    nodes = pd.read_parquet(os.path.join(cache, "nodes.parquet"))
    nodes = nodes.set_index(nodes["gene_id"].astype(str))
    shared = components.index.intersection(nodes.index)
    components, nodes = components.loc[shared], nodes.loc[shared]
    # A gene with no co-expression coordinate was in too few studies to place. Dropped and
    # counted, not filled with zero: a zero on a centred component is the middle of the map,
    # which would put every unmeasured gene in one false module.
    placed = components.notna().all(axis=1)
    print(f"{len(components):,} genes carry components; {int((~placed).sum()):,} have none and "
          f"are set aside")
    components, nodes = components[placed], nodes[placed]
    print(f"{len(components):,} genes on {components.shape[1]} co-expression components")

    assignment = modules(components)
    sizes = assignment.value_counts()
    print(f"{assignment.nunique()} modules, sizes {sizes.min()} to {sizes.max()} "
          f"(median {int(sizes.median())})")

    labels = nodes["compartment_best"].astype(str)
    # The null is matched on how much is known about a gene, not on nothing: a module of
    # well-studied genes is enriched for every compartment that is easy to assign.
    measured = nodes.notna().sum(axis=1).astype(float)
    table = name_modules(assignment, labels, measured)

    named = table[table["named"]]
    print(f"\n{len(named)} of {len(table)} modules carry a compartment above their own "
          f"matched baseline:\n")
    for row in named.itertuples():
        print(f"  module {row.module:3d}  n={row.genes:5d}  {row.compartment:26s} "
              f"{row.share:5.1%} against {row.matched_baseline:5.1%} matched  "
              f"({row.over_detection_limit:.1f}x the smallest visible difference)")
    unnamed = table[~table["named"]]
    if len(unnamed):
        print(f"\n{len(unnamed)} modules hold no compartment above their matched baseline. "
              f"Genes that\nmove together without sharing a compartment is a result, not a "
              f"failure to look harder.")

    out = os.path.join(ROOT, "results", f"coexpression_modules_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "modules.tsv"), sep="\t", index=False)
    assignment.to_csv(os.path.join(out, "gene_modules.tsv"), sep="\t")
    json.dump({"genes": int(len(components)), "modules": int(assignment.nunique()),
               "named": int(len(named)), "n_null": N_NULL,
               "matched_on": "how many categories of knowledge the gene carries",
               "table": table.to_dict("records")},
              open(os.path.join(out, "summary.json"), "w"), indent=2, default=float)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
