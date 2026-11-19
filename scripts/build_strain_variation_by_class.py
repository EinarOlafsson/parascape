"""What differs between Toxoplasma strains, and whether it is the secreted machinery.

Figure 8 of the ParaScape paper (instruction 101). The received account is that the proteins the
parasite puts into the host are the ones that differ between strains, because those are what
immunity sees and what virulence depends on. It is a testable claim and ParaScape can test it on
every gene at once.

THE MEASUREMENT. For each gene, the ratio of nonsynonymous to synonymous SNPs across the sequenced
strains: a gene changing its protein faster than its silent sites is under diversifying selection.
Classes of gene are then compared with that ratio.

THE NULL IS MATCHED ON LENGTH, AND HAS TO BE. A longer gene carries more SNPs of both kinds, and
the secreted proteins are not a random draw from the length distribution, so a comparison against
all other genes measures the length difference as much as the selection. Each class is compared
with random gene sets matched to it gene by gene on coding length, and the baseline reported is
what the same procedure gives on those.

A ratio is not dN/dS. Counts of nonsynonymous and synonymous SNPs are not normalised by the number
of sites at which each could occur, so the value is an index for comparing genes with each other
and not an estimate of selection pressure. The figure says so rather than implying otherwise.

Writes `results/strain_variation_by_class_<date>/`.
"""
from __future__ import annotations

import json
import os
from datetime import date

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: A class with fewer genes than this is not tested.
MIN_GENES = 25

#: Matched random sets behind each class.
N_NULL = 2000

#: The window a matched gene is drawn from, as a fraction of all genes, around the member it
#: stands in for in the length ordering.
MATCH_WINDOW = 0.10

#: A gene with fewer synonymous SNPs than this has too unstable a ratio to carry.
MIN_SYNONYMOUS = 1


def classes(frame: pd.DataFrame) -> dict:
    """The gene classes compared, each as a boolean over the frame."""
    compartment = frame["compartment_best"].astype(str)
    # `is_exported` is a Plasmodium target and is NOT a column of the Toxoplasma node table, so the
    # headline class of the first draft was silently all-False and dropped for being too small.
    # The Toxoplasma equivalents that DO exist are used instead and named for what they are: a
    # signal peptide, and the dense granules as hyperLOPIT assigns them.
    out = {
        "carries a signal peptide": frame["has_signal_peptide"].fillna(False).astype(bool)
        if "has_signal_peptide" in frame.columns else pd.Series(False, index=frame.index),
        "dense granules (hyperLOPIT)": frame["lopit_map"].astype(str) == "dense granules",
        "apical secretory (rhoptries, micronemes)": compartment == "apical_secretory",
        "surface and pellicle": compartment == "pellicle_PM_IMC",
        "nucleus": compartment == "nucleus",
        "mitochondrion": compartment == "mitochondrion",
        "apicoplast": compartment == "apicoplast",
        "cytosol": compartment == "cytosol",
        "required in vitro (lowest fitness tenth)":
            frame["fit_invitro_hff"] <= frame["fit_invitro_hff"].quantile(0.10)
            if "fit_invitro_hff" in frame.columns else pd.Series(False, index=frame.index),
    }
    kept = {k: v.fillna(False) for k, v in out.items() if v.fillna(False).sum() >= MIN_GENES}
    dropped = [k for k in out if k not in kept]
    if dropped:
        print("classes too small to test, named rather than omitted: " + ", ".join(dropped))
    return kept


def matched_null(rng, members, order, position, n_null: int) -> np.ndarray:
    """Random gene sets matched gene by gene on where each member sits in the length ordering."""
    n = len(order)
    half = max(1, int(round(n * MATCH_WINDOW / 2)))
    centres = np.array([position[g] for g in members], dtype=int)
    picks = np.empty((n_null, len(centres)), dtype=int)
    for i in range(n_null):
        for j, centre in enumerate(centres):
            lo, hi = max(0, centre - half), min(n, centre + half + 1)
            picks[i, j] = rng.integers(lo, hi)
    return picks


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    nodes = pd.read_parquet(os.path.join(cache, "nodes.parquet"))
    nodes = nodes.set_index(nodes["gene_id"].astype(str))

    usable = nodes[nodes["snp_synonymous"] >= MIN_SYNONYMOUS].copy()
    usable["ratio"] = usable["snp_nonsynonymous"] / usable["snp_synonymous"]
    print(f"{len(usable):,} of {len(nodes):,} genes carry at least {MIN_SYNONYMOUS} synonymous SNP "
          f"and can carry a ratio")

    order = usable["length"].sort_values().index.to_numpy()
    position = {g: i for i, g in enumerate(order)}
    ratio_by_position = usable["ratio"].reindex(order).to_numpy(dtype=float)
    rng = np.random.default_rng(0)

    rows, per_gene = [], []
    for name, mask in classes(usable).items():
        members = list(usable.index[mask])
        observed = float(usable.loc[members, "ratio"].mean())
        # The PER-GENE ratios, kept so the figure can show the distribution the class mean hides.
        # A class of 1,253 genes summarised as one bar is the thing the house rules refuse: more
        # than twenty values per group is a violin. A gene in two classes appears in both.
        per_gene.extend({"class": name, "gene_id": g, "ratio": float(v),
                         "length": float(usable.at[g, "length"])}
                        for g, v in usable.loc[members, "ratio"].items())
        picks = matched_null(rng, members, order, position, N_NULL)
        draws = np.nanmean(ratio_by_position[picks], axis=1)
        null_mean, null_sd = float(draws.mean()), float(draws.std(ddof=1))
        p = (1.0 + int((np.abs(draws - null_mean) >= abs(observed - null_mean)).sum())) \
            / (len(draws) + 1.0)
        rows.append({"class": name, "genes": len(members), "ratio": observed,
                     "length_matched_baseline": null_mean,
                     "detection_limit": 1.96 * null_sd,
                     "over_detection_limit": abs(observed - null_mean) / (1.96 * null_sd)
                     if null_sd else np.nan,
                     "p": p,
                     "median_length": float(usable.loc[members, "length"].median())})
        print(f"  {name:42s} n={len(members):4d}  ratio={observed:5.2f}  "
              f"length-matched baseline={null_mean:5.2f}  p={p:.4f}")

    table = pd.DataFrame(rows).sort_values("ratio", ascending=False)
    table["q"] = (table["p"] * len(table)).clip(upper=1.0)
    top = table.iloc[0]
    print(f"\nMost variable class: {top['class']} at {top['ratio']:.2f} against a length-matched "
          f"baseline of\n{top['length_matched_baseline']:.2f} (q = {top['q']:.4f}). "
          f"Median coding length {top['median_length']:.0f}.")
    print("\nThe ratio is an index, not dN/dS: the counts are not normalised by the number of "
          "sites\nat which each kind of change could occur.")

    out = os.path.join(ROOT, "results", f"strain_variation_by_class_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "by_class.tsv"), sep="\t", index=False)
    pd.DataFrame(per_gene).to_csv(os.path.join(out, "per_gene.tsv"), sep="\t", index=False)
    json.dump({"genes_with_a_ratio": int(len(usable)), "n_null": N_NULL,
               "match": "coding length", "classes": table.to_dict("records")},
              open(os.path.join(out, "summary.json"), "w"), indent=2, default=float)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
