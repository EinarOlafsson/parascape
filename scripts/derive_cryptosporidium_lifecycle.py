"""Cryptosporidium transcription across the life cycle, from the single-cell series GSE232438.

Instruction 41. The v4 candidate list is 11.1 per cent on target and this is the one genuine match
in it: "Transcriptional control of the Cryptosporidium lifecycle", offered against
`Cp::transcription · oocyst`.

WHAT IT IS. Single-cell matrices at five points of in vitro culture -- 24, 36, 42, 46 and 54 hours
-- and one in vivo sample, keyed by CryptoDB identifiers. 10,319 cells over 4,020 features, and
every one of the 3,923 genes of the Cryptosporidium space is in it.

WHAT THIS DERIVES, AND WHAT IT CANNOT. A pseudobulk per sample: the mean expression of each gene
over the cells of that timepoint, as log2 CPM. That is transcription at a POINT IN CULTURE, which
is what the series measures and what can be checked. It is NOT transcription in the oocyst: the
GEO supplement ships barcodes and matrices with no stage or cluster annotation, so no cell here can
be called an oocyst without redoing the paper's clustering and stage assignment. The slot it was
offered for therefore stays empty, and the reason is recorded rather than papered over by treating
the latest timepoint as if it were a stage.

Writes `results/cryptosporidium_lifecycle_<date>/`.
"""
from __future__ import annotations

import gzip
import json
import os
from datetime import date

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = "/media/carruthers/mnt1/parascape_downloads/GSE232438"
SAMPLES = ("24hrs", "36hrs", "42hrs", "46hrs", "54hrs", "in_vivo")


def features() -> list:
    path = os.path.join(DATA, "GSE232438_all_features.tsv.gz")
    return [line.split("\t")[0] for line in gzip.open(path, "rt").read().splitlines()]


def pseudobulk(sample: str, genes: list) -> pd.Series:
    """Mean expression over the cells of one sample, as log2 CPM.

    Summed then re-closed per million over the sample, so a timepoint with more cells or deeper
    cells is not simply larger. Pseudobulk and not a per-cell mean, because a per-cell mean lets a
    few deeply sequenced cells speak for the sample.
    """
    base = os.path.join(DATA, f"GSE232438_{sample}")
    triplets = pd.read_csv(base + "_matrix.mtx.gz", sep=r"\s+", comment="%", skiprows=1,
                           header=None, names=["i", "j", "v"],
                           dtype={"i": np.int32, "j": np.int32, "v": np.float32})
    n_cells = int(pd.read_csv(base + "_barcodes.tsv.gz", header=None).shape[0])
    m = sp.csr_matrix((triplets["v"], (triplets["i"] - 1, triplets["j"] - 1)),
                      shape=(len(genes), n_cells))
    total = np.asarray(m.sum(axis=1)).ravel()
    cpm = total / total.sum() * 1e6 if total.sum() else total
    return pd.Series(np.log2(cpm + 1.0), index=genes, name=sample)


def main() -> None:
    genes = features()
    cp = pd.read_parquet(os.path.join(ROOT, "results", "cryptosporidium_space_2026_11_08",
                                      "cp_gene_nodes.parquet"))
    space = set(cp["gene_id"].astype(str))
    print(f"{len(genes):,} features in the series; {len(space):,} genes in the Cp space")

    columns = {}
    for sample in SAMPLES:
        s = pseudobulk(sample, genes)
        columns[f"cp_lifecycle_{sample}"] = s
        print(f"  {sample:9s} {int((s > 0).sum()):,} genes with signal")
    frame = pd.DataFrame(columns)
    frame.index.name = "gene_id"
    kept = frame.loc[frame.index.isin(space)]
    print(f"\n{len(kept):,} of {len(frame):,} rows are genes of the Cp space "
          f"({len(kept)/len(space):.1%} of it)")

    # Does it recover the biology its own paper is about? Transcription should change most between
    # the earliest culture point and the in vivo sample, and least between adjacent timepoints.
    order = [f"cp_lifecycle_{s}" for s in SAMPLES]
    corr = kept[order].corr(method="spearman")
    adjacent = float(np.mean([corr.iloc[i, i + 1] for i in range(len(order) - 2)]))
    far = float(corr.loc["cp_lifecycle_24hrs", "cp_lifecycle_in_vivo"])
    print(f"\nadjacent in vitro timepoints agree at {adjacent:.3f}; "
          f"24 h against in vivo at {far:.3f}")
    if far < adjacent:
        print("  the time course behaves like one: neighbours resemble each other more than the "
              "ends do")

    out = os.path.join(ROOT, "results", f"cryptosporidium_lifecycle_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    kept.reset_index().to_csv(os.path.join(out, "deposit_cp_lifecycle.tsv"), sep="\t", index=False)
    corr.to_csv(os.path.join(out, "timepoint_agreement.tsv"), sep="\t")
    json.dump({"accession": "GSE232438", "samples": list(SAMPLES),
               "features_in_series": len(genes), "genes_in_space": len(space),
               "rows_kept": int(len(kept)),
               "adjacent_agreement": adjacent, "earliest_against_in_vivo": far,
               "fills": "transcription at a point in culture",
               "does_not_fill": "Cp::transcription · oocyst -- the GEO supplement carries no stage "
                                "or cluster annotation, so no cell here can be called an oocyst "
                                "without redoing the clustering and stage assignment of the paper"},
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
