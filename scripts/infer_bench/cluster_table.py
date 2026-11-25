"""One table of every E8 cluster: label, feature set, resolution, cluster, its top class and the quality metrics.

Built from the E8 gallery output (clusters.csv); adds the class base rate, fold enrichment and Hennig's
stability rating. One row per label x feature set x resolution x cluster.

    python scripts/infer_bench/cluster_table.py GALLERY_DIR OUT_STEM
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402


def rating(s: float) -> str:
    """Hennig (2007) clusterboot reading of the mean best Jaccard."""
    if not np.isfinite(s):
        return ""
    return "stable" if s > 0.75 else "pattern" if s >= 0.6 else "weak" if s >= 0.5 else "dissolved"


def main(argv=None) -> int:
    src, stem = (argv if argv is not None else sys.argv[1:])
    c = pd.read_csv(os.path.join(src, "clusters.csv"))
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    rate = {}
    for lab in c["label"].unique():
        y = pd.Series(np.asarray(L[lab][0], dtype=object)).dropna()
        rate.update({(lab, k): v for k, v in y.value_counts(normalize=True).items()})
    c["class_rate"] = [rate.get((a, b), np.nan) for a, b in zip(c["label"], c["top_class"].astype(str))]
    c["fold_enrichment"] = c["top_share"] / c["class_rate"]
    c["enriched"] = c["q"] < 0.05
    c["stability_rating"] = c["stability"].map(rating)
    out = c.rename(columns={"embedding": "feature_set", "rank": "cluster_rank", "cluster": "cluster_id",
                            "genes": "cluster_genes", "labelled": "labelled_genes",
                            "top_count": "top_class_genes", "top_share": "top_class_share",
                            "p": "enrichment_p", "q": "enrichment_q", "silhouette": "silhouette_mean"})
    cols = ["label", "feature_set", "resolution", "cluster_rank", "cluster_id", "cluster_genes", "labelled_genes",
            "top_class", "top_class_genes", "top_class_share", "class_rate", "fold_enrichment", "enrichment_p",
            "enrichment_q", "enriched", "stability", "stability_rating", "silhouette_mean", "mean_gene_confidence"]
    out = out[cols].sort_values(["label", "feature_set", "resolution", "cluster_rank"])
    out.to_csv(stem + ".csv", index=False, float_format="%.4g")
    out.to_parquet(stem + ".parquet", index=False)
    print(len(out), "clusters;", int(out["enriched"].sum()), "enriched;",
          out["stability_rating"].value_counts().to_dict())
    return 0


if __name__ == "__main__":
    sys.exit(main())
