"""Stability of the host-facing sub-clusters and confidence of every gene in them.

The secretory cluster of the all-views map (Louvain 1.0; `show_host_clusters.py`) is re-clustered at
resolution 4.0 (seed 0) as the reference. Ten redraws keep a random 80 % of the columns of every view (and every
extra view whole), rebuild the embedding, and re-cluster the same genes at resolution 4.0 with a new seed. Per
sub-cluster: Hennig stability (mean best Jaccard); per gene: confidence (mean Jaccard of its reference
sub-cluster with the redraw sub-cluster it lands in).

    python scripts/infer_bench/host_subcluster_stability.py OUT_DIR
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, HERE)
import cluster_quality as Q  # noqa: E402
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e5_clustering as E5  # noqa: E402
import e8_selected_clustering as E8  # noqa: E402

REDRAWS = 10


def main(argv=None) -> int:
    out = (argv if argv is not None else sys.argv[1:])[0]
    os.makedirs(out, exist_ok=True)
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    eff = np.asarray(L["effect on host"][0], dtype=object) == "host-facing"
    pap = np.asarray(L["host-facing (papers)"][0], dtype=object)
    allf = {v: c for v, c in views.items() if len(c) >= 2}
    E, rows = E8.embed(nodes, views, extra, allf, list(extra))
    c = E5.louvain(E, 1.0, 0)
    top = int(pd.Series(c[eff[rows]]).value_counts().index[0])
    inside = np.flatnonzero(c == top)
    genes = rows[inside]
    ref = E5.louvain(E[inside], 4.0, 0)
    labs = []
    for b in range(REDRAWS):
        rng = np.random.default_rng(500 + b)
        f2 = {v: list(rng.choice(cols, max(2, int(0.8 * len(cols))), replace=False)) for v, cols in allf.items()}
        E2, rows2 = E8.embed(nodes, views, extra, f2, list(extra))
        pos = {g: i for i, g in enumerate(rows2)}
        keep = np.array([g in pos for g in genes])
        lab = np.full(len(genes), -1)
        lab[keep] = E5.louvain(E2[[pos[g] for g in genes[keep]]], 4.0, 1 + b)
        labs.append(lab)
        print(f"redraw {b + 1}/{REDRAWS}", flush=True)
    stab, conf = Q.stability(ref, labs)
    per = []
    for k in np.unique(ref):
        m = ref == k
        g = genes[m]
        prod = nodes["product"].iloc[g].fillna("").str.lower()
        per.append({"subcluster_id": int(k), "genes": int(m.sum()), "effect_on_host": int(eff[g].sum()),
                    "host_facing_papers": int((pap[g] == "host-facing").sum()),
                    "hypothetical": int(prod.str.contains("hypothetical").sum()),
                    "stability": float(stab[k]), "mean_gene_confidence": float(conf[m].mean())})
    P = pd.DataFrame(per)
    P["effect_share"] = P["effect_on_host"] / P["genes"]
    P = P.sort_values("effect_share", ascending=False).reset_index(drop=True)
    P["figure_rank"] = np.arange(1, len(P) + 1)
    P.to_csv(os.path.join(out, "subcluster_stability.csv"), index=False)
    G = pd.DataFrame({"gene_id": nodes["gene_id"].to_numpy()[genes], "product": nodes["product"].to_numpy()[genes],
                      "subcluster_id": ref, "confidence": conf, "effect_on_host": eff[genes],
                      "host_facing_papers": np.where(pd.notna(pap[genes]), pap[genes], "")})
    G = G.merge(P[["subcluster_id", "figure_rank", "stability", "effect_share"]], on="subcluster_id")
    cand = G[(G["figure_rank"] <= 4) & ~G["effect_on_host"] & (G["host_facing_papers"] != "host-facing")]
    cand = cand.sort_values(["figure_rank", "confidence"], ascending=[True, False])
    G.to_csv(os.path.join(out, "subcluster_genes_confidence.csv"), index=False)
    cand.to_csv(os.path.join(out, "candidates.csv"), index=False)
    print(P.head(8).round(3).to_string(index=False))
    print(cand.head(40).round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
