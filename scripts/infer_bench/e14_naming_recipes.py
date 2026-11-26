"""Infer E14: do the per-label tuned recipes NAME genes better than the map? (Einar 2026-11-26: test naming first)

For one E11 label, on the untouched gene families (E7 final split; development genes supply the labels the rule
learns from), the app's naming rule (`question_templates.dominant_by_cluster`: >= 15 labelled genes, >= 2x
enrichment, >= 30 % precision) is applied to three partitions:
  - tuned: the configuration E11 chose for this label (supervised parts fitted on development genes only);
  - E10 default: whitened profile, 5-seed Louvain consensus (resolution 4);
  - map-like: the app's map approximated - UMAP (15 neighbours, min_dist 0, cosine, 10 dimensions) of the E8
    reference profile, then HDBSCAN with the app's floor rule (smallest cluster = ceil(15 / labelled fraction)).
Scored on untouched labelled genes: precision of the named label, coverage (share named).

    python scripts/infer_bench/e14_naming_recipes.py --label NAME --e11 E11_DIR --out DIR
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e7_inference as E7  # noqa: E402
import e11_optimize as E11  # noqa: E402
from parascape import question_templates as QT  # noqa: E402

FLOOR, ENRICH, PREC = 15, 2.0, 0.3


def name_score(labels_full, y, dev, fin):
    """Apply the app's rule with development labels only; score untouched labelled genes."""
    train = np.where(dev, y, None).astype(object)
    dom = QT.dominant_by_cluster(labels_full, pd.Series(train), min_label=FLOOR, enrichment=ENRICH)
    ok = dom[dom["agrees"] & (dom["share_in_cluster"] >= PREC)].set_index("cluster")["label"]
    held = np.flatnonzero(fin & pd.notna(y) & (labels_full >= 0))
    allheld = int((fin & pd.notna(y)).sum())
    named = [(g, ok[labels_full[g]]) for g in held if labels_full[g] in ok.index]
    right = sum(str(y[g]) == lbl for g, lbl in named)
    return {"held_out": allheld, "named": len(named), "right": right,
            "precision": right / len(named) if named else np.nan, "coverage": len(named) / max(allheld, 1),
            "groups": int(len(np.unique(labels_full[labels_full >= 0])))}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--label", required=True)
    ap.add_argument("--e11", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    final = None
    for f in glob.glob(os.path.join(args.e11, "L*", "final.json")):
        j = json.load(open(f))
        if j["label"] == args.label and not j["shuffled"]:
            final = j
    if final is None:
        raise SystemExit(f"no E11 result for {args.label!r}")
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    y, kind, tcols = L[args.label]
    y = np.asarray(y, dtype=object)
    drop = set(H.IM.exclusions(tcols, views))
    fam = H.families(nodes)
    fin = E7.final_split(fam)
    dev = ~fin
    D = E11.Data(nodes, views, extra, drop)
    R = E11.rankings(D, y, dev, seed=9)
    rows = []
    configs = {"tuned (E11)": final["chosen_config"],
               "E10 default": {**final["chosen_config"], "extract": "all", "seq_views": True, "pca": 0,
                               "equal_views": False, "whiten": True, "supervised": "none", "method": "louvain",
                               "knn": 15, "metric": "cosine", "resolution": 4.0, "consensus": 5}}
    for name, cfg in configs.items():
        E, rws = E11.embed(D, R, cfg, y, dev)
        c = E11.cluster(E, cfg, seed=0)
        full = np.full(D.n, -1); full[rws] = c
        rows.append({"label": args.label, "partition": name, **name_score(full, y, dev, fin)})
        print(rows[-1], flush=True)
    # map-like: UMAP + HDBSCAN of the E8 reference profile, the app's floor rule for the smallest cluster
    import umap
    from sklearn.cluster import HDBSCAN
    ref = {**configs["E10 default"], "pca": 30, "whiten": False}
    E, rws = E11.embed(D, R, ref, y, dev)
    Z = umap.UMAP(n_neighbors=15, min_dist=0.0, n_components=10, metric="cosine", random_state=0).fit_transform(E)
    frac = float(pd.notna(np.where(dev, y, None)[rws]).mean())
    floor = int(np.ceil(FLOOR / max(frac, 1e-3)))
    c = HDBSCAN(min_cluster_size=max(floor, 5), copy=True).fit_predict(Z)
    full = np.full(D.n, -1); full[rws] = c
    rows.append({"label": args.label, "partition": f"map-like (UMAP + HDBSCAN, min {max(floor, 5)})",
                 **name_score(full, y, dev, fin)})
    print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(args.out, "naming.csv"), index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
