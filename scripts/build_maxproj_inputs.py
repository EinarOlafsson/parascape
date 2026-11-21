"""Inputs of the figure 6 gene clustering, for spaCR's max-projection clustering (spacr.ops_maxproj).

Einar, 2026-11-21: the ParaScape gene clustering reuses spaCR's multi-view max-projection code
(item 992) instead of a second implementation. This script writes what that code reads, from the
same views and labels the figure 6 predictive matrix uses (scripts/build_information_matrix.py), so
the two analyses cannot drift apart:

  genes.parquet     one row per gene: gene_id and every view's columns (extra views, which are not
                    node columns, are named <view>__<i>)
  views.json        [{"name", "block", "columns", "views_ix"}], the views file format of
                    ops_maxproj.read_views plus each view's column indices in genes.parquet
                    (gene_id excluded, so index 0 is the first measurement column)
  coverage.parquet  gene x view, True where the gene has a value in that view
  pairs.csv         positive control: observed-direct parasite-parasite partner pairs,
                    gene_a,gene_b,source (the paper reference)
  markers.csv       positive controls: gene,label,label_set for the hyperLOPIT niche markers
                    (target only when 'place' is not a view) and the EC classes
  planted.csv       planted complexes, gene,complex,view,strength: per view 5 disjoint complexes
                    of 8 covered genes at each strength 1, 2, 4 and 8 (seed 20261121); the plant
                    recipe is in manifest.json, and a no-plant run uses no row of it
  labels.parquet    every figure 6 label, one column each
  manifest.json     counts, column md5s, the plant recipe and the source paths

Output: /media/carruthers/mnt2/parascape/maxproj_inputs (scratch disk).
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import build_information_matrix as IM                   # noqa: E402

PLANT_SEED = 20261121
PLANT_STRENGTHS = (1, 2, 4, 8)
PLANT_PER_STRENGTH = 5
PLANT_SIZE = 8

OUT = os.environ.get("PARASCAPE_MAXPROJ_INPUTS", "/media/carruthers/mnt2/parascape/maxproj_inputs")


def main():
    os.makedirs(OUT, exist_ok=True)
    nodes = pd.read_parquet(os.path.join(IM.CACHE, "nodes.parquet"))
    ids = nodes["gene_id"].astype(str).to_numpy()
    node_views = IM.build_views(nodes)
    extra = IM.extra_views(nodes)

    table = {"gene_id": ids}
    views, coverage = [], {}
    for name, cols in node_views.items():
        for c in cols:
            table[c] = nodes[c].to_numpy(float)
        views.append({"name": name, "block": name, "columns": list(cols)})
        coverage[name] = nodes[cols].notna().any(axis=1).to_numpy()

    for name, array in extra.items():
        cols = [f"{name}__{i}" for i in range(array.shape[1])]
        for i, c in enumerate(cols):
            table[c] = array[:, i]
        views.append({"name": name, "block": name, "columns": cols})
        coverage[name] = ~np.isnan(array).all(axis=1)
    genes = pd.DataFrame(table)
    position = {c: i - 1 for i, c in enumerate(genes.columns)}
    for view in views:
        view["views_ix"] = [position[c] for c in view["columns"]]
    genes.to_parquet(os.path.join(OUT, "genes.parquet"), index=False)
    with open(os.path.join(OUT, "views.json"), "w", encoding="utf-8") as handle:
        json.dump(views, handle, indent=1)
    cover = pd.DataFrame(coverage)
    cover.insert(0, "gene_id", ids)
    cover.to_parquet(os.path.join(OUT, "coverage.parquet"), index=False)

    pairs = pd.read_parquet(os.path.join(IM.CACHE, "information_map", "pairs.parquet"))
    pairs = pairs[(pairs.source_tier == "observed-direct") & (pairs.partner_side != "host")]
    known = set(ids)
    pairs = pairs[pairs.gene_a.isin(known) & pairs.gene_b.isin(known) & (pairs.gene_a != pairs.gene_b)]
    pairs = pd.DataFrame({"gene_a": np.minimum(pairs.gene_a, pairs.gene_b),
                          "gene_b": np.maximum(pairs.gene_a, pairs.gene_b),
                          "source": pairs.reference.astype(str)}).drop_duplicates(["gene_a", "gene_b"])
    pairs.to_csv(os.path.join(OUT, "pairs.csv"), index=False)

    imap = pd.read_parquet(os.path.join(IM.CACHE, "information_map", "values.parquet"))
    found = IM.labels(nodes, imap)
    labels = {"gene_id": ids}
    for name, (values, _kind, _sources) in found.items():
        labels[name] = values
    pd.DataFrame(labels).to_parquet(os.path.join(OUT, "labels.parquet"), index=False)
    markers = []
    for name, label_set in (("hyperLOPIT niche", "hyperLOPIT niche"),
                            ("enzymatic function (EC class)", "EC class")):
        values = pd.Series(found[name][0], index=ids)
        values = values[values.notna() & (values.astype(str) != "nan")]
        markers.append(pd.DataFrame({"gene": values.index, "label": values.astype(str).to_numpy(),
                                     "label_set": label_set}))
    markers = pd.concat(markers, ignore_index=True)
    markers.to_csv(os.path.join(OUT, "markers.csv"), index=False)

    rng = np.random.default_rng(PLANT_SEED)
    used, planted = set(), []
    for view in views:
        pool = [g for g, c in zip(ids, coverage[view["name"]]) if c and g not in used]
        for strength in PLANT_STRENGTHS:
            for k in range(PLANT_PER_STRENGTH):
                members = rng.choice(len(pool), PLANT_SIZE, replace=False)
                genes_k = [pool[i] for i in members]
                pool = [g for g in pool if g not in set(genes_k)]
                used.update(genes_k)
                name = f"plant_{view['name']}_s{strength}_{k}"
                planted += [{"gene": g, "complex": name, "view": view["name"], "strength": strength}
                            for g in genes_k]
    planted = pd.DataFrame(planted)
    planted.to_csv(os.path.join(OUT, "planted.csv"), index=False)

    manifest = {
        "genes": int(len(genes)), "views": {v["name"]: {"columns": len(v["columns"]),
                                                         "covered_genes": int(coverage[v["name"]].sum())}
                                             for v in views},
        "direct_pairs": int(len(pairs)), "labels": [k for k in labels if k != "gene_id"],
        "markers": markers.groupby("label_set").size().to_dict(),
        "planted": {"complexes": int(planted.complex.nunique()), "genes": int(len(planted)),
                    "size": PLANT_SIZE, "per_view_and_strength": PLANT_PER_STRENGTH,
                    "strengths": list(PLANT_STRENGTHS), "seed": PLANT_SEED,
                    "recipe": "for each planted complex, draw one unit direction over its view's "
                              "columns (seed = PLANT_SEED + complex index in planted.csv order) and "
                              "add strength x that column's standard deviation over covered genes x "
                              "the direction to every member's covered values in that view only; "
                              "the no-plant run adds nothing"},
        "columns_md5": hashlib.md5("\n".join(genes.columns).encode()).hexdigest(),
        "sources": {"nodes": os.path.join(IM.CACHE, "nodes.parquet"), "structure": IM.STRUCT,
                    "pairs": os.path.join(IM.CACHE, "information_map", "pairs.parquet")},
        "note": "hyperLOPIT niche is a target only in runs without the 'place' view (circularity).",
    }
    with open(os.path.join(OUT, "manifest.json"), "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=1)
    print(json.dumps(manifest, indent=1))


if __name__ == "__main__":
    main()
