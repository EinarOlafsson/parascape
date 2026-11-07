#!/usr/bin/env python3
"""Take the GSE108740 expression log once in the shipped Toxoplasma table, and re-lay out its map (instruction 64.11a).

The rna108740_* columns are ln(FPKM + 1) of GEO's GSE108740_FPKM.xlsx; expr_tachy, expr_cyst and expr_max were
log2(mean + 1) of those -- a log of a log -- and expr_max also mixed in centred-log columns of other series. This
applies exactly the corrected `build_graph` formula to the shipped `nodes.parquet` (the raw tree `build_graph`
reads is not on every machine), then recomputes the packaged layout with `build_graph.embed` and replaces only
`xyz` and `layout_metadata` in `graph.npz`; the edge layers do not use these columns and are kept byte for byte.
Writes a before/after record.

    python scripts/patch_expression_summaries.py [--out results/expression_fix_2026_11_07]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

COLUMNS = ("expr_tachy", "expr_cyst", "expr_max")


def corrected(nodes: pd.DataFrame) -> pd.DataFrame:
    """The three summaries as `build_graph.load_nodes` now computes them."""
    fpkm = [c for c in nodes.columns if c.startswith("rna108740_")]
    linear = np.expm1(nodes[fpkm].apply(pd.to_numeric, errors="coerce"))
    tach = [c for c in fpkm if "Tachyzoite" in c or "Tachyzoites" in c]
    cyst = [c for c in fpkm if "Tissue_cyst" in c or "Tissue_cysts" in c]
    return pd.DataFrame({"expr_tachy": np.log2(linear[tach].mean(axis=1) + 1),
                         "expr_cyst": np.log2(linear[cyst].mean(axis=1) + 1),
                         "expr_max": np.log2(linear[fpkm].max(axis=1) + 1)}, index=nodes.index)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "expression_fix_2026_11_07"))
    args = ap.parse_args(argv)
    from parascape import build_graph as B
    data = os.path.join(ROOT, "parascape", "data")
    nodes_path, graph_path = os.path.join(data, "nodes.parquet"), os.path.join(data, "graph.npz")
    nodes = pd.read_parquet(nodes_path)
    new = corrected(nodes)
    os.makedirs(args.out, exist_ok=True)
    record = {"before": nodes[list(COLUMNS)].describe().round(4).to_dict(),
              "after": new.describe().round(4).to_dict(),
              "spearman_before_after": {c: float(nodes[c].corr(new[c], method="spearman")) for c in COLUMNS},
              "nodes_sha256_before": hashlib.sha256(open(nodes_path, "rb").read()).hexdigest()}
    for c in COLUMNS:
        nodes[c] = new[c]
    nodes.to_parquet(nodes_path, index=False)
    with np.load(graph_path, allow_pickle=False) as g:
        flat = {k: g[k] for k in g.files}
    if not np.array_equal(flat["gene_ids"].astype(str), nodes.gene_id.astype(str).to_numpy()):
        raise ValueError("graph gene order differs from the node table")
    coords, metadata = B.embed(nodes, return_metadata=True)
    record["layout_before_sha256"] = hashlib.sha256(np.asarray(flat["xyz"]).tobytes()).hexdigest()
    flat["xyz"] = np.asarray(coords, dtype=flat["xyz"].dtype)
    # The packaged layout record ties the coordinates to the exact node file they were computed from.
    metadata["node_file_sha256"] = hashlib.sha256(open(nodes_path, "rb").read()).hexdigest()
    flat["layout_metadata"] = np.array(json.dumps(metadata))
    np.savez_compressed(graph_path, **flat)
    record.update(layout_after_sha256=hashlib.sha256(flat["xyz"].tobytes()).hexdigest(),
                  nodes_sha256_after=hashlib.sha256(open(nodes_path, "rb").read()).hexdigest())
    json.dump(record, open(os.path.join(args.out, "record.json"), "w"), indent=1)
    print(json.dumps({k: record[k] for k in ("spearman_before_after",)}, indent=1))
    print(pd.DataFrame(record["after"]).loc[["min", "50%", "max"]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
