"""E10 cluster-size balance: share of genes in the largest cluster, and clusters under 5 genes, per arm (seed 0).

A clustering that puts most genes into a few giant clusters can score well on recall-type measures (complexes
kept together, STRING recall) while being useless for inference; this table lets the E10 figure show it.

    python scripts/infer_bench/e10_sizes.py OUT_CSV [--k 100,300]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e10_ops_ledger as E  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out")
    ap.add_argument("--k", default="100,300")
    ap.add_argument("--arms", default="", help="semicolon-separated arm names (default: all)")
    args = ap.parse_args(argv)
    nodes, views, extra, _ = H.load()
    rows = []
    for arm, _d in E.ARMS:
        if args.arms and arm not in args.arms.split(";"):
            continue
        for k in (int(x) for x in args.k.split(",")):
            if arm.startswith("shuffled"):
                sn, sx = E.shuffled_nodes(nodes, views, extra, 700)
                c = E.cluster_arm(sn, views, sx, "baseline", k, 0)[0]
            else:
                c = E.cluster_arm(nodes, views, extra, arm, k, 0)[0]
            s = np.sort(np.bincount(c[c >= 0]))[::-1]
            rows.append({"arm": arm, "k": k, "largest_share": s[0] / s.sum(), "top3_share": s[:3].sum() / s.sum(),
                         "clusters_under_5": int((s < 5).sum()), "median_size": float(np.median(s))})
            print(rows[-1], flush=True)
            pd.DataFrame(rows).to_csv(args.out, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
