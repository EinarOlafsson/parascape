"""Answer every shipped question with the map method and with the benchmarked method; compare the answers.

For each question template: whether it answered, mean and best recovery F1, genes named, genes the control
corroborates, and (benchmarked) median group cohesion. Writes compare_methods.csv.

    python scripts/infer_bench/compare_infer_methods.py OUT_DIR [--limit N]
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import replace

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e0_harness as H  # noqa: E402
from parascape import questions as Q  # noqa: E402
from parascape.question_templates import run as run_template  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--part", default="", help="i/n: only questions i, i+n, ... (parallel parts)")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    nodes = H.load()[0]
    qs = Q.load()
    if args.limit:
        qs = qs[: args.limit]
    tag = ""
    if args.part:
        i_, n_ = (int(x) for x in args.part.split("/"))
        qs = qs[i_::n_]
        tag = f"_part{i_}of{n_}"
    rows = []
    path = os.path.join(args.out, f"compare_methods{tag}.csv")
    for q in qs:
        t = Q.as_template(q)
        for method in ("umap+hdbscan", "benchmarked"):
            t0 = time.time()
            try:
                r = run_template(nodes, replace(t, method=method), log=lambda *a: None)
                inf = r.inference
                row = {"question": t.question[:90], "holdout": t.holdout, "method": method, "ok": r.ok,
                       "why_not": r.stopped_because[:120], "mean_f1": r.summary.get("mean_f1", np.nan),
                       "best_f1": r.summary.get("best_f1", np.nan), "named": len(inf),
                       "corroborated": int(inf["control_agrees"].fillna(False).sum())
                       if "control_agrees" in inf else 0,
                       "groups": r.quality.get("clusters", np.nan),
                       "median_cohesion": float(r.per_cluster.drop_duplicates("cluster")["cohesion"].median())
                       if "cohesion" in r.per_cluster else np.nan}
            except Exception as e:
                row = {"question": t.question[:90], "holdout": t.holdout, "method": method, "ok": False,
                       "why_not": f"error: {e}"[:120]}
            row["seconds"] = round(time.time() - t0, 1)
            rows.append(row)
            print(row, flush=True)
            pd.DataFrame(rows).to_csv(path, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
