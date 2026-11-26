"""Infer E13: how precisely does the app's naming rule name held-out genes, per method and labelled-gene floor?

The Infer tab names an unlabelled gene when its group holds at least MIN_LABEL (15) labelled genes, the group's
commonest label is >= 2x its map-wide rate (enrichment) and holds >= 30 % of the group's labelled genes
(precision). E12's comparison found the "Tested groups" partitions recover labels far better but name fewer
genes, because finer groups rarely reach 15 labelled genes. Here the rule itself is graded: for every question
both methods answer, the labelled genes are split into 3 folds; the rule is applied with the training genes'
labels (`dominant_by_cluster`, the app's own code) and the held-out labelled genes in a naming group are scored
- precision (named label = true label) and coverage (share of held-out genes named) - at floors 15, 10 and 5.

Writes naming.csv (question x method x floor x fold) and summary.csv.

    python scripts/infer_bench/e13_naming_floor.py OUT_DIR [--part i/n]
"""
from __future__ import annotations

import argparse
import os
import sys
from dataclasses import replace

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e0_harness as H  # noqa: E402
from parascape import questions as Q  # noqa: E402
from parascape import question_templates as QT  # noqa: E402

FLOORS = (15, 10, 5)
FOLDS = 3


def grade(labels, truth, template, floor, seed=0):
    y = truth.astype("object").where(truth.notna(), None).to_numpy()
    lab = np.flatnonzero(np.array([v is not None for v in y]) & (np.asarray(labels) >= 0))
    fold = np.random.default_rng(seed).permutation(len(lab)) % FOLDS
    out = []
    for f in range(FOLDS):
        held = lab[fold == f]
        train = y.copy()
        train[held] = None
        dom = QT.dominant_by_cluster(labels, pd.Series(train), min_label=floor, enrichment=template.min_enrichment)
        ok = dom[dom["agrees"] & (dom["share_in_cluster"] >= template.min_precision)].set_index("cluster")["label"]
        named = [(g, ok[labels[g]]) for g in held if labels[g] in ok.index]
        right = sum(str(y[g]) == lbl for g, lbl in named)
        out.append({"fold": f, "held_out": len(held), "named": len(named),
                    "precision": right / len(named) if named else np.nan,
                    "coverage": len(named) / len(held) if len(held) else np.nan})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out")
    ap.add_argument("--part", default="")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    nodes = H.load()[0]
    qs = Q.load()
    tag = ""
    if args.part:
        i_, n_ = (int(x) for x in args.part.split("/"))
        qs, tag = qs[i_::n_], f"_part{i_}of{n_}"
    rows = []
    path = os.path.join(args.out, f"naming{tag}.csv")
    for q in qs:
        t = Q.as_template(q)
        for method in ("umap+hdbscan", "benchmarked"):
            r = QT.run(nodes, replace(t, method=method), log=lambda *a: None)
            if not r.ok or not len(r.truth):
                continue
            for floor in FLOORS:
                for g in grade(np.asarray(r.labels), r.truth.reset_index(drop=True), t, floor):
                    rows.append({"question": t.question[:90], "holdout": t.holdout, "method": method,
                                 "floor": floor, "groups": r.quality.get("clusters"), **g})
            print(t.holdout, method, flush=True)
            pd.DataFrame(rows).to_csv(path, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
