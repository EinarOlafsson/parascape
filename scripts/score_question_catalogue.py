#!/usr/bin/env python3
"""Score the admitted catalogue questions through the harness (instruction 73, phase 4), resumably.

Einar (2026-11-08): keep all 177 admitted questions. Each catalogue row becomes one or more yes/no questions:

* a **value** column: the extreme 20 % in the direction the question means -- the most depleted genes for a fitness
  screen (lowest values), the largest absolute change for a contrast that has both signs, otherwise the highest
  values -- against the other labelled genes;
* a **class** column: each class with at least 50 genes against the rest (a two-valued column: its True / larger
  value).

Every question runs through `harness.run_question`: the label's own dataset family and information category held
out, orthogroup-grouped nested cross-validation, calibration, the frozen test touched once and logged. One JSON per
question in the output directory; a question already there is skipped, so the run resumes. `--shard i/n` splits
the work across processes.

    python scripts/score_question_catalogue.py --shard 0/4
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

CATALOGUE = os.path.join(ROOT, "results", "question_catalogue_2026_11_08", "catalogue.tsv")
OUT = os.path.join(ROOT, "results", "question_scores_2026_11_08")
EXTREME = 0.2
MIN_CLASS = 50


def _safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(text)).strip("_")[:80]


def binary_labels(column, kind, category, series):
    """[(suffix, wording, boolean labels over the labelled genes)] for one catalogue row."""
    import numpy as np
    s = series.dropna()
    out = []
    if kind == "value":
        v = s.astype(float)
        if category == "fitness":
            cut = v.quantile(EXTREME)
            out.append(("low", f"among the most depleted 20 % ({column} <= {cut:.3g})", v <= cut))
        elif (v < 0).any() and (v > 0).any():
            cut = v.abs().quantile(1 - EXTREME)
            out.append(("change", f"among the largest 20 % of changes (|{column}| >= {cut:.3g})", v.abs() >= cut))
        else:
            cut = v.quantile(1 - EXTREME)
            out.append(("high", f"among the highest 20 % ({column} >= {cut:.3g})", v >= cut))
        return out
    values = s.astype(str)
    counts = values.value_counts()
    if len(counts) == 2:
        pos = "True" if "True" in counts.index else sorted(counts.index)[-1]
        out.append((_safe(pos), f"{column} = {pos}", values == pos))
        return out
    for cls, n in counts.items():
        if n >= MIN_CLASS and n <= len(values) - MIN_CLASS:
            out.append((_safe(cls), f"{column} = {cls}", values == cls))
    return out


def main(argv=None) -> int:
    """Score this shard's questions; skip those already scored."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args(argv)
    i, n = map(int, args.shard.split("/"))
    import numpy as np
    import pandas as pd
    from parascape import datasets as ds, harness as H, organisms as O
    from parascape.question_registry import Question
    os.makedirs(args.out, exist_ok=True)
    nodes = pd.read_parquet(os.path.join(ROOT, "parascape", "data", "nodes.parquet")).set_index("gene_id")
    cat = pd.read_csv(CATALOGUE, sep="\t")
    cat = cat[cat["admitted"]].reset_index(drop=True)

    def family_of(column):
        found = ds.provenance(column, O.TOXOPLASMA)
        return found.key if found is not None else "unregistered"
    numeric = [c for c in nodes.columns if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]
    families = pd.Series({c: family_of(c) for c in numeric}, name="family")
    groups_all = nodes["orthogroup"].astype(object)
    groups_all = groups_all.where(groups_all.notna(), pd.Series({g: "self:" + g for g in nodes.index}))
    jobs = [r for k, r in enumerate(cat.itertuples(index=False)) if k % n == i]
    print(f"shard {i}/{n}: {len(jobs)} catalogue rows", flush=True)
    for row in jobs:
        for suffix, wording, labels in binary_labels(row.column, row.kind, row.category, nodes[row.column]):
            qid = f"cat_{_safe(row.column)}_{suffix}"
            path = os.path.join(args.out, qid + ".json")
            if os.path.exists(path):
                continue
            family = family_of(row.column)
            q = Question(id=qid, text=f"{row.question} -- yes when {wording}", answer_type="boolean",
                         question_class="B", label_source=f"nodes.parquet:{row.column}", label_family=family,
                         evidence_forbidden=(family,), label_tier="observed-screen", label_category=row.category,
                         status="scored")
            entities = labels.index
            t0 = time.monotonic()
            try:
                full = nodes.loc[entities, numeric].astype("float64")
                features = H.allowed_features(q, full, families)
                r = H.run_question(q, labels=labels.astype(bool), features=features, families=families,
                                   groups=groups_all.reindex(entities),
                                   strategies=[H.Strategy.of("boosted trees", "hist_gbt", max_iter=150)], seed=0,
                                   importance_repeats=1, log_path=os.path.join(ROOT, H.HOLDOUT_LOG))
                ft = r.frozen_test or {}
                rec = dict(question=qid, column=row.column, category=row.category, wording=wording,
                           labelled=int(len(labels)), positives=int(labels.sum()), features=int(features.shape[1]),
                           cv_auc=r.cross_validated.get("roc_auc"), frozen_auc=ft.get("roc_auc"),
                           frozen_ap=ft.get("average_precision"), prevalence=float(labels.mean()),
                           precision_lower90=ft.get("precision_lower90"), may_fill=bool(r.may_fill),
                           flags=list(r.flags), refusal=r.refusal,
                           top_families=r.importance.head(3)["family"].tolist() if len(r.importance) else [],
                           seconds=round(time.monotonic() - t0, 1))
            except Exception as error:
                rec = dict(question=qid, column=row.column, category=row.category, wording=wording,
                           error=f"{type(error).__name__}: {error}"[:500], seconds=round(time.monotonic() - t0, 1))
            with open(path + ".partial", "w") as fh:
                json.dump(rec, fh, default=float)
            os.replace(path + ".partial", path)
            print(qid, rec.get("frozen_auc", rec.get("error")), rec["seconds"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
