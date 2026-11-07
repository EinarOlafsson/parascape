"""Cluster and module recovery, held out (instruction 64.14), summarised from what ships, as a notebook.

Reads the gallery label scores (maps chosen on a label-free structure score, scored against labels with a
shuffled-label chance) and the calibration self-tests of the three module mechanisms (the cluster for a label is
chosen on known genes and scored on hidden ones). Nothing is refitted.

    python scripts/build_cluster_benchmark_summary_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"cluster_benchmark_summary_{DATE}")
NB = os.path.join(REPO, "notebooks", f"cluster_benchmark_summary_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("How well do clusters recover labels they never saw?", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**Two records, read together (instruction 64.14).** (1) The UMAP gallery: every map was chosen by a "
          "label-free structure score, so no label picked the reported winner; each label is then scored on each map "
          "(category F1 against five shuffled-label runs; unclustered genes count against recall), and maps whose "
          "inputs restate the label are flagged circular. (2) The calibration self-tests of the module mechanisms: "
          "the cluster for a label is chosen on known genes and scored on hidden genes against a measured null.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, json, pandas as pd",
            "from parascape import calibration as C",
            "g = pd.read_csv(os.path.join(REPO, 'parascape/data/umap_gallery_scores.tsv'), sep='\\t')",
            "clean = g[~g.circular.astype(bool)]",
            "best = clean.sort_values('category_skill', ascending=False).groupby(['organism', 'label']).head(1)",
            "best = best[['organism', 'label', 'map', 'category_f1', 'category_f1_chance', 'category_skill', 'best_category', 'best_skill']]",
            "best.to_csv(os.path.join(OUT, 'best_noncircular_map_per_label.tsv'), sep='\\t', index=False)",
            "print(len(g), 'label x map scores;', int(g.circular.astype(bool).sum()), 'flagged circular and set aside')",
            "best.round(3)")
    nb.code("rows = []",
            "for key in ('holdout_search', 'consensus_modules', 'multiplex_modules'):",
            "    for organism in ('Tg', 'Pf'):",
            "        e = C.entry(key, organism) or {}",
            "        d = e.get('default') or {}",
            "        rows.append(dict(mechanism=key, organism=organism, metric=e.get('metric', ''), conclusive=d.get('conclusive'),",
            "                         observed=d.get('observed'), chance=d.get('chance'), skill=d.get('skill'),",
            "                         skill_low=d.get('skill_low'), skill_high=d.get('skill_high'), grade=C.grade(e) if e else None))",
            "modules = pd.DataFrame(rows)",
            "modules.to_csv(os.path.join(OUT, 'module_mechanisms_held_out.tsv'), sep='\\t', index=False)",
            "modules.round(3)")
    nb.write(NB)


if __name__ == "__main__":
    main()
