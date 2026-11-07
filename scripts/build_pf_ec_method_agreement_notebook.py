"""Paired method agreement and error dependence on the Plasmodium EC cohort (instruction 68.02), as a notebook.

The next small task 68.02 named: on the exact frozen cohort the random forest, kNN and the two baselines were
scored on (152 held-out genes, the same roles and seed), how often do the methods agree, are their errors
the same genes, and is one better than another beyond chance? Descriptive only: no new fit, no tuning.

    python scripts/build_pf_ec_method_agreement_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"pf_ec_method_agreement_{DATE}")
NB = os.path.join(REPO, "notebooks", f"pf_ec_method_agreement_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    nb = ExecutedNotebook("Do the EC methods make the same mistakes?", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**The cohort.** The frozen Plasmodium EC test genes of 68.02 (152 genes, seed 20261008): native random "
          "forest, kNN, the training-majority call and the prevalence-weighted call, each already scored. A call "
          "is correct when its EC set equals the truth set -- the definition that reproduces the recorded counts "
          "(forest 59, kNN 40, majority 43, prevalence 35).",
          "**What is measured.** For every pair: agreement of the calls, how many genes both get right, both "
          "wrong, or only one; McNemar's exact test on the discordant genes (is one better beyond chance?); and "
          "the correlation of their errors. **Limit stated first:** both models read the same source features, "
          "so shared errors partly reflect shared inputs; this says nothing about biological accuracy.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code(
        "import itertools, json, os",
        "import numpy as np, pandas as pd",
        "from scipy.stats import binomtest",
        "R = os.path.join(REPO, 'results')",
        "rf = pd.read_parquet(f'{R}/pf_functional_ec_forest_2026_10_08_v2/rows.parquet').set_index('entity')",
        "kn = pd.read_parquet(f'{R}/pf_functional_ec_knn_2026_10_08_v2/rows.parquet').set_index('entity')",
        "ct = pd.read_parquet(f'{R}/pf_functional_ec_controls_2026_10_08/control_rows.parquet').set_index('entity')",
        "assert set(rf.index) == set(kn.index) == set(ct.index), 'not the same cohort'",
        "S = lambda v: frozenset(json.loads(v)) if isinstance(v, str) else frozenset()",
        "truth = ct['truth'].map(S)",
        "calls = pd.DataFrame({'random forest': rf['prediction'].map(S), 'kNN': kn['prediction'].reindex(ct.index).map(S),",
        "                      'majority': ct['majority'].map(S), 'prevalence': ct['prevalence_call'].map(S)}).reindex(ct.index)",
        "rf_ = rf.reindex(ct.index)",
        "right = calls.apply(lambda c: c == truth)",
        "right.sum().rename('correct of 152')")
    nb.code(
        "rows = []",
        "for a, b in itertools.combinations(calls.columns, 2):",
        "    ra, rb = right[a], right[b]",
        "    only_a, only_b = int((ra & ~rb).sum()), int((rb & ~ra).sum())",
        "    p = binomtest(only_a, only_a + only_b, 0.5).pvalue if only_a + only_b else float('nan')",
        "    ea, eb = (~ra).astype(float), (~rb).astype(float)",
        "    phi = float(np.corrcoef(ea, eb)[0, 1]) if ea.std() and eb.std() else float('nan')",
        "    rows.append(dict(method_a=a, method_b=b, same_call=round(float((calls[a] == calls[b]).mean()), 3),",
        "                     both_right=int((ra & rb).sum()), both_wrong=int((~ra & ~rb).sum()),",
        "                     only_a_right=only_a, only_b_right=only_b, mcnemar_p=round(float(p), 4),",
        "                     error_correlation=round(phi, 3)))",
        "pairs = pd.DataFrame(rows)",
        "pairs.to_csv(os.path.join(OUT, 'pairs.tsv'), sep='\\t', index=False)",
        "pairs")
    nb.md("### Would the two models together do better than either?",
          "An upper bound, not a method: the share of genes at least one model gets right (an oracle that always "
          "picks the right one). If it is barely above the forest alone, kNN adds nothing the forest lacks.")
    nb.code(
        "either = int((right['random forest'] | right['kNN']).sum())",
        "summary = pd.DataFrame([dict(genes=len(right), forest=int(right['random forest'].sum()), knn=int(right['kNN'].sum()),",
        "                             either_oracle=either, both=int((right['random forest'] & right['kNN']).sum()))])",
        "summary.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\\t', index=False)",
        "summary")
    nb.write(NB)


if __name__ == "__main__":
    main()
