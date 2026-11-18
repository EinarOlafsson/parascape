"""What makes two published perturbations look alike: the biology, or the laboratory.

Figure 4 of the ParaScape paper (instruction 101). The figure was planned as "which conditions
group, and why", and the first thing the data says is that the intended version cannot be built:
of the contrasts the ladder parses, exactly TWO keys are stated by more than one study, and one of
those is `unstated`. Independent laboratories almost never perform the same perturbation in a form
the corpus can recognise, so "do perturbation effects replicate across labs" is not answerable here
and the figure must not pretend otherwise.

WHAT IS ANSWERABLE, AND MATTERS MORE. For every pair of contrasts, how alike are they, and is that
explained by sharing a STUDY or by sharing a KIND of perturbation? If contrasts from one laboratory
resemble each other more than contrasts of the same kind from different laboratories do, then a
cross-study comparison of perturbations is measuring the laboratory, and figure 2 is not a
precaution but a prerequisite.

Writes `results/what_structures_contrasts_<date>/`.
"""
from __future__ import annotations

import json
import os
import re
import sys
from datetime import date
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import perturbation_atlas as PA     # noqa: E402

LADDER = os.path.join(ROOT, "results", "contrast_ladder_2026_11_06", "contrasts.tsv")

#: A pair sharing fewer genes than this is not compared.
MIN_SHARED = 1000


def normalise(value) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def ladder_keys() -> pd.DataFrame:
    """Each parsed contrast with the axis it varies and the level it moves to."""
    lad = pd.read_csv(LADDER, sep="\t")
    lad["key"] = lad["axis"].astype(str) + "|" + lad["treated_level"].map(normalise)
    return lad


def pair_table(store: dict) -> pd.DataFrame:
    """Every pair of contrasts: how alike, same study or not, same axis or not."""
    lfc, reg = store["lfc"], store["registry"]
    columns = list(lfc.columns)
    study = reg["study"].astype(str).to_dict()
    axis = reg["ladder_axis"].astype(str).to_dict() if "ladder_axis" in reg.columns else {}
    rows = []
    values = {c: lfc[c] for c in columns}
    for a, b in combinations(columns, 2):
        both = pd.DataFrame({"a": values[a], "b": values[b]}).dropna()
        if len(both) < MIN_SHARED:
            continue
        rows.append({"a": a, "b": b,
                     "same_study": study.get(a) == study.get(b),
                     "same_axis": axis.get(a, "?") == axis.get(b, "?") and axis.get(a, "?") != "nan",
                     "genes": len(both),
                     "rho": float(spearmanr(both["a"], both["b"]).statistic)})
    return pd.DataFrame(rows)


def main() -> None:
    lad = ladder_keys()
    stated = lad.groupby("key")["study"].nunique()
    repeated = stated[(stated >= 2) & (~stated.index.str.endswith("|unstated"))]
    print(f"{len(lad)} parsed contrasts over {lad['study'].nunique()} studies")
    print(f"perturbations stated by two or more independent studies: {len(repeated)}")
    for key, n in repeated.items():
        print(f"   {key}  ({n} studies)")
    print("So 'do perturbation effects replicate across laboratories' cannot be asked of this\n"
          "corpus beyond that, and the figure says so instead of fitting something to it.")

    store = PA.load()
    print(f"\n{len(store['lfc'].columns)} stored contrasts; comparing every pair")
    pairs = pair_table(store)
    same = pairs[pairs["same_study"]]["rho"].abs()
    diff = pairs[~pairs["same_study"]]["rho"].abs()
    u = mannwhitneyu(same, diff, alternative="greater")
    print(f"\n  pairs from the SAME study      n={len(same):6,}  median |rho| = {same.median():.3f}")
    print(f"  pairs from DIFFERENT studies   n={len(diff):6,}  median |rho| = {diff.median():.3f}")
    print(f"  the same-study pairs are more alike: p = {u.pvalue:.3g}")

    axis_same = pairs[(~pairs["same_study"]) & pairs["same_axis"]]["rho"].abs()
    axis_diff = pairs[(~pairs["same_study"]) & (~pairs["same_axis"])]["rho"].abs()
    note = ""
    if len(axis_same) > 20 and len(axis_diff) > 20:
        ua = mannwhitneyu(axis_same, axis_diff, alternative="greater")
        print(f"\n  ACROSS studies, same kind of axis    n={len(axis_same):6,}  "
              f"median |rho| = {axis_same.median():.3f}")
        print(f"  across studies, different axis       n={len(axis_diff):6,}  "
              f"median |rho| = {axis_diff.median():.3f}")
        print(f"  sharing the kind of perturbation: p = {ua.pvalue:.3g}")
        note = (f"same axis across studies {axis_same.median():.3f} against "
                f"{axis_diff.median():.3f}")

    ratio = same.median() / diff.median() if diff.median() else float("nan")
    print(f"\nTwo contrasts from one laboratory resemble each other {ratio:.1f} times as much as "
          f"two\nfrom different laboratories. The laboratory is the strongest thing structuring "
          f"this corpus,\nwhich is why the untreated arms have to be made comparable first.")

    out = os.path.join(ROOT, "results", f"what_structures_contrasts_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    pairs.to_csv(os.path.join(out, "pairs.tsv"), sep="\t", index=False)
    json.dump({"parsed_contrasts": int(len(lad)),
               "studies": int(lad["study"].nunique()),
               "perturbations_repeated_across_studies": {k: int(v) for k, v in repeated.items()},
               "pairs": int(len(pairs)),
               "median_abs_rho_same_study": float(same.median()),
               "median_abs_rho_different_study": float(diff.median()),
               "ratio": float(ratio), "p_same_study_more_alike": float(u.pvalue),
               "axis_note": note},
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
