"""Which kind of measurement explains which other kind, across every category of knowledge.

Figure 6 of the ParaScape paper (instruction 101), the panel that generalises the two already
built. The localisation panel asks whether transcription places a protein; the protein-abundance
panel asks what sets how much of it there is. This asks the same question for every pair of
categories at once: if I have measurements of kind A, how much of a measurement of kind B can I
recover?

WHY THAT IS A BIOLOGICAL QUESTION AND NOT A BOOKKEEPING ONE. A pair that predicts well is a pair
the parasite couples: transcription places the organelles it builds in waves because their genes
share a schedule. A pair that does not is a place where the cell has a degree of freedom. The
answer is about the organism, and only incidentally about which experiment to fund.

EVERY SCORE IS AGAINST A SHUFFLED BASELINE of the same procedure on the same genes, because a
cross-validated R2 on a few hundred genes with many columns is not zero by default.

A PREDICTOR NEVER INCLUDES ITS OWN TARGET, nor any column of the target's own category: otherwise
`place` predicts `place` and the figure reports that a measurement equals itself.

AND A PAIR CAN BE CIRCULAR WITHOUT SHARING A COLUMN. Translation efficiency is ribosome footprints
divided by RNA, so a predictor set holding both the footprints and the efficiency of one study
determines that study's RNA arithmetically. Predicting `transcript` from `translation` therefore
recovers a quantity that is already inside the predictors, and it scores 0.81 for that reason as
much as for any coupling in the parasite. Such pairs are marked STRUCTURAL and left out of the
headline, rather than quietly reported as the strongest result.

Writes `results/which_data_explains_which_<date>/`.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import information_categories as IC    # noqa: E402

NOT_KNOWLEDGE = ("identity", "attention")

#: One well-measured continuous target per category, chosen because it is the quantity a reader
#: would name for that category rather than because it scores well.
TARGETS = {
    "transcript": "expr_tachy",
    "protein_abundance": "protein_ibaq_log2",
    "fitness": "fit_invitro_hff",
    "translation": "te302107_tachy_r1",
    "modification": "n_phosphosites",
    "conservation": "snp_nonsynonymous",
}

#: Pairs where the predictors determine the target by arithmetic rather than by biology.
#: (translation -> transcript): translation efficiency is footprints over RNA, so footprints and
#: efficiency together give RNA back.
STRUCTURAL = {("translation", "transcript")}

#: A category needs this many usable columns to be a predictor.
MIN_PREDICTORS = 3

#: Genes needed before a pair is scored at all.
MIN_GENES = 200

N_FOLDS, N_REPEATS, N_SHUFFLES = 5, 4, 10

#: Columns are capped per category so one wide category does not win on width alone.
MAX_COLUMNS = 40


def usable_columns(nodes: pd.DataFrame) -> dict:
    """Category -> its measured numeric columns, capped, excluding derivations and flags."""
    derived = set(getattr(IC, "DERIVED_FROM", {}))
    out = {}
    for column in nodes.columns:
        series = nodes[column]
        if column in derived or not series.isna().any():
            continue
        if not pd.api.types.is_numeric_dtype(series):
            continue
        for category in IC.categories_of(column):
            if category not in NOT_KNOWLEDGE:
                out.setdefault(category, []).append(column)
    return {k: sorted(v, key=lambda c: -nodes[c].notna().sum())[:MAX_COLUMNS]
            for k, v in out.items() if len(v) >= MIN_PREDICTORS}


def scored(x: np.ndarray, y: np.ndarray, seed: int = 0) -> float:
    """Variance of `y` explained out of sample, averaged over repeated folds."""
    out = []
    for repeat in range(N_REPEATS):
        folds = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed + repeat)
        predicted = np.empty_like(y)
        for train, test in folds.split(x):
            scaler = StandardScaler().fit(x[train])
            model = RidgeCV(alphas=np.logspace(-3, 4, 16)).fit(scaler.transform(x[train]),
                                                               y[train])
            predicted[test] = model.predict(scaler.transform(x[test]))
        out.append(1.0 - ((y - predicted) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    return float(np.mean(out))


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    nodes = pd.read_parquet(os.path.join(cache, "nodes.parquet"))
    nodes = nodes.set_index(nodes["gene_id"].astype(str))
    columns = usable_columns(nodes)
    print(f"{len(columns)} categories usable as predictors: "
          + ", ".join(f"{k} ({len(v)})" for k, v in sorted(columns.items())))

    rng = np.random.default_rng(0)
    rows = []
    for target_category, target in TARGETS.items():
        if target not in nodes.columns:
            print(f"  target {target!r} absent; {target_category} not asked")
            continue
        y_all = pd.to_numeric(nodes[target], errors="coerce")
        own = set(columns.get(target_category, [])) | {target}
        for predictor, cols in sorted(columns.items()):
            if predictor == target_category:
                continue
            use = [c for c in cols if c not in own]
            if len(use) < MIN_PREDICTORS:
                continue
            frame = nodes[use].apply(pd.to_numeric, errors="coerce")
            frame = frame.assign(_y=y_all).dropna()
            if len(frame) < MIN_GENES:
                continue
            y = frame.pop("_y").to_numpy(dtype=float)
            x = frame.to_numpy(dtype=float)
            real = scored(x, y)
            null = np.mean([scored(x, rng.permutation(y), seed=50 + i)
                            for i in range(N_SHUFFLES)])
            structural = (predictor, target_category) in STRUCTURAL
            rows.append({"target": target, "target_category": target_category,
                         "predictor_category": predictor, "columns": len(use),
                         "genes": int(len(frame)), "r2": real, "shuffled_r2": null,
                         "lift": real - null, "structural": structural})
            print(f"  {target_category:18s} <- {predictor:18s} n={len(frame):5,} "
                  f"R2={real:+.3f} (shuffled {null:+.3f})")

    table = pd.DataFrame(rows).sort_values("lift", ascending=False)
    structural = table[table["structural"]]
    table_real = table[~table["structural"]]
    best = table_real[table_real["lift"] > 0.05]
    if len(structural):
        for row in structural.itertuples():
            print(f"  STRUCTURAL, left out of the headline: {row.predictor_category} -> "
                  f"{row.target_category} at R2 = {row.r2:.3f}; the predictors determine the "
                  f"target by arithmetic.")
    print(f"\n{len(best)} of {len(table_real)} non-structural pairs recover more than 5 points "
          f"above their own shuffled baseline.\n")
    for row in best.head(10).itertuples():
        print(f"  {row.predictor_category:18s} explains {row.target_category:18s} "
              f"R2 = {row.r2:.3f} on {row.genes:,} genes")

    out = os.path.join(ROOT, "results", f"which_data_explains_which_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "pairs.tsv"), sep="\t", index=False)
    json.dump({"targets": TARGETS, "pairs": int(len(table)),
               "pairs_above_5_points": int(len(best)),
               "structural_pairs_excluded": sorted(f"{a}->{b}" for a, b in STRUCTURAL),
               "max_columns_per_category": MAX_COLUMNS,
               "table": table.to_dict("records")},
              open(os.path.join(out, "summary.json"), "w"), indent=2, default=float)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
