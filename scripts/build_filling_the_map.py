"""How often a call is right on genes the method was never shown, and how many genes it reaches.

Figure 11 of the ParaScape paper (instruction 101). Without this the paper is a description of a
database; with it, it is the inference engine ParaScape exists to be.

WHAT IS REPORTED, PER KIND OF FACT. Three numbers, and all three are needed for any of them to
mean anything:

* how often the call is right **when the method answers**, on held-out genes;
* how often it answers at all, because a method that abstains on everything is right on nothing;
* what the same accuracy would be by always guessing the commonest answer, because a target where
  nine genes in ten are cytosol is 90 per cent right by saying cytosol.

A method that beats the commonest-answer baseline is doing work. One that does not is reporting
the shape of the label.

THE RECORD IS NOT RE-RUN HERE. `track_record.parquet` holds a million held-out predictions made by
the harness under frozen partitions. Recomputing them in a figure script would be a second copy of
the evaluation, free to drift from the one the application reports.

Writes `results/filling_the_map_<date>/`.
"""
from __future__ import annotations

import json
import os
from datetime import date

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: A target with fewer held-out calls than this is not summarised: a rate on a handful of calls
#: is a number with no precision.
MIN_CALLS = 50

COLUMNS = ["organism", "strategy", "target", "truth", "prediction", "correct", "abstained"]


def summarise(record: pd.DataFrame) -> pd.DataFrame:
    """Per organism, target and strategy: answered, right when answered, and the baseline."""
    rows = []
    for (organism, target, strategy), group in record.groupby(
            ["organism", "target", "strategy"], observed=True):
        answered = group[~group["abstained"].astype(bool)]
        if len(group) < MIN_CALLS:
            continue
        truth = group["truth"].astype(str)
        commonest = float(truth.value_counts(normalize=True).max()) if len(truth) else float("nan")
        rows.append({
            "organism": organism, "target": target, "strategy": strategy,
            "held_out_calls": int(len(group)),
            "answered": int(len(answered)),
            "answer_rate": float(len(answered) / len(group)),
            "right_when_answered": float(answered["correct"].astype(bool).mean())
            if len(answered) else float("nan"),
            "commonest_answer_baseline": commonest,
        })
    table = pd.DataFrame(rows)
    table["lift_over_baseline"] = table["right_when_answered"] - table["commonest_answer_baseline"]
    return table.sort_values(["target", "lift_over_baseline"], ascending=[True, False])


def best_per_target(table: pd.DataFrame) -> pd.DataFrame:
    """The strategy that does best on each target, which is what a user would be given."""
    answered = table[table["answered"] >= MIN_CALLS]
    if answered.empty:
        return answered
    return answered.loc[answered.groupby("target")["lift_over_baseline"].idxmax()]


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    record = pd.read_parquet(os.path.join(cache, "track_record.parquet"), columns=COLUMNS)
    print(f"{len(record):,} held-out predictions, "
          f"{record['target'].nunique()} targets, {record['strategy'].nunique()} strategies")

    table = summarise(record)
    best = best_per_target(table)
    print(f"\nbest strategy per target (of {len(table)} target-strategy pairs):\n")
    print(f"  {'target':26s} {'calls':>7s} {'answers':>8s} {'right':>7s} {'baseline':>9s} "
          f"{'lift':>7s}  strategy")
    for row in best.sort_values("lift_over_baseline", ascending=False).itertuples():
        print(f"  {row.target[:26]:26s} {row.held_out_calls:7,} {row.answer_rate:7.1%} "
              f"{row.right_when_answered:7.1%} {row.commonest_answer_baseline:9.1%} "
              f"{row.lift_over_baseline:+7.1%}  {row.strategy}")

    beats = best[best["lift_over_baseline"] > 0]
    print(f"\n{len(beats)} of {len(best)} targets are called better than by guessing the "
          f"commonest answer.")
    if len(best) > len(beats):
        failing = ", ".join(best[best["lift_over_baseline"] <= 0]["target"].astype(str))
        print(f"Not beaten, and reported rather than dropped: {failing}.")

    out = os.path.join(ROOT, "results", f"filling_the_map_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "by_target_and_strategy.tsv"), sep="\t", index=False)
    best.to_csv(os.path.join(out, "best_per_target.tsv"), sep="\t", index=False)
    json.dump({"held_out_predictions": int(len(record)),
               "targets": int(record["target"].nunique()),
               "strategies": int(record["strategy"].nunique()),
               "targets_summarised": int(len(best)),
               "targets_beating_baseline": int(len(beats)),
               "min_calls": MIN_CALLS,
               "best": best.to_dict("records")},
              open(os.path.join(out, "summary.json"), "w"), indent=2, default=float)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
