"""Instruction 62 stage 2: does the pooled hold-out track record agree with the strategy scorecard?

The track record (`parascape/data/track_record*.parquet`) holds one row per hidden gene per strategy;
the scorecard (`strategy_calibration.json`, from `results/calibration_2026-09-26b/runs.jsonl`) holds
one number per calibration run. Pooled to the same strategy, target and setting, and read with the
same definition, the two should give the same answer. This notebook measures whether they do.

    python scripts/build_track_record_check_notebook.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_06"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"track_record_check_{DATE}")
NB = os.path.join(REPO, "notebooks", f"track_record_check_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Does the hold-out track record agree with the scorecard?",
                          namespace={"REPO": REPO, "OUT": OUT})
    nb.md(
        "**The check.** Instruction 62 built a per-gene hold-out record for every label strategy and "
        "required that, pooled to the strategy level, it agree with the scorecard the strategy cards "
        "already show. That agreement was never measured; this measures it.",
        "**Like with like.** The scorecard metric of most label strategies is *correct calls per "
        "hidden gene* (an abstention counts as not correct); three are scored by *precision of calls* "
        "(abstentions left out). The track record is reduced to the same quantity for each, at the "
        "same setting, per organism and target. The two use different folds and seeds, so they should "
        "agree within sampling error, not to the digit.")
    nb.code(
        "import json, os",
        "import numpy as np, pandas as pd",
        "from scipy.stats import spearmanr",
        "from parascape import track_record as T",
        "runs = pd.read_json(os.path.join(REPO, 'results/calibration_2026-09-26b/runs.jsonl'), lines=True)",
        "ledger = pd.concat([T.shipped(o) for o in ('Tg', 'Pf')], ignore_index=True)",
        "ledger = ledger[ledger['mode'].astype(str) == 'together']",
        "for c in ('organism', 'strategy', 'target', 'setting_key'):",
        "    ledger[c] = ledger[c].astype(str)",
        "print(len(runs), 'calibration runs;', len(ledger), 'hold-out rows over',",
        "      ledger.strategy.nunique(), 'strategies and', ledger.target.nunique(), 'targets')")

    nb.md("## 1 · The same quantity from both",
          "A calibration run records only the parameters it varied, and writes floats unrounded "
          "(`min_share=0.30000000000000004`); the record writes the full setting. A calibration run "
          "is matched to a record setting when every parameter it names has the same value there, "
          "after rounding. A record setting the calibration never tried has no scorecard cell, and is "
          "counted as such rather than compared with a neighbour.")
    nb.code(
        "def norm(v):",
        "    try:",
        "        return str(round(float(v), 6)).rstrip('0').rstrip('.') if not isinstance(v, bool) else str(v)",
        "    except (TypeError, ValueError):",
        "        return str(v)",
        "def parse(key):",
        "    return {k.strip(): norm(v) for k, v in (part.split('=', 1) for part in key.split(', ') if '=' in part)}",
        "lab = runs[runs.strategy.isin(set(ledger.strategy)) & runs.verdict.isin(['PASS', 'FAIL'])].copy()",
        "lab['target'] = lab.settings.map(lambda d: str(d.get('target', '')))",
        "lab['named'] = lab.settings.map(lambda d: {k: norm(v) for k, v in d.items() if k != 'target'})",
        "record_keys = ledger[['organism', 'strategy', 'setting_key']].drop_duplicates()",
        "def match(row):",
        "    keys = record_keys[(record_keys.organism == row.organism) & (record_keys.strategy == row.strategy)].setting_key",
        "    hits = [k for k in keys if all(parse(k).get(p) == v for p, v in row.named.items())]",
        "    return hits[0] if len(hits) == 1 else ''",
        "lab['setting_key'] = lab.apply(match, axis=1)",
        "unmatched = lab[lab.setting_key == ''].groupby(['organism', 'strategy']).size().rename('runs at settings the record does not hold')",
        "lab = lab[lab.setting_key != '']",
        "PRECISION = set(lab.loc[lab.metric.str.startswith('precision'), 'strategy'])",
        "card = (lab.groupby(['organism', 'strategy', 'target', 'setting_key', 'metric'])",
        "        ['observed'].agg(['mean', 'std', 'size']).rename(columns={'mean': 'scorecard', 'std': 'scorecard_sd', 'size': 'runs'})",
        "        .reset_index())",
        "never = sorted(set(map(tuple, record_keys[['organism', 'strategy']].values)) - set(map(tuple, card[['organism', 'strategy']].values)))",
        "print('record strategies whose recorded setting was never calibrated:', never)",
        "def record(part):",
        "    n = len(part); answered = int((~part.abstained.astype(bool)).sum())",
        "    right = int(part.correct.fillna(False).astype(bool).sum())",
        "    return pd.Series(dict(genes=n, answered=answered, right=right))",
        "rec = ledger.groupby(['organism', 'strategy', 'target', 'setting_key']).apply(record).reset_index()",
        "both = card.merge(rec, on=['organism', 'strategy', 'target', 'setting_key'], how='inner')",
        "both['track_record'] = np.where(both.strategy.isin(PRECISION), both.right / both.answered.replace(0, np.nan),",
        "                                both.right / both.genes)",
        "both['difference'] = both.track_record - both.scorecard",
        "# Sampling error of the record itself (binomial) and of the scorecard (spread over its runs).",
        "n = np.where(both.strategy.isin(PRECISION), both.answered, both.genes)",
        "se = np.sqrt(both.track_record * (1 - both.track_record) / np.clip(n, 1, None) + (both.scorecard_sd.fillna(0) ** 2) / both.runs)",
        "both['z'] = both.difference / np.where(se > 0, se, np.nan)",
        "print(len(both), 'organism x strategy x target x setting cells joined;',",
        "      int(card.shape[0] - both.shape[0]), 'scorecard cells have no matching record (other settings)')",
        "both.round(3).head(12)")

    nb.md("## 2 · Do they agree?",
          "Agreement is read three ways: the rank correlation across cells (does the record order the "
          "strategies the way the scorecard does), the typical absolute difference, and how many cells "
          "differ by more than sampling error explains (|z| > 3).")
    nb.code(
        "rho = spearmanr(both.track_record, both.scorecard, nan_policy='omit').statistic",
        "mad = float(both.difference.abs().median())",
        "outliers = both[both.z.abs() > 3].sort_values('z', key=abs, ascending=False)",
        "check = pd.DataFrame([dict(cells=len(both), spearman=round(float(rho), 3),",
        "                           median_abs_difference=round(mad, 3),",
        "                           mean_difference=round(float(both.difference.mean()), 3),",
        "                           beyond_sampling_error=len(outliers),",
        "                           agrees=bool(rho >= 0.9 and mad <= 0.05))])",
        "both.to_csv(os.path.join(OUT, 'cells.tsv'), sep='\\t', index=False)",
        "check.to_csv(os.path.join(OUT, 'check.tsv'), sep='\\t', index=False)",
        "check")
    nb.code(
        "outliers[['organism', 'strategy', 'target', 'setting_key', 'scorecard', 'track_record', 'difference', 'z', 'runs', 'genes']].round(3).head(20)")

    nb.md("### The one systematic difference: map neighbours counts only the genes it placed",
          "The recorded map-neighbours setting draws its map from a sample of 4,000 genes, so a hidden "
          "gene outside the sample is never placed and gets no answer. The record counts it as a hidden "
          "gene without a correct call; the scorecard counted only the genes on the map. Over the genes "
          "it placed, the record agrees with the scorecard. Both are true; the scorecard line for map "
          "neighbours should say *over the genes on the sampled map*, because it overstates how many "
          "hidden genes the strategy answers.")
    nb.code(
        "mn = both[both.strategy == 'map_neighbours'].copy()",
        "mn['placed_share'] = mn.answered / mn.genes",
        "mn['record_over_placed'] = mn.right / mn.answered",
        "aligned = both.copy()",
        "sel = aligned.strategy == 'map_neighbours'",
        "aligned.loc[sel, 'track_record'] = aligned.loc[sel, 'right'] / aligned.loc[sel, 'answered']",
        "aligned['difference'] = aligned.track_record - aligned.scorecard",
        "rho2 = spearmanr(aligned.track_record, aligned.scorecard).statistic",
        "mad2 = float(aligned.difference.abs().median())",
        "check2 = pd.DataFrame([dict(cells=len(aligned), spearman=round(float(rho2), 3), median_abs_difference=round(mad2, 3),",
        "                            mean_difference=round(float(aligned.difference.mean()), 3),",
        "                            within_005=int((aligned.difference.abs() <= 0.05).sum()),",
        "                            agrees=bool(rho2 >= 0.9 and mad2 <= 0.05))])",
        "check_all = pd.concat([check.assign(reading='as recorded'), check2.assign(reading='map neighbours over the genes it placed')], ignore_index=True)",
        "check_all.to_csv(os.path.join(OUT, 'check.tsv'), sep='\\t', index=False)",
        "display_cols = ['organism', 'target', 'scorecard', 'track_record', 'placed_share', 'record_over_placed']",
        "print(mn[display_cols].round(3).to_string(index=False))",
        "check_all")
    nb.md("### By strategy",
          "Pooled to the strategy level -- the level the requirement names -- with the cells it rests on.")
    nb.code(
        "by_strategy = both.groupby(['organism', 'strategy']).agg(cells=('target', 'size'),",
        "    scorecard=('scorecard', 'mean'), track_record=('track_record', 'mean')).reset_index()",
        "by_strategy['difference'] = by_strategy.track_record - by_strategy.scorecard",
        "by_strategy.to_csv(os.path.join(OUT, 'by_strategy.tsv'), sep='\\t', index=False)",
        "by_strategy.round(3)")
    nb.write(NB)


if __name__ == "__main__":
    main()
