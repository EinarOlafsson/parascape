"""The borrowed profile-inference strategies (instruction 77), scored by the harness, as an executed notebook.

Every strategy runs through `harness.run_question` unchanged -- nested cross-validation, guard_fit,
calibration, the baselines, permutation importance by family -- and touches a frozen partition once.
Because the holdout budget lets only one strategy definition look at a partition, each strategy gets
its own seed (its own fresh partition), recorded in `results/harness_holdout_log.jsonl`. Boosted
trees over the same evidence are run the same way as the reference.

    python scripts/build_borrowed_strategies_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_06"
#: DRYRUN=1 shuffles every label and logs to a scratch file: checks the code without spending a partition.
DRYRUN = os.environ.get("DRYRUN") == "1"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"borrowed_strategies_{DATE}") if os.environ.get("DRYRUN") != "1" else os.environ["DRYRUN_OUT"]
NB = os.path.join(REPO, "notebooks", f"borrowed_strategies_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("The borrowed inference strategies, each scored by the harness",
                          namespace={"REPO": REPO, "OUT": OUT, "DRYRUN": DRYRUN,
                                     "LOG": os.environ.get("DRYRUN_LOG", "") if DRYRUN else ""})
    nb.md(
        "**What this is.** Six inference techniques borrowed from large-scale profile inference "
        "(instruction 77), each written as a model the validation harness can score, and each scored "
        "on the two registry questions that have enough known answers for a frozen test.",
        "**The rule.** A strategy ships when its frozen-test result beats the best baseline; it may "
        "fill a cell only when its frozen-test lower bound clears 0.90. A cross-validated number alone "
        "never fills anything. A strategy that could not be tested is listed with the reason.",
        "**Why each has its own partition.** The harness lets one strategy definition look at a "
        "frozen partition once; a second, different strategy on the same partition would turn it into "
        "a tuning set and is refused. So every strategy below gets a fresh seed. The comparison "
        "between strategies therefore carries the noise of different partitions, which the "
        "cross-validated column (many folds) shows more steadily.",
        f"Built {_dt.date.today().isoformat()}.")

    nb.code(
        "import json, os, time",
        "from pathlib import Path",
        "import numpy as np, pandas as pd",
        "from parascape import harness as H, borrowed_models as B, question_registry as Q",
        "from parascape import datasets as ds, organisms as O, information_categories as IC",
        "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet')).set_index('gene_id')",
        "imap = pd.read_parquet(os.path.join(REPO, 'results/information_map_v0_2026_11_04/values.parquet'))",
        "def family_of(column):",
        "    found = ds.provenance(column, O.TOXOPLASMA)",
        "    return found.key if found is not None else 'unregistered'",
        "numeric = [c for c in nodes.columns if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]",
        "families = pd.Series({c: family_of(c) for c in numeric}, name='family')",
        "categories = {c: IC.categories_of(c) for c in numeric}",
        "B.set_families(families, categories)",
        "print(len(numeric), 'numeric columns in', families.nunique(), 'families')")

    nb.md("## 1 · The two questions with a frozen test",
          "* `required_without_glucose`: measured for 7,393 genes by the glucose-withdrawal screen "
          "(`fit_no_glucose <= -2`); every pooled fitness screen is held out as evidence.",
          "* `host_facing`: does the protein reach beyond the parasite plasma membrane into the host "
          "cell -- 585 genes with a curated answer in the information map (65 yes). Every place column "
          "and every proximity screen is held out.")
    nb.code(
        "QUESTIONS = {}",
        "q = Q.get('required_without_glucose')",
        "y = (nodes['fit_no_glucose'] <= -2.0)[nodes['fit_no_glucose'].notna()]",
        "QUESTIONS[q.id] = (q, y)",
        "q = Q.get('host_facing')",
        "hf = imap['host_facing'].dropna().astype(bool)",
        "QUESTIONS[q.id] = (q, hf[hf.index.isin(nodes.index)])",
        "INPUTS = {}",
        "for qid, (question, labels) in QUESTIONS.items():",
        "    full = nodes.loc[labels.index, numeric].astype('float64')",
        "    features = H.allowed_features(question, full, families)",
        "    groups = nodes['orthogroup'].reindex(labels.index)",
        "    groups = groups.where(groups.notna(), pd.Series({g: 'self:' + g for g in labels.index}))",
        "    if DRYRUN:",
        "        labels = pd.Series(np.random.default_rng(0).permutation(labels.to_numpy()), index=labels.index)",
        "    INPUTS[qid] = (question, labels, features, groups)",
        "    print(f'{qid}: {len(labels)} genes, {int(labels.sum())} yes; {features.shape[1]} allowed columns '",
        "          f'in {families.reindex(features.columns).nunique()} families; frozen test possible: '",
        "          f'{question.frozen_test_possible(len(labels))}')")

    nb.md("## 2 · The strategies",
          "Boosted trees are the reference. The six borrowed strategies are `parascape.borrowed_models`; "
          "the docstring of each says what it does. Selective prediction is not a separate model: it is "
          "how the harness grades every one of them (coverage at the 0.90 bar).")
    nb.code(
        "STRATEGIES = [H.Strategy.of('boosted trees (reference)', 'hist_gbt', max_iter=150),",
        "              H.Strategy.of('positive-unlabelled bagging', 'pu_bagging', n_bags=10),",
        "              H.Strategy.of('inverse-propensity weighting', 'propensity_weighted'),",
        "              H.Strategy.of('low-rank completion', 'low_rank_completion', rank=30),",
        "              H.Strategy.of('per-family stacking', 'family_stacking'),",
        "              H.Strategy.of('lookalike ranking', 'lookalike', k=10),",
        "              H.Strategy.of('network propagation', 'graph_propagation', k=15)]",
        "SEED0 = 10      # seeds 0-2 of required_without_glucose are spent; one fresh seed per strategy",
        "pd.DataFrame([dict(strategy=s.name, model=s.model, seed=SEED0 + i, hash=s.identity[:12]) for i, s in enumerate(STRATEGIES)])")

    nb.md("## 3 · The runs",
          "One `run_question` per strategy and question: nested 3 x 3 cross-validation inside train + "
          "tune, refit, isotonic calibration on the calibration partition, then the frozen partition "
          "once. Every use is appended to the holdout log.")
    nb.code(
        "rows, importance = [], []",
        "for qid, (question, labels, features, groups) in INPUTS.items():",
        "    for i, strategy in enumerate(STRATEGIES):",
        "        t0 = time.time()",
        "        r = H.run_question(question, labels=labels, features=features, families=families,",
        "                           groups=groups, strategies=[strategy], seed=SEED0 + i,",
        "                           outer_folds=3, inner_folds=3, importance_repeats=2,",
        "                           log_path=LOG or os.path.join(REPO, H.HOLDOUT_LOG))",
        "        cv, ft = r.cross_validated, r.frozen_test or {}",
        "        base = r.baselines[r.baselines.partition == 'frozen_test'] if len(r.baselines) else r.baselines",
        "        best = base.sort_values('precision_lower90', ascending=False).iloc[0] if len(base) else None",
        "        rows.append(dict(question=qid, strategy=strategy.name, seed=SEED0 + i,",
        "            cv_auc=round(cv.get('roc_auc', np.nan), 3), cv_ap=round(cv.get('average_precision', np.nan), 3),",
        "            cv_coverage_at_bar=round(cv.get('coverage_at_bar', np.nan), 4),",
        "            test_n=ft.get('n'), test_auc=round(ft.get('roc_auc', np.nan), 3),",
        "            test_ap=round(ft.get('average_precision', np.nan), 3),",
        "            test_calls=ft.get('calls'), test_precision=ft.get('precision'),",
        "            test_lower90=round(ft.get('precision_lower90', np.nan), 3) if ft else np.nan,",
        "            test_coverage_at_bar=round(ft.get('coverage_at_bar', np.nan), 4) if ft else np.nan,",
        "            best_baseline=None if best is None else best['baseline'] if 'baseline' in best else best.iloc[0],",
        "            best_baseline_lower90=np.nan if best is None else round(float(best['precision_lower90']), 3),",
        "            calibration_error=round(r.calibration.get('calibration_error', np.nan), 3),",
        "            may_fill=H.may_fill(r), flags=';'.join(r.flags), minutes=round((time.time() - t0) / 60, 1)))",
        "        if r.importance is not None and len(r.importance):",
        "            importance.append(r.importance.head(5).assign(question=qid, strategy=strategy.name))",
        "        print(qid, '|', strategy.name, '|', rows[-1]['cv_auc'], rows[-1]['test_auc'], rows[-1]['minutes'], 'min')",
        "scores = pd.DataFrame(rows)",
        "scores.to_csv(os.path.join(OUT, 'scores.tsv'), sep='\\t', index=False)",
        "pd.concat(importance).to_csv(os.path.join(OUT, 'importance_top5.tsv'), sep='\\t', index=False)",
        "scores")

    nb.md("## 4 · Which ship",
          "A strategy **ships** on a question when its frozen-test average precision beats the best "
          "baseline there and is not worse than the reference by more than 0.02 on the frozen "
          "test **and** on the cross-validated estimate. It **may fill** only when the frozen-test lower "
          "bound clears 0.90 -- `may_fill` above. Anything else is listed, with its numbers, as not "
          "better than what is already there.")
    nb.code(
        "ref = scores[scores.strategy == 'boosted trees (reference)'].set_index('question')",
        "def verdict(row):",
        "    if not np.isfinite(row.test_ap):",
        "        return 'untested: no frozen test'",
        "    beats_baseline = row.test_lower90 > row.best_baseline_lower90 if np.isfinite(row.best_baseline_lower90) else True",
        "    r = ref.loc[row.question]",
        "    close = row.test_ap >= r.test_ap - 0.02 and row.cv_ap >= r.cv_ap - 0.02",
        "    if row.strategy.startswith('boosted'):",
        "        return 'reference'",
        "    return 'ships' if beats_baseline and close else ('beats baseline, worse than reference' if beats_baseline else 'does not beat the baseline')",
        "scores['verdict'] = scores.apply(verdict, axis=1)",
        "untested = pd.DataFrame([dict(strategy=k, verdict='untested', reason=v) for k, v in B.UNTESTED.items()])",
        "scores.to_csv(os.path.join(OUT, 'scores.tsv'), sep='\\t', index=False)",
        "untested.to_csv(os.path.join(OUT, 'untested.tsv'), sep='\\t', index=False)",
        "scores[['question', 'strategy', 'cv_ap', 'test_ap', 'test_lower90', 'best_baseline_lower90', 'may_fill', 'verdict']]")
    nb.code("untested")
    if not DRYRUN:
        nb.write(NB)


if __name__ == "__main__":
    main()
