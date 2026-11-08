"""The first filled cells of the information map (instruction 72), as an executed notebook.

The questions of the catalogue that MAY fill -- their frozen-test lower bound on precision reached 0.90 -- and whose
labels are clean biology (three may-fill questions with bookkeeping labels are left out and named) are re-run with
the identical strategy and split (the harness logs a repeat, not a second look with anything changed), and their
deployed, calibrated fit is applied to the genes with no label. A cell is written only where the answer is yes and
the calibrated probability reaches the threshold chosen on train + tune. Each filled column holds only fills, with
the calibrated probability and the source `filled:hist_gbt@2026-11-08 <question>`; evidence rows say the same.
Filled cells are never training labels.

    python scripts/build_information_map_fills_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"information_map_fills_{DATE}")
NB = os.path.join(REPO, "notebooks", f"information_map_fills_{DATE}.ipynb")
MAP = os.path.join(REPO, "parascape", "data", "information_map")
NOT_BIOLOGY = ("fit_density_dependence_log10padj", "crispr_gra17_candidate", "n_holes")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Information map: the first filled cells",
                          namespace={"REPO": REPO, "OUT": OUT, "MAP": MAP, "NOT_BIOLOGY": NOT_BIOLOGY})
    nb.md("**What this writes (instruction 72).** Filled cells, only from questions whose frozen test cleared the "
          "90 % precision bar at its lower bound, only where the calibrated probability of yes reaches the threshold "
          "that produced that precision. Everything else the methods say stays a ranking, outside the map.",
          f"Left out although they may fill, because the label is not biology: {', '.join(NOT_BIOLOGY)}.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, numpy as np, pandas as pd",
            "from parascape import datasets as ds, harness as H, harness_answer as HA, information_categories as IC, organisms as O",
            "from parascape.question_registry import Question",
            "S = pd.read_csv(os.path.join(REPO, 'results/question_scoreboard_2026_11_08/scores.tsv'), sep='\\t')",
            "chosen = S[S.may_fill & ~S.column.isin(NOT_BIOLOGY)]",
            "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet')).set_index('gene_id')",
            "numeric = [c for c in nodes.columns if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]",
            "fam = lambda c: (ds.provenance(c, O.TOXOPLASMA).key if ds.provenance(c, O.TOXOPLASMA) else 'unregistered')",
            "families = pd.Series({c: fam(c) for c in numeric})",
            "groups = nodes['orthogroup'].astype(object)",
            "groups = groups.where(groups.notna(), pd.Series({g: 'self:' + g for g in nodes.index}))",
            "cat = pd.read_csv(os.path.join(REPO, 'results/question_catalogue_2026_11_08/catalogue.tsv'), sep='\\t').set_index('column')",
            "chosen[['question', 'column', 'category', 'frozen_auc', 'precision_lower90']]")
    nb.md("## 1 · Does the bar hold on other partitions?",
          "Einar (2026-11-08): a question fills only if its frozen-test lower bound reaches 0.90 on **at least 2 "
          "of 3** partitions -- the original (seed 0) and two fresh ones (seeds 4 and 5, unused by any other run "
          "of these questions), the same strategy each time. One frozen test at the edge of the bar can pass by "
          "the luck of the split: the nucleus question cleared it on seeds 0 and 3 and missed on seed 2 (0.885).")
    nb.code("SEEDS = (0, 4, 5)",
            "def as_question(r):",
            "    row = cat.loc[r.column]",
            "    for suffix, wording, labels in HA.binary_labels(r.column, row.kind, row.category, nodes[r.column]):",
            "        qid = f'cat_{HA._safe(r.column)}_{suffix}'",
            "        if qid == r.question:",
            "            family = fam(r.column)",
            "            q = Question(id=qid, text=f'{row.question} -- yes when {wording}', answer_type='boolean', question_class='B',",
            "                         label_source=f'nodes.parquet:{r.column}', label_family=family, evidence_forbidden=(family,),",
            "                         label_tier='observed-screen', label_category=row.category, status='scored')",
            "            return q, labels",
            "    raise KeyError(r.question)",
            "stab = []",
            "for r in chosen.itertuples(index=False):",
            "    q, labels = as_question(r)",
            "    feats = H.allowed_features(q, nodes.loc[labels.index, numeric].astype('float64'), families)",
            "    for seed in SEEDS:",
            "        res = H.run_question(q, labels=labels.astype(bool), features=feats, families=families,",
            "                             groups=groups.reindex(labels.index),",
            "                             strategies=[H.Strategy.of('boosted trees', 'hist_gbt', max_iter=150)], seed=seed,",
            "                             importance_repeats=1, log_path=H.HOLDOUT_LOG)",
            "        stab.append(dict(question=r.question, seed=seed, frozen_auc=round(res.frozen_test['roc_auc'], 3),",
            "                         lower90=round(res.frozen_test['precision_lower90'], 3), may_fill=bool(res.may_fill)))",
            "St = pd.DataFrame(stab)",
            "St.to_csv(os.path.join(OUT, 'partition_stability.tsv'), sep='\\t', index=False)",
            "passes = St.groupby('question').may_fill.sum()",
            "stable = set(passes[passes >= 2].index)",
            "print('fill on >= 2 of 3 partitions:', sorted(stable))",
            "print('withheld (bar not met on 2 of 3):', sorted(set(passes.index) - stable))",
            "chosen = chosen[chosen.question.isin(stable)]",
            "St.pivot(index='question', columns='seed', values='lower90')")
    nb.md("## 2 · Re-run, deploy, keep the fills")
    nb.code("fills, summary = [], []",
            "TAG = 'filled:hist_gbt@2026-11-08'",
            "for r in chosen.itertuples(index=False):",
            "    row = cat.loc[r.column]",
            "    for suffix, wording, labels in HA.binary_labels(r.column, row.kind, row.category, nodes[r.column]):",
            "        qid = f'cat_{HA._safe(r.column)}_{suffix}'",
            "        if qid != r.question:",
            "            continue",
            "        family = fam(r.column)",
            "        q = Question(id=qid, text=f'{row.question} -- yes when {wording}', answer_type='boolean', question_class='B',",
            "                     label_source=f'nodes.parquet:{r.column}', label_family=family, evidence_forbidden=(family,),",
            "                     label_tier='observed-screen', label_category=row.category, status='scored')",
            "        full = nodes.loc[labels.index, numeric].astype('float64')",
            "        feats = H.allowed_features(q, full, families)",
            "        res = H.run_question(q, labels=labels.astype(bool), features=feats, families=families,",
            "                             groups=groups.reindex(labels.index),",
            "                             strategies=[H.Strategy.of('boosted trees', 'hist_gbt', max_iter=150)], seed=0,",
            "                             importance_repeats=1, log_path=H.HOLDOUT_LOG)",
            "        unl = nodes.index[nodes[r.column].isna()]",
            "        if not len(unl):   # every gene is measured: nothing to fill",
            "            summary.append(dict(column='filled_' + qid[len('cat_'):], question=q.text, unlabelled=0, called=0, filled=0,",
            "                                may_fill=res.may_fill, frozen_lower90=round(res.frozen_test['precision_lower90'], 3), flags=', '.join(res.flags)))",
            "            continue",
            "        dep = H.deploy(res, nodes.loc[unl, feats.columns].astype('float64'))",
            "        hit = dep[dep.fills]",
            "        # Applicability: the 90 % precision was measured on labelled genes. If the genes without a label are",
            "        # filled far more often than the label occurs (> 3x its prevalence), they are a different population",
            "        # (missing not at random) and the measured precision does not carry over: no fills, reason recorded.",
            "        prevalence = float(labels.mean()); rate = len(hit) / len(unl)",
            "        applicable = rate <= 3 * prevalence",
            "        if not applicable:",
            "            hit = hit.iloc[:0]",
            "        col = 'filled_' + qid[len('cat_'):]",
            "        for g, p in hit.calibrated.items():",
            "            fills.append(dict(gene_id=g, trait=col, value='True', probability=round(float(p), 4),",
            "                              source_tier=f'{TAG} {qid}', reference='results/question_scoreboard_2026_11_08',",
            "                              dataset='question catalogue (instruction 73)'))",
            "        summary.append(dict(column=col, question=q.text, unlabelled=len(unl), called=int(dep.called.sum()),",
            "                            filled=len(hit), fill_rate=round(rate, 3), prevalence=round(prevalence, 3), applicable=applicable,",
            "                            may_fill=res.may_fill, frozen_lower90=round(res.frozen_test['precision_lower90'], 3),",
            "                            flags=', '.join(res.flags)))",
            "F = pd.DataFrame(fills); Sm = pd.DataFrame(summary)",
            "F.to_csv(os.path.join(OUT, 'fills.tsv'), sep='\\t', index=False); Sm.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\\t', index=False)",
            "Sm")
    nb.md("## 3 · Into the published map",
          "One column per question, holding only its fills; probability = calibrated probability; source = the tag "
          "and the question. Evidence rows appended in the same form. Observed columns are not touched. Every "
          "earlier `filled_` column is removed first, so a question withheld by the partition rule leaves no cells.")
    nb.code("values = pd.read_parquet(os.path.join(MAP, 'values.parquet')); prob = pd.read_parquet(os.path.join(MAP, 'probability.parquet'))",
            "source = pd.read_parquet(os.path.join(MAP, 'source.parquet')); evidence = pd.read_parquet(os.path.join(MAP, 'evidence.parquet'))",
            "before = {k: v.shape for k, v in dict(values=values, probability=prob, source=source, evidence=evidence).items()}",
            "for frame in (values, prob, source):",
            "    frame.drop(columns=[c for c in frame.columns if c.startswith('filled_')], inplace=True)",
            "for col in Sm[Sm.filled > 0].column:",
            "    rows = F[F.trait == col].set_index('gene_id')",
            "    values[col] = rows.value.map(lambda v: v == 'True').reindex(values.index).astype('object')",
            "    prob[col] = rows.probability.reindex(values.index)",
            "    source[col] = rows.source_tier.reindex(values.index)",
            "evidence = evidence[~evidence.trait.astype(str).str.startswith('filled_')]",
            "evidence = pd.concat([evidence, F[evidence.columns.intersection(F.columns)]], ignore_index=True)",
            "values.to_parquet(os.path.join(MAP, 'values.parquet')); prob.to_parquet(os.path.join(MAP, 'probability.parquet'))",
            "source.to_parquet(os.path.join(MAP, 'source.parquet')); evidence.to_parquet(os.path.join(MAP, 'evidence.parquet'), index=False)",
            "desc_path = os.path.join(MAP, 'column_descriptions.json')",
            "desc = json.load(open(desc_path)) if os.path.exists(desc_path) else {}",
            "desc = {k: v for k, v in desc.items() if not k.startswith('filled_')}",
            "desc.update({r.column: f'FILLED (predicted, never observed): {r.question}' for r in Sm[Sm.filled > 0].itertuples(index=False)})",
            "json.dump(desc, open(desc_path, 'w'), indent=1)",
            "after = {k: v.shape for k, v in dict(values=values, probability=prob, source=source, evidence=evidence).items()}",
            "from parascape import information_map as IM",
            "m = IM.load_published(MAP)",
            "print(before, '->', after)",
            "print(m.banner())",
            "m.predicted_cells().groupby('column').size().to_frame('filled cells')")
    nb.write(NB)


if __name__ == "__main__":
    main()
