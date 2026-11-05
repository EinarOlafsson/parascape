"""Instruction 73 phase 1: the question registry and the validation harness, demonstrated.

Builds notebooks/harness_demo_2026_11_05.ipynb (executed) and the tables under
results/harness_demo_2026_11_05/. It runs ONE real question end to end through
`parascape.harness` -- `required_without_glucose`, class A, 7,393 labelled genes -- and reports the
nested cross-validated estimate and the frozen-test result side by side, which is the pair of
numbers the three earlier scoring runs did not have. A second question, `screen_egress_phenotype`
(318 genes), is run to show the other branch: below its own label threshold it keeps nested
cross-validation only, is flagged `no_frozen_test`, and may not fill cells.

The host-facing, function and partner questions are deliberately NOT re-run here; that is a later
phase of instruction 73.
"""
import os
import sys
from pathlib import Path
from textwrap import dedent

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT))
from notebook_runner import ExecutedNotebook  # noqa: E402

nb = ExecutedNotebook('The question registry and the validation harness, demonstrated on one real question')


def code(src: str):
    """Run one dedented code cell in the notebook's namespace."""
    return nb.code(dedent(src).strip('\n'))


nb.md(
    '# Phase 1 of instruction 73 — the registry, the harness, and both numbers',
    'Einar, 2026-11-05: *"Whenever predictions are made there should always be a holdout set with '
    'known answers that the model predicts, and the quality of the strategy measured that way. Did '
    'you bake in a mechanism like this, or cross validation? Preferably both."*',
    '**The honest answer that prompted this work: cross-validation yes, a frozen holdout no.** All '
    'three scoring runs so far (host-facing v1 and v2, function, partners) used 5-fold '
    'cross-validation grouped by orthogroup. Out-of-fold predictions on known labels are a held-out '
    'measurement, so those numbers are not circular — but **no partition was held back from '
    'everything**. The same folds that produced the reported number also chose the feature set (14 '
    'variants for host-facing, 10 for function). Picking the best by cross-validation and then '
    'quoting its cross-validated score is optimistic, so 0.727, 0.617 and 0.891 are **upper bounds** '
    'on what a frozen test would give. And `parascape/splits.py` already provided exactly the right '
    'primitive — `ROLES`, `make_split`, `guard_fit`, `FIT_ACCESS` — which the runs did not call.',
    'This notebook is the architecture that fixes it, and the gate on phase 1: one question run end '
    'to end through the harness, with the frozen-test number beside the cross-validated one.',
    '| What the harness enforces | How |\n|---|---|\n'
    '| Four partitions, whole orthogroups to one side, seed recorded | `splits.make_split`, '
    '`splits.cold_node_pairs` for pair questions |\n'
    '| Nested cross-validation inside train + tune, where all searching happens | `GroupKFold` '
    'outer and inner, strategies compared per outer fold |\n'
    '| A leak **raises** instead of being noticed in review | `StageGuard` calls '
    '`split.guard_fit` at feature selection, model fit, setting selection and calibration, and '
    'records every call |\n'
    '| The frozen partition is touched **once** per strategy | `results/harness_holdout_log.jsonl`, '
    'keyed by question id, split identity and a hash of the strategy definition |\n'
    '| Both numbers, never one | `HarnessResult.summary()` |\n'
    '| A small-label question says so rather than being judged per run | '
    '`Question.min_labels_for_frozen_test` → the `no_frozen_test` flag |\n'
    '| Filling a cell requires the **frozen-test** lower bound to clear 0.90 | `harness.may_fill` |')

# --------------------------------------------------------------------------- 1 the registry
nb.md(
    '## 1 · The question registry',
    '`parascape/question_registry.py` — named apart from the existing `parascape/questions.py`, '
    'which is instruction 45\'s catalogue of what the *table* can be asked. One row per question: '
    'its id, the question as a biologist would ask '
    'it, the answer type, the label source and its tier, which evidence families are allowed and '
    'which must be **held out**, the class A–D of instruction 73, the label count below which a '
    'frozen partition is unaffordable, and the status. The three admission tests — answerable, '
    'non-redundant, useful — are fields too, so the catalogue can be audited rather than argued '
    'about.',
    'The registry has one structural guarantee, and it is the lesson of the first host-facing run: '
    '**a question must hold out the family its own labels came from.** Registering a question whose '
    '`label_family` is not in `evidence_forbidden` raises.')
code("""
    import os
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[v] = '1'
    import json, warnings
    from pathlib import Path
    import numpy as np, pandas as pd
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    import sys
    sys.path.insert(0, os.getcwd())
    from parascape import question_registry as Q, harness as H, splits as S, datasets as ds
    from parascape import organisms as O

    ROOT = Path('.').resolve()
    OUT = ROOT / 'results/harness_demo_2026_11_05'; OUT.mkdir(parents=True, exist_ok=True)
    catalogue = Q.catalogue()
    catalogue.to_csv(OUT / 'question_catalogue.tsv', sep='\\t', index=False)
    print('registered questions:', len(catalogue))
    print(catalogue[['id', 'question_class', 'answer_type', 'label_tier',
                     'min_labels_for_frozen_test', 'status', 'admitted']].to_string(index=False))
""")
code("""
    print('What each question may NOT be answered from (the label family is always in the list):')
    print()
    for row in catalogue.itertuples():
        print('%-34s class %s  labels: %s' % (row.id, row.question_class, row.label_source))
        print('%-34s forbidden: %s' % ('', row.evidence_forbidden))
    print()
    print('A question that does not hold out its own label source is refused:')
    try:
        Q.Question(id='bad', text='Is it so?', answer_type='boolean', question_class='B',
                   label_source='screen', label_family='splitcas9_imaging_screen',
                   evidence_forbidden=('pubmed',))
    except ValueError as error:
        print('  ValueError:', error)
    print('An unknown answer type or class is refused:')
    for kwargs in (dict(answer_type='probability', question_class='B'),
                   dict(answer_type='boolean', question_class='E')):
        try:
            Q.Question(id='bad', text='?', label_source='s', label_family='f',
                       evidence_forbidden=('f',), **kwargs)
        except ValueError as error:
            print('  ValueError:', error)
""")

# --------------------------------------------------------------------------- 2 the question
nb.md(
    '## 2 · The question this is demonstrated on, and the caveat that ships with it',
    '`required_without_glucose` — *"Which genes does the parasite need in order to grow when glucose '
    'is withdrawn?"* Class **A**: the answer is **measured** for 7,393 genes '
    '(`fit_no_glucose`, the glucose-limitation CRISPR screen), which is why class A exists as a '
    'separate class at all.',
    '**Two caveats, stated before the numbers.**',
    '1. **Rung 0 of the strategy ladder carries the measurement.** For a gene the screen measured, '
    'the right answer is the screen\'s own score with the screen\'s own uncertainty — no model. The '
    'prediction problem is only the ~750 genes the screen missed. What is run below is therefore the '
    'harness being demonstrated on a question with enough labels to show a frozen test, not a claim '
    'that these genes needed predicting.\n'
    '2. **"Important for X" is not "involved in X."** A gene required for growth without glucose may '
    'be a glycolytic enzyme or simply essential everywhere. The conditional form of this question '
    'needs the screen\'s published per-gene statistics, which `nodes.parquet` does not carry — the '
    'same blocker recorded for the Wang screen. The label used here is the absolute one, '
    '`fit_no_glucose <= -2.0`, and it is named in the registry as such.',
    'Every CRISPR fitness family is held out, not only the column the labels came from: '
    '`fit_complete_medium_2025` would otherwise answer `fit_no_glucose` and nothing would have been '
    'learnt.')
code("""
    question = Q.get('required_without_glucose')
    nodes = pd.read_parquet(ROOT / 'parascape/data/nodes.parquet').set_index('gene_id')

    labels = (nodes['fit_no_glucose'] <= -2.0)[nodes['fit_no_glucose'].notna()]
    groups = nodes['orthogroup'].reindex(labels.index)
    groups = groups.where(groups.notna(), pd.Series({g: 'self:' + g for g in labels.index}))
    print(question.text)
    print()
    print('label source :', question.label_source)
    print('class        :', question.question_class, '| answer type:', question.answer_type,
          '| tier:', question.label_tier)
    print('labelled genes:', len(labels), '| positives:', int(labels.sum()),
          '(%.1f %%)' % (100 * labels.mean()))
    print('genes the screen did not measure (the real prediction problem):',
          int(nodes['fit_no_glucose'].isna().sum()))
    print('distinct orthogroups among the labelled genes:', groups.nunique())
    print('frozen test affordable? ', question.frozen_test_possible(len(labels)),
          '(threshold', question.min_labels_for_frozen_test, 'labels)')
""")

# --------------------------------------------------------------------------- 3 evidence
nb.md(
    '## 3 · Evidence, by family, with the forbidden families removed — and the harness refusing them',
    'Features are the numeric and boolean columns of `nodes.parquet`, each tagged with the dataset '
    'family that produced it (`datasets.provenance`). A family is held out as a unit, because 36 '
    'expression columns from one experiment are one piece of evidence and not 36.',
    '`harness.run_question` **raises** if a forbidden family is present rather than dropping it '
    'quietly: a leak should be a failure, not a silent correction. `harness.allowed_features` is how '
    'a caller asks for the filtering on purpose, which is what is done here.',
    'The forbidden list holds out **every** pooled CRISPR fitness screen, enumerated by a rule: any '
    'family contributing a column named `crispr*` or `fit_*`, plus the curated enteric fitness set. '
    'Four of those families do not have "crispr" in their name, and the harness found each of them '
    'itself, over two runs of this very notebook. Section 6 is that story and what the holdout '
    'budget charged for it.')
code("""
    def family_of(column):
        found = ds.provenance(column, O.TOXOPLASMA)
        return found.key if found is not None else 'unregistered'

    numeric = [c for c in nodes.columns
               if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]
    families = pd.Series({c: family_of(c) for c in numeric}, name='family')
    full = nodes.loc[labels.index, numeric].astype('float64')

    features = H.allowed_features(question, full, families)
    dropped = sorted(set(full.columns) - set(features.columns))
    print('numeric/boolean columns in nodes.parquet:', len(numeric))
    print('HELD OUT as forbidden families (%d columns):' % len(dropped))
    for fam in sorted(question.evidence_forbidden):
        members = [c for c in dropped if families[c] == fam]
        if members:
            print('   %-30s %s' % (fam, ', '.join(members)))
    print()
    print('features used:', features.shape[1], 'in', families.reindex(features.columns).nunique(),
          'families')
    print(families.reindex(features.columns).value_counts().head(12).to_string())
""")
code("""
    print('The harness refuses the forbidden families rather than dropping them quietly:')
    try:
        H.run_question(question, labels=labels, features=full, families=families, groups=groups,
                       strategies=[H.Strategy.of('all', 'hist_gbt')])
    except ValueError as error:
        print('  ValueError:', str(error)[:220], '...')
    print()
    print('And it refuses a strategy that asks for one:')
    try:
        H.run_question(question, labels=labels, features=features, families=families, groups=groups,
                       strategies=[H.Strategy.of('leaky', 'hist_gbt',
                                                 families=('crispr_invitro',))])
    except ValueError as error:
        print('  ValueError:', error)
""")

# --------------------------------------------------------------------------- 4 partitions
nb.md(
    '## 4 · Four partitions, and the guard that makes a leak raise',
    '`splits.make_split`, grouped by OrthoMCL orthogroup so paralogs never straddle a boundary, seed '
    '0, fractions 0.55 / 0.15 / 0.15 / 0.15. `FIT_ACCESS` is the permission table the guard enforces: '
    'model fitting and feature selection may see only `train`; setting selection and refitting may '
    'see `train` + `tune`; calibration may see only `calibration`; **nothing may see `test`**.',
    'Nested cross-validation inside train + tune therefore fits under the `refit` stage rather than '
    '`model_fit` — the distinction `splits.FIT_ACCESS` already draws, and the one place the existing '
    'API had to be read carefully rather than guessed at.')
code("""
    SEED = 2          # the seed-0 and seed-1 partitions are already spent -- see section 6
    split, pair_partition = H.partition_question(question, labels.index, groups.loc[labels.index],
                                                 seed=SEED)
    sizes = {role: len(split.entities(role)) for role in S.ROLES}
    print('split identity:', split.identity)
    print('seed', SEED, '| grouping:', split.group_kind, '| access regime:', split.feature_access)
    print('partition sizes:', sizes, '| total', sum(sizes.values()))
    print('positives per partition:',
          {role: int(labels.loc[list(split.entities(role))].sum()) for role in S.ROLES})
    crossing = sum(1 for g, sub in groups.groupby(groups)
                   if len({r.role for r in split.assignments if r.entity in set(sub.index)}) > 1)
    print('orthogroups straddling a partition boundary:', crossing)

    guard = H.StageGuard(split)
    test = list(split.entities('test'))
    for stage in ('model_fit', 'feature_selection', 'setting_selection', 'refit', 'calibration'):
        try:
            guard.guard(stage, test)
            print('  LEAK NOT CAUGHT:', stage)
        except ValueError as error:
            print('  %-18s reaching into `test` raises: %s' % (stage, str(error)[:60]))
    print('refused calls recorded:', len(guard.calls), '(a refused call is not a recorded one)')
""")

# --------------------------------------------------------------------------- 5 the run
nb.md(
    '## 5 · The run: nested cross-validation, then the frozen partition once',
    'Three strategies are compared, and the comparison happens **inside train + tune**: boosted '
    'trees over every allowed family, a logistic model over the same (the linear comparison), and '
    'boosted trees restricted to the expression and proteome families. The outer loop gives the '
    'cross-validated estimate; the strategy chosen in the majority of outer folds is refitted on '
    'train + tune, calibrated on `calibration` with `claims.CertaintyModel` (isotonic, cross-fitted '
    'error), and scored **once** on `test`.',
    'The fill threshold is chosen on the train + tune out-of-fold predictions — the lowest calibrated '
    'probability whose Clopper–Pearson lower 90 % bound reaches the bar — and then applied unchanged '
    'to the frozen partition. Choosing it on the test partition is the exact mistake the whole '
    'mechanism exists to prevent.')
code("""
    expression_like = sorted({f for f in families.reindex(features.columns).unique()
                              if f.startswith('gse') or 'proteome' in f or 'expression' in f})
    strategies = [H.Strategy.of('boosted trees, every allowed family', 'hist_gbt', max_iter=150),
                  H.Strategy.of('logistic, every allowed family', 'logistic'),
                  H.Strategy.of('boosted trees, expression and proteome only', 'hist_gbt',
                                families=expression_like, max_iter=150)]
    for s in strategies:
        print('%-46s %s  %s' % (s.name, s.identity[:12], json.dumps(s.definition['settings'])))
    print('expression/proteome families in the restricted strategy:', len(expression_like))
""")
code("""
    result = H.run_question(question, labels=labels, features=features, families=families,
                            groups=groups, strategies=strategies, seed=SEED,
                            outer_folds=3, inner_folds=3)
    print('strategy chosen per outer fold:', result.chosen_per_outer_fold)
    print('deployed strategy :', result.strategy, '| hash', result.strategy_hash[:12])
    print('fill threshold    : %.4f (calibrated probability, chosen on train+tune)' % result.threshold)
    print('flags             :', result.flags or '(none)')
    print()
    print('inner-fold average precision per strategy, per outer fold:')
    print(pd.DataFrame(result.inner_scores).round(4).to_string())
""")

nb.md(
    '### The two numbers, side by side',
    'This is what the harness exists to produce. The left column is the nested cross-validated '
    'estimate — how the strategy compares, over many folds. The right column is the frozen-test '
    'result — what it achieves, once. **Filling a cell depends on the right column only.**')
code("""
    summary = result.summary()
    summary.to_csv(OUT / 'summary.tsv', sep='\\t', index=False)
    print(summary.set_index('partition').T.to_string())
    print()
    print('may this question fill cells?', H.may_fill(result),
          '(frozen-test lower bound %.4f against a bar of %.2f)'
          % (result.frozen_test['precision_lower90'], result.bar))
    summary.set_index('partition').T
""")

nb.md(
    '### Calibration, baselines and permutation importance by family',
    'Calibration is fitted on the `calibration` partition **only**, reusing '
    '`claims.CertaintyModel` and its isotonic fit rather than a fourth copy of it, with the error '
    'cross-fitted inside that partition. Nothing is reported except as improvement over the best '
    'baseline. Permutation importance is measured by **family** on the calibration partition — a '
    'partition the model was not fitted on, and not the frozen one, which is spent.')
code("""
    print('calibration:', result.calibration)
    print()
    print('BASELINES (nothing is reported except as improvement over the best of these):')
    print(result.baselines.round(4).to_string(index=False))
    result.baselines.to_csv(OUT / 'baselines.tsv', sep='\\t', index=False)
""")
code("""
    result.importance.to_csv(OUT / 'importance_by_family.tsv', sep='\\t', index=False)
    print('PERMUTATION IMPORTANCE BY FAMILY (average precision lost when the whole family is shuffled)')
    print('baseline average precision on the calibration partition: %.4f'
          % result.importance['baseline_average_precision'].iloc[0])
    print(result.importance.head(15)[['family', 'n_features', 'importance_mean',
                                      'importance_sd']].round(4).to_string(index=False))
    result.importance.head(10)[['family', 'n_features', 'importance_mean', 'importance_sd']].round(4)
""")

# --------------------------------------------------------------------------- 6 the budget
nb.md(
    '## 6 · The holdout budget, on disk — and the first thing it refused',
    '`results/harness_holdout_log.jsonl` is persistent and is the record of every look at a frozen '
    'partition. The key is the question id, the split identity, a hash of the strategy definition '
    '**and a hash of the evidence**: the question\'s own forbidden-family list plus the exact feature '
    'set. A second look by a changed strategy *or* changed evidence is refused — it would make the '
    'frozen test a third cross-validation — and the question is marked as needing a fresh partition. '
    'Re-running the identical strategy on the identical evidence is allowed and recorded as a repeat.',
    '**This has already happened twice, which is why the run above used seed 2.** The first run of '
    'this question (seed 0, logged below) held out the seven families whose names begin with '
    '`crispr` plus the curated fitness sets. The harness\'s own permutation importance then put '
    '`gra17_synthlethal` at the top and its best-single-feature baseline chose '
    '`crispr_invivo_platform_lfc` — two further **pooled CRISPR screens** under family names that do '
    'not contain "crispr", answering the question with the same kind of measurement. They were added '
    'to `evidence_forbidden`, the question moved to a fresh partition (seed 1), and the same two '
    'diagnostics then found a *third*: `gra12`, whose columns are `crispr_gra12s1_l2fc_invitro` and '
    'its siblings.',
    'The forbidden list is now enumerated by a rule instead of by eye — every family contributing a '
    'column named `crispr*` or `fit_*`, plus the curated enteric measurements — and the question is '
    'on its third partition (seed 2). Two partitions are spent; the honest cost of a corrected '
    'definition is a fresh partition, and this is what paying it looks like. It is also exactly why '
    'the budget hashes the **evidence** as well as the strategy: the model, its settings and its '
    'version never changed, so the strategy hash alone would have called all three runs a repeat.')
code("""
    log = H.holdout_log(root=ROOT)
    mine = log[log.question_id == question.id]
    print('every recorded look at a frozen partition of this question:')
    print(mine[['seed', 'strategy', 'strategy_hash', 'evidence_hash', 'n_test', 'calls',
                'precision', 'precision_lower90', 'when']].to_string(index=False))
    print()
    for row in mine.itertuples():
        print('seed %s -> split %s' % (row.seed, row.split_identity[:16]))
    print()
    print('this run:  strategy %s | evidence %s' % (result.strategy_hash[:12],
                                                    result.evidence_hash[:12]))
""")
code("""
    print('Asking again for a SPENT partition (seed 1), with the corrected forbidden list:')
    spent = H.run_question(question, labels=labels, features=features, families=families,
                           groups=groups, strategies=strategies, seed=1,
                           outer_folds=3, inner_folds=3, record=False)
    print('  flags      :', spent.flags)
    print('  frozen test:', spent.frozen_test)
    print('  refusal    :', spent.refusal)
    print('  may fill   :', H.may_fill(spent))
    print('  the cross-validated estimate is still produced, and is still only an estimate: '
          'lower bound %.4f' % spent.cross_validated['precision_lower90'])
""")
code("""
    print('And a CHANGED STRATEGY asking for the fresh seed-2 partition (400 rounds, not 150):')
    changed = H.Strategy.of('boosted trees, every allowed family', 'hist_gbt', max_iter=400)
    refused = H.run_question(question, labels=labels, features=features, families=families,
                             groups=groups, strategies=[changed], seed=SEED,
                             outer_folds=3, inner_folds=3, record=False)
    print('  strategy hash:', changed.identity[:12], 'against', result.strategy_hash[:12])
    print('  flags        :', refused.flags)
    print('  frozen test  :', refused.frozen_test)
    print('  refusal      :', refused.refusal)
    print('  may fill     :', H.may_fill(refused))
""")

# --------------------------------------------------------------------------- 7 the other branch
nb.md(
    '## 7 · The other branch: a question that cannot afford a frozen partition',
    '`screen_egress_phenotype` has 318 labelled genes against a threshold of 400. A held-back '
    'partition costs labels a 318-gene question cannot spare, so the harness keeps nested '
    'cross-validation only, flags the question `no_frozen_test`, and **refuses to let it fill '
    'cells** — it can rank and prioritise, which is what those questions are for. The threshold is '
    'a property of the question in the registry, not a judgement made per run.')
code("""
    egress = Q.get('screen_egress_phenotype')
    y_egress = nodes['screen_egress_phenotype'].dropna().astype(bool)
    g_egress = nodes['orthogroup'].reindex(y_egress.index)
    g_egress = g_egress.where(g_egress.notna(), pd.Series({g: 'self:' + g for g in y_egress.index}))
    x_egress = H.allowed_features(egress, nodes.loc[y_egress.index, numeric].astype('float64'),
                                  families)
    print('labelled genes:', len(y_egress), '| positives:', int(y_egress.sum()))
    print('held out:', sorted(set(nodes[numeric].columns) - set(x_egress.columns)))
    egress_result = H.run_question(egress, labels=y_egress, features=x_egress, families=families,
                                   groups=g_egress,
                                   strategies=[H.Strategy.of('boosted trees', 'hist_gbt', max_iter=150),
                                               H.Strategy.of('logistic', 'logistic')],
                                   seed=0, outer_folds=3, inner_folds=3)
    print()
    print('flags      :', egress_result.flags)
    print('refusal    :', egress_result.refusal)
    print('frozen test:', egress_result.frozen_test)
    print('may fill   :', H.may_fill(egress_result))
    print()
    print('cross-validated estimate (ranking is what this question is for):')
    print(egress_result.summary().set_index('partition').T.to_string())
    egress_result.summary().to_csv(OUT / 'summary_egress.tsv', sep='\\t', index=False)
""")

# --------------------------------------------------------------------------- 8 what it means
nb.md('## 8 · What this run does and does not establish')
code("""
    lines = []
    cv, ft = result.cross_validated, result.frozen_test
    lines.append('QUESTION: %s' % question.text)
    lines.append('%d labels, %d positives, %d features in %d families, seed %d, split %s'
                 % (result.n_labels, int(labels.sum()), features.shape[1],
                    families.reindex(features.columns).nunique(), result.seed,
                    result.split_identity[:12]))
    lines.append('')
    lines.append('%-28s %12s %12s' % ('', 'cross-validated', 'frozen test'))
    for key in ('n', 'calls', 'precision', 'precision_lower90', 'coverage_at_bar',
                'calibration_error', 'roc_auc', 'average_precision'):
        lines.append('%-28s %12s %12s' % (key, round(cv[key], 4), round(ft[key], 4)))
    lines.append('')
    best = result.baselines.query("partition == 'frozen_test'")
    lines.append('best baseline on the frozen partition: %s, precision lower bound %.4f'
                 % (best.loc[best.precision_lower90.idxmax(), 'baseline'],
                    best.precision_lower90.max()))
    lines.append('fill bar %.2f | may fill: %s' % (result.bar, H.may_fill(result)))
    (OUT / 'verdict.txt').write_text('\\n'.join(lines) + '\\n')
    print('\\n'.join(lines))
""")
code("""
    json.dump({'required_without_glucose': result.as_dict(),
               'screen_egress_phenotype': egress_result.as_dict()},
              (OUT / 'results.json').open('w'), indent=2, default=str)
    sorted(p.name for p in OUT.iterdir())
""")
nb.md(
    '**What is established.** The architecture runs: four partitions grouped by orthogroup with the '
    'seed recorded, a guard that raises at every stage rather than being audited afterwards, nested '
    'cross-validation inside train + tune where the strategy and the threshold are chosen, isotonic '
    'calibration on the calibration partition alone, one scored look at the frozen partition, that '
    'look recorded in a persistent log, a refusal when a changed strategy asks for a second look, '
    'and both numbers reported together. A question below its own label threshold is flagged and '
    'cannot fill.',
    '**What is not established.** Nothing about glucose biology. The class A answer for these genes '
    'is the screen\'s own measurement (rung 0), not this model, and the conditional form of the '
    'question — *specifically* important without glucose rather than essential everywhere — still '
    'needs the screen\'s published per-gene statistics. The frozen-test number here is a property of '
    'the harness working, not a biological claim; and the gap between it and the cross-validated '
    'estimate is the size of the optimism the three earlier scoreboards carry.',
    '**Next in instruction 73:** rung 0 across class A (every genome-wide screen carried as a '
    'calibrated answer with its own uncertainty), then the catalogue with its redundancy report, '
    'then class B through this harness — host-facing, function and partners re-run with a frozen '
    'partition, which is the rest of phase 1\'s gate.')

path = nb.write(str(ROOT / 'notebooks/harness_demo_2026_11_05.ipynb'))
print('wrote', path)
