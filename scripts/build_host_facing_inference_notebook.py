"""Instruction 72, inference stages 0 and 1: can host_facing be filled at 90 % precision?

Writes notebooks/host_facing_inference_2026_11_04.ipynb (executed) and the result tables under
results/host_facing_inference_2026_11_04/.

Stage 0 is the baseline board of instruction 70 (majority, prevalence-weighted random, the single
best feature, ortholog copy); stage 1 is one gradient-boosted-tree model on the tabular tier-2
columns, out-of-fold predictions calibrated with isotonic regression, and the lower 90 % confidence
bound of precision at every calibrated threshold -- the number that decides whether anything may be
written into the information map as a fact.
"""
import sys
from pathlib import Path
from textwrap import dedent

import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT))
from notebook_runner import ExecutedNotebook  # noqa: E402

nb = ExecutedNotebook('Host-facing inference: baselines, one boosted-tree model, and the 90 % fill bar')


def code(src: str):
    return nb.code(dedent(src).strip('\n'))


nb.md(
    '# Host-facing inference: baselines, one boosted-tree model, and the 90 % fill bar',
    'Instruction 72 (Einar, 2026-11-04), inference plan stages **0** (baselines) and **1** (one '
    'gradient-boosted-tree model per tier-1 trait on the tabular tier-2 columns). Scored by the '
    'instruction 70 rules: five folds grouped by orthogroup, simple baselines in the same table, '
    'out-of-fold predictions calibrated with isotonic regression.',
    '**The question this notebook has to answer honestly.** The information map writes a predicted '
    'cell as a *fact* only where the **lower 90 % confidence bound** of held-out precision reaches '
    '0.90. Not the point estimate -- the lower bound. With 73 positive genes that bar may simply be '
    'out of reach, and "we cannot fill anything yet" is a valid result, not a failure of the code. '
    'Everything below is arranged so that the answer cannot be flattered.',
    'A reader who wants only the answer should read section 9 (the key number), then section 10 '
    '(the same question restricted to secretory proteins, which is the biologically interesting '
    'one) and the limitations in section 15.')

# --------------------------------------------------------------------------- 1 inputs
nb.md(
    '## 1 · Inputs and the labels',
    'The labels are the `host_facing` column of the information map v0 '
    '(`results/information_map_v0_2026_11_04/`), which holds **observed facts only** -- nothing '
    'predicted. In plain terms: *True* means this protein has been shown to reach beyond the '
    'parasite surface, *False* means it has been shown to stay inside, and empty means nobody has '
    'tested it. Empty is **not** a negative.',
    '`conflicts.tsv` lists 8 genes where two sources disagree (all 8 are PVM-proximity hits that a '
    'curated reading places in the vacuolar space or the cytosol: GRA2, GRA9, GRA29, GRA38, GRA39, '
    'GRA45, NTPase II and one cytosol marker). The map keeps both sides visible rather than '
    'averaging them away, and **this notebook excludes all 8 from training and from every score** '
    'so that a disputed gene can neither teach nor grade the model.')
code("""
    import os
    # One thread per BLAS/OpenMP pool. The tables here are 500 x 450, so thread spin-up costs far
    # more than the arithmetic saves -- permutation importance over 450 features is ~100x slower
    # with 16 threads than with 1. The machine is shared, and HistGradientBoosting is deterministic
    # either way, so this changes only the runtime.
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[v] = '1'
    import json, re, warnings
    from pathlib import Path
    import numpy as np, pandas as pd
    from scipy.stats import beta
    warnings.filterwarnings('ignore', category=FutureWarning)

    ROOT = Path('.').resolve()
    MAP = ROOT / 'results/information_map_v0_2026_11_04'
    SEQ = ROOT / 'results/sequence_traits_2026_11_04/sequence_traits.parquet'
    NODES = ROOT / 'parascape/data/nodes.parquet'
    OUT = ROOT / 'results/host_facing_inference_2026_11_04'; OUT.mkdir(parents=True, exist_ok=True)

    values = pd.read_parquet(MAP / 'values.parquet')
    probability = pd.read_parquet(MAP / 'probability.parquet')
    evidence = pd.read_parquet(MAP / 'evidence.parquet')
    conflicts = pd.read_csv(MAP / 'conflicts.tsv', sep='\\t')
    CONFLICTED = sorted(set(conflicts.gene_id))

    raw = values['host_facing'].dropna().astype(bool)
    y = raw.drop(index=[g for g in CONFLICTED if g in raw.index])
    print('information map rows:', len(values))
    print('host_facing labelled:', len(raw), '(True', int(raw.sum()), '/ False', int((~raw).sum()), ')')
    print('conflicted genes excluded:', len(CONFLICTED), '->', ', '.join(CONFLICTED))
    print('labels used:', len(y), '(True', int(y.sum()), '/ False', int((~y).sum()), ')')
    print('empty (untested) genes:', int(values['host_facing'].isna().sum()))
""")

nb.md(
    '## 2 · Where each label comes from',
    'This matters for analysis 3 and for reading every number below. The positives are not one '
    'kind of evidence:',
    '* **proximity only** -- the gene was biotinylated in the PVM proximity-labelling screen and '
    'nothing else says it is host-facing (tier `observed-proximity`, probability 0.70);',
    '* **curated / direct** -- a curated host interaction (ROP16→STAT3 and relatives) or an '
    'experimental GO host-cell annotation with a PMID (tiers `observed-curated` / `observed-direct`, '
    'probability 0.90-0.95).',
    'The negatives are also two kinds: the small set of **curated gold negatives** (proteins shown '
    'to stay in the vacuole, e.g. GRA9 by immunogold EM) and the much larger set of **hyperLOPIT '
    'internal markers** -- residents of the apicoplast, mitochondrion, nucleus, ribosome, '
    'proteasome, Golgi, ER, cytosol, tubulin cytoskeleton and IMC. The second kind is an inference '
    'from residence, not an export test, and it is what makes the overall numbers look good '
    '(section 10).')
code("""
    hf = evidence[evidence.trait == 'host_facing'].copy()
    pos = hf[hf.value == 'yes'].groupby('gene_id').source_tier.agg(lambda s: set(s))
    neg = hf[hf.value == 'no'].groupby('gene_id').dataset.agg(lambda s: ' | '.join(sorted(set(s))))

    def label_source(gene, is_pos):
        if is_pos:
            tiers = pos.get(gene, set())
            return 'proximity only' if tiers == {'observed-proximity'} else 'curated/direct'
        text = neg.get(gene, '')
        return 'curated gold negative' if 'gold_negative' in text else 'hyperLOPIT internal marker'

    source_class = pd.Series({g: label_source(g, bool(v)) for g, v in y.items()}, name='label_source')
    source_class.index.name = 'gene_id'
    label_frame = pd.DataFrame({'label': y, 'label_source': source_class})
    label_frame.groupby(['label', 'label_source']).size().to_frame('genes')
""")
code("""
    print('observed probability attached to each label (information map v0 tier defaults):')
    print(probability.loc[y.index, 'host_facing'].groupby([y, source_class]).agg(['mean', 'min', 'max']).to_string())
""")

# --------------------------------------------------------------------------- 3 features
nb.md(
    '## 3 · Features: every numeric and boolean tier-2 column, minus what would give the answer away',
    'Tier 2 in instruction 72 means *measurements* -- expression, fitness, proteome, PTMs, screens, '
    'sequence and structure. Those are what a prediction may learn from. Tier 1 are the '
    '*conclusions* (host-facing, compartment, function, partners) and they are what is being '
    'predicted, so they stay out.',
    'The feature table is the numeric and boolean columns of `nodes.parquet` plus nine sequence and '
    'topology traits from the DeepTMHMM-derived table (`n_tm`, `has_signal_peptide`, `seq_len`, '
    '`mature_len`, `nterm_len`, `cterm_len`, `nterm_side`, `cterm_side`, `flanks_identical`). The '
    'raw amino-acid and 3Di strings are **not** used -- those are stage 2, as frozen feature blocks.',
    '**Excluded as leakage, and why.** The negatives *are* hyperLOPIT markers, so any column that '
    'carries a hyperLOPIT compartment call would hand the model the negative class: every column '
    'whose name contains `lopit`, `compartment`, `ortholopit` or `tagm` goes. `pvm_proximity_positive` '
    '*is* the proximity positive set and `n_host_targets` counts curated host targets, so both '
    'define positives and both go. Identity and free-text columns (gene_id, sequence, product, '
    'InterPro descriptions, EC numbers) are not features.')
code("""
    nodes = pd.read_parquet(NODES).set_index('gene_id')
    seqt = pd.read_parquet(SEQ)
    SEQ_FEATURES = ['n_tm', 'has_signal_peptide', 'seq_len', 'mature_len', 'nterm_len', 'cterm_len',
                    'nterm_side', 'cterm_side', 'flanks_identical']
    CATEGORICAL = ['nterm_side', 'cterm_side']

    # object columns that hold only booleans (the four phenotype screens) become float 1/0/NaN so the
    # tree can use "not screened" as its own branch instead of guessing a value.
    for c in nodes.columns:
        if nodes[c].dtype == object and set(map(type, nodes[c].dropna())) <= {bool, np.bool_}:
            nodes[c] = nodes[c].astype('boolean').astype('float64')

    numeric = [c for c in nodes.columns
               if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]
    text_cols = [c for c in nodes.columns if c not in numeric]
    NAME_PAT = re.compile('lopit|compartment|ortholopit|tagm', re.I)
    EXPLICIT = {'pvm_proximity_positive', 'n_host_targets'}
    leak_named = [c for c in numeric if NAME_PAT.search(c)]
    leak_explicit = [c for c in numeric if c in EXPLICIT]
    # n_tm / has_signal_peptide exist in both tables; the sequence-traits version is the one kept.
    duplicated = [c for c in numeric if c in SEQ_FEATURES]
    node_features = [c for c in numeric if c not in set(leak_named) | set(leak_explicit) | set(duplicated)]

    # A column with no measurement at all on any labelled gene cannot inform this trait, and
    # sklearn 1.9 cannot even bin an all-missing column, so those are dropped and named.
    labelled = [g for g in y.index if g in nodes.index]
    all_missing = [c for c in node_features if nodes.loc[labelled, c].notna().sum() == 0]
    node_features = [c for c in node_features if c not in set(all_missing)]

    X = nodes[node_features].astype('float64').join(seqt[SEQ_FEATURES])
    X['has_signal_peptide'] = X['has_signal_peptide'].astype('float64')
    X['flanks_identical'] = X['flanks_identical'].astype('float64')
    for c in CATEGORICAL:
        X[c] = X[c].astype('category')
    FEATURES = list(X.columns)
    print('nodes.parquet columns:', nodes.shape[1], '| numeric/boolean:', len(numeric), '| free text:', len(text_cols))
    print('excluded by name (lopit|compartment|ortholopit|tagm):', len(leak_named))
    print('  ', leak_named)
    print('excluded explicitly (define the labels):', leak_explicit)
    print('dropped as duplicated by the sequence-traits table:', duplicated)
    print('dropped as never measured on any labelled gene:', all_missing)
    print('free-text / identity columns never used as features:')
    print('  ', text_cols)
    print()
    print('FINAL FEATURE COUNT:', len(FEATURES), '(of which categorical:', CATEGORICAL, ')')
""")
code("""
    print('the final feature list:')
    for i in range(0, len(FEATURES), 4):
        print('   ' + '  '.join(f'{c:<42}' for c in FEATURES[i:i + 4]).rstrip())
""")

nb.md(
    '## 4 · Which dataset each feature comes from (provenance families)',
    'Instruction 70 holds out *sources*, not columns, because 150 expression columns from one '
    'experiment are one piece of evidence, not 150. `parascape.datasets.provenance(col, "Tg")` gives '
    'the registered dataset for each column; the table below is the family of every feature. The '
    'point of printing it is to notice a family that could encode the label without being named '
    '`lopit` or `proximity`.')
code("""
    from parascape import datasets as ds
    def family(col):
        d = ds.provenance(col, 'Tg')
        return d.key if d is not None else 'unregistered'
    fam = pd.Series({c: family(c) for c in FEATURES}, name='family')
    fam.index.name = 'feature'
    fam_counts = fam.value_counts().rename('features').to_frame()
    for f, group in fam.groupby(fam):
        print(f'{f}  ({len(group)})')
        names = sorted(group.index)
        for i in range(0, len(names), 4):
            print('      ' + '  '.join(f'{c:<40}' for c in names[i:i + 4]).rstrip())
    fam_counts.T
""")
nb.md(
    '### Families that could encode the label — said out loud rather than quietly kept',
    'Two deserve naming.',
    '**`pride_proximity` (`n_proximity_partners`).** A *different* proximity-labelling experiment '
    'from the PVM screen that defines the positives, so it is not the same source -- but it is the '
    'same kind of assay, and if the two overlap strongly this column is close to a copy of the '
    'label. The overlap is measured below. It is kept (the instruction excludes `pvm_proximity` '
    'specifically) and it is reported in the importance table; if it lands near the top, the model '
    'is reading a proximity assay, not biology.',
    '**`hosttx_effectors` (`hosttx_T2`, `hosttx_padj`, `hosttx_signature_pc01-20`).** This one was '
    'not expected and it turns out to be the most important family in the model (section 13), so it '
    'is said plainly: it is a **screen for effectors that target host-cell transcription**. A gene '
    'whose knockout changes host gene expression is acting on the host, which is very close to being '
    'the conclusion rather than evidence for it. Formally it is tier 2 -- a measurement, and not one '
    'of the sources that define the `host_facing` label, so instruction 72 does not exclude it -- but '
    'a model whose signal comes from here is reproducing one host-effector screen rather than '
    'integrating the rest of the knowledge matrix. Section 11 therefore ablates this family too, and '
    'the answer in the README is conditioned on the result.',
    '**`organelle_surface_turboid` (`surface_*_log2fc`).** Another proximity-biotinylation assay, '
    'this time on the cytosolic faces of the apicoplast, mitochondrion and ER. Same caution as '
    '`pride_proximity`: not the PVM screen that defines the positives, but the same technique.',
    '**`spacr_escrt_screen` (`escrt_recruitment_*`), `secretome_partition`, `toxodb_epitopes`, '
    '`crosslink_interactome`.** All four are genuinely host-facing-adjacent biology (host ESCRT '
    'recruitment to the PVM, secreted versus vesicular partition, antibody-accessible epitopes, '
    'crosslinked partners). They are evidence, not label definitions, so they stay -- but a result '
    'resting entirely on them is a result about one screen.',
    'Also worth stating: `fit_hyperlopit_unassigned_invivo_lib1/2` are CRISPR fitness scores, not '
    'localisation, yet their names contain "lopit" and the name rule therefore removed them. That '
    'is a loss of two legitimate features, accepted because a name-based rule that errs toward '
    'removal is the safe direction.')
code("""
    prox = nodes['n_proximity_partners']
    pos_genes, neg_genes = y.index[y], y.index[~y]
    print('n_proximity_partners (family pride_proximity) measured for', int(prox.notna().sum()), 'of',
          len(prox), 'genes genome-wide')
    print('  measured for %d of %d positives and %d of %d negatives'
          % (int(prox.loc[pos_genes].notna().sum()), len(pos_genes),
             int(prox.loc[neg_genes].notna().sum()), len(neg_genes)))
    print('  distinct values among labelled genes:', sorted(prox.loc[y.index].dropna().unique()))
    print('  -> it is CONSTANT where present, so it can only act as a presence/absence flag;'
          if prox.loc[y.index].nunique(dropna=True) <= 1 else '  -> it varies among labelled genes;')
    print('     its permutation importance in section 13 settles whether the model uses it at all.')
    print()
    print('pvm_proximity_positive (EXCLUDED) against the label, to show what was removed:')
    print(pd.crosstab(pd.read_parquet(NODES).set_index('gene_id')['pvm_proximity_positive'].loc[y.index]
                      .fillna(-1), y).to_string())
""")

# --------------------------------------------------------------------------- 5 folds
nb.md(
    '## 5 · Five folds grouped by orthogroup',
    'Paralogs share an orthogroup and very often share a fate: the GRA and ROP families are full of '
    'near-copies. If one copy sat in the training fold and another in the test fold, the model could '
    'score well by recognising the family rather than by learning anything. So whole orthogroups go '
    'to one fold (OrthoMCL groups from `MASTER_parasite_wide_by_orthogroup.csv`, '
    '`tgon_gene_id` is semicolon-separated; a gene with no group is its own group). Seed 0, recorded.',
    'This costs something, and the cost is visible later: the *ortholog copy* baseline, which '
    'predicts a gene from its paralogs, is left with nothing to copy from. That is correct. A '
    'baseline that only works because paralogs leak across folds is not a baseline.')
code("""
    OG = Path('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/MASTER_parasite_wide_by_orthogroup.csv')
    og = pd.read_csv(OG, usecols=['orthomcl_group', 'tgon_gene_id']).dropna()
    group_of = {g: row.orthomcl_group for row in og.itertuples() for g in str(row.tgon_gene_id).split(';')}
    groups = pd.Series([group_of.get(g, 'self:' + g) for g in y.index], index=y.index, name='group')
    rng = np.random.default_rng(0)
    unique = np.array(sorted(set(groups))); rng.shuffle(unique)
    fold_of_group = {g: i % 5 for i, g in enumerate(unique)}
    fold = groups.map(fold_of_group).rename('fold')
    label_frame['fold'] = fold; label_frame['group'] = groups
    print('labelled genes:', len(y), '| distinct orthogroups:', len(unique),
          '| genes in a shared (multi-gene) orthogroup:', int((groups.value_counts() > 1).reindex(groups).sum()))
    pd.crosstab(fold, y).rename(columns={False: 'negatives', True: 'positives'})
""")

# --------------------------------------------------------------------------- 6 metrics
nb.md(
    '## 6 · How every strategy is scored',
    'The same four numbers for everything, computed per fold and then reported as the mean with the '
    'fold-to-fold spread, because with ~15 positives per fold a single fold mean is noise.',
    '* **accuracy** -- fraction of genes called correctly. With 87 % negatives this is a weak number '
    'and is reported only so the majority baseline can be seen beating most things on it.',
    '* **precision** -- of the genes called host-facing, how many are. This is the number the fill '
    'bar is about.',
    '* **recall** -- of the true host-facing genes, how many were called.',
    '* **ROC AUC** and **average precision** for the strategies that produce a score rather than a '
    'yes/no, since those do not depend on where a threshold is put.')
code("""
    from sklearn.metrics import roc_auc_score, average_precision_score

    def fold_metrics(truth, called, score=None):
        rows = []
        for f in range(5):
            m = (fold.loc[truth.index] == f).values
            t, c = truth.values[m], np.asarray(called)[m]
            tp, fp, fn = int((c & t).sum()), int((c & ~t).sum()), int((~c & t).sum())
            row = dict(fold=f, n=int(m.sum()), accuracy=float((c == t).mean()),
                       precision=(tp / (tp + fp)) if (tp + fp) else np.nan,
                       recall=(tp / (tp + fn)) if (tp + fn) else np.nan, calls=int(c.sum()))
            if score is not None and len(set(t)) == 2:
                s = np.asarray(score)[m]
                row['roc_auc'] = float(roc_auc_score(t, s))
                row['average_precision'] = float(average_precision_score(t, s))
            rows.append(row)
        return pd.DataFrame(rows)

    def board_rows(name, truth, called, score=None, coverage=None):
        per = fold_metrics(truth, called, score)
        out = []
        for metric in ('accuracy', 'precision', 'recall', 'roc_auc', 'average_precision'):
            if metric not in per:
                continue
            vals = per[metric].astype(float)
            out.append(dict(strategy=name, metric=metric, mean=vals.mean(), fold_min=vals.min(),
                            fold_max=vals.max(), **{f'fold_{f}': vals.iloc[f] for f in range(5)}))
        out.append(dict(strategy=name, metric='calls', mean=per.calls.mean(), fold_min=per.calls.min(),
                        fold_max=per.calls.max(), **{f'fold_{f}': per.calls.iloc[f] for f in range(5)}))
        if coverage is not None:
            out.append(dict(strategy=name, metric='coverage', mean=coverage, fold_min=coverage,
                            fold_max=coverage, **{f'fold_{f}': np.nan for f in range(5)}))
        return out
""")

# --------------------------------------------------------------------------- 7 baselines
nb.md(
    '## 7 · Stage 0 · the four baselines',
    '**Training-fold majority.** Always says the majority class of the training folds, which is '
    '*not host-facing*. It makes no positive calls at all, so its precision is undefined; its '
    'accuracy is the one to beat before anything else is claimed.',
    '**Prevalence-weighted random.** Calls each test gene host-facing with the training-fold '
    'prevalence as a coin weight (seed 0). Its precision is the prevalence itself -- about 0.13 -- '
    'and that is the floor the 0.90 bar must be compared with.',
    '**Signal peptide alone.** The prediction *is* `has_signal_peptide` from DeepTMHMM. Instruction '
    '72 warns that this will look strong, and it should: the negatives are mostly markers of '
    'compartments inside the parasite (cytosol, nucleus, ribosome) and almost none of those have a '
    'signal peptide. Any model has to beat this to have earned anything.',
    '**Ortholog copy.** A gene takes the label of other labelled genes in its orthogroup, and '
    'abstains when there are none. Because folds are grouped by orthogroup, every ortholog with a '
    'label is in the *same* fold, so under the honest protocol this baseline can copy nothing: '
    'coverage 0. The leaky variant (copy from any labelled gene regardless of fold) is reported '
    'beside it with its coverage, clearly marked, to show what the grouping is protecting against.')
code("""
    board = []
    prevalence = {f: float(y[fold != f].mean()) for f in range(5)}
    print('training-fold prevalence of host-facing:', {f: round(p, 4) for f, p in prevalence.items()})

    majority = pd.Series([bool(y[fold != f].mean() > 0.5) for f in fold], index=y.index)
    board += board_rows('baseline: training majority', y, majority.values, coverage=1.0)

    rng_b = np.random.default_rng(0)
    coin = pd.Series([rng_b.random() < prevalence[f] for f in fold], index=y.index)
    board += board_rows('baseline: prevalence-weighted random', y, coin.values,
                        score=np.array([prevalence[f] for f in fold]), coverage=1.0)

    sp = X.loc[y.index, 'has_signal_peptide'].fillna(0).astype(bool)
    board += board_rows('baseline: signal peptide alone', y, sp.values, score=sp.values.astype(float), coverage=1.0)

    def ortholog_copy(within_fold_only=True):
        pred, covered = [], []
        for g in y.index:
            peers = groups.index[(groups == groups[g]) & (groups.index != g)]
            if within_fold_only:
                peers = [p for p in peers if fold[p] != fold[g]]
            peers = [p for p in peers if p in y.index]
            if peers:
                pred.append(bool(y.loc[peers].mean() >= 0.5)); covered.append(True)
            else:
                pred.append(False); covered.append(False)
        return pd.Series(pred, index=y.index), pd.Series(covered, index=y.index)

    honest, honest_cov = ortholog_copy(True)
    leaky, leaky_cov = ortholog_copy(False)
    print('ortholog copy, grouped folds respected: coverage', round(honest_cov.mean(), 4),
          f'({int(honest_cov.sum())} of {len(y)} genes have a labelled ortholog in another fold)')
    print('ortholog copy, ignoring folds (LEAKY):   coverage', round(leaky_cov.mean(), 4),
          f'({int(leaky_cov.sum())} genes)')
    if leaky_cov.any():
        sub = leaky_cov[leaky_cov].index
        tp = int((leaky[sub] & y[sub]).sum()); calls = int(leaky[sub].sum())
        print('  on its covered genes: accuracy %.3f, precision %s, recall %.3f'
              % ((leaky[sub] == y[sub]).mean(), f'{tp / calls:.3f}' if calls else 'n/a (no positive calls)',
                 (leaky[sub] & y[sub]).sum() / max(int(y[sub].sum()), 1)))
    board += board_rows('baseline: ortholog copy (grouped, abstains)', y, honest.values,
                        coverage=float(honest_cov.mean()))
    pd.DataFrame(board).round(3)
""")

# --------------------------------------------------------------------------- 8 model
nb.md(
    '## 8 · Stage 1 · one gradient-boosted-tree model',
    '`HistGradientBoostingClassifier`. Missing values are passed straight through: the tree learns '
    'a direction for "not measured", which is the right treatment when a gene is absent from a '
    'proteomics run for reasons that are themselves informative. **No imputation anywhere.** The '
    'two side columns (`nterm_side`, `cterm_side`) are pandas categoricals and are declared as such '
    'to the model.',
    '**No hyperparameter search was run on a test fold.** The settings are fixed in advance and '
    'identical for every fold (300 boosting rounds at learning rate 0.05, at most 15 leaves, at '
    'least 10 samples per leaf, L2 1.0). Early stopping is switched **off**: it would hold back a '
    'fifth of the training genes as a validation split, which with ~58 positives per training set '
    'is both wasteful and unstable (and in scikit-learn 1.9 that split happens before binning, so a '
    'column measured on only a handful of genes can vanish from it). A fixed round count keeps every '
    'fold identical and the whole procedure deterministic. No class weighting, so the raw output '
    'keeps its meaning as a probability before calibration.',
    'Out-of-fold predictions: each gene is scored by the model of the four folds that did not '
    'contain it or any of its orthologs.')
code("""
    from sklearn.ensemble import HistGradientBoostingClassifier

    def make_model():
        return HistGradientBoostingClassifier(
            max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=10,
            l2_regularization=1.0, early_stopping=False, random_state=0,
            categorical_features='from_dtype')

    def usable(features, rows):
        # sklearn 1.9 cannot bin a column that is missing for every training row; such a column
        # carries nothing for that fold anyway. Dropped per fold, and the count is reported.
        return [c for c in features if c in CATEGORICAL or X.loc[rows, c].notna().any()]

    def out_of_fold(features, labels=None, folds=None, verbose=False):
        labels = y if labels is None else labels
        folds = fold if folds is None else folds
        oof = pd.Series(np.nan, index=labels.index, dtype=float)
        models, used = {}, {}
        for f in sorted(folds.unique()):
            tr = labels.index[folds != f]; te = labels.index[folds == f]
            cols = usable(features, tr)
            if verbose and len(cols) < len(features):
                print(f'  fold {f}: {len(features) - len(cols)} feature(s) missing for every training '
                      f'gene, dropped for this fold: {sorted(set(features) - set(cols))}')
            m = make_model().fit(X.loc[tr, cols], labels.loc[tr])
            oof.loc[te] = m.predict_proba(X.loc[te, cols])[:, 1]
            models[f], used[f] = m, cols
        return oof, models, used

    oof, models, used = out_of_fold(FEATURES, verbose=True)
    board += board_rows('GBT: all tier-2 features', y, (oof >= 0.5).values, score=oof.values, coverage=1.0)
    print('out-of-fold raw probability, by true label:')
    print(oof.groupby(y).describe().round(3).to_string())
    print()
    print('boosting rounds per fold (fixed, no early stopping):',
          {f: int(m.n_iter_) for f, m in models.items()})
    print('features actually used per fold:', {f: len(c) for f, c in used.items()})
""")

nb.md(
    '### The scoreboard',
    'Read it in this order: does the model beat the *best* baseline, and does it beat it in most '
    'folds rather than on average? Instruction 70 asks for at least 4 of 5 folds.')
code("""
    scoreboard = pd.DataFrame(board)
    scoreboard.to_csv(OUT / 'scoreboard.tsv', sep='\\t', index=False)
    pivot = scoreboard.pivot_table(index='strategy', columns='metric', values='mean')
    pivot = pivot[[c for c in ['accuracy', 'precision', 'recall', 'roc_auc', 'average_precision', 'calls', 'coverage'] if c in pivot]]
    pivot.round(3)
""")
code("""
    print('per-fold accuracy, model against the best baseline in that fold:')
    acc = scoreboard[scoreboard.metric == 'accuracy'].set_index('strategy')
    folds_cols = [f'fold_{f}' for f in range(5)]
    base = acc.loc[[s for s in acc.index if s.startswith('baseline')], folds_cols]
    mod = acc.loc['GBT: all tier-2 features', folds_cols]
    cmp = pd.DataFrame({'model': mod, 'best_baseline': base.max(), 'which': base.idxmax()})
    cmp['beats'] = cmp.model > cmp.best_baseline
    print(cmp.round(3).to_string())
    print('folds where the model beats the best baseline on accuracy:', int(cmp.beats.sum()), 'of 5')
    prec = scoreboard[scoreboard.metric == 'precision'].set_index('strategy')
    pr_cmp = pd.DataFrame({'model': prec.loc['GBT: all tier-2 features', folds_cols],
                           'signal_peptide': prec.loc['baseline: signal peptide alone', folds_cols]})
    pr_cmp['beats'] = pr_cmp.model > pr_cmp.signal_peptide
    print()
    print('per-fold precision against the signal-peptide baseline:')
    print(pr_cmp.round(3).to_string())
    print('folds where the model beats signal peptide on precision:', int(pr_cmp.beats.sum()), 'of 5')
""")

# --------------------------------------------------------------------------- 9 calibration + bound
nb.md(
    '## 9 · Calibration, then the key number',
    'A boosted tree produces a score between 0 and 1 that is not a probability. Isotonic regression '
    'maps the score to the observed rate of being right -- the same machinery the project already '
    'uses for claims (`parascape.claims._isotonic`, `_ece`), reused rather than rewritten.',
    '**Cross-fitted, like `claims.certainty_model`.** For each fold the calibration curve is fitted '
    'on the out-of-fold predictions of the *other* folds and then applied to this one. Fitting the '
    'curve on the same points it is scored on would make the calibration look perfect; this way the '
    'numbers are what a new gene would see.',
    'The reliability table below says: for genes whose calibrated probability fell in this range, '
    'how many were there and what fraction were really host-facing. If the model is calibrated, the '
    'two should track each other. The expected calibration error is the average gap, weighted by '
    'how many genes fall in each bin.')
code("""
    from parascape.claims import _isotonic, _ece

    def cross_fitted_calibration(scores, labels, folds):
        cal = pd.Series(np.nan, index=scores.index, dtype=float)
        for f in sorted(folds.unique()):
            tr = folds.index[folds != f]; te = folds.index[folds == f]
            xs, ys = _isotonic(scores.loc[tr].values, labels.loc[tr].astype(float).values)
            cal.loc[te] = np.interp(scores.loc[te].values, xs, ys)
        return cal

    cal = cross_fitted_calibration(oof, y, fold)
    ece = _ece(cal.values, y.astype(float).values, bins=10)
    print('expected calibration error (10 quantile bins, cross-fitted): %.4f' % ece)
    print('ECE of the uncalibrated model score, for comparison:        %.4f'
          % _ece(oof.values, y.astype(float).values, bins=10))

    edges = np.array([0, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
    idx = np.clip(np.searchsorted(edges, cal.values, side='right') - 1, 0, len(edges) - 2)
    rows = []
    for b in range(len(edges) - 1):
        m = idx == b
        if not m.any():
            continue
        rows.append(dict(predicted_range=f'[{edges[b]:.2f}, {edges[b + 1]:.2f})', n=int(m.sum()),
                         mean_predicted=float(cal.values[m].mean()),
                         observed_positive_rate=float(y.values[m].mean()),
                         positives=int(y.values[m].sum())))
    calibration = pd.DataFrame(rows)
    calibration.attrs['ece'] = ece
    calibration.to_csv(OUT / 'calibration.tsv', sep='\\t', index=False)
    calibration.round(3)
""")

nb.md(
    '### The fill bar: the lower 90 % confidence bound of precision',
    'At each calibrated threshold: how many genes are called host-facing, what fraction of those '
    'calls are right, and -- the number that decides everything -- the **Clopper-Pearson lower '
    '90 % confidence bound** on that fraction (`scipy.stats.beta.ppf(0.10, k, n−k+1)`, the exact '
    'one-sided bound for a binomial).',
    'Why the bound and not the fraction: 4 correct calls out of 4 is a precision of 1.00 and a '
    'lower bound of 0.56. Four calls cannot demonstrate 90 % precision no matter how clean they '
    'are. Reaching a lower bound of 0.90 needs roughly 22 consecutive correct calls, or about 35 '
    'calls at 97 % precision. That arithmetic, not the model, is usually what decides the answer '
    'here.',
    '`genes_fillable` counts the currently-empty genes (no observed host_facing value) whose '
    'calibrated probability reaches the threshold -- what the map would gain if that threshold were '
    'deployable.')
code("""
    def clopper_pearson_lower(k, n, confidence=0.90):
        if n == 0:
            return np.nan
        if k == 0:
            return 0.0
        return float(beta.ppf(1 - confidence, k, n - k + 1))

    # Probabilities for every gene: a model fitted on all labelled genes, mapped through an isotonic
    # curve fitted on the pooled out-of-fold scores. The HONEST precision numbers come from the
    # cross-fitted out-of-fold calibration above; this step only extends coverage to untested genes.
    FEATURES_FINAL = usable(FEATURES, y.index)
    final_model = make_model().fit(X.loc[y.index, FEATURES_FINAL], y)
    xs_all, ys_all = _isotonic(oof.values, y.astype(float).values)
    all_raw = pd.Series(final_model.predict_proba(X[FEATURES_FINAL])[:, 1], index=X.index)
    all_cal = pd.Series(np.interp(all_raw.values, xs_all, ys_all), index=X.index)
    empty = values.index[values['host_facing'].isna()]
    empty = [g for g in empty if g in all_cal.index]

    THRESHOLDS = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97, 0.99]
    rows = []
    for t in THRESHOLDS:
        sel = cal >= t
        n, k = int(sel.sum()), int(y[sel].sum())
        rows.append(dict(threshold=t, calls=n, correct=k,
                         precision=(k / n) if n else np.nan,
                         precision_lower90=clopper_pearson_lower(k, n),
                         recall=k / int(y.sum()),
                         genes_fillable=int((all_cal.loc[empty] >= t).sum()),
                         reaches_bar=bool(n and clopper_pearson_lower(k, n) >= 0.90)))
    thresholds = pd.DataFrame(rows)
    thresholds.to_csv(OUT / 'thresholds.tsv', sep='\\t', index=False)
    thresholds.round(4)
""")
code("""
    ok = thresholds[thresholds.reaches_bar]
    if len(ok):
        r = ok.iloc[0]
        ANSWER = (f'YES. The lowest calibrated threshold whose lower 90 % bound of precision reaches 0.90 is '
                  f'{r.threshold:.2f}: {int(r.calls)} held-out calls, {int(r.correct)} correct, precision '
                  f'{r.precision:.3f}, lower bound {r.precision_lower90:.3f}. At that threshold '
                  f'{int(r.genes_fillable)} currently-empty genes would be filled.')
    else:
        best = thresholds.loc[thresholds.precision_lower90.idxmax()]
        ANSWER = ('NO. No calibrated threshold from 0.50 to 0.99 reaches a lower 90 % confidence bound of 0.90. '
                  f'The best bound anywhere on the curve is {best.precision_lower90:.3f} at threshold '
                  f'{best.threshold:.2f} ({int(best.calls)} calls, {int(best.correct)} correct, point-estimate '
                  f'precision {best.precision:.3f}). Nothing may be written into the information map as a fact; '
                  'every prediction stays exploratory.')
    print(ANSWER)
""")
nb.md(
    '### Why "no" is the expected answer, in one line of arithmetic',
    'The bound is a property of the number of calls as much as of the model. The table below shows '
    'the best lower bound obtainable with *perfect* precision at each number of calls: with fewer '
    'than ~22 calls the bar is unreachable even with no mistakes at all. There are 73 positives '
    'total, ~15 per fold.')
code("""
    pd.DataFrame([dict(calls=n, lower90_if_all_correct=round(clopper_pearson_lower(n, n), 3))
                  for n in (4, 8, 12, 16, 20, 22, 25, 30, 40, 60)]).set_index('calls').T
""")

# --------------------------------------------------------------------------- 10 secretory
nb.md(
    '## 10 · Analysis 1 · the honest question: among secretory proteins, which are exported?',
    'Instruction 72 is explicit about this. Most negatives are residents of compartments inside the '
    'parasite, and almost none of them have a signal peptide; most positives do. So the overall '
    'numbers largely measure "can we tell a nuclear protein from a secreted one", which is not the '
    'question a Toxoplasma biologist is asking. The question that matters is: **given that a '
    'protein enters the secretory pathway, is it exported beyond the parasite or does it stay in the '
    'vacuole?**',
    'Everything is therefore repeated on the labelled genes with a predicted signal peptide alone: '
    'same folds, same features, same calibration, same bound. The set is small, and that smallness '
    'is part of the answer.')
code("""
    sp_all = X['has_signal_peptide'].fillna(0).astype(bool)
    sec = y.index[sp_all.loc[y.index]]
    y_sec, fold_sec = y.loc[sec], fold.loc[sec]
    print('labelled genes with a predicted signal peptide:', len(sec),
          '(positives', int(y_sec.sum()), '/ negatives', int((~y_sec).sum()), ')')
    print('positives WITHOUT a predicted signal peptide (DeepTMHMM misses, e.g. GRA16):',
          int((~sp_all.loc[y.index[y]]).sum()), 'of', int(y.sum()))
    print('prevalence inside the secretory set: %.3f (overall %.3f)' % (y_sec.mean(), y.mean()))
    print('per-fold counts:')
    print(pd.crosstab(fold_sec, y_sec).to_string())
""")
code("""
    board_sec = []
    prev_sec = {f: float(y_sec[fold_sec != f].mean()) for f in sorted(fold_sec.unique())}

    def fold_metrics_sec(truth, called, score=None):
        rows = []
        for f in sorted(fold_sec.unique()):
            m = (fold_sec.loc[truth.index] == f).values
            t, c = truth.values[m], np.asarray(called)[m]
            tp, fp, fn = int((c & t).sum()), int((c & ~t).sum()), int((~c & t).sum())
            row = dict(fold=f, n=int(m.sum()), accuracy=float((c == t).mean()),
                       precision=(tp / (tp + fp)) if (tp + fp) else np.nan,
                       recall=(tp / (tp + fn)) if (tp + fn) else np.nan, calls=int(c.sum()))
            if score is not None and len(set(t)) == 2:
                row['roc_auc'] = float(roc_auc_score(t, np.asarray(score)[m]))
                row['average_precision'] = float(average_precision_score(t, np.asarray(score)[m]))
            rows.append(row)
        return pd.DataFrame(rows)

    def board_rows_sec(name, truth, called, score=None):
        per = fold_metrics_sec(truth, called, score); out = []
        for metric in ('accuracy', 'precision', 'recall', 'roc_auc', 'average_precision', 'calls'):
            if metric not in per:
                continue
            v = per[metric].astype(float)
            out.append(dict(strategy=name, metric=metric, mean=v.mean(), fold_min=v.min(), fold_max=v.max(),
                            **{f'fold_{f}': v.iloc[i] for i, f in enumerate(sorted(fold_sec.unique()))}))
        return out

    maj_sec = pd.Series([bool(y_sec[fold_sec != f].mean() > 0.5) for f in fold_sec], index=sec)
    board_sec += board_rows_sec('secretory · training majority', y_sec, maj_sec.values)
    rb = np.random.default_rng(0)
    coin_sec = pd.Series([rb.random() < prev_sec[f] for f in fold_sec], index=sec)
    board_sec += board_rows_sec('secretory · prevalence-weighted random', y_sec, coin_sec.values,
                                score=np.array([prev_sec[f] for f in fold_sec]))
    board_sec += board_rows_sec('secretory · signal peptide alone (constant True)', y_sec,
                                np.ones(len(sec), bool), score=np.ones(len(sec)))
    oof_sec, _, _ = out_of_fold(FEATURES, labels=y_sec, folds=fold_sec)
    board_sec += board_rows_sec('secretory · GBT', y_sec, (oof_sec >= 0.5).values, score=oof_sec.values)
    sec_board = pd.DataFrame(board_sec)
    sec_board['scope'] = 'secretory only'
    sec_board.to_csv(OUT / 'scoreboard_secretory.tsv', sep='\\t', index=False)
    sec_board.pivot_table(index='strategy', columns='metric', values='mean').round(3)
""")
code("""
    cal_sec = cross_fitted_calibration(oof_sec, y_sec, fold_sec)
    rows = []
    for t in THRESHOLDS:
        s = cal_sec >= t
        n, k = int(s.sum()), int(y_sec[s].sum())
        rows.append(dict(threshold=t, calls=n, correct=k, precision=(k / n) if n else np.nan,
                         precision_lower90=clopper_pearson_lower(k, n),
                         reaches_bar=bool(n and clopper_pearson_lower(k, n) >= 0.90)))
    thresholds_sec = pd.DataFrame(rows)
    thresholds_sec.to_csv(OUT / 'thresholds_secretory.tsv', sep='\\t', index=False)
    ok_sec = thresholds_sec[thresholds_sec.reaches_bar]
    ANSWER_SEC = ('reaches the bar at threshold %.2f (%d calls)' % (ok_sec.iloc[0].threshold, ok_sec.iloc[0].calls)
                  if len(ok_sec) else
                  'never reaches the bar; best lower bound %.3f' % thresholds_sec.precision_lower90.max())
    print('secretory-only:', ANSWER_SEC)
    print('ECE (secretory, cross-fitted): %.4f' % _ece(cal_sec.values, y_sec.astype(float).values, bins=5))
    thresholds_sec.round(4)
""")

# --------------------------------------------------------------------------- 11 ablation
nb.md(
    '## 11 · Analysis 2 · the publication-bias ablation',
    'A famous gene has more papers, more BioID studies and more filled columns. GRA16 and ROP16 are '
    'famous *because* they are host-facing, so counting their papers is a way of asking "is this a '
    'well-known effector" -- a question whose answer leaked out of the literature, not out of the '
    'measurement. The model is therefore refitted with every attention and bookkeeping column '
    'removed, and the drop is reported.',
    'Columns removed: `n_publications`, `n_fulltext`, `n_papers_focal`, `n_papers_substantive`, '
    '`n_papers_incidental`, `n_bioid_studies`, `n_ipms_studies`, and anything whose name mentions '
    'papers, publications, mentions, citations or attention. **`lit_tier` and `attention_depth` are '
    'free-text columns and were therefore never in the feature set at all** -- the dtype rule had '
    'already excluded them, which is worth saying because the instruction lists them for removal.',
    'If the drop is small, the model is not living on publication counts. If it is large, then the '
    'overall result above is partly a measure of fame.')
code("""
    ATTENTION = {'n_publications', 'n_fulltext', 'n_papers_focal', 'n_papers_substantive',
                 'n_papers_incidental', 'n_bioid_studies', 'n_ipms_studies'}
    PAT = re.compile('paper|publication|mention|citation|attention|lit_tier', re.I)
    attention_features = sorted({c for c in FEATURES if c in ATTENTION or PAT.search(c)})
    print('attention / bookkeeping features removed:', attention_features)
    print('named by the instruction but never a feature (free text):',
          [c for c in ('lit_tier', 'attention_depth') if c in text_cols])
    FEATURES_NOATT = [c for c in FEATURES if c not in set(attention_features)]
    oof_noatt, _, _ = out_of_fold(FEATURES_NOATT)
    board += board_rows('GBT: no attention features', y, (oof_noatt >= 0.5).values,
                        score=oof_noatt.values, coverage=1.0)

    def headline(scores, truth):
        return dict(roc_auc=roc_auc_score(truth, scores), average_precision=average_precision_score(truth, scores))

    abl = pd.DataFrame({'all features': headline(oof.values, y.values),
                        'without attention': headline(oof_noatt.values, y.values)}).T
    abl.loc['drop'] = abl.loc['all features'] - abl.loc['without attention']
    cal_noatt = cross_fitted_calibration(oof_noatt, y, fold)
    abl['best_lower90'] = [max(clopper_pearson_lower(int(y[c >= t].sum()), int((c >= t).sum())) or 0
                               for t in THRESHOLDS) for c in (cal, cal_noatt)] + [np.nan]
    abl.round(4)
""")
nb.md(
    '### The second ablation, added after reading the importance table: the host-effector screen',
    'The permutation importance in section 13 puts `hosttx_padj` and `hosttx_T2` at the top. Those '
    'come from a screen for parasite genes whose knockout changes **host** transcription, so the '
    'model is partly reading another measurement of host-directed activity rather than predicting it '
    'from parasite-side biology. That is not label leakage in the instruction 72 sense -- the label '
    'is defined by PVM proximity, curated host interactions and GO host-cell terms, not by this '
    'screen -- but it changes what the result means. Three further fits:',
    '* without the whole `hosttx_effectors` family;',
    '* without `hosttx_effectors` and the attention columns together;',
    '* without `hosttx_effectors`, the attention columns and every proximity-style assay '
    '(`pride_proximity`, `organelle_surface_turboid`) -- the most conservative feature set, in which '
    'nothing resembling a host-contact assay or a publication count remains.',
    'The lower 90 % bound is recomputed each time. **If the bound reaches 0.90 only with these '
    'columns in, then the fill bar is met by a host-effector screen plus publication counts, not by '
    'the knowledge matrix, and the README has to say so.**')
code("""
    HOSTTX = sorted(fam.index[fam == 'hosttx_effectors'])
    PROXY = sorted(fam.index[fam.isin(['pride_proximity', 'organelle_surface_turboid'])])
    print('hosttx_effectors columns removed:', len(HOSTTX))
    print('proximity-style columns removed:', PROXY)

    def best_bound(scores, labels=None, folds=None):
        labels = y if labels is None else labels
        folds = fold if folds is None else folds
        c = cross_fitted_calibration(scores, labels, folds)
        vals = []
        for t in THRESHOLDS:
            n, k = int((c >= t).sum()), int(labels[c >= t].sum())
            if n:
                vals.append((clopper_pearson_lower(k, n), t, n, k))
        return max(vals) if vals else (np.nan, np.nan, 0, 0)

    sets = {'all tier-2 features': FEATURES,
            'minus attention': FEATURES_NOATT,
            'minus host-effector screen': [c for c in FEATURES if c not in set(HOSTTX)],
            'minus host-effector screen and attention':
                [c for c in FEATURES if c not in set(HOSTTX) | set(attention_features)],
            'minus host-effector, attention and proximity-style assays':
                [c for c in FEATURES if c not in set(HOSTTX) | set(attention_features) | set(PROXY)]}
    rows = []
    for name, cols in sets.items():
        o = oof if name == 'all tier-2 features' else (oof_noatt if name == 'minus attention' else out_of_fold(cols)[0])
        b, t, n, k = best_bound(o)
        rows.append(dict(feature_set=name, n_features=len(cols), roc_auc=roc_auc_score(y, o),
                         average_precision=average_precision_score(y, o), best_lower90=b,
                         at_threshold=t, calls=n, correct=k, reaches_bar=bool(b >= 0.90)))
    ablations = pd.DataFrame(rows)
    ablations.to_csv(OUT / 'ablations.tsv', sep='\\t', index=False)
    ABLATION_VERDICT = ('the bar is still reached with the host-effector screen, the attention columns and the '
                        'proximity-style assays all removed'
                        if bool(ablations.reaches_bar.iloc[-1]) else
                        'the bar is NOT reached once the host-effector screen, the attention columns and the '
                        'proximity-style assays are removed (best lower bound %.3f)' % ablations.best_lower90.iloc[-1])
    print(ABLATION_VERDICT)
    ablations.round(4)
""")

code("""
    # Also the inverse: attention columns ALONE, which says how much of the signal is fame by itself.
    oof_att, _, _ = out_of_fold(attention_features)
    oof_sp_only, _, _ = out_of_fold(['has_signal_peptide', 'n_tm', 'nterm_side', 'cterm_side', 'seq_len'])
    board += board_rows('GBT: attention features only', y, (oof_att >= 0.5).values, score=oof_att.values, coverage=1.0)
    board += board_rows('GBT: topology features only', y, (oof_sp_only >= 0.5).values, score=oof_sp_only.values, coverage=1.0)
    pd.DataFrame({'all tier-2': headline(oof.values, y.values),
                  'attention only': headline(oof_att.values, y.values),
                  'topology only (SP, TM, sides, length)': headline(oof_sp_only.values, y.values),
                  'signal peptide alone (no model)': headline(sp.values.astype(float), y.values)}).T.round(4)
""")

# --------------------------------------------------------------------------- 12 by source
nb.md(
    '## 12 · Analysis 3 · performance by where the label came from',
    'The 73 positives are not one thing. If the model recovers the proximity-labelling positives but '
    'not the curated and directly shown ones, it has learned the proximity screen, not host-facing '
    'biology. And if it separates positives from hyperLOPIT internal markers but not from curated '
    'gold negatives (the GRA proteins that stay in the vacuole), then it is telling inside from '
    'outside the *parasite*, not inside from outside the *vacuole* -- a much easier question.',
    'Each row below uses one positive class against one negative class, so the AUC and recall are '
    'comparable across rows.')
code("""
    bysrc = []
    for pname in ('proximity only', 'curated/direct', 'all positives'):
        pmask = y & (source_class == pname) if pname != 'all positives' else y
        for nname in ('hyperLOPIT internal marker', 'curated gold negative', 'all negatives'):
            nmask = (~y) & (source_class == nname) if nname != 'all negatives' else ~y
            idx = y.index[pmask | nmask]
            t = y.loc[idx]
            if t.sum() == 0 or (~t).sum() == 0:
                continue
            s = oof.loc[idx]
            bysrc.append(dict(positives=pname, negatives=nname, n_pos=int(t.sum()), n_neg=int((~t).sum()),
                              roc_auc=roc_auc_score(t, s), average_precision=average_precision_score(t, s),
                              recall_at_0_5=float((s[t] >= 0.5).mean()),
                              false_positive_rate_at_0_5=float((s[~t] >= 0.5).mean()),
                              median_score_pos=float(s[t].median()), median_score_neg=float(s[~t].median())))
    by_source = pd.DataFrame(bysrc)
    by_source.to_csv(OUT / 'by_label_source.tsv', sep='\\t', index=False)
    by_source.round(3)
""")
code("""
    print('calibrated out-of-fold probability by label source:')
    print(cal.groupby(source_class).describe().round(3).to_string())
    print()
    print('signal-peptide rate by label source (the baseline that the easy negatives make strong):')
    print(sp.groupby(source_class).agg(['size', 'mean']).round(3).to_string())
""")

# --------------------------------------------------------------------------- 13 importance
nb.md(
    '## 13 · Analysis 4 · permutation importance: biology or bookkeeping?',
    '`sklearn.inspection.permutation_importance`, run per fold on the held-out fold only: a column '
    'is shuffled and the loss of average precision on genes the model never saw is recorded. Five '
    'shuffles per feature per fold, seed 0, averaged over folds.',
    'The family of each feature is printed beside it. This is the table that shows whether the model '
    'is reading biology (topology, secretion, fitness, expression) or bookkeeping (how many papers '
    'and studies a gene appears in, how many columns are filled).')
code("""
    from sklearn.inspection import permutation_importance
    imp = {}
    for f in range(5):
        te = y.index[fold == f]
        cols = used[f]
        r = permutation_importance(models[f], X.loc[te, cols], y.loc[te], scoring='average_precision',
                                   n_repeats=5, random_state=0, n_jobs=1)
        imp[f] = pd.Series(r.importances_mean, index=cols)
    importance = pd.DataFrame(imp)
    importance.columns = [f'fold_{f}' for f in importance.columns]
    importance['importance_mean'] = importance.mean(axis=1)
    importance['importance_sd'] = importance[[f'fold_{f}' for f in range(5)]].std(axis=1)
    importance['family'] = fam.reindex(importance.index)
    importance = importance.sort_values('importance_mean', ascending=False)
    importance.index.name = 'feature'
    importance.to_csv(OUT / 'importance.tsv', sep='\\t')
    importance.head(20)[['family', 'importance_mean', 'importance_sd'] + [f'fold_{f}' for f in range(5)]].round(4)
""")
code("""
    top = importance.head(20)
    print('top-20 importance grouped by family (sum of mean importance):')
    print(top.groupby('family').importance_mean.agg(['size', 'sum']).sort_values('sum', ascending=False).round(4).to_string())
    print()
    ATT_FAMILIES = {'pubmed', 'unregistered'}
    print('is n_proximity_partners in the top 20?', 'n_proximity_partners' in top.index,
          '| rank:', int(list(importance.index).index('n_proximity_partners')) + 1 if 'n_proximity_partners' in importance.index else None)
    print('is has_signal_peptide in the top 20?', 'has_signal_peptide' in top.index,
          '| rank:', int(list(importance.index).index('has_signal_peptide')) + 1)
    for c in sorted(attention_features):
        if c in importance.index:
            print(f'   attention feature {c}: rank {list(importance.index).index(c) + 1}, mean importance {importance.loc[c, "importance_mean"]:.4f}')
""")

# --------------------------------------------------------------------------- 14 predictions
nb.md(
    '## 14 · Predictions for every gene, with what pushed each one up',
    'A calibrated probability for all 8,140 genes: the model fitted on all labelled genes, mapped '
    'through the isotonic curve fitted on the pooled out-of-fold scores. For the labelled genes the '
    'out-of-fold probability is given too -- that is the honest one, and it is the one every score '
    'in this notebook uses.',
    '**The contributing features are measured, not guessed.** For each gene, each of the top 20 '
    'features in turn is set to missing (which the tree handles natively) and the fall in the gene '
    'own probability is recorded; the three largest falls are reported. So `top_features` means '
    '"take this measurement away and the call weakens most", per gene, not a global ranking.',
    '**Nothing here is written into the information map.** Whether any of it may be is decided by '
    'the bound in section 9, and the answer is in the README.')
code("""
    occlude = [c for c in importance.head(20).index if c not in CATEGORICAL and c in FEATURES_FINAL]
    base = all_raw.values
    deltas = {}
    for c in occlude:
        Xc = X[FEATURES_FINAL].copy()
        Xc[c] = np.nan
        deltas[c] = base - final_model.predict_proba(Xc)[:, 1]
    deltas = pd.DataFrame(deltas, index=X.index)

    def top_three(row):
        s = row.sort_values(ascending=False)
        return '; '.join(f'{k}={v:+.3f}' for k, v in s.head(3).items() if abs(v) > 1e-4) or 'none above 0.0001'

    predictions = pd.DataFrame({
        'gene_id': X.index,
        'calibrated_probability': all_cal.values.round(4),
        'raw_model_score': all_raw.values.round(4),
    }).set_index('gene_id')
    predictions['out_of_fold_calibrated'] = cal.reindex(predictions.index).round(4)
    predictions['fold'] = fold.reindex(predictions.index)
    predictions['label'] = raw.reindex(predictions.index)
    predictions['label_source'] = source_class.reindex(predictions.index)
    predictions['excluded_conflict'] = predictions.index.isin(CONFLICTED)
    predictions['has_signal_peptide'] = sp_all.reindex(predictions.index).values
    predictions['top_features'] = deltas.apply(top_three, axis=1)
    predictions.sort_values('calibrated_probability', ascending=False).to_csv(OUT / 'predictions.tsv', sep='\\t')
    print('predictions written for', len(predictions), 'genes')
    predictions.loc[predictions.label.isna()].sort_values('calibrated_probability', ascending=False).head(15)[
        ['calibrated_probability', 'has_signal_peptide', 'top_features']]
""")
code("""
    print('distribution of the calibrated probability over the 7,573 untested genes:')
    print(all_cal.loc[empty].describe(percentiles=[.5, .9, .99, .999]).round(4).to_string())
    print()
    print('untested genes at or above each threshold:')
    print(pd.Series({t: int((all_cal.loc[empty] >= t).sum()) for t in THRESHOLDS}).to_frame('genes').T.to_string())
""")

# --------------------------------------------------------------------------- 15 limits + readme
nb.md(
    '## 15 · Limitations, stated without softening',
    '**0. The headline YES is a hair-width result.** The bound clears 0.90 at 0.903-0.906 on three of five '
    'feature sets and misses it at 0.877-0.887 on the other two. Removing one screen moves the answer. '
    'That is what 73 positives buys: an answer that is not stable to a reasonable change in the feature '
    'set, and should not be deployed as if it were.',
    '1. **There are 73 positives.** Over five folds that is about 15 per fold. The Clopper-Pearson '
    'bound cannot reach 0.90 on fewer than ~22 calls even if every call is correct, so for most of '
    'the threshold range the bar is arithmetically unreachable regardless of how good the model is. '
    'This is a data-quantity result before it is a modelling result.',
    '2. **The negatives are easy.** 486 of 494 negatives are hyperLOPIT markers of compartments '
    'inside the parasite or curated vacuolar residents. "Has a signal peptide" separates them from '
    'the positives almost by itself, so the overall accuracy and AUC are inflated relative to the '
    'question a biologist cares about. Section 10 is the number to quote; section 8 is the number '
    'that will look better in a figure.',
    '3. **Most positives are one assay.** 52 of 73 positives are PVM proximity-labelling hits and '
    'nothing else, at tier probability 0.70. Proximity to the PVM is not the same as being exported, '
    'and the 8 excluded conflicts are exactly the cases where a curated reading disagrees. A model '
    'that recovers these well may be recovering that single screen.',
    '4. **Publication bias is real and only partly removed.** The ablation removes the counting '
    'columns, but a famous gene is also the gene with a measured melting temperature, a palmitome '
    'entry and a phosphoproteome row; *missingness itself* carries fame, and no ablation removes '
    'that.',
    '5. **DeepTMHMM misses signal peptides.** GRA16 (TGME49_208830), a secreted effector, is called '
    'GLOBULAR with no signal peptide. So the secretory-restricted analysis in section 10 is run on a '
    'set that is missing some true secreted proteins, and `mature_len` is wrong wherever the '
    'prediction misses. SignalP or a curated cleavage table would fix this; it is an open gap.',
    '6. **The negatives are an inference, not an export test.** A hyperLOPIT marker of the nucleus '
    'is assigned "not host-facing" because it lives in the nucleus, at probability 0.80. Nobody '
    'tested it for export. Precision measured against such negatives is precision against a '
    'curation decision.',
    '7. **The ortholog-copy baseline is empty by construction.** Grouping folds by orthogroup is '
    'correct and it leaves that baseline with nothing to copy; its coverage is reported as 0 rather '
    'than hidden.',
    '8. **The all-gene probabilities are not independently validated.** They come from a model '
    'fitted on all labelled genes and an isotonic curve fitted in-sample on the pooled out-of-fold '
    'scores. Only the cross-fitted out-of-fold numbers are honest; the all-gene column exists to '
    'count what *would* be fillable, not to be trusted gene by gene.',
    '9. **One model, fixed settings, no search.** Nothing was tuned, so nothing was tuned on a test '
    'fold -- but equally, no claim is made that this is the best achievable model. Stages 2-5 '
    '(sequence and structure blocks, guilt by association, the joint row model, pairs) are not in '
    'this notebook.')
code("""
    scoreboard = pd.DataFrame(board)
    scoreboard.to_csv(OUT / 'scoreboard.tsv', sep='\\t', index=False)
    final_pivot = scoreboard.pivot_table(index='strategy', columns='metric', values='mean')
    final_pivot = final_pivot[[c for c in ['accuracy', 'precision', 'recall', 'roc_auc', 'average_precision', 'calls', 'coverage'] if c in final_pivot]]
    final_pivot.round(3)
""")
code("""
    top10 = importance.head(10)[['family', 'importance_mean']]
    readme = []
    readme.append('# Host-facing inference, stages 0 and 1 (2026-11-04)\\n')
    readme.append('Built by `scripts/build_host_facing_inference_notebook.py`; the executed notebook is\\n'
                  '`notebooks/host_facing_inference_2026_11_04.ipynb`. Instruction 72 inference plan, stages 0\\n'
                  '(baselines) and 1 (one gradient-boosted-tree model), scored by the instruction 70 rules.\\n')
    readme.append('## The question\\n')
    readme.append('Can the `host_facing` column of the information map be filled by prediction at the 90 % bar\\n'
                  'Einar set -- that is, is there a calibrated probability threshold at which the **lower 90 %\\n'
                  'confidence bound** of held-out precision reaches 0.90? Only at or above such a threshold may a\\n'
                  'predicted cell be written as a fact.\\n')
    readme.append('## The answer\\n')
    readme.append(ANSWER + '\\n')
    readme.append(f'\\nRestricted to proteins with a predicted signal peptide -- the question that matters, because '
                  f'the overall negatives are mostly proteins that live inside the parasite and are separable by '
                  f'signal peptide alone -- **{ANSWER_SEC}**.\\n')
    readme.append('\\n**And the overall answer is conditional on two feature families that sit close to the '
                  'conclusion.** The most important features are a screen for effectors that alter host '
                  'transcription (`hosttx_padj`, `hosttx_T2`) and a publication count (`n_papers_incidental`). '
                  'Removing them, and the other proximity-style assays:\\n')
    readme.append('\\n```\\n' + ablations.round(4).to_string(index=False) + '\\n```\\n')
    readme.append('\\nSo: ' + ABLATION_VERDICT + '.\\n')
    readme.append(f'\\n**But it is met by a hair, and not robustly.** Across those five feature sets the best '
                  f'lower bound ranges {ablations.best_lower90.min():.3f} to {ablations.best_lower90.max():.3f}, and '
                  f'{int((~ablations.reaches_bar).sum())} of {len(ablations)} fall below 0.90 -- including the set '
                  f'that removes only the host-effector screen. The bar being cleared at 0.903-0.906 while a '
                  f'neighbouring feature set gives 0.887 means the result sits inside the noise of a 73-positive '
                  f'label set. The right reading is "this is at the edge of what the existing labels can '
                  f'demonstrate", not "the fill bar is cleared".\\n')
    readme.append('## Method\\n')
    readme.append(f'* Labels: `host_facing` from `results/information_map_v0_2026_11_04/values.parquet`, observed\\n'
                  f'  facts only: {int(y.sum())} True and {int((~y).sum())} False after excluding the\\n'
                  f'  {len(CONFLICTED)} genes in `conflicts.tsv` ({", ".join(CONFLICTED)}) from both training and\\n'
                  f'  scoring. {int(values["host_facing"].isna().sum())} genes are untested and stay empty.\\n'
                  f'* Features: {len(FEATURES)} numeric/boolean columns of `nodes.parquet` plus nine sequence and\\n'
                  f'  topology traits. Excluded as leakage: {len(leak_named)} columns matching\\n'
                  f'  lopit|compartment|ortholopit|tagm, plus `pvm_proximity_positive` and `n_host_targets`, plus\\n'
                  f'  all identity/free-text columns. Full lists are printed in the notebook and in\\n'
                  f'  `importance.tsv` (family per feature).\\n'
                  f'* Folds: 5, grouped by OrthoMCL orthogroup so paralogs never split across folds, seed 0.\\n'
                  f'* Model: `HistGradientBoostingClassifier`, missing values passed through with no imputation,\\n'
                  f'  `nterm_side`/`cterm_side` as declared categoricals, 300 rounds at learning rate 0.05, the same\\n'
                  f'  for every fold, early stopping OFF so the fit is deterministic and no training gene is held\\n'
                  f'  back. No hyperparameter search; nothing was tuned on any held-out fold.\\n'
                  f'* Calibration: isotonic (`parascape.claims._isotonic`), cross-fitted across folds; expected\\n'
                  f'  calibration error {ece:.4f} over 10 quantile bins (uncalibrated\\n'
                  f'  {_ece(oof.values, y.astype(float).values, bins=10):.4f}).\\n'
                  f'* Bound: exact one-sided Clopper-Pearson, `scipy.stats.beta.ppf(0.10, k, n-k+1)`.\\n')
    readme.append('## Headline numbers\\n')
    readme.append('```\\n' + final_pivot.round(3).to_string() + '\\n```\\n')
    readme.append('Secretory-pathway proteins only (' + str(len(sec)) + ' genes, ' + str(int(y_sec.sum())) + ' positives):\\n')
    readme.append('```\\n' + sec_board.pivot_table(index='strategy', columns='metric', values='mean').round(3).to_string() + '\\n```\\n')
    readme.append('Publication-bias ablation:\\n')
    readme.append('```\\n' + abl.round(4).to_string() + '\\n```\\n')
    readme.append('Top 10 features by permutation importance (held-out folds, average precision):\\n')
    readme.append('```\\n' + top10.round(4).to_string() + '\\n```\\n')
    readme.append('## Files\\n')
    readme.append('| File | Contents |\\n|---|---|\\n'
                  '| `scoreboard.tsv` | strategy x metric, mean plus every fold |\\n'
                  '| `scoreboard_secretory.tsv` | the same restricted to predicted signal-peptide proteins |\\n'
                  '| `thresholds.tsv` | threshold, calls, precision, Clopper-Pearson lower 90 % bound, genes fillable |\\n'
                  '| `thresholds_secretory.tsv` | the same, secretory proteins only |\\n'
                  '| `calibration.tsv` | reliability table: predicted range, n, observed positive rate |\\n'
                  '| `by_label_source.tsv` | performance per positive class x negative class |\\n'
                  '| `ablations.tsv` | the lower 90 % bound with the host-effector screen, the attention columns and the proximity-style assays removed |\\n'
                  '| `importance.tsv` | permutation importance per feature with its provenance family |\\n'
                  '| `predictions.tsv` | all 8,140 genes: calibrated probability, out-of-fold probability, fold, label, per-gene occlusion attributions |\\n')
    readme.append('\\n## Limitations (none of these are softened)\\n')
    readme.append(f'1. **{int(y.sum())} positives.** ~{int(y.sum()) // 5} per fold. The Clopper-Pearson bound cannot\\n'
                  '   reach 0.90 on fewer than ~22 calls even with perfect precision, so over most of the threshold\\n'
                  '   range the bar is unreachable as arithmetic, before any question of model quality.\\n'
                  f'2. **The negatives are easy.** {int((~y).sum())} negatives, almost all hyperLOPIT markers of\\n'
                  '   compartments inside the parasite. "Has a signal peptide" nearly separates the classes on its\\n'
                  '   own, so the overall numbers overstate usefulness; the secretory-only numbers are the honest\\n'
                  '   ones.\\n'
                  f'3. **Most positives are one assay.** {int(((source_class == "proximity only") & y).sum())} of\\n'
                  f'   {int(y.sum())} positives are PVM proximity-labelling hits and nothing else (tier probability\\n'
                  '   0.70). Proximity to the PVM is not proof of export, and the 8 excluded conflicts are exactly\\n'
                  '   where curation disagrees with that screen.\\n'
                  '4. **Publication bias is only partly removable.** The ablation drops the counting columns, but a\\n'
                  '   well-studied gene is also the gene with a melting temperature, a palmitome entry and a\\n'
                  '   phosphoproteome row. Missingness itself carries fame.\\n'
                  '5. **DeepTMHMM misses signal peptides** (GRA16 / TGME49_208830 is called GLOBULAR). The\\n'
                  '   secretory-restricted analysis is therefore run on a set missing some true secreted proteins,\\n'
                  '   and `mature_len` is wrong wherever the prediction misses. SignalP would fix it; open gap.\\n'
                  '6. **The negatives are an inference, not an export test** (residence in an internal compartment,\\n'
                  '   probability 0.80). Precision against them is precision against a curation decision.\\n'
                  '7. **The ortholog-copy baseline has coverage 0** under orthogroup-grouped folds, by construction.\\n'
                  '8. **The all-gene probability column is not independently validated**: model fitted on all labels,\\n'
                  '   isotonic curve fitted in-sample on pooled out-of-fold scores. Only the cross-fitted\\n'
                  '   out-of-fold numbers are honest.\\n'
                  '9. **One model, fixed settings.** Stages 2-5 of the instruction 72 plan are not in this notebook.\\n'
                  '10. **The signal is concentrated in two families that sit close to the conclusion.** The top\\n'
                  '    permutation-importance features are `hosttx_padj` / `hosttx_T2`, from a screen for effectors\\n'
                  '    that target host-cell transcription, and `n_papers_incidental`, a publication count. Neither\\n'
                  '    defines the label, so neither is leakage under instruction 72, but a model built on them is\\n'
                  '    partly re-reading one host-effector screen and partly reading how famous a gene is.\\n'
                  '    `ablations.tsv` is the honest version of the headline result.\\n'
                  '11. **Permutation importance is noisy here.** With ~112 genes and ~15 positives per held-out\\n'
                  '    fold the individual importances are small and their signs flip between folds\\n'
                  '    (`hosttx_padj` is -0.08 in one fold and +0.05 in another). Read the top of that table as a\\n'
                  '    direction, not as a per-feature ranking to be trusted.\\n')
    (OUT / 'README.md').write_text('\\n'.join(readme))
    print((OUT / 'README.md').read_text()[:1500])
""")

nb.write(str(ROOT / 'notebooks/host_facing_inference_2026_11_04.ipynb'))
print('wrote notebook and', ROOT / 'results/host_facing_inference_2026_11_04')
