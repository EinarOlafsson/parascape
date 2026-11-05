"""Instruction 72, target trait 2 (FUNCTION): scoring function prediction on the information map.

The labels are `results/function_truth_2026_11_04/` -- 242 genes, 15 coarse classes, 57 % of the
genes carrying more than one. The truth-set README already argues that no single class can be
scored on its own (largest non-screen class 37 genes; a Clopper-Pearson lower 90 % bound needs
~65 held-out calls even at a true precision of 0.95), so the headline here is the POOLED
multi-label question: over all (gene, class) pairs, how precisely can a class be assigned?

Three lessons from the two host-facing runs shape this notebook and are stated in it:
  * a wrong label set flatters a model far more than a leaky feature does;
  * a feature must be reported against the hard comparison, never against "all genes";
  * "cannot yet be scored" is a valid and valuable answer.

Writes notebooks/function_inference_2026_11_04.ipynb (executed) and the tables under
results/function_inference_2026_11_04/.
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

nb = ExecutedNotebook('Function inference: the pooled multi-label question on the 242-gene truth set')


def code(src: str):
    return nb.code(dedent(src).strip('\n'))


# --------------------------------------------------------------------------- 0 intro
nb.md(
    '# Function inference (2026-11-04) — the pooled multi-label question',
    'Instruction 72, target trait **2 of 3** (host-facing was scored twice; binding partners come '
    'later). Labels: `results/function_truth_2026_11_04/` — **242 genes, 15 coarse classes, 138 of '
    'the genes (57 %) multi-label**, at three tiers (32 experimental GO, 113 curated Swiss-Prot '
    'FUNCTION comments, 97 arrayed-imaging-screen phenotypes).',
    '**Why the headline is pooled and not per class.** The truth-set README does the arithmetic: a '
    'Clopper–Pearson lower 90 % bound reaches 0.90 only with about **65 held-out calls** even at a '
    'true precision of 0.95, and the largest non-screen class has **37** genes — before splitting '
    'into folds. So no class can be scored alone. What the set can support is the question *over all '
    '(gene, class) pairs at once*: **how precisely can a class be assigned?** Per-class numbers are '
    'reported anyway, each with its call count, so it is visible that the classes are under-powered '
    'rather than silently pooled away.',
    '**Three lessons this notebook is built around**, learned the expensive way in the two '
    'host-facing runs:',
    '1. **A wrong label set flatters a model far more than a leaky feature does.** Host-facing run 1 '
    'reported a lower bound of 0.906 and it came from counting 22 vacuole-resident dense granule '
    'proteins as positives; putting the leaky features back on the corrected labels was worth only '
    '+0.013. No feature audit catches that — only human review. Here the labels are curator text and '
    'experimental GO, so the analogous risk is the **class rules**, and §1 shows what they contain.',
    '2. **A feature must be reported against the hard comparison, never against "all genes".** '
    '`domain_novel` looked spectacular against all 8,140 genes and vanished against the negatives '
    'that answered the question. For function the hard comparisons are (a) the **CAFA no-knowledge** '
    'genes — no prior annotation of any kind in the feature set — and (b) the **non-screen tiers** '
    'scored apart from the 97 screen-only genes.',
    '3. **"Cannot yet be scored" is a valid, valuable answer.** The fill bar is a **lower 90 % '
    'confidence bound** of held-out precision at 0.90, not a point estimate. Nothing is written into '
    'the information map unless the bound clears it.',
    'A reader who wants only the answer should read §8 (pooled headline), §9 (per class), §10 '
    '(no-knowledge), §11 (tier split) and §15 (what is needed).')

# --------------------------------------------------------------------------- 1 labels
nb.md(
    '## 1 · Labels: 242 genes × 15 classes, as a long table of (gene, class) pairs',
    'The unit of scoring is a **(gene, class) pair**. A pair is positive when the truth set assigns '
    'that class to that gene, and negative when the gene is labelled but that class is not among its '
    'classes. **Unlabelled genes are not used as negatives anywhere** — they are the blanks the map '
    'is trying to fill, and treating them as negatives is the error instruction 70 forbids.',
    'That leaves one honest caveat which is stated here and repeated in the limitations: a negative '
    'pair means *"the curator text and the experimental GO annotations of this gene do not mention '
    'this class"*, not *"this gene does not do this"*. Function annotation is open-world and '
    'incomplete; a model that correctly predicts an unrecorded function is scored as wrong. This '
    'biases precision **downward**, so a bound that clears the bar is still trustworthy, while a '
    'bound that fails may be pessimistic. CAFA has the same property.',
    'Every pair keeps its **tier** (`observed-direct` 0.95 / `observed-curated` 0.90 / '
    '`observed-screen` 0.80) and its source, because the tiers are not equally hard and §11 scores '
    'them separately.')
code("""
    import os
    # One thread per BLAS/OpenMP pool: the tables are ~242 x 700 and thread spin-up costs more than
    # the arithmetic saves. The machine is shared.
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[v] = '1'
    import json, re, warnings
    from pathlib import Path
    import numpy as np, pandas as pd
    from scipy.stats import beta
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    ROOT = Path('.').resolve()
    TRUTH = ROOT / 'results/function_truth_2026_11_04'
    SEQ = ROOT / 'results/sequence_traits_2026_11_04/sequence_traits.parquet'
    DOM = ROOT / 'results/domain_traits_2026_11_04/domain_traits.parquet'
    CX = ROOT / 'results/coexpression_localization_2026_11_04/coexpression_neighbours.parquet'
    NODES = ROOT / 'parascape/data/nodes.parquet'
    OUT = ROOT / 'results/function_inference_2026_11_04'; OUT.mkdir(parents=True, exist_ok=True)

    truth = pd.read_csv(TRUTH / 'function_truth.tsv', sep='\\t').set_index('gene_id')
    long = pd.read_csv(TRUTH / 'function_classes_long.tsv', sep='\\t')
    ic = pd.read_csv(TRUTH / 'label_information_content.tsv', sep='\\t')
    coverage = pd.read_csv(TRUTH / 'coverage.tsv', sep='\\t')

    CLASSES = sorted(coverage.cls)
    GENES = list(truth.index)
    Y = pd.DataFrame(False, index=GENES, columns=CLASSES)
    for r in long.itertuples():
        if r.gene_id in Y.index and r.cls in Y.columns:
            Y.loc[r.gene_id, r.cls] = True
    IC = ic[ic.kind == 'coarse_class'].set_index('label').ic_corpus.reindex(CLASSES)

    print('labelled genes: %d | classes: %d | (gene, class) pairs: %d'
          % (len(GENES), len(CLASSES), Y.size))
    print('positive pairs: %d (prevalence %.4f) | multi-label genes: %d (%.0f %%)'
          % (int(Y.values.sum()), Y.values.mean(), int((Y.sum(axis=1) > 1).sum()),
             100 * (Y.sum(axis=1) > 1).mean()))
    print()
    print('per-gene tier (one tier per gene, the deciding one):')
    print(truth.tier.value_counts().to_string())
    print()
    tab = pd.DataFrame({'genes': Y.sum(), 'ic_corpus': IC}).sort_values('genes', ascending=False)
    tab['prevalence'] = (tab.genes / len(GENES)).round(4)
    print('classes, their size and their information content (a vague class must not score like a specific one):')
    print(tab.to_string())
""")
code("""
    # Which tiers carry which class: the screen supports exactly four process classes and nothing else.
    src_tier = long.groupby(['cls', 'tier']).gene_id.nunique().unstack(fill_value=0)
    print('class x tier (genes):')
    print(src_tier.to_string())
    print()
    SCREEN_GENES = list(truth.index[truth.tier == 'observed-screen'])
    NONSCREEN_GENES = [g for g in GENES if g not in set(SCREEN_GENES)]
    print('screen-only genes: %d | non-screen (experimental GO + curated Swiss-Prot): %d | overlap: %d'
          % (len(SCREEN_GENES), len(NONSCREEN_GENES),
             len(set(SCREEN_GENES) & set(NONSCREEN_GENES))))
    SCREEN_CLASSES = sorted(long.loc[long.tier == 'observed-screen', 'cls'].unique())
    print('classes the screen can support:', SCREEN_CLASSES)
    print('positive pairs at screen-only genes: %d of %d (%.0f %% of all positives)'
          % (int(Y.loc[SCREEN_GENES].values.sum()), int(Y.values.sum()),
             100 * Y.loc[SCREEN_GENES].values.sum() / Y.values.sum()))
    print()
    print('classes per gene by tier:')
    print(pd.DataFrame({'classes_per_gene': Y.sum(axis=1)}).join(truth.tier)
          .groupby('tier').classes_per_gene.agg(['size', 'mean', 'max']).round(2).to_string())
""")
nb.md(
    '### The label risk here is the class rules, not the gene list',
    'Lesson 1 says the labels are where a flattering result comes from. For function the labels are '
    'not a human ruling gene by gene; they are **regular expressions over curator text plus GO '
    'descendant sets** (`class_rules.md`, printed verbatim in the truth set). Two consequences are '
    'recorded rather than discovered later:',
    '* **A class can be a restatement of a name.** `metabolism` matches `\\bsynthase\\b`, '
    '`\\btransferase\\b`, `\\bdehydrogenase\\b`; `signalling` matches `\\bkinase\\b`. The curator text '
    'and the ToxoDB product name and the InterPro domain of the same protein all tend to carry the '
    'same word. So a domain column is **not** an independent predictor of such a class — it is close '
    'to the same annotation act. That is why §12 reports the bound again with every annotation-derived '
    'family removed, and why the no-knowledge split in §10 is the number that matters for filling '
    'blanks.',
    '* **The screen tier is a different kind of label.** `observed-screen` genes carry only the four '
    'phenotype classes, assigned from four boolean columns of one imaging screen. Those columns are '
    'the label source and are dropped from the features throughout (§2). But the screen also chose '
    '**which 319 genes to look at**, so "was in the screen" is itself informative; §11 scores the 97 '
    'screen-only genes apart from the 145 others for exactly that reason.')

# --------------------------------------------------------------------------- 2 features
nb.md(
    '## 2 · Features, and everything held out',
    'Tier 2 of instruction 72 means **measurements**: expression, fitness, proteome, PTMs, screens, '
    'sequence, structure, domains. Tier 1 are the **conclusions** — function, localization, '
    'host-facing, partners — and they are what is predicted, so they are not features.',
    '**Held out, non-negotiable:**',
    '* the whole `splitcas9_imaging_screen` family — `screen_egress_phenotype`, '
    '`screen_actin_phenotype`, `screen_apicoplast_phenotype`, `screen_replication_phenotype`, '
    '`screen_any_phenotype`, `screen_scorers_agree`. These four phenotype columns **are** the source '
    'of the 97 screen-tier labels. **They are dropped throughout**, for every gene and every '
    'question, not only when screen genes are scored. Dropping them everywhere rather than '
    'conditionally is the simpler and stricter of the two options instruction 72 allows: it means no '
    'result in this notebook, including the non-screen ones, can be read as circular through them.\n'
    '* every `lopit*`, `compartment*`, `ortholopit*`, `tagm*` column — localization is a **tier-1 '
    'conclusion**, and it is also a direct restatement of some classes (a hyperLOPIT apicoplast '
    'marker is an `apicoplast_mitochondrion` gene almost by definition). Keeping them would be the '
    'function analogue of host-facing run 1.\n'
    '* `hosttx_effectors` (24 columns), `pvm_proximity_positive`, `n_host_targets` — a screen for '
    'parasite genes whose loss changes host transcription, and the PVM proximity screen. Both are '
    'nearly the `host_cell_modulation` class itself.\n'
    '* every free-text and identity column: `sequence`, `product`, `interpro_id`, `interpro_desc`, '
    '`pfam_id`, `ec_number`, `orthogroup`, `gene_names`, the resistance/compound strings. The class '
    'rules match curator text, and `product` is the same kind of text; using it would score agreement '
    'between two annotations of one protein.',
    '**Kept, in five families**, and every one of them is printed below with its family: the numeric '
    'and boolean `nodes.parquet` columns that survive the filter; the DeepTMHMM sequence/topology '
    'traits (`sequence_traits`); the **313 domain columns** (`domain`); and the **co-expression '
    'neighbour-label features** (`coexpression`), built in §3.')
code("""
    nodes = pd.read_parquet(NODES).set_index('gene_id')
    seqt = pd.read_parquet(SEQ)
    dom = pd.read_parquet(DOM)

    # object columns holding only booleans (the phenotype screens) -> float 1/0/NaN so a tree can
    # branch on "not screened" instead of guessing a value.
    for c in nodes.columns:
        if nodes[c].dtype == object and set(map(type, nodes[c].dropna())) <= {bool, np.bool_}:
            nodes[c] = nodes[c].astype('boolean').astype('float64')

    numeric = [c for c in nodes.columns
               if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]
    text_cols = [c for c in nodes.columns if c not in numeric]

    from parascape import datasets as ds
    def family_of(col):
        d = ds.provenance(col, 'Tg')
        return d.key if d is not None else 'unregistered'
    node_family = pd.Series({c: family_of(c) for c in nodes.columns})

    LABEL_SOURCE_FAMILIES = ['splitcas9_imaging_screen']          # the 97 screen-tier labels
    NEAR_CONCLUSION_FAMILIES = ['hosttx_effectors', 'pvm_proximity']
    EXPLICIT = ['n_host_targets', 'pvm_proximity_positive']
    NAME_PAT = re.compile('lopit|compartment|ortholopit|tagm|^screen_', re.I)

    held_label_src = sorted(c for c in nodes.columns if node_family.get(c) in LABEL_SOURCE_FAMILIES)
    held_near = sorted(c for c in nodes.columns if node_family.get(c) in NEAR_CONCLUSION_FAMILIES)
    held_explicit = sorted(c for c in nodes.columns if c in EXPLICIT)
    held_named = sorted(c for c in nodes.columns if NAME_PAT.search(c))
    HELD = set(held_label_src) | set(held_near) | set(held_explicit) | set(held_named)

    print('nodes.parquet: %d columns, %d numeric/boolean, %d free text/identity' %
          (nodes.shape[1], len(numeric), len(text_cols)))
    for name, cols in [('LABEL SOURCE - the splitCas9 imaging screen that defines the 97 screen-tier labels', held_label_src),
                       ('near-conclusion families (hosttx_effectors, pvm_proximity)', held_near),
                       ('columns that restate a conclusion', held_explicit),
                       ('matching lopit|compartment|ortholopit|tagm|^screen_ (tier-1 localization + the screen)', held_named)]:
        print('HELD OUT - %s (%d):' % (name, len(cols)))
        for i in range(0, len(cols), 3):
            print('      ' + '  '.join(f'{c:<38}' for c in cols[i:i + 3]).rstrip())
    print()
    print('held out as free text / identity (never features):')
    for i in range(0, len(text_cols), 3):
        print('      ' + '  '.join(f'{c:<38}' for c in text_cols[i:i + 3]).rstrip())
    print()
    print('TOTAL held out: %d columns (%d of them numeric)' % (len(HELD), len(HELD & set(numeric))))
""")
code("""
    SEQ_FEATURES = ['n_tm', 'has_signal_peptide', 'seq_len', 'mature_len', 'nterm_len', 'cterm_len',
                    'nterm_side', 'cterm_side', 'flanks_identical',
                    'tdi_coverage', 'tdi_n_fragments', 'tdi_align_identity', 'tdi_align_coverage',
                    'tdi_source']
    CATEGORICAL = ['nterm_side', 'cterm_side', 'tdi_source']
    DOMAIN_FEATURES = list(dom.columns)

    node_features = [c for c in numeric if c not in HELD and c not in SEQ_FEATURES]
    BASE = (nodes[node_features].astype('float64')
            .join(seqt[SEQ_FEATURES])
            .join(dom.astype('float64')))
    for c in ('has_signal_peptide', 'flanks_identical'):
        BASE[c] = BASE[c].astype('float64')
    for c in CATEGORICAL:
        BASE[c] = BASE[c].astype('category')

    all_missing = [c for c in BASE.columns
                   if c not in CATEGORICAL and BASE.loc[GENES, c].notna().sum() == 0]
    BASE = BASE.drop(columns=all_missing)
    print('dropped as never measured on any labelled gene (%d):' % len(all_missing), all_missing[:12],
          '...' if len(all_missing) > 12 else '')
    print('base feature block: %d columns over %d genes' % (BASE.shape[1], BASE.shape[0]))
""")

# --------------------------------------------------------------------------- 3 coexpression
nb.md(
    '## 3 · Co-expression neighbour features — guilt by association, as its own family',
    'Instruction 72 stage 3 is guilt by association. `coexpression_neighbours.parquet` gives, per '
    'gene, its strongest co-expression partners with a `combined_percentile` across studies and a '
    '`study_support` count. The features built from it are the **classes of a gene’s strongest '
    'partners**: for each class, the similarity-weighted fraction of a gene’s labelled partners that '
    'carry it, and the similarity of the best partner that carries it; plus four label-free summaries '
    '(number of partners, best and mean similarity, mean study support). Top **20** partners per '
    'gene.',
    '**These features use labels, so they are built inside the fold loop.** For fold *f* the partner '
    'labels come only from genes in the **training** folds; a test gene never sees the label of a '
    'partner that is in its own fold. Getting this wrong would be the most direct leak available in '
    'the whole notebook — a gene’s paralog or operon-mate handing over the answer — which is also why '
    'the folds are grouped by orthogroup (§4).',
    'They are kept as the separate family **`coexpression`** so §12 can remove them as a block, '
    'because guilt by association is a different claim from a measurement: it says *this gene '
    'resembles genes that do X*, which is evidence about neighbours, not about the gene.')
code("""
    cx = pd.read_parquet(CX)
    cx = cx[cx.gene.isin(BASE.index) & cx.partner.isin(BASE.index)]
    K = 20
    top = (cx.sort_values('combined_percentile', ascending=False)
             .groupby('gene', sort=False).head(K).copy())
    print('co-expression rows: %d over %d genes (%.1f partners per gene, capped at K=%d)'
          % (len(top), top.gene.nunique(), len(top) / top.gene.nunique(), K))
    print('combined_percentile of the kept partners: median %.3f, min %.3f'
          % (top.combined_percentile.median(), top.combined_percentile.min()))
    print('labelled genes with at least one co-expression partner: %d of %d'
          % (len(set(top.gene) & set(GENES)), len(GENES)))

    CX_LABEL_FREE = ['cx_n_partners', 'cx_top_sim', 'cx_mean_sim', 'cx_mean_support']
    CX_FEATURES = CX_LABEL_FREE + ['cx_n_labelled'] + \\
        ['cx_frac_' + c for c in CLASSES] + ['cx_best_' + c for c in CLASSES]

    base_cx = top.groupby('gene').agg(cx_n_partners=('partner', 'size'),
                                      cx_top_sim=('combined_percentile', 'max'),
                                      cx_mean_sim=('combined_percentile', 'mean'),
                                      cx_mean_support=('study_support', 'mean'))

    def coexpression_block(label_genes):
        '''Neighbour-label features using ONLY the labels of `label_genes`.'''
        lab = Y.loc[[g for g in label_genes if g in Y.index]]
        m = top.merge(lab, left_on='partner', right_index=True, how='inner')
        out = base_cx.copy()
        if not len(m):
            for c in CLASSES:
                out['cx_frac_' + c] = np.nan
                out['cx_best_' + c] = np.nan
            out['cx_n_labelled'] = 0
            return out.reindex(BASE.index)[CX_FEATURES]
        w = m.combined_percentile
        wsum = w.groupby(m.gene).sum()
        out['cx_n_labelled'] = m.groupby('gene').size()
        for c in CLASSES:
            out['cx_frac_' + c] = (w * m[c].astype(float)).groupby(m.gene).sum() / wsum
            hit = m[m[c]]
            out['cx_best_' + c] = hit.groupby('gene').combined_percentile.max() if len(hit) else np.nan
        out['cx_n_labelled'] = out['cx_n_labelled'].fillna(0)
        return out.reindex(BASE.index)[CX_FEATURES]

    FEATURES = list(BASE.columns) + CX_FEATURES
    fam = pd.Series({c: ('coexpression' if c in CX_FEATURES
                         else 'domain' if c in DOMAIN_FEATURES
                         else 'sequence_traits' if c in SEQ_FEATURES
                         else node_family.get(c, 'unregistered')) for c in FEATURES}, name='family')
    fam.index.name = 'feature'
    print()
    print('co-expression features: %d (%d label-free, %d neighbour-label)'
          % (len(CX_FEATURES), len(CX_LABEL_FREE), len(CX_FEATURES) - len(CX_LABEL_FREE)))
    print('TOTAL FEATURE COUNT: %d | categorical: %s' % (len(FEATURES), CATEGORICAL))
""")
code("""
    print('THE FINAL FEATURE LIST, each with its provenance family.')
    print('A family is held out as a unit: 36 expression columns from one experiment are one piece of')
    print('evidence, not 36. The 313 domain columns are printed as a count, then in full.')
    print()
    for f, group in sorted(fam.groupby(fam), key=lambda kv: (-len(kv[1]), kv[0])):
        names = sorted(group.index)
        print('%s  (%d)' % (f, len(names)))
        for i in range(0, len(names), 4):
            print('      ' + '  '.join(f'{c:<34}' for c in names[i:i + 4]).rstrip())
    features_table = fam.to_frame()
    features_table['tier'] = 2
    features_table.to_csv(OUT / 'features.tsv', sep='\\t')
    fam.value_counts().rename('features').to_frame().T
""")
nb.md(
    '### Which family that is still in could encode a label? — asked in advance, not after the fact',
    'This is the question host-facing run 1 answered too late. Said now, for every kept family that '
    'stands close to a class:',
    '* **`domain` (313 columns)** — InterPro and Pfam family flags. The truth set refuses them as '
    '*labels* because they are predictions from sequence, and instruction 72 asks for them as '
    '*features*. Both are right, but they are **annotation-derived**: `dom_IPR000719` (protein '
    'kinase) and the `signalling` class both ultimately come from somebody recognising a kinase. This '
    'is not leakage of the label *values* — it is correlation between two annotation acts on the same '
    'protein. It is the single biggest reason §10 (no-knowledge genes) and §12 (annotation-derived '
    'families removed) exist, and a result that lives only in this family is a result about '
    'annotation agreement, not about biology.\n'
    '* **`toxodb_ec_numbers` (`has_ec_number`)** — an EC number is an enzyme classification from '
    'sequence similarity. For `metabolism` it is nearly the class. Kept, flagged, and removed in §12 '
    'with the domains.\n'
    '* **`coexpression`** — uses the labels of neighbours. Fold-safe by construction (§3), but it is '
    'guilt by association and gets its own family so it can be taken out.\n'
    '* **`pubmed`, `bioid_corpus_membership`, `ipms_corpus_membership`, `ipms_baits`, `pmc_oa`, '
    '`lit_tier`** — attention columns. Fame predicts *whether* a gene is labelled, which is why the '
    'host-facing runs held them out. Here every gene in the scoring set is already labelled, so fame '
    'is largely conditioned away and the columns stay; §12 removes them to check.\n'
    '* **`oxidative_stress_screen`, `differentiation_screen`, `brady_subtypes`, `sexual_stages`, '
    '`gse132248_stress`** — phenotype and condition screens that sit very close to '
    '`stress_differentiation`. Evidence, not label definitions: none of them assigned a class. Kept '
    'and watched in the importance table.\n'
    '* **`crispr_*` fitness** — the truth set is explicit that essentiality is not a function and '
    'assigns no class from it, so these are clean features.',
    'Nothing kept here is the *source* of a label. The one family that was — the splitCas9 imaging '
    'screen — is dropped throughout.')
code("""
    print('sanity check: are the four label-source screen columns really gone?')
    leaked = [c for c in FEATURES if c.startswith('screen_')]
    print('  feature names starting with screen_:', leaked if leaked else 'none')
    print('  features in the splitcas9_imaging_screen family:',
          [c for c in FEATURES if node_family.get(c) == 'splitcas9_imaging_screen'] or 'none')
    print()
    print('and what was removed, against the labels, so the size of the drop is visible:')
    scr = pd.read_parquet(NODES).set_index('gene_id')
    for col, cls in [('screen_egress_phenotype', 'egress'),
                     ('screen_actin_phenotype', 'cytoskeleton'),
                     ('screen_apicoplast_phenotype', 'apicoplast_mitochondrion'),
                     ('screen_replication_phenotype', 'cell_cycle_division')]:
        v = scr[col].reindex(GENES)
        v = v.map(lambda x: bool(x) if x is not None and x is not np.nan else np.nan)
        agree = (v.astype('float64') == Y[cls].astype(float)).mean()
        print('  %-32s -> %-26s identical on %.0f %% of labelled genes'
              % (col, cls, 100 * agree))
    print()
    print('localization (held out) against the apicoplast class, which is why it is held out:')
    comp = scr['compartment'].reindex(GENES)
    print(pd.crosstab(comp.isin(['apicoplast', 'mitochondrion']),
                      Y['apicoplast_mitochondrion']).to_string())
""")

# --------------------------------------------------------------------------- 4 folds
nb.md(
    '## 4 · Five folds grouped by orthogroup, seed 0 — the same construction as the host-facing runs',
    'Paralogs share an orthogroup and very often share a function; a paralog in the training fold and '
    'its twin in the test fold would let the model score by recognising the family. So whole '
    'orthogroups go to one fold (OrthoMCL groups from `MASTER_parasite_wide_by_orthogroup.csv`; a '
    'gene with no group is its own group). Seed 0, assigned once over all 242 labelled genes.',
    'The same fold assignment serves every class: one binary model per class on the same folds '
    '(one-vs-rest), so the held-out predictions of all 15 classes can be pooled over identical '
    'splits. The cost shows up in the ortholog-copy baseline, which is left with nothing to copy from '
    '(coverage near 0). That is correct: a baseline that works only because paralogs leak is not a '
    'baseline.')
code("""
    OG = Path('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/MASTER_parasite_wide_by_orthogroup.csv')
    og = pd.read_csv(OG, usecols=['orthomcl_group', 'tgon_gene_id']).dropna()
    group_of = {g: row.orthomcl_group for row in og.itertuples() for g in str(row.tgon_gene_id).split(';')}
    groups = pd.Series([group_of.get(g, 'self:' + g) for g in GENES], index=GENES, name='group')
    rng = np.random.default_rng(0)
    unique = np.array(sorted(set(groups))); rng.shuffle(unique)
    fold = groups.map({g: i % 5 for i, g in enumerate(unique)}).rename('fold')
    FOLDS = sorted(fold.unique())
    print('242 labelled genes in %d orthogroups; genes sharing an orthogroup with another labelled gene: %d'
          % (len(unique), int((groups.value_counts() > 1).reindex(groups).sum())))
    print()
    print('genes per fold, and positive pairs per fold:')
    print(pd.DataFrame({'genes': fold.value_counts().sort_index(),
                        'positive_pairs': Y.groupby(fold).sum().sum(axis=1),
                        'screen_only': pd.Series(fold.loc[SCREEN_GENES]).value_counts().sort_index(),
                        'non_screen': pd.Series(fold.loc[NONSCREEN_GENES]).value_counts().sort_index()}
                       ).to_string())
    print()
    print('positives per class per fold - this is the under-powering, visible:')
    print(Y.groupby(fold).sum().T.assign(total=Y.sum()).sort_values('total', ascending=False).to_string())
""")

# --------------------------------------------------------------------------- 5 machinery
nb.md(
    '## 5 · The machinery: model, pooling, the bound, and information-content weighting',
    '`HistGradientBoostingClassifier`, one per (class, fold). Missing values pass straight through — '
    'absence from a proteomics run is itself informative, so **no imputation anywhere**. Fixed '
    'settings for every class and fold (300 rounds, learning rate 0.05, 15 leaves, 10 samples per '
    'leaf, L2 1.0), early stopping off so the fit is deterministic. **No hyperparameter search**, so '
    'nothing was tuned on a held-out fold.',
    '**Pooling.** The held-out prediction for every (gene, class) pair is collected into one 242 × 15 '
    'matrix and flattened; precision, recall, ROC AUC and average precision are computed over the '
    'pooled pairs. That is the headline question: *over all (gene, class) pairs, how precisely can a '
    'class be assigned?*',
    '**The bound.** Exact one-sided **Clopper–Pearson** lower 90 % bound on precision, '
    '`beta.ppf(0.10, k, n−k+1)`. The fill rule is the bound, not the point estimate.',
    '**Information-content weighting.** Instruction 72 requires that a vague class not count as much '
    'as a specific one, so every figure is reported twice: unweighted, and with each call credited by '
    'the class corpus information content (`apicoplast_mitochondrion` 1.60 up to '
    '`transcription_chromatin` 4.92). The weighted precision is Σ w·correct / Σ w. A weighted count '
    'has no exact binomial bound, so its bound uses the **Kish effective sample size** '
    '(Σw)²/Σw² — an approximation, labelled as one wherever it appears.',
    '**Which is the headline?** The **unweighted** Clopper–Pearson bound is the number the fill '
    'decision turns on, because the fill rule in instruction 72 is a statement about *cells* — each '
    'written cell must be right with probability 0.90 — and only the unweighted bound is exact for '
    'counting cells. The IC-weighted precision is the headline for **how much the filled cells are '
    'worth**: it says whether a passing model passes by calling specific things or by calling '
    '`apicoplast_mitochondrion` on everything. Both are reported side by side at every threshold and '
    'the notebook never quotes one alone.')
code("""
    from sklearn.metrics import roc_auc_score, average_precision_score
    from sklearn.ensemble import HistGradientBoostingClassifier
    from parascape.claims import _isotonic, _ece

    def make_model():
        return HistGradientBoostingClassifier(
            max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=10,
            l2_regularization=1.0, early_stopping=False, random_state=0,
            categorical_features='from_dtype')

    def usable(frame, features, rows):
        '''Columns sklearn can bin: categorical, or at least two distinct measured values.'''
        out = []
        for c in features:
            if c in CATEGORICAL:
                out.append(c)
            elif frame.loc[rows, c].nunique(dropna=True) >= 2:
                out.append(c)
        return out

    # One feature frame per fold: the fixed block plus the co-expression block built from the
    # TRAINING labels of that fold only.
    XF, CXB = {}, {}
    for f in FOLDS:
        tr = [g for g in GENES if fold[g] != f]
        CXB[f] = coexpression_block(tr)
        XF[f] = BASE.join(CXB[f])
    CXB['all'] = coexpression_block(GENES)
    XF['all'] = BASE.join(CXB['all'])
    print('per-fold feature frames built: %s' % {k: XF[k].shape for k in XF})
    print('cx_n_labelled among labelled genes, per fold (partners whose label the model may use):')
    print(pd.DataFrame({f: CXB[f].loc[GENES, 'cx_n_labelled'].describe()[['mean', 'max']]
                        for f in FOLDS}).round(2).to_string())

    def clopper_pearson_lower(k, n, confidence=0.90):
        if n == 0:
            return np.nan
        if k == 0:
            return 0.0
        if k >= n:
            return float(beta.ppf(1 - confidence, k, 1))
        return float(beta.ppf(1 - confidence, k, n - k + 1))

    def weighted_precision_and_bound(correct, weights):
        '''IC-weighted precision, with a Kish-effective-n Clopper-Pearson bound (an approximation).'''
        w = np.asarray(weights, dtype=float)
        c = np.asarray(correct, dtype=float)
        if not len(w) or w.sum() == 0:
            return np.nan, np.nan, 0.0
        p = float((w * c).sum() / w.sum())
        n_eff = float(w.sum() ** 2 / (w ** 2).sum())
        k_eff = int(round(p * n_eff))
        return p, clopper_pearson_lower(k_eff, int(round(n_eff))), n_eff

    THRESHOLDS = [0.20, 0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97, 0.99]
    print()
    print('the arithmetic that decides most of this notebook - calls needed before the bound CAN reach 0.90:')
    print(pd.DataFrame([dict(calls=n, lower90_if_every_call_correct=round(clopper_pearson_lower(n, n), 3))
                        for n in (4, 8, 12, 16, 20, 22, 25, 30, 40, 60, 100)]).set_index('calls').T.to_string())
""")
nb.md(
    'With fewer than about **22 calls** the lower 90 % bound cannot reach 0.90 even if every single '
    'call is correct. The pooled question has 3,630 pairs and ~470 positives, which is the first '
    'reason it is the headline: unlike any single class, it can produce enough calls for the bound to '
    'be reachable at all. Whether the precision is there is a different matter, and that is what §8 '
    'measures.')
code("""
    def fit_oof(features, train_genes, score_genes=None, verbose=False):
        '''One-vs-rest out-of-fold scores: 15 classes x 5 folds on the shared folds.'''
        train_genes = [g for g in train_genes]
        target = list(score_genes) if score_genes is not None else train_genes
        oof = pd.DataFrame(np.nan, index=target, columns=CLASSES, dtype=float)
        models, used = {}, {}
        for f in FOLDS:
            tr = [g for g in train_genes if fold[g] != f]
            te = [g for g in target if fold[g] == f]
            if not tr or not te:
                continue
            Xf = XF[f]
            cols = usable(Xf, features, tr)
            used[f] = cols
            Xtr = Xf.loc[tr, cols]
            for c in CLASSES:
                ytr = Y.loc[tr, c]
                if ytr.nunique() < 2:
                    oof.loc[te, c] = float(ytr.mean())
                    continue
                m = make_model().fit(Xtr, ytr)
                models[(c, f)] = m
                oof.loc[te, c] = m.predict_proba(Xf.loc[te, cols])[:, 1]
        if verbose:
            print('   trained on %d genes, scored %d; columns per fold: %s'
                  % (len(train_genes), len(target), {f: len(c) for f, c in used.items()}))
        return oof, models, used

    def pairs_frame(score_df, genes=None, classes=None):
        '''Flatten (gene, class) pairs with score, truth, fold, tier and information content.'''
        genes = list(score_df.index) if genes is None else [g for g in genes if g in score_df.index]
        classes = CLASSES if classes is None else classes
        s = score_df.loc[genes, classes]
        out = s.stack().rename('score').reset_index()
        out.columns = ['gene_id', 'cls', 'score']
        out['truth'] = [bool(Y.loc[g, c]) for g, c in zip(out.gene_id, out.cls)]
        out['fold'] = out.gene_id.map(fold)
        out['tier'] = out.gene_id.map(truth.tier)
        out['ic'] = out.cls.map(IC)
        return out

    def pooled_metrics(pf, threshold=0.5):
        called = pf.score >= threshold
        tp = int((called & pf.truth).sum()); fp = int((called & ~pf.truth).sum())
        fn = int((~called & pf.truth).sum())
        wp, wlo, neff = weighted_precision_and_bound(pf.truth[called].values, pf.ic[called].values)
        d = dict(pairs=len(pf), positives=int(pf.truth.sum()), calls=tp + fp, correct=tp,
                 precision=(tp / (tp + fp)) if (tp + fp) else np.nan,
                 recall=(tp / (tp + fn)) if (tp + fn) else np.nan,
                 precision_lower90=clopper_pearson_lower(tp, tp + fp),
                 precision_ic_weighted=wp, precision_ic_lower90_approx=wlo, ic_effective_n=neff)
        if pf.truth.nunique() == 2:
            d['roc_auc'] = roc_auc_score(pf.truth, pf.score)
            d['average_precision'] = average_precision_score(pf.truth, pf.score)
        return d

    def best_bound(pf):
        '''The best lower 90 % bound anywhere on the threshold curve, and where it is.'''
        best = (np.nan, np.nan, 0, 0)
        for t in THRESHOLDS:
            sel = pf.score >= t
            n = int(sel.sum())
            if n == 0:
                continue
            k = int(pf.truth[sel].sum())
            lo = clopper_pearson_lower(k, n)
            if not (best[0] == best[0]) or lo > best[0]:
                best = (lo, t, n, k)
        return best

    print('machinery ready')
""")

# --------------------------------------------------------------------------- 6 baselines
nb.md(
    '## 6 · Baselines on the same folds',
    'Instruction 70: every strategy is reported as **improvement over the best baseline**, not over '
    'shuffled labels.',
    '* **class prior** — predict the single most common class of the training folds '
    '(`apicoplast_mitochondrion`) for every gene and nothing else. Its pooled precision is that '
    'class’s prevalence; its recall is low because it makes one call per gene.\n'
    '* **prevalence-weighted random** — for each class, call each test gene with the training-fold '
    'prevalence of that class as a coin weight, seed 0. Its precision is the prevalence, and that is '
    'the floor.\n'
    '* **best single domain column per class** — for each class and fold, the one domain column with '
    'the highest average precision **on the training folds** is chosen and its value is the '
    'prediction. This is the strongest simple annotation-only strategy, and the one the model has to '
    'beat to be worth anything: it is what "look up whether this protein has a kinase domain" '
    'achieves.\n'
    '* **ortholog copy** — a gene takes the classes of other labelled genes in its orthogroup, '
    'abstaining when there are none. Under orthogroup-grouped folds every labelled ortholog is in the '
    'same fold, so coverage is ~0 by construction; reported rather than hidden.')
code("""
    baseline_scores = {}

    # 1 class prior
    pri = pd.DataFrame(0.0, index=GENES, columns=CLASSES)
    for f in FOLDS:
        tr = [g for g in GENES if fold[g] != f]
        te = [g for g in GENES if fold[g] == f]
        best = Y.loc[tr].sum().idxmax()
        pri.loc[te, best] = 1.0
    baseline_scores['baseline: class prior (most common class)'] = pri

    # 2 prevalence-weighted random
    rb = np.random.default_rng(0)
    rnd = pd.DataFrame(0.0, index=GENES, columns=CLASSES)
    for f in FOLDS:
        tr = [g for g in GENES if fold[g] != f]
        te = [g for g in GENES if fold[g] == f]
        prev = Y.loc[tr].mean()
        for c in CLASSES:
            rnd.loc[te, c] = (rb.random(len(te)) < prev[c]).astype(float)
    baseline_scores['baseline: prevalence-weighted random'] = rnd

    # 3 best single domain column per class, chosen on the training folds
    domcols = [c for c in DOMAIN_FEATURES if c in BASE.columns]
    dmat = BASE.loc[GENES, domcols].astype(float)
    bestdom = pd.DataFrame(0.0, index=GENES, columns=CLASSES)
    chosen = []
    for f in FOLDS:
        tr = [g for g in GENES if fold[g] != f]
        te = [g for g in GENES if fold[g] == f]
        for c in CLASSES:
            ytr = Y.loc[tr, c].values
            if ytr.sum() == 0:
                continue
            aps = {}
            for col in domcols:
                v = dmat.loc[tr, col].fillna(0).values
                if v.std() == 0:
                    continue
                aps[col] = average_precision_score(ytr, v)
            if not aps:
                continue
            col = max(aps, key=aps.get)
            chosen.append(dict(fold=f, cls=c, column=col, train_ap=round(aps[col], 3)))
            bestdom.loc[te, c] = dmat.loc[te, col].fillna(0).values
    baseline_scores['baseline: best single domain column per class'] = bestdom
    chosen_df = pd.DataFrame(chosen)

    # 4 ortholog copy
    ortho = pd.DataFrame(0.0, index=GENES, columns=CLASSES)
    cov = []
    for g in GENES:
        peers = [p for p in GENES if p != g and groups[p] == groups[g] and fold[p] != fold[g]]
        cov.append(bool(peers))
        if peers:
            ortho.loc[g] = (Y.loc[peers].mean() >= 0.5).astype(float).values
    ORTHO_COVERAGE = float(np.mean(cov))
    baseline_scores['baseline: ortholog copy (grouped, abstains)'] = ortho

    print('ortholog-copy coverage under orthogroup-grouped folds: %.3f' % ORTHO_COVERAGE)
    print()
    print('the domain column chosen per class (most frequent choice over the five folds):')
    print(chosen_df.groupby(['cls', 'column']).size().rename('folds').reset_index()
          .sort_values(['cls', 'folds'], ascending=[True, False]).groupby('cls').head(1).to_string(index=False))
""")

# --------------------------------------------------------------------------- 7 run
nb.md(
    '## 7 · The model, in two training variants',
    '* **`all-labels`** — one binary model per class fitted on all 242 labelled genes of the training '
    'folds. More data; the screen-only genes help place the boundary for the four phenotype classes.\n'
    '* **`non-screen`** — fitted on the 145 experimental-GO and curated-Swiss-Prot genes only, and '
    'scored only on them. This is the variant that cannot have learned anything from the screen’s '
    'choice of which genes to look at, and §11 reads it against the other.',
    'Both are out-of-fold throughout: no test gene, and no paralog of one, is ever in its own '
    'training set.')
code("""
    import time
    t0 = time.time()
    print('fitting the all-labels variant (15 classes x 5 folds)')
    OOF, MODELS, USED = fit_oof(FEATURES, GENES, verbose=True)
    print('   %.0f s' % (time.time() - t0))
    t0 = time.time()
    print('fitting the non-screen variant (145 genes)')
    OOF_NS, MODELS_NS, USED_NS = fit_oof(FEATURES, NONSCREEN_GENES, verbose=True)
    print('   %.0f s' % (time.time() - t0))

    PF = pairs_frame(OOF)
    PF_NS = pairs_frame(OOF_NS)
    print()
    print('pooled held-out pairs: %d (%d positive)' % (len(PF), int(PF.truth.sum())))
    print('pooled ROC AUC %.3f | average precision %.3f (prevalence %.3f)'
          % (roc_auc_score(PF.truth, PF.score), average_precision_score(PF.truth, PF.score),
             PF.truth.mean()))
""")
code("""
    board = []
    def board_rows(strategy, scores, genes=None, variant='all-labels', coverage=1.0):
        pf = pairs_frame(scores, genes=genes)
        rows = []
        per_fold = {}
        for f in FOLDS:
            sub = pf[pf.fold == f]
            if not len(sub):
                continue
            per_fold[f] = pooled_metrics(sub)
        agg = pooled_metrics(pf)
        for metric in ('precision', 'recall', 'precision_lower90', 'precision_ic_weighted',
                       'roc_auc', 'average_precision', 'calls'):
            vals = [per_fold[f].get(metric, np.nan) for f in per_fold]
            row = dict(strategy=strategy, variant=variant, metric=metric,
                       pooled=agg.get(metric, np.nan),
                       fold_mean=float(np.nanmean(vals)) if len(vals) else np.nan,
                       fold_min=float(np.nanmin(vals)) if len(vals) else np.nan,
                       fold_max=float(np.nanmax(vals)) if len(vals) else np.nan)
            row.update({f'fold_{f}': per_fold[f].get(metric, np.nan) for f in per_fold})
            rows.append(row)
        rows.append(dict(strategy=strategy, variant=variant, metric='coverage', pooled=coverage,
                         fold_mean=coverage, fold_min=coverage, fold_max=coverage))
        return rows

    for name, sc in baseline_scores.items():
        board += board_rows(name, sc,
                            coverage=ORTHO_COVERAGE if 'ortholog' in name else 1.0)
    board += board_rows('GBT one-vs-rest (all labels)', OOF)
    board += board_rows('GBT one-vs-rest (non-screen genes only)', OOF_NS,
                        genes=NONSCREEN_GENES, variant='non-screen')
    scoreboard = pd.DataFrame(board)
    scoreboard.to_csv(OUT / 'scoreboard.tsv', sep='\\t', index=False)

    piv = scoreboard.pivot_table(index=['variant', 'strategy'], columns='metric', values='pooled')
    cols = ['precision', 'precision_lower90', 'precision_ic_weighted', 'recall', 'roc_auc',
            'average_precision', 'calls', 'coverage']
    print('POOLED SCOREBOARD at the 0.5 cut (precision_lower90 = Clopper-Pearson lower 90 %):')
    print(piv[[c for c in cols if c in piv]].round(3).to_string())
""")
code("""
    print('per-fold spread, and whether the model beats the best baseline in at least 4 of 5 folds')
    print('(instruction 70 deployment rule):')
    print()
    fcols = [f'fold_{f}' for f in FOLDS]
    for metric in ('average_precision', 'precision'):
        acc = scoreboard[(scoreboard.metric == metric) & (scoreboard.variant == 'all-labels')]
        acc = acc.set_index('strategy')
        base = acc.reindex(columns=fcols).loc[[s for s in acc.index if s.startswith('baseline')]]
        for strat in [s for s in acc.index if s.startswith('GBT')]:
            mod = acc.reindex(columns=fcols).loc[strat].astype(float)
            beats = int((mod.values > base.max().values).sum())
            print('%-20s %-42s folds beaten: %d/5' % (metric, strat, beats))
            print('%-20s   model per fold: %s' % ('', np.round(mod.values, 3)))
            print('%-20s   best baseline:  %s' % ('', np.round(base.max().values.astype(float), 3)))
    BEST_BASELINE = (scoreboard[(scoreboard.metric == 'average_precision') &
                                (scoreboard.variant == 'all-labels') &
                                scoreboard.strategy.str.startswith('baseline')]
                     .sort_values('pooled', ascending=False).iloc[0])
    MODEL_AP = float(scoreboard[(scoreboard.metric == 'average_precision') &
                                (scoreboard.strategy == 'GBT one-vs-rest (all labels)')].pooled.iloc[0])
    print()
    print('best baseline: %s, pooled average precision %.3f; model %.3f (improvement %+.3f)'
          % (BEST_BASELINE.strategy, BEST_BASELINE.pooled, MODEL_AP, MODEL_AP - BEST_BASELINE.pooled))
""")

# --------------------------------------------------------------------------- 8 calibration + thresholds
nb.md(
    '## 8 · Calibration, then the pooled threshold table — the headline',
    'A boosted tree emits a score that is not a probability. Isotonic regression maps it to the '
    'observed rate of being right (`parascape.claims._isotonic`), fitted **across the pooled pairs** '
    'and **cross-fitted**: for each fold the curve comes from the pooled out-of-fold pairs of the '
    '*other* folds and is applied to this one, so no point helps calibrate itself.',
    'Then, at each calibrated threshold: how many (gene, class) calls are made, how many are right, '
    'the precision, the **Clopper–Pearson lower 90 % bound**, the IC-weighted precision with its '
    'effective-n bound, and how many currently-unlabelled cells would be filled.')
code("""
    def cross_fitted_calibration(score_df, train_genes):
        cal = pd.DataFrame(np.nan, index=score_df.index, columns=CLASSES, dtype=float)
        for f in FOLDS:
            tr = [g for g in train_genes if fold[g] != f]
            te = [g for g in score_df.index if fold[g] == f]
            if not tr or not te:
                continue
            s = score_df.loc[tr].values.ravel()
            t = Y.loc[tr, CLASSES].values.ravel().astype(float)
            ok = ~np.isnan(s)
            xs, ys = _isotonic(s[ok], t[ok])
            cal.loc[te] = np.interp(score_df.loc[te].values, xs, ys)
        return cal

    CAL = cross_fitted_calibration(OOF, GENES)
    CAL_NS = cross_fitted_calibration(OOF_NS.loc[NONSCREEN_GENES], NONSCREEN_GENES)
    PFC = pairs_frame(CAL)
    PFC_NS = pairs_frame(CAL_NS)

    rows = []
    edges = np.array([0, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
    for name, pf, raw in [('all-labels', PFC, PF), ('non-screen', PFC_NS, PF_NS)]:
        s, t = pf.score.values, pf.truth.values.astype(float)
        ece = _ece(s, t, bins=10)
        ece_raw = _ece(raw.score.values, raw.truth.values.astype(float), bins=10)
        idx = np.clip(np.searchsorted(edges, s, side='right') - 1, 0, len(edges) - 2)
        for b in range(len(edges) - 1):
            m = idx == b
            if not m.any():
                continue
            rows.append(dict(variant=name, ece=round(ece, 4), ece_uncalibrated=round(ece_raw, 4),
                             predicted_range=f'[{edges[b]:.2f}, {edges[b + 1]:.2f})', n_pairs=int(m.sum()),
                             mean_predicted=round(float(s[m].mean()), 4),
                             observed_positive_rate=round(float(t[m].mean()), 4),
                             positives=int(t[m].sum())))
    calibration = pd.DataFrame(rows)
    calibration.to_csv(OUT / 'calibration.tsv', sep='\\t', index=False)
    print('expected calibration error over the pooled pairs:')
    print(calibration.groupby('variant')[['ece', 'ece_uncalibrated']].first().to_string())
    print()
    print(calibration[calibration.variant == 'all-labels'].to_string(index=False))
""")
code("""
    # genome-wide models, for counting what would be fillable. Fitted on all 242 labels, scored on
    # all 8,140 genes, mapped through the isotonic curve fitted in-sample on the pooled out-of-fold
    # pairs. Not independently validated - it exists to count cells.
    Xall = XF['all']
    cols_all = usable(Xall, FEATURES, GENES)
    FULL = pd.DataFrame(np.nan, index=BASE.index, columns=CLASSES, dtype=float)
    FINAL_MODELS = {}
    for c in CLASSES:
        m = make_model().fit(Xall.loc[GENES, cols_all], Y.loc[GENES, c])
        FINAL_MODELS[c] = m
        FULL[c] = m.predict_proba(Xall[cols_all])[:, 1]
    s = OOF.values.ravel(); t = Y[CLASSES].values.ravel().astype(float)
    xs, ys = _isotonic(s[~np.isnan(s)], t[~np.isnan(s)])
    FULL_CAL = pd.DataFrame(np.interp(FULL.values, xs, ys), index=FULL.index, columns=CLASSES)
    UNLABELLED = [g for g in BASE.index if g not in set(GENES)]
    print('genome-wide calibrated matrix: %d genes x %d classes; unlabelled genes %d (%d cells)'
          % (FULL_CAL.shape[0], FULL_CAL.shape[1], len(UNLABELLED), len(UNLABELLED) * len(CLASSES)))

    def threshold_table(pf, variant, label='pooled'):
        rows = []
        for th in THRESHOLDS:
            d = pooled_metrics(pf, th)
            sub = FULL_CAL.loc[UNLABELLED] >= th
            d.update(question=label, variant=variant, threshold=th,
                     fillable_cells=int(sub.values.sum()),
                     fillable_genes=int((sub.sum(axis=1) > 0).sum()),
                     reaches_bar=bool(d['calls'] and d['precision_lower90'] >= 0.90))
            rows.append(d)
        return pd.DataFrame(rows)

    thr_pooled = pd.concat([threshold_table(PFC, 'all-labels'),
                            threshold_table(PFC_NS, 'non-screen', 'pooled, non-screen genes only')],
                           ignore_index=True)
    order = ['question', 'variant', 'threshold', 'pairs', 'positives', 'calls', 'correct', 'precision',
             'precision_lower90', 'precision_ic_weighted', 'precision_ic_lower90_approx',
             'ic_effective_n', 'recall', 'fillable_cells', 'fillable_genes', 'reaches_bar']
    thr_pooled = thr_pooled[[c for c in order if c in thr_pooled]]
    thr_pooled.to_csv(OUT / 'thresholds_pooled.tsv', sep='\\t', index=False)
    print()
    print(thr_pooled[thr_pooled.variant == 'all-labels'].round(4).to_string(index=False))
""")
code("""
    def verdict(tab):
        ok = tab[tab.reaches_bar]
        if len(ok):
            r = ok.sort_values('threshold').iloc[0]
            best = tab.loc[tab.precision_lower90.idxmax()]
            return dict(reaches_bar=True, threshold=float(r.threshold), calls=int(r.calls),
                        correct=int(r.correct), precision=float(r.precision),
                        lower90=float(r.precision_lower90),
                        ic_weighted=float(r.precision_ic_weighted),
                        fillable_cells=int(r.fillable_cells), fillable_genes=int(r.fillable_genes),
                        best_lower90=float(best.precision_lower90))
        b = tab.loc[tab.precision_lower90.idxmax()]
        return dict(reaches_bar=False, threshold=float(b.threshold), calls=int(b.calls),
                    correct=int(b.correct), precision=float(b.precision),
                    lower90=float(b.precision_lower90), ic_weighted=float(b.precision_ic_weighted),
                    fillable_cells=int(b.fillable_cells), fillable_genes=int(b.fillable_genes),
                    best_lower90=float(b.precision_lower90))

    V = verdict(thr_pooled[thr_pooled.variant == 'all-labels'])
    V_NS = verdict(thr_pooled[thr_pooled.variant == 'non-screen'])
    if V['reaches_bar']:
        POOLED_ANSWER = ('YES for the POOLED multi-label question. The lowest calibrated threshold whose '
                         'Clopper-Pearson lower 90 %% bound of precision reaches 0.90 is %.2f: %d held-out '
                         '(gene, class) calls, %d correct, precision %.3f, lower bound %.4f, '
                         'IC-weighted precision %.3f. At that threshold %d currently-empty cells over %d '
                         'genes would be filled.'
                         % (V['threshold'], V['calls'], V['correct'], V['precision'], V['lower90'],
                            V['ic_weighted'], V['fillable_cells'], V['fillable_genes']))
    else:
        POOLED_ANSWER = ('NO for the POOLED multi-label question. No calibrated threshold from 0.20 to 0.99 '
                         'reaches a lower 90 %% confidence bound of 0.90. The best bound anywhere on the '
                         'curve is %.4f (threshold %.2f: %d calls, %d correct, point-estimate precision '
                         '%.3f, IC-weighted %.3f). Nothing may be written into the information map as a '
                         'fact for the function trait.'
                         % (V['best_lower90'], V['threshold'], V['calls'], V['correct'],
                            V['precision'], V['ic_weighted']))
    print(POOLED_ANSWER)
    print()
    print('non-screen-only variant:', 'reaches the bar' if V_NS['reaches_bar'] else 'does not reach the bar',
          '(best bound %.4f)' % V_NS['best_lower90'])
    pd.DataFrame({'all-labels': V, 'non-screen': V_NS}).T
""")
nb.md(
    '### Not a calibration artefact: precision at the top of the pooled ranking',
    'A "no" in the threshold table could in principle mean only that the isotonic curve refuses ever '
    'to say 0.90, rather than that the ranking is bad. This removes the doubt: the pooled pairs are '
    'sorted by raw out-of-fold score and the precision of the top *k* is computed directly with its '
    'bound, for every *k* from 10 up. **Any** threshold rule is a choice of *k*, so the best bound '
    'over all *k* is the ceiling for the whole approach.')
code("""
    rows = []
    for name, pf in [('all-labels', PF), ('non-screen', PF_NS)]:
        srt = pf.sort_values('score', ascending=False)
        tv, icv = srt.truth.values, srt.ic.values
        for k in range(10, len(srt) + 1, 1 if len(srt) < 2000 else 5):
            correct = int(tv[:k].sum())
            wp, wlo, _ = weighted_precision_and_bound(tv[:k], icv[:k])
            rows.append(dict(variant=name, top_k=k, correct=correct, precision=correct / k,
                             precision_lower90=clopper_pearson_lower(correct, k),
                             precision_ic_weighted=wp,
                             reaches_bar=bool(clopper_pearson_lower(correct, k) >= 0.90)))
    at_k = pd.DataFrame(rows)
    at_k.to_csv(OUT / 'precision_at_k.tsv', sep='\\t', index=False)
    best_k = at_k.loc[at_k.groupby('variant').precision_lower90.idxmax()]
    ATK_NOTE = '; '.join('%s: top %d, %d correct, precision %.3f, lower bound %.4f'
                         % (r.variant, r.top_k, r.correct, r.precision, r.precision_lower90)
                         for r in best_k.itertuples())
    print('best bound anywhere on the pooled ranking -', ATK_NOTE)
    print('does ANY k clear the bar?', bool(at_k.reaches_bar.any()))
    best_k.set_index('variant').round(4)
""")

# --------------------------------------------------------------------------- 9 per class
nb.md(
    '## 9 · Per class, with every call count — so the under-powering is visible',
    'The pooled headline could hide a model that is excellent on two big classes and useless on '
    'thirteen. This table is the pooled threshold table repeated per class at the 0.5 cut and at the '
    'pooled operating threshold, with each class’s **call count** beside its bound, so a reader can '
    'see which classes are under-powered rather than finding them averaged away.',
    'The truth-set README predicted what this table shows: with 65 held-out calls needed even at a '
    'true precision of 0.95, and the largest non-screen class holding 37 genes, **no class can clear '
    'the bar on its own**. The table is here to demonstrate that rather than assert it.')
code("""
    OP = V['threshold'] if V['reaches_bar'] else 0.5
    rows = []
    for c in CLASSES:
        for th_name, th in [('0.50', 0.5), ('pooled operating (%.2f)' % OP, OP)]:
            pf = PFC[PFC.cls == c]
            d = pooled_metrics(pf, th)
            d.update(cls=c, threshold_name=th_name, threshold=th, class_genes=int(Y[c].sum()),
                     ic_corpus=float(IC[c]),
                     screen_supported=c in SCREEN_CLASSES,
                     reaches_bar=bool(d['calls'] and d['precision_lower90'] >= 0.90),
                     calls_needed_for_bar=22)
            rows.append(d)
    per_class = pd.DataFrame(rows)
    order = ['cls', 'class_genes', 'ic_corpus', 'screen_supported', 'threshold_name', 'threshold',
             'calls', 'correct', 'precision', 'precision_lower90', 'recall', 'roc_auc',
             'average_precision', 'reaches_bar']
    per_class = per_class[[c for c in order if c in per_class]]
    per_class.to_csv(OUT / 'per_class.tsv', sep='\\t', index=False)
    print('PER CLASS at the 0.5 cut (sorted by class size):')
    show = per_class[per_class.threshold_name == '0.50'].sort_values('class_genes', ascending=False)
    print(show.round(3).to_string(index=False))
    print()
    print('classes whose held-out call count could even in principle support the bound (>= 22 calls):',
          int((show.calls >= 22).sum()), 'of', len(CLASSES))
    print('classes that reach the bar:', int(show.reaches_bar.sum()))
""")

# --------------------------------------------------------------------------- 10 no knowledge
nb.md(
    '## 10 · The CAFA no-knowledge split — the number that matters for filling blanks',
    '**The definition, stated explicitly.** A gene is *no-knowledge* when it carries **no prior '
    'functional annotation of any kind inside the feature set**, which here means both of:',
    '1. **no domain annotation** — `domain_novel` is True: no InterPro entry and no Pfam entry in the '
    'ToxoDB gene-attribute download, so all 313 domain columns are False for it; and\n'
    '2. **no informative name** — its ToxoDB `product` is "hypothetical protein" or an unspecified '
    'product, so nothing in the annotation names a function. (`product` is held out as a feature, but '
    'a gene whose name is "guanylate cyclase" is not a blank, and the split is about blanks.)',
    'A looser second definition — **domain-free only**, ignoring the product name — is reported '
    'beside it, because it keeps more genes and is the weaker claim.',
    '**Why this is the honest number.** Lesson 2: report a feature against the hard comparison. The '
    'genes the map actually needs filled are the ones nobody has annotated; a model that scores well '
    'overall by agreeing with existing domain annotations fills nothing new. If the bound holds on '
    'the no-knowledge genes it is a real gap-filling capability. If it only holds overall, the result '
    'is annotation agreement.',
    '**Expect this split to be small, and say so.** A gene reaches the function truth set because '
    'somebody wrote about it, and those genes mostly have domains and names.')
code("""
    raw_nodes = pd.read_parquet(NODES).set_index('gene_id')
    product = raw_nodes['product'].fillna('')
    UNINFORMATIVE = re.compile(r'hypothetical|unspecified product|unknown function|putative uncharacteri',
                               re.I)
    no_name = product.str.contains(UNINFORMATIVE)
    no_domain = dom['domain_novel'].reindex(raw_nodes.index).fillna(False).astype(bool)
    NK_STRICT = (no_domain & no_name).reindex(GENES).fillna(False)
    NK_LOOSE = no_domain.reindex(GENES).fillna(False)
    print('no-knowledge definitions over the 242 labelled genes:')
    print('  strict (no InterPro/Pfam domain AND uninformative product name): %d genes'
          % int(NK_STRICT.sum()))
    print('  loose  (no InterPro/Pfam domain, any name):                      %d genes'
          % int(NK_LOOSE.sum()))
    print('  genome-wide: strict %d of %d genes, loose %d'
          % (int((no_domain & no_name).sum()), len(raw_nodes), int(no_domain.sum())))
    print()
    print('tier of the no-knowledge labelled genes (strict):')
    print(truth.loc[NK_STRICT.index[NK_STRICT], 'tier'].value_counts().to_string())
    print()
    rows = []
    for name, mask in [('no-knowledge strict (no domain, no name)', NK_STRICT),
                       ('no-knowledge loose (no domain)', NK_LOOSE),
                       ('has prior annotation (the complement of strict)', ~NK_STRICT),
                       ('all labelled genes', pd.Series(True, index=GENES))]:
        genes = list(mask.index[mask])
        if not genes:
            continue
        pf = PFC[PFC.gene_id.isin(genes)]
        for th in (0.5, 0.7, 0.9):
            d = pooled_metrics(pf, th)
            d.update(split=name, genes=len(genes), threshold=th,
                     reaches_bar=bool(d['calls'] and d['precision_lower90'] >= 0.90))
            rows.append(d)
        d = pooled_metrics(pf, 0.5)
        best = best_bound(pf)
        rows.append(dict(split=name, genes=len(genes), threshold=np.nan,
                         pairs=len(pf), positives=int(pf.truth.sum()),
                         calls=best[2], correct=best[3], precision=np.nan,
                         precision_lower90=best[0], roc_auc=d.get('roc_auc'),
                         average_precision=d.get('average_precision'),
                         reaches_bar=bool(best[0] >= 0.90) if best[0] == best[0] else False))
    no_knowledge = pd.DataFrame(rows)
    cols = ['split', 'genes', 'threshold', 'pairs', 'positives', 'calls', 'correct', 'precision',
            'precision_lower90', 'precision_ic_weighted', 'recall', 'roc_auc', 'average_precision',
            'reaches_bar']
    no_knowledge = no_knowledge[[c for c in cols if c in no_knowledge]]
    no_knowledge.to_csv(OUT / 'no_knowledge.tsv', sep='\\t', index=False)
    print('the rows with threshold = NaN hold the BEST bound over all thresholds for that split:')
    print(no_knowledge.round(3).to_string(index=False))
""")
code("""
    nk_pf = PFC[PFC.gene_id.isin(NK_STRICT.index[NK_STRICT])]
    kn_pf = PFC[PFC.gene_id.isin(NK_STRICT.index[~NK_STRICT])]
    NK_NOTE = ('no-knowledge (strict) %d genes, %d pooled pairs, %d positives; pooled AUC %.3f '
               'against %.3f on the %d annotated genes. Calls at 0.5: %d, correct %d.'
               % (int(NK_STRICT.sum()), len(nk_pf), int(nk_pf.truth.sum()),
                  roc_auc_score(nk_pf.truth, nk_pf.score) if nk_pf.truth.nunique() == 2 else np.nan,
                  roc_auc_score(kn_pf.truth, kn_pf.score), int((~NK_STRICT).sum()),
                  int((nk_pf.score >= 0.5).sum()), int(nk_pf.truth[nk_pf.score >= 0.5].sum())))
    print(NK_NOTE)
    print()
    print('the strict no-knowledge genes, their classes and their best-scoring class:')
    show = []
    for g in NK_STRICT.index[NK_STRICT]:
        s = CAL.loc[g].sort_values(ascending=False)
        show.append(dict(gene_id=g, tier=truth.loc[g, 'tier'],
                         truth_classes='|'.join(sorted(Y.columns[Y.loc[g]])),
                         top_prediction=s.index[0], top_prob=round(float(s.iloc[0]), 3),
                         top_is_correct=bool(Y.loc[g, s.index[0]])))
    print(pd.DataFrame(show).to_string(index=False))
""")

# --------------------------------------------------------------------------- 11 tiers
nb.md(
    '## 11 · The tier split — and the function analogue of the easy contrast',
    'Lesson 2 again. The 145 genes at the two non-screen tiers (32 experimental GO, 113 curated '
    'Swiss-Prot) and the 97 screen-only genes are **different questions**, for two reasons that have '
    'nothing to do with each other:',
    '* the screen genes carry only the **four phenotype classes**, which are the four largest and '
    'vaguest in the set (`apicoplast_mitochondrion` IC 1.60, `cell_cycle_division` 1.96, '
    '`cytoskeleton` 2.06, `egress` 2.53), so a model that calls them well scores a lot of pooled '
    'pairs cheaply;\n'
    '* the screen chose **which 319 genes to look at**, so the screen-only genes are a non-random '
    'slice of the genome even with the four label columns dropped.',
    '**If the model does well only on the screen-derived classes, that is the function analogue of '
    'the host-facing easy contrast** — a flattering overall number carried by the easy half. The '
    'table below reports: pooled metrics on the non-screen genes and on the screen-only genes, for '
    'both training variants, and pooled metrics split by **class group** (the four screen-supported '
    'classes against the other eleven).')
code("""
    rows = []
    SPLITS = [('non-screen genes (experimental GO + curated Swiss-Prot)', NONSCREEN_GENES, None),
              ('screen-only genes', SCREEN_GENES, None),
              ('all genes, the 4 screen-supported classes', GENES, SCREEN_CLASSES),
              ('all genes, the other 11 classes', GENES, [c for c in CLASSES if c not in SCREEN_CLASSES]),
              ('non-screen genes, the other 11 classes', NONSCREEN_GENES,
               [c for c in CLASSES if c not in SCREEN_CLASSES])]
    for variant, cal in [('all-labels', CAL), ('non-screen', CAL_NS)]:
        for name, genes, classes in SPLITS:
            genes = [g for g in genes if g in cal.index]
            if not genes:
                continue
            pf = pairs_frame(cal, genes=genes, classes=classes)
            if pf.truth.nunique() < 2:
                continue
            best = best_bound(pf)
            d = pooled_metrics(pf, 0.5)
            d.update(variant=variant, split=name, genes=len(genes),
                     classes=len(classes) if classes else len(CLASSES),
                     best_lower90=best[0], best_at_threshold=best[1], best_calls=best[2],
                     best_correct=best[3], reaches_bar=bool(best[0] >= 0.90))
            rows.append(d)
    by_tier = pd.DataFrame(rows)
    cols = ['variant', 'split', 'genes', 'classes', 'pairs', 'positives', 'calls', 'correct',
            'precision', 'precision_lower90', 'precision_ic_weighted', 'recall', 'roc_auc',
            'average_precision', 'best_lower90', 'best_at_threshold', 'best_calls', 'best_correct',
            'reaches_bar']
    by_tier = by_tier[[c for c in cols if c in by_tier]]
    by_tier.to_csv(OUT / 'by_tier.tsv', sep='\\t', index=False)
    print(by_tier.round(3).to_string(index=False))
""")
code("""
    def grab(variant, split, col):
        m = (by_tier.variant == variant) & (by_tier.split == split)
        return float(by_tier[m][col].iloc[0]) if m.any() else np.nan

    TIER_NOTE = ('screen-only genes: pooled AUC %.3f, best bound %.4f over %d calls. Non-screen genes: '
                 'pooled AUC %.3f, best bound %.4f. The four screen-supported classes: AUC %.3f; the '
                 'other eleven: AUC %.3f. Non-screen genes restricted to the other eleven classes '
                 '(the hardest honest cut): AUC %.3f, best bound %.4f.'
                 % (grab('all-labels', 'screen-only genes', 'roc_auc'),
                    grab('all-labels', 'screen-only genes', 'best_lower90'),
                    grab('all-labels', 'screen-only genes', 'best_calls'),
                    grab('all-labels', 'non-screen genes (experimental GO + curated Swiss-Prot)', 'roc_auc'),
                    grab('all-labels', 'non-screen genes (experimental GO + curated Swiss-Prot)', 'best_lower90'),
                    grab('all-labels', 'all genes, the 4 screen-supported classes', 'roc_auc'),
                    grab('all-labels', 'all genes, the other 11 classes', 'roc_auc'),
                    grab('all-labels', 'non-screen genes, the other 11 classes', 'roc_auc'),
                    grab('all-labels', 'non-screen genes, the other 11 classes', 'best_lower90')))
    print(TIER_NOTE)
""")

# --------------------------------------------------------------------------- 12 importance
nb.md(
    '## 12 · Permutation importance by family, on held-out folds only',
    'A feature (or a whole family) is shuffled **in the held-out fold**, every class model is '
    're-scored, the pooled pairs are rebuilt and the fall in pooled average precision is recorded. '
    'Families are permuted as a block with one shared row permutation, so the correlation inside a '
    'family survives and the number means "this piece of evidence, removed".',
    'Families first, because that is the unit of evidence: 36 expression columns from one experiment '
    'are one piece of evidence, not 36. Then individual features, for the top-10 list.',
    'Read the family column rather than the feature names. With ~48 genes held out per fold the '
    'individual numbers are small and noisy.')
code("""
    import time
    def pooled_ap_from(scores):
        s = scores[CLASSES].values.ravel()
        t = Y[CLASSES].values.ravel().astype(float)
        return float(average_precision_score(t, s))

    BASE_AP = pooled_ap_from(OOF)
    print('baseline pooled average precision: %.4f' % BASE_AP)

    def permuted_ap(cols_to_permute, seed):
        rng = np.random.default_rng(seed)
        out = OOF.copy()
        for f in FOLDS:
            te = [g for g in GENES if fold[g] == f]
            cols = USED[f]
            touch = [c for c in cols_to_permute if c in cols]
            if not touch or not te:
                continue
            Xte = XF[f].loc[te, cols].copy()
            perm = rng.permutation(len(te))
            for c in touch:
                Xte[c] = Xte[c].values[perm]
            for c in CLASSES:
                m = MODELS.get((c, f))
                if m is not None:
                    out.loc[te, c] = m.predict_proba(Xte)[:, 1]
        return pooled_ap_from(out)

    t0 = time.time()
    fam_rows = []
    for family, group in fam.groupby(fam):
        cols = list(group.index)
        vals = [BASE_AP - permuted_ap(cols, 10_000 + 7 * r) for r in range(5)]
        fam_rows.append(dict(family=family, n_features=len(cols),
                             importance_mean=float(np.mean(vals)), importance_sd=float(np.std(vals))))
    family_importance = pd.DataFrame(fam_rows).sort_values('importance_mean', ascending=False)
    print('family permutation importance (%.0f s):' % (time.time() - t0))
    print(family_importance.head(25).round(4).to_string(index=False))
""")
code("""
    t0 = time.time()
    rows = []
    for feat in FEATURES:
        if not any(feat in USED[f] for f in FOLDS):
            continue
        vals = [BASE_AP - permuted_ap([feat], 50_000 + 13 * r) for r in range(3)]
        rows.append(dict(feature=feat, family=fam[feat], importance_mean=float(np.mean(vals)),
                         importance_sd=float(np.std(vals))))
    importance = pd.DataFrame(rows).sort_values('importance_mean', ascending=False)
    importance['scope'] = 'feature'
    fi = family_importance.rename(columns={'family': 'feature'}).assign(scope='family')
    fi['family'] = fi.feature
    pd.concat([importance, fi[['feature', 'family', 'importance_mean', 'importance_sd', 'scope',
                               'n_features']]], ignore_index=True).to_csv(
        OUT / 'importance.tsv', sep='\\t', index=False)
    TOP10 = importance.head(10)[['feature', 'family', 'importance_mean', 'importance_sd']]
    print('per-feature permutation importance over %d features (%.0f s). Top 20:'
          % (len(importance), time.time() - t0))
    print(importance.head(20).round(4).to_string(index=False))
    TOP10.round(4)
""")
nb.md(
    '### Robustness: does the verdict survive removing a family?',
    'Host-facing run 1’s sharpest lesson was that its answer moved when one family was removed. So '
    'the pooled question is refitted with each worrying block taken out and the best bound reported '
    'each time. **A verdict that depends on one family is not a verdict.**',
    'The set that matters most is **minus every annotation-derived family** — the 313 domain columns '
    'and `has_ec_number`. Those are the ones that could be reading the same annotation act the class '
    'rules read, and §10 is the same question asked by gene instead of by feature.')
code("""
    ANNOTATION_DERIVED = [c for c in FEATURES if fam[c] in ('domain', 'interpro', 'toxodb_ec_numbers')]
    ATTENTION = [c for c in FEATURES if fam[c] in ('pubmed', 'bioid_corpus_membership',
                                                   'ipms_corpus_membership', 'ipms_baits', 'pmc_oa')]
    STRESS_SCREENS = [c for c in FEATURES if fam[c] in ('oxidative_stress_screen', 'differentiation_screen',
                                                        'brady_subtypes', 'sexual_stages',
                                                        'gse132248_stress')]
    SETS = {
        'all features': FEATURES,
        'minus coexpression (guilt by association)': [c for c in FEATURES if fam[c] != 'coexpression'],
        'minus the 313 domain columns': [c for c in FEATURES if fam[c] != 'domain'],
        'minus ALL annotation-derived (domains + EC number)': [c for c in FEATURES if c not in set(ANNOTATION_DERIVED)],
        'minus attention columns': [c for c in FEATURES if c not in set(ATTENTION)],
        'minus stress/differentiation screens': [c for c in FEATURES if c not in set(STRESS_SCREENS)],
        'minus sequence/topology traits': [c for c in FEATURES if fam[c] != 'sequence_traits'],
        'ONLY domains': [c for c in FEATURES if fam[c] == 'domain'],
        'ONLY coexpression': [c for c in FEATURES if fam[c] == 'coexpression'],
        'ONLY expression + fitness + proteome (no annotation, no coexpression)':
            [c for c in FEATURES if fam[c] not in ('domain', 'interpro', 'toxodb_ec_numbers',
                                                   'coexpression', 'sequence_traits')],
    }
    rows = []
    for name, cols in SETS.items():
        if not cols:
            continue
        o, _m, _u = fit_oof(cols, GENES)
        c = cross_fitted_calibration(o, GENES)
        pf = pairs_frame(c)
        best = best_bound(pf)
        praw = pairs_frame(o)
        rows.append(dict(feature_set=name, n_features=len(cols),
                         roc_auc=roc_auc_score(praw.truth, praw.score),
                         average_precision=average_precision_score(praw.truth, praw.score),
                         best_lower90=best[0], at_threshold=best[1], calls=best[2], correct=best[3],
                         reaches_bar=bool(best[0] >= 0.90)))
    robust = pd.DataFrame(rows)
    robust.to_csv(OUT / 'robustness.tsv', sep='\\t', index=False)
    ROBUST_NOTE = ('%d of %d feature sets reach the bar; best bound ranges %.3f to %.3f; pooled AUC ranges '
                   '%.3f to %.3f' % (int(robust.reaches_bar.sum()), len(robust), robust.best_lower90.min(),
                                     robust.best_lower90.max(), robust.roc_auc.min(), robust.roc_auc.max()))
    print(ROBUST_NOTE)
    robust.round(4)
""")

# --------------------------------------------------------------------------- 13 predictions
nb.md(
    '## 13 · Predictions for every gene and every class',
    'A calibrated probability for all 8,140 genes × 15 classes, from models fitted on all 242 labels '
    'and mapped through the isotonic curve fitted on the pooled out-of-fold pairs. For labelled genes '
    'the **out-of-fold** probability is given too — that is the honest one, and it is what every '
    'score in this notebook uses. The all-gene columns exist to count what *would* be fillable.',
    '**Nothing here is written into the information map.** Whether any of it may be is decided by the '
    'bound in §8.')
code("""
    pred = pd.DataFrame(index=BASE.index)
    pred.index.name = 'gene_id'
    for c in CLASSES:
        pred['prob_' + c] = FULL_CAL[c].round(4)
    for c in CLASSES:
        pred['oof_' + c] = CAL[c].reindex(BASE.index).round(4)
    for c in CLASSES:
        pred['label_' + c] = Y[c].reindex(BASE.index)
    pred['tier'] = truth.tier.reindex(BASE.index)
    pred['fold'] = fold.reindex(BASE.index)
    pred['labelled'] = pred.index.isin(GENES)
    pred['no_knowledge_strict'] = (no_domain & no_name).reindex(BASE.index).fillna(False)
    pred['no_knowledge_loose'] = no_domain.reindex(BASE.index).fillna(False)
    pred['top_class'] = FULL_CAL.idxmax(axis=1)
    pred['top_prob'] = FULL_CAL.max(axis=1).round(4)
    pred['n_classes_above_0_5'] = (FULL_CAL >= 0.5).sum(axis=1)
    pred.sort_values('top_prob', ascending=False).to_csv(OUT / 'predictions.tsv', sep='\\t')
    print('predictions written for %d genes x %d classes' % (len(pred), len(CLASSES)))
    print()
    print('unlabelled cells at or above each threshold (cells / genes with at least one):')
    print(pd.DataFrame({t: [int((FULL_CAL.loc[UNLABELLED] >= t).values.sum()),
                            int(((FULL_CAL.loc[UNLABELLED] >= t).sum(axis=1) > 0).sum())]
                        for t in THRESHOLDS}, index=['cells', 'genes']).to_string())
    print()
    print('highest-scoring unlabelled genes (exploratory, NOT facts):')
    print(pred.loc[UNLABELLED].sort_values('top_prob', ascending=False)
          .head(15)[['top_class', 'top_prob', 'n_classes_above_0_5', 'no_knowledge_strict']].to_string())
""")

# --------------------------------------------------------------------------- 14 limits
nb.md(
    '## 14 · Limitations, none softened',
    '1. **Negative pairs are "not recorded", not "not true".** Function annotation is open-world: a '
    'gene that really does X but whose curator text never says so scores a correct prediction as a '
    'false positive. Precision is therefore biased **downward**, which makes a passing bound '
    'trustworthy and a failing bound possibly pessimistic. CAFA has the same property and no fix '
    'exists short of experiments.\n'
    '2. **The classes are coarse and of very unequal specificity** (IC 1.60 to 4.92). The pooled '
    'number is dominated by the big vague classes, which is exactly why the IC-weighted figure is '
    'reported beside it at every threshold.\n'
    '3. **The labels come from text rules.** `class_rules.md` matches regular expressions in curator '
    'text and GO descendant sets. A protein whose FUNCTION comment names a domain gets the class the '
    'domain implies, and the domain columns see the same annotation — not label leakage, but '
    'correlation between two annotation acts. §10 and §12 are the two ways this notebook measures it; '
    'neither removes it.\n'
    '4. **The screen tier is a different experiment with a different label space** (four process '
    'classes) and a non-random gene list. Its four source columns are dropped throughout, but §11 is '
    'needed because the gene selection cannot be dropped.\n'
    '5. **242 genes is small, and 145 at the non-screen tiers is smaller.** Per class it is too small '
    'to score at all (§9), which is what the truth-set README said in advance.\n'
    '6. **One-vs-rest ignores the structure between classes.** Nothing stops the model calling '
    '`translation` and `apicoplast_mitochondrion` with no regard for how often they co-occur; a '
    'joint model (instruction 72 stage 4) is not in this notebook.\n'
    '7. **The co-expression features are fold-safe but not study-held-out.** Partner labels come from '
    'training folds only, yet the co-expression matrix itself is built from the same expression '
    'studies as the expression features, so the two families are not independent evidence.\n'
    '8. **The all-gene probability columns are not independently validated** — the isotonic curve is '
    'fitted in-sample on the pooled out-of-fold pairs. Only the cross-fitted out-of-fold numbers are '
    'honest.\n'
    '9. **The ortholog-copy baseline has coverage near 0** under orthogroup-grouped folds, by '
    'construction, and is reported rather than hidden.\n'
    '10. **One model, fixed settings, no hyperparameter search.** Stages 2–5 of the instruction 72 '
    'plan (frozen sequence/structure blocks, label propagation, the joint row model, pairs) are not '
    'here.\n'
    '11. **Permutation importance is noisy at this size** (~48 genes held out per fold); signs flip '
    'between folds. Read families, not features.')
nb.md(
    '## 15 · Labels or evidence? — the arithmetic, not a preference',
    'Two different things can block the bar and they need different fixes:',
    '* **If the true precision at the operating point is already above 0.90**, the only obstacle is '
    'sample size. The Clopper–Pearson bound converges upward to the true precision, so enough calls '
    'will demonstrate it, and the table below says how many — which is directly actionable, because '
    'Einar is reviewing a queue of 254 genes.\n'
    '* **If the true precision is below 0.90**, the bound converges to that precision and **no number '
    'of labels will ever clear the bar.** More labels then only measure the failure more precisely, '
    'and the fix has to be better evidence or a better model. That is what the second host-facing run '
    'found for host-facing.',
    'The truth-set README argued in advance that for function the constraint is **labels**, not '
    'evidence. The cells below check that claim three ways instead of assuming it: the measured '
    'precision against the bar, the **learning curve** (does the model still improve as labels are '
    'added?), and the call count on the hardest honest cut.')
nb.md(
    '### The learning curve — the only honest way to say whether more labels would help',
    '"The precision is 0.68, so no number of labels can reach 0.90" is true **only if the model does '
    'not get better with more training data**, and that is an assumption, not a measurement. So the '
    'whole pipeline is refitted on 40 %, 60 %, 80 % and 100 % of the training labels — dropped by '
    '**orthogroup**, so the subsample respects the same grouping as the folds — and scored on the '
    'same held-out genes every time. Two independent draws per fraction, seed 1.',
    'If pooled average precision and precision are still climbing at 100 %, labels are the binding '
    'constraint and Einar’s review queue is the fix. If the curve has flattened, the evidence is the '
    'ceiling and more curation will not move it.')
code("""
    import time
    t0 = time.time()
    group_list = sorted(set(groups))
    rng_lc = np.random.default_rng(1)
    lc_rows = []
    for frac in (0.4, 0.6, 0.8, 1.0):
        reps = 1 if frac == 1.0 else 2
        for rep in range(reps):
            if frac == 1.0:
                sub = list(GENES)
            else:
                keep = set(rng_lc.choice(group_list, size=int(round(len(group_list) * frac)),
                                         replace=False))
                sub = [g for g in GENES if groups[g] in keep]
            o, _m, _u = fit_oof(FEATURES, sub, score_genes=GENES)
            c = cross_fitted_calibration(o, sub)
            pf_raw, pf_cal = pairs_frame(o), pairs_frame(c)
            best = best_bound(pf_cal)
            srt = pf_raw.sort_values('score', ascending=False)
            tv = srt.truth.values
            prec_at_k = {k: tv[:k].sum() / k for k in (10, 20, 25, 50)}
            lc_rows.append(dict(fraction=frac, rep=rep, train_genes=len(sub),
                                roc_auc=roc_auc_score(pf_raw.truth, pf_raw.score),
                                average_precision=average_precision_score(pf_raw.truth, pf_raw.score),
                                precision_top25=prec_at_k[25], precision_top50=prec_at_k[50],
                                best_lower90=best[0], best_calls=best[2], best_correct=best[3]))
    learning_curve = pd.DataFrame(lc_rows)
    learning_curve.to_csv(OUT / 'learning_curve.tsv', sep='\\t', index=False)
    lc = learning_curve.groupby('fraction')[['train_genes', 'roc_auc', 'average_precision',
                                             'precision_top25', 'precision_top50',
                                             'best_lower90']].mean()
    print('learning curve (%.0f s), mean over draws, trained on a fraction of the label set,'
          % (time.time() - t0))
    print('scored on the same held-out genes throughout:')
    print(lc.round(4).to_string())
    slope_ap = float(lc.average_precision.iloc[-1] - lc.average_precision.iloc[-2])
    slope_p25 = float(lc.precision_top25.iloc[-1] - lc.precision_top25.iloc[-2])
    LC_NOTE = ('from 80 %% to 100 %% of the labels pooled average precision moves %+.4f and precision of '
               'the top 25 calls moves %+.4f; over the whole curve average precision goes %.3f -> %.3f'
               % (slope_ap, slope_p25, lc.average_precision.iloc[0], lc.average_precision.iloc[-1]))
    print()
    print(LC_NOTE)
    print('still climbing at the top of the curve:', bool(slope_ap > 0.01 or slope_p25 > 0.02))
""")
code("""
    rows = []
    for tp in (1.00, 0.99, 0.98, 0.97, 0.95, 0.93, 0.92, 0.91):
        n = 1
        while n < 50000 and clopper_pearson_lower(int(round(n * tp)), n) < 0.90:
            n += 1
        rows.append(dict(true_precision=tp, calls_needed=n if n < 50000 else np.nan))
    needed = pd.DataFrame(rows)
    print('calls needed for the lower 90 % bound to reach 0.90, by TRUE precision:')
    print(needed.to_string(index=False))
    print('(at a true precision of 0.90 or below the bound never reaches 0.90, at any number of calls)')
    print()
    best_thr = thr_pooled[thr_pooled.variant == 'all-labels'].loc[
        thr_pooled[thr_pooled.variant == 'all-labels'].precision_lower90.idxmax()]
    atk_best = at_k[at_k.variant == 'all-labels'].sort_values('precision_lower90').iloc[-1]
    print('pooled, best calibrated threshold: %.2f, %d calls, %d correct, precision %.3f, bound %.4f'
          % (best_thr.threshold, best_thr.calls, best_thr.correct, best_thr.precision,
             best_thr.precision_lower90))
    print('pooled, best point on the raw ranking: top %d, %d correct, precision %.3f, bound %.4f'
          % (atk_best.top_k, atk_best.correct, atk_best.precision, atk_best.precision_lower90))
    obs_p = float(max(best_thr.precision, at_k[at_k.variant == 'all-labels'].precision.max()))
    hardest = by_tier[(by_tier.variant == 'all-labels') &
                      (by_tier.split == 'non-screen genes, the other 11 classes')]
    h_calls = int(hardest.best_calls.iloc[0]); h_correct = int(hardest.best_correct.iloc[0])
    h_bound = float(hardest.best_lower90.iloc[0])
    call_rate = float(best_thr.calls) / len(PFC)
    still_climbing = bool(slope_ap > 0.01 or slope_p25 > 0.02)
    if obs_p < 0.90:
        WHAT_IS_NEEDED = (
            'WHICH CASE IS THIS? BOTH, at different places on the curve, and the two have to be '
            'separated.\\n\\n'
            '(a) AT THE POOLED OPERATING POINT IT IS THE EVIDENCE. The best precision this model reaches '
            'anywhere - at any calibrated threshold and at any point on the pooled ranking - is %.3f, '
            'which is BELOW the 0.90 bar itself. The lower bound converges upward to the TRUE precision, '
            'so if this model is what we have, no number of labels brings the bound to 0.90: more labels '
            'would measure %.3f more precisely.\\n\\n'
            '(b) BUT THE MODEL IS NOT FIXED, AND THE LEARNING CURVE SAYS LABELS STILL BUY ACCURACY: %s. '
            'So (a) is a statement about the model trained on 242 genes, not a ceiling on the problem. '
            'This is where the truth-set README is right: it argued the constraint for function is the '
            'label count, and %s.\\n\\n'
            '(c) ON THE HARDEST HONEST CUT THE BOUND IS CAPPED BY ARITHMETIC, NOT BY ERRORS. Restricted to '
            'the 145 non-screen genes and the 11 non-screen-supported classes, the best threshold makes '
            '%d calls and %d are correct - precision %.3f - and the lower 90 %% bound is %.3f ONLY because '
            '%d calls cannot support more. That cut is label-limited in the plainest possible way: the '
            'bound needs about 22 calls before it can reach 0.90 at all, and 65 at a true precision of '
            '0.95.\\n\\n'
            'ACTIONABLE NUMBER. At the current call rate (%.1f %% of held-out pairs called at the best '
            'threshold) a 65-call demonstration needs about %d labelled genes against the %d there are '
            'now - roughly %d more. Einar\\'s 254-gene review queue is larger than that, and its top 100 '
            'alone would take the labelled set from %d to about %d genes and the non-screen set from %d to '
            '%d. That is worth doing. What it will NOT do on its own is lift a model whose precision is '
            '%.3f to 0.90: the queue buys the measurement and some accuracy, and the rest has to come '
            'from evidence that distinguishes the classes better than expression, fitness and domain '
            'flags do.'
            % (obs_p, obs_p, LC_NOTE,
               'the curve has not flattened' if still_climbing else 'the curve HAS flattened, which '
               'weakens that argument',
               h_calls, h_correct, (h_correct / h_calls) if h_calls else float('nan'), h_bound, h_calls,
               100 * call_rate,
               int(round(65 / max(call_rate, 1e-9) / len(CLASSES))), len(GENES),
               max(0, int(round(65 / max(call_rate, 1e-9) / len(CLASSES))) - len(GENES)),
               len(GENES), len(GENES) + 100, len(NONSCREEN_GENES), len(NONSCREEN_GENES) + 100,
               obs_p))
    else:
        need = int(needed.set_index('true_precision').loc[0.95, 'calls_needed'])
        WHAT_IS_NEEDED = (
            'WHICH CASE IS THIS? The LABEL case, as the truth-set README argued. The observed precision is '
            '%.3f, above the bar, so the block is sample size: about %d held-out calls are needed to '
            'demonstrate it at a true precision of 0.95. The model calls %.1f %% of held-out (gene, class) '
            'pairs at its best threshold, so that needs roughly %d labelled pairs, i.e. about %d labelled '
            'genes against the %d there are now - roughly %d more reviewed genes. Einar\\'s 254-gene review '
            'queue covers that, and the top 100 alone would take the labelled set from %d to %d genes.'
            % (obs_p, need, 100 * call_rate, int(round(need / max(call_rate, 1e-9))),
               int(round(need / max(call_rate, 1e-9) / len(CLASSES))), len(GENES),
               max(0, int(round(need / max(call_rate, 1e-9) / len(CLASSES))) - len(GENES)),
               len(GENES), len(GENES) + 100))
    print()
    print(WHAT_IS_NEEDED)
""")
code("""
    # per-class: how many more labelled genes would each class need before it could be scored alone?
    rows = []
    for c in CLASSES:
        pf = PFC[PFC.cls == c]
        calls = int((pf.score >= 0.5).sum())
        rate = calls / len(pf) if len(pf) else 0
        rows.append(dict(cls=c, genes_now=int(Y[c].sum()), held_out_calls_now=calls,
                         call_rate=round(rate, 3),
                         genes_for_22_calls=int(np.ceil(22 / rate)) if rate > 0 else np.nan,
                         genes_for_65_calls=int(np.ceil(65 / rate)) if rate > 0 else np.nan))
    per_class_need = pd.DataFrame(rows).sort_values('genes_now', ascending=False)
    print('what each class would need before its OWN bound could be demonstrated')
    print('(22 calls = the bound can reach 0.90 only with perfect precision; 65 calls = at a true')
    print(' precision of 0.95, which is the truth-set README number):')
    print(per_class_need.to_string(index=False))
    per_class_need.to_csv(OUT / 'per_class_label_need.tsv', sep='\\t', index=False)
""")

# --------------------------------------------------------------------------- 16 README
code("""
    def fmt(df, index=True):
        return '```\\n' + df.to_string(index=index) + '\\n```\\n'

    piv = scoreboard.pivot_table(index=['variant', 'strategy'], columns='metric', values='pooled')
    scols = [c for c in ['precision', 'precision_lower90', 'precision_ic_weighted', 'recall',
                         'roc_auc', 'average_precision', 'calls', 'coverage'] if c in piv]
    R = []
    R.append('# Function inference on the 242-gene truth set (2026-11-04)\\n')
    R.append('Built by `scripts/build_function_inference_notebook.py`; the executed notebook is\\n'
             '`notebooks/function_inference_2026_11_04.ipynb`. Instruction 72, target trait 2 of 3.\\n'
             'Labels: `results/function_truth_2026_11_04/`.\\n')
    R.append('## The answer\\n')
    R.append('**' + POOLED_ANSWER + '**\\n')
    R.append('\\nNon-screen-tier variant: %s (best bound %.4f).\\n'
             % ('reaches the bar' if V_NS['reaches_bar'] else 'does not reach the bar', V_NS['best_lower90']))
    R.append('\\nNot a calibration artefact: ranking the pooled pairs by raw out-of-fold score and taking\\n'
             'the top k directly, for every k from 10 up - which is what any threshold rule amounts to -\\n'
             'the bar is cleared by no k. ' + ATK_NOTE + '. (`precision_at_k.tsv`)\\n')
    R.append('\\nRobustness across feature sets: ' + ROBUST_NOTE + '. (`robustness.tsv`)\\n')
    R.append('\\n## Labels or evidence?\\n')
    R.append(WHAT_IS_NEEDED + '\\n')
    R.append('\\n### Learning curve (mean over draws; training labels dropped by orthogroup)\\n')
    R.append(fmt(lc.round(4)))
    R.append('## Why the headline is pooled\\n')
    R.append('The truth-set README does the arithmetic: a Clopper-Pearson lower 90 %% bound needs ~65\\n'
             'held-out calls even at a true precision of 0.95, and the largest non-screen class has 37\\n'
             'genes. So no class can be scored alone, and the headline is the pooled multi-label\\n'
             'question over all %d (gene, class) pairs (%d positive). Per-class numbers are in\\n'
             '`per_class.tsv` with every call count, so the under-powering is visible rather than\\n'
             'averaged away: %d of %d classes have >= 22 held-out calls at the 0.5 cut, which is the\\n'
             'minimum at which the bound could reach 0.90 even with perfect precision.\\n'
             % (len(PFC), int(PFC.truth.sum()),
                int((per_class[per_class.threshold_name == '0.50'].calls >= 22).sum()), len(CLASSES)))
    R.append('## Unweighted or information-content weighted?\\n')
    R.append('Both, at every threshold. The **unweighted** Clopper-Pearson bound is the number the fill\\n'
             'decision turns on: the fill rule is a statement about cells, and only the unweighted bound\\n'
             'is exact for counting cells. The **IC-weighted** precision is the headline for what the\\n'
             'filled cells are worth - it says whether a model passes by calling specific things or by\\n'
             'calling `apicoplast_mitochondrion` (IC 1.60) on everything. The weighted bound uses a Kish\\n'
             'effective sample size and is an approximation, labelled as such in every table.\\n')
    R.append('## Pooled scoreboard (pooled over held-out pairs, 0.5 cut)\\n')
    R.append(fmt(piv[scols].round(3)))
    R.append('Best baseline: **%s** (pooled average precision %.3f). Model %.3f, improvement %+.3f.\\n'
             % (BEST_BASELINE.strategy, BEST_BASELINE.pooled, MODEL_AP, MODEL_AP - BEST_BASELINE.pooled))
    R.append('\\nPer-fold values for every strategy and metric are in `scoreboard.tsv`.\\n')
    R.append('## The pooled threshold table (the headline)\\n')
    R.append(fmt(thr_pooled[thr_pooled.variant == 'all-labels'][
        ['threshold', 'calls', 'correct', 'precision', 'precision_lower90', 'precision_ic_weighted',
         'recall', 'fillable_cells', 'fillable_genes', 'reaches_bar']].round(4), index=False))
    R.append('## Per class\\n')
    R.append(fmt(per_class[per_class.threshold_name == '0.50'].sort_values('class_genes', ascending=False)[
        ['cls', 'class_genes', 'ic_corpus', 'screen_supported', 'calls', 'correct', 'precision',
         'precision_lower90', 'recall', 'roc_auc', 'average_precision', 'reaches_bar']].round(3), index=False))
    R.append('## CAFA no-knowledge split\\n')
    R.append('**Definition.** A gene is *no-knowledge* when it has no prior functional annotation of any\\n'
             'kind in the feature set: (1) `domain_novel` is True - no InterPro and no Pfam entry, so all\\n'
             '313 domain columns are False - and (2) its ToxoDB `product` is "hypothetical protein" or an\\n'
             'unspecified product. A looser definition (domain-free only, any name) is reported beside it.\\n'
             'Strict: %d of the 242 labelled genes. Loose: %d. Genome-wide the strict definition covers\\n'
             '%d of %d genes, which is the population the map most needs filled.\\n'
             % (int(NK_STRICT.sum()), int(NK_LOOSE.sum()), int((no_domain & no_name).sum()), len(raw_nodes)))
    R.append('\\n' + NK_NOTE + '\\n')
    R.append(fmt(no_knowledge.round(3), index=False))
    R.append('## Tier split, and the function analogue of the easy contrast\\n')
    R.append(TIER_NOTE + '\\n')
    R.append(fmt(by_tier[['variant', 'split', 'genes', 'classes', 'pairs', 'positives', 'calls',
                          'correct', 'precision', 'precision_lower90', 'roc_auc', 'average_precision',
                          'best_lower90', 'reaches_bar']].round(3), index=False))
    R.append('## Method\\n')
    R.append('* **Labels**: %d genes x %d classes = %d (gene, class) pairs, %d positive. Unlabelled genes\\n'
             '  are never negatives. A negative pair means "not recorded", not "not true".\\n'
             '* **Features**: %d, in %d families - numeric/boolean `nodes.parquet` columns, the DeepTMHMM\\n'
             '  sequence/topology traits, the 313 domain columns, and %d co-expression neighbour-label\\n'
             '  features. Every feature with its family is in `features.tsv` and printed in the notebook.\\n'
             '* **Held out**: the whole `splitcas9_imaging_screen` family (the 6 `screen_*` columns),\\n'
             '  dropped THROUGHOUT rather than conditionally, because those four phenotype columns are the\\n'
             '  source of the 97 screen-tier labels; every `lopit|compartment|ortholopit|tagm` column\\n'
             '  (localization is a tier-1 conclusion and restates some classes); `hosttx_effectors`,\\n'
             '  `pvm_proximity_positive`, `n_host_targets`; and every free-text or identity column\\n'
             '  (`sequence`, `product`, `interpro_id`, `interpro_desc`, `pfam_id`, `ec_number`,\\n'
             '  `orthogroup`, ...). %d columns held out in total.\\n'
             '* **Folds**: 5, grouped by OrthoMCL orthogroup so paralogs never split, seed 0, assigned once\\n'
             '  over the 242 genes and shared by all 15 class models.\\n'
             '* **Model**: one `HistGradientBoostingClassifier` per class per fold (one-vs-rest), missing\\n'
             '  values passed through with no imputation, 300 rounds at learning rate 0.05, 15 leaves,\\n'
             '  early stopping off, no hyperparameter search. Two variants: `all-labels` (242 genes) and\\n'
             '  `non-screen` (145 genes, scored on 145).\\n'
             '* **Co-expression features are built inside the fold loop**: partner labels come only from\\n'
             '  the training folds, so no test gene sees a same-fold partner label.\\n'
             '* **Calibration**: isotonic (`parascape.claims._isotonic`), cross-fitted over the pooled\\n'
             '  pairs. ECE in `calibration.tsv`.\\n'
             '* **Bound**: exact one-sided Clopper-Pearson, `scipy.stats.beta.ppf(0.10, k, n-k+1)`.\\n'
             % (len(GENES), len(CLASSES), len(PFC), int(PFC.truth.sum()), len(FEATURES),
                fam.nunique(), len(CX_FEATURES), len(HELD) + len(text_cols)))
    R.append('## Calibration\\n')
    R.append(fmt(calibration.groupby('variant')[['ece', 'ece_uncalibrated']].first().round(4)))
    R.append('## Top 10 features by permutation importance (held-out folds, pooled average precision)\\n')
    R.append(fmt(TOP10.round(4), index=False))
    R.append('\\nBy family, which is the unit of evidence:\\n')
    R.append(fmt(family_importance.round(4), index=False))
    R.append('## Files\\n')
    R.append('| File | Contents |\\n|---|---|\\n'
             '| `scoreboard.tsv` | strategy x variant x metric, pooled and per fold |\\n'
             '| `thresholds_pooled.tsv` | the headline: threshold, calls, correct, precision, lower 90 % bound, IC-weighted precision, cells fillable, reaches 0.90 |\\n'
             '| `per_class.tsv` | the same per class with its call count, at 0.5 and at the pooled operating threshold |\\n'
             '| `per_class_label_need.tsv` | labelled genes each class would need before its own bound could be demonstrated |\\n'
             '| `no_knowledge.tsv` | the CAFA no-knowledge split, both definitions, against the annotated complement |\\n'
             '| `by_tier.tsv` | non-screen vs screen-only genes, and screen-supported vs other classes |\\n'
             '| `calibration.tsv` | reliability table and ECE over the pooled pairs |\\n'
             '| `importance.tsv` | permutation importance per feature and per family |\\n'
             '| `features.tsv` | every feature with its family |\\n'
             '| `robustness.tsv` | the bound with each feature family removed, and with each one alone |\\n'
             '| `learning_curve.tsv` | pooled metrics when trained on 40/60/80/100 % of the labels, dropped by orthogroup |\\n'
             '| `precision_at_k.tsv` | precision and bound of the top k pooled pairs, every k from 10 up |\\n'
             '| `predictions.tsv` | all 8,140 genes x 15 classes: calibrated probability, out-of-fold probability for labelled genes, label, tier, fold, no-knowledge flags |\\n')
    R.append('\\n## Limitations (not softened)\\n')
    R.append('1. **Negative pairs are "not recorded", not "not true".** Open-world annotation: a correct\\n'
             '   prediction of an unrecorded function scores as a false positive, so precision is biased\\n'
             '   downward. A passing bound is trustworthy; a failing bound may be pessimistic.\\n'
             '2. **The classes are coarse and of unequal specificity** (IC 1.60 to 4.92); the pooled number\\n'
             '   is dominated by the big vague classes, which is why the IC-weighted figure sits beside it.\\n'
             '3. **The labels come from regular expressions over curator text and GO descendant sets.** The\\n'
             '   domain columns read the same annotation act. Not label leakage, but not independent\\n'
             '   evidence either; the no-knowledge split and the annotation-free feature set measure it.\\n'
             '4. **The screen tier is a different experiment** with four process classes and a non-random\\n'
             '   319-gene list. Its four source columns are dropped throughout; the gene selection cannot be.\\n'
             '5. **242 genes, 145 at the non-screen tiers.** Per class this is too small to score at all.\\n'
             '6. **One-vs-rest ignores the structure between classes**; no joint model here.\\n'
             '7. **Co-expression features are fold-safe but not study-held-out** - the same expression\\n'
             '   studies feed the expression family.\\n'
             '8. **The all-gene probability columns are not independently validated** (isotonic curve fitted\\n'
             '   in-sample on pooled out-of-fold pairs). Only the out-of-fold numbers are honest.\\n'
             '9. **The ortholog-copy baseline has coverage %.2f** under orthogroup-grouped folds, by\\n'
             '   construction.\\n'
             '10. **One model, fixed settings, no search.** Stages 2-5 of the instruction 72 plan are absent.\\n'
             '11. **Permutation importance is noisy** at ~48 held-out genes per fold; read families.\\n'
             % ORTHO_COVERAGE)
    (OUT / 'README.md').write_text('\\n'.join(R))
    print((OUT / 'README.md').read_text()[:3000])
""")

nb.write(str(ROOT / 'notebooks/function_inference_2026_11_04.ipynb'))
print('wrote notebook and', ROOT / 'results/function_inference_2026_11_04')
