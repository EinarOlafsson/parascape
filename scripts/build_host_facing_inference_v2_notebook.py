"""Instruction 72, host-facing inference rerun (v2) on Einar's reviewed labels.

Run 1 (`results/host_facing_inference_2026_11_04/`) looked excellent and was not trustworthy: it
separated "inside the parasite" from "secreted" rather than "stays in the vacuole" from "exported
into the host", its top features were a host-transcription effector screen and a publication count,
and against the 9 curated vacuole-resident negatives its AUC collapsed to 0.855.

This rerun uses the labels after Einar reviewed all 140 published claims (43 vacuole-resident or
curated negatives instead of 9), holds out every declared dependency of `host_facing`, and makes the
HARD question the headline: among proteins that enter the secretory pathway, which are exported
beyond the parasite?

Writes notebooks/host_facing_inference_v2_2026_11_04.ipynb (executed) and the tables under
results/host_facing_inference_v2_2026_11_04/. Run 1's outputs are not touched.
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

nb = ExecutedNotebook('Host-facing inference v2: the hard question on Einar’s reviewed labels')


def code(src: str):
    return nb.code(dedent(src).strip('\n'))


nb.md(
    '# Host-facing inference v2 — the hard question, on reviewed labels',
    'Instruction 72 (Einar, 2026-11-04). This is a **rerun**, and the reason for the rerun is the '
    'honest part of the first one.',
    '**What run 1 found, and why its headline was not usable.** The first scoring run '
    '(`results/host_facing_inference_2026_11_04/`) reported ROC-AUC 0.995 and a lower 90 % precision '
    'bound of 0.906 — apparently clearing the fill bar. It was not a usable result:',
    '* the negatives were 477 hyperLOPIT markers of compartments *inside the parasite* (nucleus, '
    'ribosome, apicoplast, cytosol) against only **9** curated vacuole-resident negatives. So the '
    'model was asked "is this protein inside the parasite or secreted?", which a signal peptide '
    'alone nearly answers, and not "does this secreted protein stay in the vacuole or reach the host '
    'cell?", which is the biological question;',
    '* against those 9 curated vacuole residents its AUC fell to **0.855 with a 33 % false-positive '
    'rate**;',
    '* its two strongest features were a screen for parasite genes whose loss changes **host** '
    'transcription (nearly the conclusion itself) and a **publication count**.',
    '**What changed.** Einar reviewed all 140 published candidate claims over 98 genes. He rejected '
    '22 dense granule proteins that the PVM proximity screen had called host-facing (GRA1, GRA4, '
    'GRA22, GRA36, GRA37, GRA40, GRA43, GRA46, GRA48, GRA49, GRA51–GRA59, GRA12B, ROP34) while '
    'accepting 46 other proximity claims. The map now holds **43 vacuole-resident or curated '
    'negatives**. Those are the hard negatives the question needs: proteins that *are* secreted and '
    'that *do not* reach the host.',
    '**Three questions are scored, and the hard one is the headline.**',
    '| # | Question | Set |\n|---|---|---|\n'
    '| 1 | **HARD (headline)** — among proteins that enter the secretory pathway, which are exported '
    'beyond the parasite? | predicted signal peptide **or** a dense-granule / rhoptry / microneme '
    'hyperLOPIT marker |\n'
    '| 2 | **HARDEST** — can a human ruling on published evidence be predicted? | the 98 genes Einar '
    'reviewed |\n'
    '| 3 | *easy, for comparison only* — host-facing against everything labelled, internal markers '
    'included | all 585 labelled genes |',
    'The fill bar, unchanged: a predicted cell may be written as a **fact** only at a calibrated '
    'probability whose held-out precision has a **Clopper–Pearson lower 90 % confidence bound** of '
    '0.90 or better. Not the point estimate. "Nothing may be filled yet" is a result.',
    'A reader who wants only the answer should read §9 (hard), §10 (hardest), §11 (easy) and the '
    'limitations in §15.')

# --------------------------------------------------------------------------- 1 labels
nb.md(
    '## 1 · Labels, and where each one comes from',
    'Labels are the `host_facing` column of the information map v0 — **observed facts only**. *True* '
    'means the protein has been shown to reach beyond the parasite; *False* means it has been shown '
    'to stay inside; empty means untested. Empty is never a negative.',
    'Every label carries its **provenance**, and it travels with the gene through every table below, '
    'because the whole point of this rerun is that the provenances are not equally hard:',
    '* `einar_direct` (0.95) — a direct experimental claim Einar accepted (or rejected);\n'
    '* `einar_curated` (0.90) — a proximity or inferred claim he ruled on;\n'
    '* `gold_negative` (0.90–0.95) — the curated vacuole-resident negatives (e.g. GRA2 shaping the '
    'IVN);\n'
    '* `hyperlopit_internal` (0.80) — a marker of a compartment physically inside the parasite. An '
    'inference from residence, not an export test, and the reason run 1 looked good;\n'
    '* `host_interaction` / `go_hostcell` — a curated host interaction or an experimental GO '
    'host-cell term.',
    '**Conflicts.** `conflicts.tsv` lists 33 genes. All 33 have Einar’s review as one of their '
    'sources, so under his rule every one of them is **resolved, not open**: his ruling wins and the '
    'gene is kept. Run 1 dropped 8 genes here; this run drops none, and the check below shows that '
    'no conflicted gene lacks his ruling.')
code("""
    import os
    # One thread per BLAS/OpenMP pool. These tables are ~600 x 470; thread spin-up costs far more
    # than the arithmetic saves, and run 1 found permutation importance takes ~8 CPU-hours with 16
    # threads against minutes with 1. The machine is shared.
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[v] = '1'
    import json, re, warnings
    from pathlib import Path
    import numpy as np, pandas as pd
    from scipy.stats import beta
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    ROOT = Path('.').resolve()
    MAP = ROOT / 'results/information_map_v0_2026_11_04'
    TRUTH = ROOT / 'results/host_facing_truth_2026_11_04'
    SEQ = ROOT / 'results/sequence_traits_2026_11_04/sequence_traits.parquet'
    NODES = ROOT / 'parascape/data/nodes.parquet'
    OUT = ROOT / 'results/host_facing_inference_v2_2026_11_04'; OUT.mkdir(parents=True, exist_ok=True)

    values = pd.read_parquet(MAP / 'values.parquet')
    probability = pd.read_parquet(MAP / 'probability.parquet')
    evidence = pd.read_parquet(MAP / 'evidence.parquet')
    conflicts = pd.read_csv(MAP / 'conflicts.tsv', sep='\\t')
    truth = pd.read_csv(TRUTH / 'host_facing_truth.tsv', sep='\\t').set_index('gene_id')

    y = values['host_facing'].dropna().astype(bool)
    print('information map rows:', len(values))
    print('host_facing labelled:', len(y), '(True', int(y.sum()), '/ False', int((~y).sum()), ')')
    print('empty (untested) genes:', int(values['host_facing'].isna().sum()))
    print()
    conflicts['einar'] = conflicts.dataset.str.contains('Einar review', case=False, na=False)
    by_gene = conflicts.groupby('gene_id').einar.any()
    print('conflicts.tsv: %d rows over %d genes' % (len(conflicts), conflicts.gene_id.nunique()))
    print('  conflicted genes whose sources INCLUDE Einar\\'s review (resolved, kept):', int(by_gene.sum()))
    print('  conflicted genes without his review (would be dropped):', int((~by_gene).sum()),
          sorted(by_gene.index[~by_gene]))
    print('genes dropped by this run:', int((~by_gene).sum()))
""")
code("""
    hf_ev = evidence[evidence.trait == 'host_facing'].copy()

    def provenance_class(row):
        d = str(row.dataset)
        if 'Einar review' in d:
            return 'einar_direct' if row.source_tier == 'observed-direct' else 'einar_curated'
        if 'gold_negative' in d:
            return 'gold_negative'
        if 'hyperLOPIT marker' in d:
            return 'hyperlopit_internal'
        if 'pvm_proximity' in d:
            return 'pvm_proximity'
        if 'host_interactions' in d:
            return 'host_interaction'
        return 'go_hostcell'

    hf_ev['prov'] = hf_ev.apply(provenance_class, axis=1)
    hf_ev['as_bool'] = hf_ev.value.eq('yes')
    PRECEDENCE = ['einar_direct', 'einar_curated', 'gold_negative', 'host_interaction',
                  'go_hostcell', 'hyperlopit_internal', 'pvm_proximity']
    rank = {p: i for i, p in enumerate(PRECEDENCE)}
    by_gene_ev = {g: df for g, df in hf_ev.groupby('gene_id')}

    def deciding_provenance(gene, label):
        df = by_gene_ev.get(gene)
        agree = df[df.as_bool == label]
        if not len(agree):
            return 'unknown'
        return sorted(agree.prov, key=lambda p: rank[p])[0]

    prov = pd.Series({g: deciding_provenance(g, bool(v)) for g, v in y.items()}, name='provenance')
    prov.index.name = 'gene_id'
    labels = pd.DataFrame({'label': y, 'provenance': prov})
    labels['observed_probability'] = probability.loc[y.index, 'host_facing']
    labels['reviewed_by_einar'] = labels.index.isin(truth.index)
    print('every label with the source that decided it:')
    print(pd.crosstab(labels.provenance, labels.label).rename(
        columns={False: 'not host-facing', True: 'host-facing'}).to_string())
    print()
    HARDNEG = ['einar_direct', 'einar_curated', 'gold_negative']
    print('the hard negatives (the point of this rerun):',
          int(((~labels.label) & labels.provenance.isin(HARDNEG)).sum()))
    print()
    print('observed probability attached to each provenance:')
    print(labels.groupby(['label', 'provenance']).observed_probability.agg(['size', 'min', 'max']).to_string())
""")

nb.md(
    '### Two genes where the map and Einar disagree — said out loud',
    'Einar’s review accepted the PVM-proximity claim for **GRA9 (TGME49_251540)** and **NTPase II '
    '/ NTP1 (TGME49_277270)** as host-facing. The pre-existing curated gold-negative set says the '
    'opposite for both (immunogold EM places GRA9 in the dense granules and then on the IVN *inside* '
    'the vacuole; NTPase II is secreted in bulk into the PV lumen). Both sources carry probability '
    '0.90, and the map’s collapse resolved the tie to **not host-facing**.',
    'This matters enough to state rather than to average away:',
    '* questions 1 (hard) and 3 (easy) use the **map labels as they stand** (65 True / 520 False), '
    'which is the column the instruction names;',
    '* question 2 (hardest) uses **Einar’s own decisions** (66 yes / 32 no), which is what that '
    'question is about — so in that question these two genes are positives;',
    '* §12 reports the hard question again with the two genes flipped, so the reader can see whether '
    'the verdict depends on them.',
    'Either way this is a **label-quality finding, not a modelling one**: two of the 43 hard '
    'negatives are contested, and the contested direction is exactly the one that decides the '
    'question.')
code("""
    disagree = truth.join(y.rename('map_label'))
    disagree = disagree[(disagree.decision == 'yes') != (disagree.map_label == True)]
    print('reviewed genes whose map label differs from Einar\\'s decision:', len(disagree))
    print(disagree[['protein', 'decision', 'deciding_strength', 'deciding_pmid', 'map_label', 'basis']].to_string())
    EINAR_FLIPS = list(disagree.index)
""")

# --------------------------------------------------------------------------- 2 features
nb.md(
    '## 2 · Features, and everything held out as a declared dependency',
    'Tier 2 of instruction 72 means *measurements* — expression, fitness, proteome, PTMs, screens, '
    'sequence and structure. Those a prediction may learn from. Tier 1 are the *conclusions* and they '
    'are what is being predicted.',
    '**Held out, non-negotiable (instruction 72 `host_facing.depends_on`).** Run 1 is the reason this '
    'list is longer than last time:',
    '* every column of `hosttx_effectors` — a screen for parasite genes whose knockout changes host '
    'transcription. A gene that rewires the host is acting on the host; this is nearly the conclusion '
    'and it was run 1’s top feature;\n'
    '* every column of `organelle_surface_turboid` — proximity biotinylation on organelle surfaces, '
    'the same technique that defines the positives;\n'
    '* every column of `pubmed`, and the attention/bookkeeping columns `n_papers_*`, '
    '`n_publications`, `n_fulltext`, `n_bioid_studies`, `n_ipms_studies`, `lit_tier`, '
    '`attention_depth` — run 1 reached AUC 0.902 on these alone, which is a measure of fame;\n'
    '* `pvm_proximity_positive` and `n_host_targets` — these *define* the positives;\n'
    '* anything matching `lopit|compartment|ortholopit|tagm` — the negatives are hyperLOPIT markers, '
    'so a compartment call hands over the negative class.',
    'Kept: the numeric and boolean columns of `nodes.parquet` that survive that filter, the sequence '
    'and topology traits from DeepTMHMM (`n_tm`, `has_signal_peptide`, lengths, flank sides), the 3Di '
    'quality columns that instruction 72 requires to travel with any 3Di use (`tdi_coverage`, '
    '`tdi_align_identity`, `tdi_align_coverage`, `tdi_n_fragments`, `tdi_source`), and a new family '
    '**`foldseek_3di`**: the composition of the Foldseek 3Di structural alphabet, which is a crude '
    'shape summary and is kept as its own family so it can be held out on its own.')
code("""
    nodes = pd.read_parquet(NODES).set_index('gene_id')
    seqt = pd.read_parquet(SEQ)

    # object columns holding only booleans (the phenotype screens) -> float 1/0/NaN, so the tree can
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

    HELD_FAMILIES = ['hosttx_effectors', 'organelle_surface_turboid', 'pubmed', 'pvm_proximity']
    ATTENTION = ['n_papers_focal', 'n_papers_substantive', 'n_papers_incidental', 'n_publications',
                 'n_fulltext', 'n_bioid_studies', 'n_ipms_studies', 'lit_tier', 'attention_depth']
    EXPLICIT = ['pvm_proximity_positive', 'n_host_targets']
    NAME_PAT = re.compile('lopit|compartment|ortholopit|tagm', re.I)

    held_family = sorted(c for c in numeric if node_family.get(c) in HELD_FAMILIES)
    held_attention = sorted(c for c in numeric if c in ATTENTION)
    held_explicit = sorted(c for c in numeric if c in EXPLICIT)
    held_named = sorted(c for c in numeric if NAME_PAT.search(c))
    HELD = set(held_family) | set(held_attention) | set(held_explicit) | set(held_named)

    print('nodes.parquet: %d columns, %d numeric/boolean, %d free text' % (nodes.shape[1], len(numeric), len(text_cols)))
    for name, cols in [('declared dependency families (hosttx_effectors, organelle_surface_turboid, pubmed, pvm_proximity)', held_family),
                       ('attention / bookkeeping columns', held_attention),
                       ('columns that define the positives', held_explicit),
                       ('matching lopit|compartment|ortholopit|tagm', held_named)]:
        print('HELD OUT - %s (%d):' % (name, len(cols)))
        for i in range(0, len(cols), 4):
            print('      ' + '  '.join(f'{c:<34}' for c in cols[i:i + 4]).rstrip())
    print()
    print('attention columns named by the instruction that were free text and never features anyway:',
          [c for c in ATTENTION if c in text_cols])
    print('TOTAL held out:', len(HELD), 'of', len(numeric), 'numeric columns')
""")
nb.md(
    '### The 3Di composition family',
    'The Foldseek 3Di alphabet writes each residue as a letter describing its local tertiary '
    'neighbourhood, so the letter frequencies of a protein are a cheap summary of its shape — '
    'helix-rich, strand-rich, disordered. 20 frequency features over `tdi_full` plus the fraction of '
    'positions with no structural letter (`X`). Not a substitute for structure comparison; it is the '
    'coarsest usable form, kept as its own family `foldseek_3di` so §13 can hold it out.',
    'Instruction 72 is explicit that 3Di quality must travel with 3Di use: 520 of these strings were '
    'filled by aligning a near-identical structure and copying its letters, so `tdi_source`, '
    '`tdi_coverage`, `tdi_align_identity` and `tdi_align_coverage` are features too. An aligned 3Di '
    'is weaker evidence than a measured one, and 637 genes have none at all (missing, not zero).')
code("""
    TDI_ALPHABET = list('ACDEFGHIKLMNPQRSTVWY')
    tdi = seqt['tdi_full'].fillna('')
    comp = {}
    total = tdi.str.len().replace(0, np.nan)
    for letter in TDI_ALPHABET:
        comp['tdi3_' + letter] = tdi.str.count(letter) / total
    comp['tdi3_frac_X'] = tdi.str.count('X') / total
    tdi_comp = pd.DataFrame(comp, index=seqt.index)
    tdi_comp.loc[tdi.str.len() == 0, :] = np.nan   # no structure is missing, not a composition of zeros
    FOLDSEEK_3DI = list(tdi_comp.columns)

    SEQ_FEATURES = ['n_tm', 'has_signal_peptide', 'seq_len', 'mature_len', 'nterm_len', 'cterm_len',
                    'nterm_side', 'cterm_side', 'flanks_identical',
                    'tdi_coverage', 'tdi_n_fragments', 'tdi_align_identity', 'tdi_align_coverage',
                    'tdi_source']
    CATEGORICAL = ['nterm_side', 'cterm_side', 'tdi_source']

    node_features = [c for c in numeric if c not in HELD and c not in SEQ_FEATURES]
    X = nodes[node_features].astype('float64').join(seqt[SEQ_FEATURES]).join(tdi_comp)
    for c in ('has_signal_peptide', 'flanks_identical'):
        X[c] = X[c].astype('float64')
    for c in CATEGORICAL:
        X[c] = X[c].astype('category')

    # A column never measured on any labelled gene carries nothing here (and sklearn cannot bin it).
    all_missing = [c for c in X.columns if c not in CATEGORICAL and X.loc[y.index, c].notna().sum() == 0]
    X = X.drop(columns=all_missing)
    FEATURES = list(X.columns)
    fam = pd.Series({c: ('foldseek_3di' if c in FOLDSEEK_3DI
                         else 'sequence_traits' if c in SEQ_FEATURES
                         else node_family.get(c, 'unregistered')) for c in FEATURES}, name='family')
    fam.index.name = 'feature'
    print('dropped as never measured on a labelled gene:', all_missing)
    print('3Di composition features:', len(FOLDSEEK_3DI), '| genes with no 3Di string at all:',
          int((tdi.str.len() == 0).sum()))
    print('FINAL FEATURE COUNT:', len(FEATURES), '| categorical:', CATEGORICAL)
""")
code("""
    print('THE FINAL FEATURE LIST, each with its provenance family')
    print('(a family is held out as a unit, because 36 expression columns from one experiment are one')
    print(' piece of evidence, not 36):')
    print()
    for f, group in sorted(fam.groupby(fam), key=lambda kv: (-len(kv[1]), kv[0])):
        print('%s  (%d)' % (f, len(group)))
        names = sorted(group.index)
        for i in range(0, len(names), 4):
            print('      ' + '  '.join(f'{c:<36}' for c in names[i:i + 4]).rstrip())
    fam.value_counts().rename('features').to_frame().T
""")
nb.md(
    '### Does any family that is still in could encode the label? — asked, not assumed',
    'Run 1’s mistake was discovering this after the fact. Said in advance, for every family that '
    'is close to the conclusion and is **kept**:',
    '* **`pride_proximity` (`n_proximity_partners`)** — a different proximity-labelling corpus from '
    'the PVM screen. Same *technique* as the label source. Kept because it is not that experiment; '
    'checked in §13 and in the importance table.\n'
    '* **`bioid_corpus_membership`, `ipms_corpus_membership`, `ipms_baits`, `pmc_oa`** — whether the '
    'gene appears in a BioID or IP-MS corpus at all. These are **literature-presence columns in '
    'disguise**: a gene is a BioID bait because somebody found it interesting, and effectors are '
    'interesting. They are not in the declared hold-out list, so they stay, but they are '
    'attention-adjacent and §13 removes them as a set.\n'
    '* **`secretome_partition`, `toxodb_epitopes`, `iedb_bcell`, `crosslink_interactome`, '
    '`cyst_wall_interactome`, `spacr_escrt_screen`** — genuinely host-facing-adjacent biology '
    '(secreted-versus-vesicular partition, antibody-accessible epitopes, host ESCRT recruitment to '
    'the PVM). Evidence, not label definitions. A result resting entirely on one of them is a result '
    'about one screen, so the importance table decides.\n'
    '* **`gra12`, `gra17_synthlethal`** — fitness screens named after host-facing effectors. They '
    'measure fitness in a knockout background, not localisation; kept and watched.\n'
    '* **`foldseek_3di`, `sequence_traits`, `alphafold`, `foldseek_struct`, `has_structure`** — '
    'structure *availability* correlates with fame (a famous protein gets solved and deposited). The '
    '3Di quality columns are required by the instruction and are kept, with that caveat recorded.',
    '* **`has_signal_peptide`** is the one feature that is deliberately *not* removed, because it is '
    'also the strongest baseline, and the hard question is defined by conditioning on it.',
    'Two legitimate features are lost to the name rule: '
    '`fit_hyperlopit_unassigned_invivo_lib1/2` are CRISPR fitness scores, not localisation, but their '
    'names contain "lopit". Accepted — a name rule that errs toward removal errs in the safe '
    'direction.')
code("""
    prox = nodes['n_proximity_partners']
    pos_genes, neg_genes = y.index[y], y.index[~y]
    print('n_proximity_partners (pride_proximity): measured for %d of %d positives, %d of %d negatives'
          % (int(prox.loc[pos_genes].notna().sum()), len(pos_genes),
             int(prox.loc[neg_genes].notna().sum()), len(neg_genes)))
    print('  distinct values among labelled genes:', sorted(prox.loc[y.index].dropna().unique())[:10])
    print()
    print('what was removed: pvm_proximity_positive against the label (HELD OUT as a feature):')
    print(pd.crosstab(nodes['pvm_proximity_positive'].loc[y.index].fillna(-1), y).to_string())
    print()
    print('and the attention columns against the label, to show what fame looks like:')
    att = [c for c in ATTENTION if c in nodes.columns and pd.api.types.is_numeric_dtype(nodes[c])]
    print(nodes.loc[y.index, att].groupby(y).median().to_string())
""")

# --------------------------------------------------------------------------- 3 questions
nb.md(
    '## 3 · The three questions, defined by their gene sets',
    '**1 · HARD (the headline).** A gene is in the secretory-pathway set if **either** DeepTMHMM '
    'predicts a signal peptide **or** it is a hyperLOPIT marker of dense granules, rhoptries 1, '
    'rhoptries 2 or micronemes. The marker arm matters: DeepTMHMM misses signal peptides (GRA16, a '
    'secreted effector, is called GLOBULAR), and the 22 dense granule proteins Einar rejected are '
    'exactly the genes this question needs as negatives. Conditioning this way removes the easy '
    'contrast — the nuclear, ribosomal and apicoplast markers go — and leaves the question a '
    'biologist asks: **this protein is secreted; does it stay in the vacuole or reach the host?**',
    '**2 · HARDEST.** The 98 genes Einar ruled on, with **his** decisions as labels (66 yes, 32 no). '
    'Every label here is a human reading of published evidence; there is no marker inference '
    'anywhere in it, and the negatives are the vacuole residents.',
    '**3 · easy, for comparison only.** All 585 labelled genes, hyperLOPIT internal markers included. '
    'This is the number run 1 reported, and it is kept only so the inflation is visible.')
code("""
    SECRETORY_MARKERS = {'dense granules', 'rhoptries 1', 'rhoptries 2', 'micronemes'}
    sp_all = X['has_signal_peptide'].fillna(0).astype(bool)
    marker = values['compartment_hyperlopit']
    in_secretory = (sp_all.reindex(values.index).fillna(False) | marker.isin(SECRETORY_MARKERS))
    in_secretory.name = 'in_secretory'

    y_einar = truth.decision.eq('yes')

    QUESTIONS = {
        'hard': dict(genes=[g for g in y.index if in_secretory.get(g, False)], labels=y,
                     title='HARD (headline): among secretory-pathway proteins, which are exported?'),
        'reviewed': dict(genes=list(truth.index), labels=y_einar,
                         title='HARDEST: the 98 genes Einar reviewed, his rulings as labels'),
        'easy': dict(genes=list(y.index), labels=y,
                     title='easy (comparison only): all labels, internal markers included'),
    }
    for k, q in QUESTIONS.items():
        lab = q['labels'].loc[q['genes']]
        q['y'] = lab
        print('%-9s %4d genes | %3d host-facing / %3d not | prevalence %.3f'
              % (k, len(lab), int(lab.sum()), int((~lab).sum()), lab.mean()))
    print()
    print('secretory set genome-wide: %d of %d genes (%d labelled, %d untested)'
          % (int(in_secretory.sum()), len(in_secretory), int(in_secretory.loc[y.index].sum()),
             int(in_secretory.sum()) - int(in_secretory.loc[y.index].sum())))
    print('  by arm: signal peptide %d, marker-only %d'
          % (int(sp_all.reindex(values.index).fillna(False).sum()),
             int((~sp_all.reindex(values.index).fillna(False) & marker.isin(SECRETORY_MARKERS)).sum())))
    print('host-facing positives EXCLUDED from the hard set (no signal peptide, no secretory marker):',
          int((~in_secretory.loc[y.index[y]]).sum()), 'of', int(y.sum()))
    print()
    print('provenance of the hard set - this is the table that says the question got harder:')
    print(pd.crosstab(labels.loc[QUESTIONS['hard']['genes'], 'provenance'],
                      labels.loc[QUESTIONS['hard']['genes'], 'label']).to_string())
""")

# --------------------------------------------------------------------------- 4 folds
nb.md(
    '## 4 · Five folds grouped by orthogroup, seed 0 — the same construction as run 1',
    'Paralogs share an orthogroup and very often share a fate; the GRA and ROP families are full of '
    'near-copies. A paralog in the training fold and its twin in the test fold would let the model '
    'score by recognising the family. So whole orthogroups go to one fold (OrthoMCL groups from '
    '`MASTER_parasite_wide_by_orthogroup.csv`; a gene with no group is its own group). Seed 0.',
    'The fold assignment is made **once, over all labelled genes**, and the three questions inherit '
    'it, so the hard set and the easy set are not scored on differently-shuffled data.',
    'The cost shows up in the baselines: the *ortholog copy* baseline is left with nothing to copy '
    'from, coverage 0. That is correct. A baseline that only works because paralogs leak is not a '
    'baseline.')
code("""
    OG = Path('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/MASTER_parasite_wide_by_orthogroup.csv')
    og = pd.read_csv(OG, usecols=['orthomcl_group', 'tgon_gene_id']).dropna()
    group_of = {g: row.orthomcl_group for row in og.itertuples() for g in str(row.tgon_gene_id).split(';')}
    ALL_LABELLED = sorted(set(y.index) | set(truth.index))
    groups = pd.Series([group_of.get(g, 'self:' + g) for g in ALL_LABELLED], index=ALL_LABELLED, name='group')
    rng = np.random.default_rng(0)
    unique = np.array(sorted(set(groups))); rng.shuffle(unique)
    fold_all = groups.map({g: i % 5 for i, g in enumerate(unique)}).rename('fold')
    labels['fold'] = fold_all.reindex(labels.index)
    labels['group'] = groups.reindex(labels.index)
    labels['in_secretory'] = in_secretory.reindex(labels.index)
    print('labelled genes:', len(ALL_LABELLED), '| distinct orthogroups:', len(unique),
          '| genes sharing an orthogroup with another labelled gene:',
          int((groups.value_counts() > 1).reindex(groups).sum()))
    for k, q in QUESTIONS.items():
        q['fold'] = fold_all.loc[q['genes']]
        print()
        print(k, 'per-fold counts (negatives / positives):')
        print(pd.crosstab(q['fold'], q['y']).to_string())
""")

# --------------------------------------------------------------------------- 5 machinery
nb.md(
    '## 5 · How every strategy is scored, and the model',
    'The same numbers for everything, per fold, reported as the mean with the fold-to-fold spread, '
    'because with ~12 positives per fold a single mean is noise: **accuracy**, **precision** (the '
    'number the fill bar is about), **recall**, **ROC AUC** and **average precision**.',
    '`HistGradientBoostingClassifier`. Missing values pass straight through — the tree learns a '
    'direction for "not measured", which is the right treatment when absence from a proteomics run is '
    'itself informative. **No imputation anywhere.** Fixed settings, identical for every fold and '
    'every question (300 rounds, learning rate 0.05, 15 leaves, 10 samples per leaf, L2 1.0), early '
    'stopping **off** so the fit is deterministic and no training gene is held back. **No '
    'hyperparameter search**, so nothing was tuned on a held-out fold.',
    '**Two training variants per question, both honest, both reported.**',
    '* `within` — the model is fitted only on that question’s own genes. For the hard question '
    'that means it never sees a nuclear or ribosomal marker, so it cannot learn the easy contrast; '
    'it also has only ~84 training genes per fold.\n'
    '* `all-labels` — the model is fitted on all 585 map labels and then **scored only on that '
    'question’s genes**. More training data; still out-of-fold, so no test gene or paralog of one '
    'is ever in its own training set. The easy negatives can help it place the decision boundary but '
    'they cannot inflate the score, because the score is computed only on hard genes.\n'
    'Neither is a cheat and they answer slightly different deployment questions, so the verdict is '
    'stated for both and the headline takes the better one, named.')
code("""
    from sklearn.metrics import roc_auc_score, average_precision_score
    from sklearn.ensemble import HistGradientBoostingClassifier
    from parascape.claims import _isotonic, _ece

    def make_model():
        return HistGradientBoostingClassifier(
            max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=10,
            l2_regularization=1.0, early_stopping=False, random_state=0,
            categorical_features='from_dtype')

    def usable(features, rows):
        return [c for c in features if c in CATEGORICAL or X.loc[rows, c].notna().any()]

    def out_of_fold(features, lab, folds, score_on=None, verbose=False):
        '''Fit per fold on lab/folds; return out-of-fold scores for score_on (default: lab's genes).'''
        target = lab.index if score_on is None else pd.Index(score_on)
        target_fold = fold_all.reindex(target)
        oof = pd.Series(np.nan, index=target, dtype=float)
        models, used = {}, {}
        for f in sorted(set(folds.unique()) | set(target_fold.dropna().unique())):
            tr = lab.index[folds != f]
            te = target[target_fold == f]
            if not len(tr):
                continue
            cols = usable(features, tr)
            m = make_model().fit(X.loc[tr, cols], lab.loc[tr])
            if len(te):
                oof.loc[te] = m.predict_proba(X.loc[te, cols])[:, 1]
            models[f], used[f] = m, cols
        if verbose:
            print('   trained on %d genes, scored %d; features used per fold: %s'
                  % (len(lab), int(oof.notna().sum()), {f: len(c) for f, c in used.items()}))
        return oof, models, used

    def fold_metrics(truth_s, called, folds, score=None):
        rows = []
        for f in sorted(folds.unique()):
            m = (folds.loc[truth_s.index] == f).values
            if not m.any():
                continue
            t, c = truth_s.values[m], np.asarray(called)[m]
            tp, fp, fn = int((c & t).sum()), int((c & ~t).sum()), int((~c & t).sum())
            row = dict(fold=int(f), n=int(m.sum()), accuracy=float((c == t).mean()),
                       precision=(tp / (tp + fp)) if (tp + fp) else np.nan,
                       recall=(tp / (tp + fn)) if (tp + fn) else np.nan, calls=int(c.sum()))
            if score is not None and len(set(t)) == 2:
                s = np.asarray(score)[m]
                row['roc_auc'] = float(roc_auc_score(t, s))
                row['average_precision'] = float(average_precision_score(t, s))
            rows.append(row)
        return pd.DataFrame(rows)

    def board_rows(question, name, truth_s, called, folds, score=None, coverage=None):
        per = fold_metrics(truth_s, called, folds, score)
        out = []
        for metric in ('accuracy', 'precision', 'recall', 'roc_auc', 'average_precision', 'calls'):
            if metric not in per:
                continue
            v = per[metric].astype(float)
            row = dict(question=question, strategy=name, metric=metric, mean=v.mean(),
                       fold_min=v.min(), fold_max=v.max())
            row.update({f'fold_{int(f)}': val for f, val in zip(per.fold, v)})
            out.append(row)
        if coverage is not None:
            out.append(dict(question=question, strategy=name, metric='coverage', mean=coverage,
                            fold_min=coverage, fold_max=coverage))
        return out

    def cross_fitted_calibration(scores, lab, folds):
        cal = pd.Series(np.nan, index=scores.index, dtype=float)
        for f in sorted(folds.unique()):
            tr = lab.index[folds != f]
            te = scores.index[folds.reindex(scores.index) == f]
            if not len(tr) or not len(te):
                continue
            xs, ys = _isotonic(scores.loc[tr].dropna().values,
                               lab.loc[scores.loc[tr].dropna().index].astype(float).values)
            cal.loc[te] = np.interp(scores.loc[te].values, xs, ys)
        return cal

    def clopper_pearson_lower(k, n, confidence=0.90):
        if n == 0:
            return np.nan
        if k == 0:
            return 0.0
        return float(beta.ppf(1 - confidence, k, n - k + 1))

    THRESHOLDS = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97, 0.99]
    print('machinery ready; bound = exact one-sided Clopper-Pearson, beta.ppf(0.10, k, n-k+1)')
    pd.DataFrame([dict(calls=n, lower90_if_every_call_correct=round(clopper_pearson_lower(n, n), 3))
                  for n in (4, 8, 12, 16, 20, 22, 25, 30, 40, 60)]).set_index('calls').T
""")
nb.md(
    'That last table is the arithmetic that decides most of this notebook: **with fewer than about '
    '22 calls the lower 90 % bound cannot reach 0.90 even if every single call is correct.** The hard '
    'question has 58 positives in total, about 12 per fold. Reaching the bar needs either a very '
    'confident model making many calls or more labels — nothing else.')

# --------------------------------------------------------------------------- 6 baselines + run
nb.md(
    '## 6 · Baselines and the model, run once per question',
    '**Training-fold majority** — always the majority class of the training folds; makes no positive '
    'calls where negatives dominate, so its precision is undefined. In the hard question the majority '
    'is *host-facing*, which makes it a much stronger baseline there (prevalence 0.55).',
    '**Prevalence-weighted random** — calls each test gene host-facing with the training-fold '
    'prevalence as a coin weight, seed 0. Its precision is the prevalence, and that is the floor.',
    '**Signal peptide alone** — the prediction *is* `has_signal_peptide`. In the easy question this '
    'is the baseline to beat (run 1: precision 0.75). In the hard question nearly every gene has one '
    'by construction, so it degenerates to "call everything host-facing" and its precision is the '
    'prevalence — which is the point of conditioning.',
    '**Ortholog copy** — a gene takes the label of other labelled genes in its orthogroup, abstaining '
    'when there are none. Under orthogroup-grouped folds every labelled ortholog is in the same fold, '
    'so coverage is 0 by construction.')
code("""
    board, cals, oofs, thr_tables = [], {}, {}, {}

    for key, q in QUESTIONS.items():
        yq, fq, genes = q['y'], q['fold'], q['genes']
        prev = {f: float(yq[fq != f].mean()) for f in sorted(fq.unique())}
        maj = pd.Series([bool(yq[fq != f].mean() > 0.5) for f in fq], index=yq.index)
        board += board_rows(key, 'baseline: training majority', yq, maj.values, fq, coverage=1.0)
        rb = np.random.default_rng(0)
        coin = pd.Series([rb.random() < prev[f] for f in fq], index=yq.index)
        board += board_rows(key, 'baseline: prevalence-weighted random', yq, coin.values, fq,
                            score=np.array([prev[f] for f in fq]), coverage=1.0)
        spq = X.loc[genes, 'has_signal_peptide'].fillna(0).astype(bool)
        board += board_rows(key, 'baseline: signal peptide alone', yq, spq.values, fq,
                            score=spq.values.astype(float), coverage=1.0)
        pred, cov = [], []
        for g in genes:
            peers = [p for p in groups.index[(groups == groups[g]) & (groups.index != g)]
                     if p in yq.index and fold_all[p] != fold_all[g]]
            pred.append(bool(yq.loc[peers].mean() >= 0.5) if peers else False)
            cov.append(bool(peers))
        board += board_rows(key, 'baseline: ortholog copy (grouped, abstains)', yq,
                            np.array(pred), fq, coverage=float(np.mean(cov)))

        # model, both training variants
        print(key, '- GBT trained within the question set')
        o_within, m_within, u_within = out_of_fold(FEATURES, yq, fq, verbose=True)
        board += board_rows(key, 'GBT (trained within the question set)', yq, (o_within >= 0.5).values,
                            fq, score=o_within.values, coverage=1.0)
        oofs[(key, 'within')] = (o_within, m_within, u_within, yq, fq)
        if key == 'easy':
            o_all = o_within
            oofs[(key, 'all-labels')] = oofs[(key, 'within')]
        else:
            print(key, '- GBT trained on all 585 map labels, scored on this question only')
            o_all, m_all, u_all = out_of_fold(FEATURES, y, fold_all.loc[y.index],
                                              score_on=genes, verbose=True)
            board += board_rows(key, 'GBT (trained on all labels, scored here)', yq,
                                (o_all >= 0.5).values, fq, score=o_all.values, coverage=1.0)
            oofs[(key, 'all-labels')] = (o_all, m_all, u_all, yq, fq)
        print()
    scoreboard = pd.DataFrame(board)
    scoreboard.to_csv(OUT / 'scoreboard.tsv', sep='\\t', index=False)
    print('scoreboard rows:', len(scoreboard))
    scoreboard.pivot_table(index=['question', 'strategy'], columns='metric', values='mean')[
        ['accuracy', 'precision', 'recall', 'roc_auc', 'average_precision', 'calls', 'coverage']].round(3)
""")
code("""
    print('per-fold spread of the model against the best baseline in that fold')
    print('(instruction 70 asks the model to beat the best baseline in at least 4 of 5 folds)')
    fcols = [f'fold_{f}' for f in range(5)]
    for key in QUESTIONS:
        sub = scoreboard[(scoreboard.question == key)]
        for metric in ('average_precision', 'precision'):
            acc = sub[sub.metric == metric].set_index('strategy')
            base = acc.reindex(columns=fcols).loc[[s for s in acc.index if s.startswith('baseline')]]
            for strat in [s for s in acc.index if s.startswith('GBT')]:
                mod = acc.reindex(columns=fcols).loc[strat]
                beats = int((mod.values > base.max().values).sum())
                print('%-9s %-18s %-42s folds beaten: %d/5   model %s' %
                      (key, metric, strat, beats,
                       np.round(mod.values.astype(float), 3)))
        print()
""")

# --------------------------------------------------------------------------- 7 calibration
nb.md(
    '## 7 · Calibration, cross-fitted',
    'A boosted tree emits a score between 0 and 1 that is not a probability. Isotonic regression maps '
    'it to the observed rate of being right, using the project’s own machinery '
    '(`parascape.claims._isotonic`, `_ece`). **Cross-fitted**: for each fold the curve is fitted on '
    'the out-of-fold scores of the *other* folds and applied to this one, so the calibration is not '
    'scored on the points that fitted it.',
    'The reliability table says: of the genes whose calibrated probability fell in this range, how '
    'many were there and what fraction really were host-facing. The expected calibration error is the '
    'average gap, weighted by how many genes land in each bin. A small set cannot be well calibrated '
    'in the tails — the hard question has ~105 genes, so read its ECE as indicative.')
code("""
    rows = []
    edges = np.array([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
    for (key, variant), (o, _m, _u, yq, fq) in oofs.items():
        if key == 'easy' and variant == 'all-labels':
            continue
        c = cross_fitted_calibration(o, yq, fq)
        cals[(key, variant)] = c
        ece = _ece(c.dropna().values, yq.loc[c.dropna().index].astype(float).values,
                   bins=min(10, max(2, len(c.dropna()) // 12)))
        ece_raw = _ece(o.dropna().values, yq.loc[o.dropna().index].astype(float).values,
                       bins=min(10, max(2, len(o.dropna()) // 12)))
        cc = c.dropna(); yy = yq.loc[cc.index]
        idx = np.clip(np.searchsorted(edges, cc.values, side='right') - 1, 0, len(edges) - 2)
        for b in range(len(edges) - 1):
            m = idx == b
            if not m.any():
                continue
            rows.append(dict(question=key, model=variant, ece=round(ece, 4), ece_uncalibrated=round(ece_raw, 4),
                             predicted_range=f'[{edges[b]:.1f}, {edges[b + 1]:.1f})', n=int(m.sum()),
                             mean_predicted=round(float(cc.values[m].mean()), 4),
                             observed_positive_rate=round(float(yy.values[m].mean()), 4),
                             positives=int(yy.values[m].sum())))
    calibration = pd.DataFrame(rows)
    calibration.to_csv(OUT / 'calibration.tsv', sep='\\t', index=False)
    print('expected calibration error per question/model:')
    print(calibration.groupby(['question', 'model'])[['ece', 'ece_uncalibrated']].first().to_string())
    calibration[calibration.question == 'hard']
""")

# --------------------------------------------------------------------------- 8 thresholds
nb.md(
    '## 8 · The threshold tables: where the fill bar is, or is not, reached',
    'At each calibrated threshold: how many genes are called host-facing, how many of those calls are '
    'right, the precision, and the **Clopper–Pearson lower 90 % bound** on it. Then what the map '
    'would gain: `fillable_all` counts currently-empty genes above the threshold genome-wide, and '
    '`fillable_secretory` counts them inside the secretory set.',
    '`fillable_secretory` is the honest one for the hard question — a model fitted on secretory genes '
    'may only be applied to secretory genes. For the reviewed question both are extrapolations and '
    'are labelled as such.')
code("""
    all_cal, all_raw_scores = {}, {}
    empty_genes = [g for g in values.index[values['host_facing'].isna()] if g in X.index]
    empty_sec = [g for g in empty_genes if in_secretory.get(g, False)]
    print('currently-empty genes: %d total, %d of them secretory' % (len(empty_genes), len(empty_sec)))

    def threshold_table(key, variant):
        o, _m, _u, yq, fq = oofs[(key, variant)]
        c = cals[(key, variant)]
        train_lab = yq if variant == 'within' else y
        cols = usable(FEATURES, train_lab.index)
        final = make_model().fit(X.loc[train_lab.index, cols], train_lab)
        raw = pd.Series(final.predict_proba(X[cols])[:, 1], index=X.index)
        xs, ys = _isotonic(o.dropna().values, yq.loc[o.dropna().index].astype(float).values)
        cal_all = pd.Series(np.interp(raw.values, xs, ys), index=X.index)
        all_cal[(key, variant)] = cal_all
        all_raw_scores[(key, variant)] = (raw, final, cols)
        rows = []
        for t in THRESHOLDS:
            sel = c >= t
            n, k = int(sel.sum()), int(yq[sel].sum())
            lo = clopper_pearson_lower(k, n)
            rows.append(dict(question=key, model=variant, threshold=t, calls=n, correct=k,
                             precision=(k / n) if n else np.nan, precision_lower90=lo,
                             recall=k / int(yq.sum()),
                             fillable_all=int((cal_all.loc[empty_genes] >= t).sum()),
                             fillable_secretory=int((cal_all.loc[empty_sec] >= t).sum()),
                             reaches_bar=bool(n and lo >= 0.90)))
        return pd.DataFrame(rows)

    for key in QUESTIONS:
        variants = ['within'] if key == 'easy' else ['within', 'all-labels']
        thr_tables[key] = pd.concat([threshold_table(key, v) for v in variants], ignore_index=True)
        thr_tables[key].to_csv(OUT / f'thresholds_{key}.tsv', sep='\\t', index=False)
    print('threshold tables written for', list(thr_tables))

    def verdict(key):
        t = thr_tables[key]
        ok = t[t.reaches_bar]
        if len(ok):
            r = ok.sort_values(['threshold', 'model']).iloc[0]
            best = ok.loc[ok.precision_lower90.idxmax()]
            return dict(question=key, reaches_bar=True, model=r.model, threshold=r.threshold,
                        calls=int(r.calls), correct=int(r.correct), precision=r.precision,
                        lower90=r.precision_lower90, fillable_all=int(r.fillable_all),
                        fillable_secretory=int(r.fillable_secretory), best_lower90=best.precision_lower90)
        b = t.loc[t.precision_lower90.idxmax()]
        return dict(question=key, reaches_bar=False, model=b.model, threshold=b.threshold,
                    calls=int(b.calls), correct=int(b.correct), precision=b.precision,
                    lower90=b.precision_lower90, fillable_all=int(b.fillable_all),
                    fillable_secretory=int(b.fillable_secretory), best_lower90=b.precision_lower90)

    VERDICTS = pd.DataFrame([verdict(k) for k in QUESTIONS]).set_index('question')
    VERDICTS.round(4)
""")

nb.md(
    '## 9 · The HARD question — the headline',
    '**Among proteins that enter the secretory pathway, which are exported beyond the parasite?** '
    'This is the question instruction 72 names as the one that matters, and the one run 1 failed '
    '(best lower bound 0.868 on its weaker version of the set).')
code("""
    t = thr_tables['hard']
    print(t.round(4).to_string(index=False))
    print()
    v = VERDICTS.loc['hard']
    if v.reaches_bar:
        HARD_ANSWER = ('YES for the HARD question. The lowest calibrated threshold whose lower 90 %% bound of '
                       'precision reaches 0.90 is %.2f (%s model): %d held-out calls, %d correct, precision '
                       '%.3f, lower bound %.4f. At that threshold %d currently-empty secretory genes would be '
                       'filled (%d genome-wide).'
                       % (v.threshold, v.model, v.calls, v.correct, v.precision, v.lower90,
                          v.fillable_secretory, v.fillable_all))
    else:
        HARD_ANSWER = ('NO for the HARD question. No calibrated threshold from 0.50 to 0.99, under either '
                       'training variant, reaches a lower 90 %% confidence bound of 0.90. The best bound '
                       'anywhere on the curve is %.4f (%s model, threshold %.2f: %d calls, %d correct, '
                       'point-estimate precision %.3f). Among secretory-pathway proteins we still cannot say '
                       'at 90 %% which are exported, so nothing may be written into the information map as a '
                       'fact for this question.'
                       % (v.best_lower90, v.model, v.threshold, v.calls, v.correct, v.precision))
    print(HARD_ANSWER)
""")
code("""
    print('what the hard question is made of, and what the model does with it:')
    hg = QUESTIONS['hard']
    c_hard = cals[('hard', 'all-labels')] if ('hard', 'all-labels') in cals else cals[('hard', 'within')]
    for variant in ('within', 'all-labels'):
        o = oofs[('hard', variant)][0]
        print('  %-11s ROC-AUC %.3f   average precision %.3f'
              % (variant, roc_auc_score(hg['y'], o), average_precision_score(hg['y'], o)))
    print()
    print('calibrated out-of-fold probability by provenance (all-labels model):')
    pv = labels.provenance.reindex(hg['genes'])
    print(pd.DataFrame({'prob': c_hard, 'provenance': pv, 'label': hg['y']}
                       ).groupby(['label', 'provenance']).prob.agg(['size', 'median', 'min', 'max']).round(3).to_string())
""")

nb.md(
    '### Being fair to the model: precision at the top of the ranking',
    'The threshold table lives in **calibrated probability** space, so a "no" there could in '
    'principle mean only that the calibration refuses to ever say 0.90 — not that the ranking is '
    'bad. This check removes that doubt. For each question, the genes are sorted by their raw '
    'out-of-fold score and the precision of the top *k* is computed directly, with its '
    'Clopper–Pearson lower bound, for every *k* that could possibly clear the bar.',
    'If the best lower bound over all *k* is below 0.90, then **no threshold rule on this model '
    'could clear the bar**, whatever the calibration does. That is the strongest form of the '
    'negative answer.')
code("""
    rows = []
    for (key, variant), (o, _m, _u, yq, fq) in oofs.items():
        if key == 'easy' and variant == 'all-labels':
            continue
        order = o.sort_values(ascending=False).index
        truth_sorted = yq.loc[order].values
        for k in range(10, len(order) + 1):
            correct = int(truth_sorted[:k].sum())
            lo = clopper_pearson_lower(correct, k)
            rows.append(dict(question=key, model=variant, top_k=k, correct=correct,
                             precision=correct / k, precision_lower90=lo,
                             reaches_bar=bool(lo >= 0.90)))
    at_k = pd.DataFrame(rows)
    at_k.to_csv(OUT / 'precision_at_k.tsv', sep='\\t', index=False)
    best_k = at_k.loc[at_k.groupby(['question', 'model']).precision_lower90.idxmax()]
    ATK_NOTE = ('best precision anywhere on the ranking: ' +
                '; '.join('%s/%s k=%d %d/%d precision %.3f lower bound %.4f'
                          % (r.question, r.model, r.top_k, r.correct, r.top_k, r.precision,
                             r.precision_lower90) for r in best_k.itertuples()))
    print('for every k >= 10, the precision of the top k genes by raw out-of-fold score.')
    print('ANY threshold rule on this model is a choice of k, so the best bound over all k is the')
    print('ceiling for the whole approach:')
    print()
    print(ATK_NOTE.replace('; ', '\\n  '))
    print()
    print('does ANY k, on ANY question, clear the bar?', bool(at_k.reaches_bar.any()))
    best_k.set_index(['question', 'model']).round(4)
""")

nb.md(
    '## 10 · The HARDEST question — the 98 genes Einar ruled on',
    'Every label here is a human reading of published evidence: 66 host-facing, 32 not, with the '
    'negatives being vacuole residents rather than nuclear or ribosomal markers. No marker inference '
    'is involved anywhere. This is the set where a correct call means the most and where there is the '
    'least data to make one.')
code("""
    print(thr_tables['reviewed'].round(4).to_string(index=False))
    print()
    v = VERDICTS.loc['reviewed']
    if v.reaches_bar:
        REVIEWED_ANSWER = ('YES for the REVIEWED question: threshold %.2f (%s), %d calls, %d correct, precision '
                           '%.3f, lower bound %.4f; %d empty genes fillable genome-wide (%d secretory).'
                           % (v.threshold, v.model, v.calls, v.correct, v.precision, v.lower90,
                              v.fillable_all, v.fillable_secretory))
    else:
        REVIEWED_ANSWER = ('NO for the REVIEWED question. Best lower bound anywhere on the curve %.4f (%s, '
                           'threshold %.2f: %d calls, %d correct, precision %.3f).'
                           % (v.best_lower90, v.model, v.threshold, v.calls, v.correct, v.precision))
    print(REVIEWED_ANSWER)
""")

nb.md(
    '## 11 · The EASY question — for comparison only',
    '**This is the number run 1 reported and it should not be quoted on its own.** All 585 labels, '
    'including the 477 hyperLOPIT markers of compartments inside the parasite. It measures mostly '
    '"can we tell a nuclear protein from a secreted one", which a signal peptide nearly answers by '
    'itself. It is here so that the gap between it and §9 is visible.')
code("""
    print(thr_tables['easy'].round(4).to_string(index=False))
    print()
    v = VERDICTS.loc['easy']
    if v.reaches_bar:
        EASY_ANSWER = ('YES for the EASY question: threshold %.2f, %d calls, %d correct, precision %.3f, lower '
                       'bound %.4f, %d empty genes fillable. This is the flattering number.'
                       % (v.threshold, v.calls, v.correct, v.precision, v.lower90, v.fillable_all))
    else:
        EASY_ANSWER = ('NO even for the EASY question. Best lower bound %.4f at threshold %.2f (%d calls, %d '
                       'correct).' % (v.best_lower90, v.threshold, v.calls, v.correct))
    print(EASY_ANSWER)
    print()
    print('side by side, the three questions:')
    VERDICTS[['reaches_bar', 'model', 'threshold', 'calls', 'correct', 'precision', 'lower90',
              'best_lower90', 'fillable_all', 'fillable_secretory']].round(4)
""")

# --------------------------------------------------------------------------- 12 provenance
nb.md(
    '### Where did run 1’s 0.906 go? — putting the held-out families back',
    'Two things changed between run 1 and this run: the **labels** (22 dense granule proteins moved '
    'from positive to negative, 43 hard negatives instead of 9) and the **features** (the declared '
    'dependencies are now held out). It is worth knowing which of the two destroyed the apparent '
    'result, so the question is asked directly: refit the easy question **with** '
    '`hosttx_effectors`, `organelle_surface_turboid`, `pubmed`, `pvm_proximity_positive`, '
    '`n_host_targets` and the attention columns put back, on the *new* labels.',
    '**This is a leakage demonstration, not a result.** Nothing from this cell may be used to fill '
    'anything: those columns are held out by instruction 72 precisely because they stand too close '
    'to the conclusion. It is run only to attribute the change.')
code("""
    LEAKY = sorted(HELD)
    leaky_cols = [c for c in LEAKY if c in nodes.columns and
                  (pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c]))]
    X_leaky = X.join(nodes[leaky_cols].astype('float64'), rsuffix='_held')
    rows = []
    for key in ('easy', 'hard'):
        yq, fq = QUESTIONS[key]['y'], QUESTIONS[key]['fold']
        # At least 10 measured values: a column observed on one or two genes has too few distinct values for
        # sklearn's binning and raises "window shape cannot be larger than input array shape".
        cols = [c for c in X_leaky.columns if c in CATEGORICAL or X_leaky.loc[yq.index, c].notna().sum() >= 10]
        oof_leak = pd.Series(np.nan, index=yq.index, dtype=float)
        for f in sorted(fq.unique()):
            tr, te = yq.index[fq != f], yq.index[fq == f]
            m = make_model().fit(X_leaky.loc[tr, cols], yq.loc[tr])
            oof_leak.loc[te] = m.predict_proba(X_leaky.loc[te, cols])[:, 1]
        c = cross_fitted_calibration(oof_leak, yq, fq)
        best = max(((clopper_pearson_lower(int(yq[c >= t].sum()), int((c >= t).sum())), t,
                     int((c >= t).sum()), int(yq[c >= t].sum()))
                    for t in THRESHOLDS if int((c >= t).sum()) > 0), default=(np.nan,) * 4)
        for name, o in (('declared dependencies HELD OUT (this run)', oofs[(key, 'within')][0]),
                        ('declared dependencies PUT BACK (leaky, run 1 style)', oof_leak)):
            bb = best if 'PUT BACK' in name else None
            if bb is None:
                cc = cals[(key, 'within')]
                bb = max(((clopper_pearson_lower(int(yq[cc >= t].sum()), int((cc >= t).sum())), t,
                           int((cc >= t).sum()), int(yq[cc >= t].sum()))
                          for t in THRESHOLDS if int((cc >= t).sum()) > 0), default=(np.nan,) * 4)
            rows.append(dict(question=key, feature_set=name, n_features=len(cols) if 'PUT BACK' in name else len(FEATURES),
                             roc_auc=roc_auc_score(yq, o), average_precision=average_precision_score(yq, o),
                             best_lower90=bb[0], at_threshold=bb[1], calls=bb[2], correct=bb[3],
                             reaches_bar=bool(bb[0] >= 0.90)))
    leak_test = pd.DataFrame(rows)
    leak_test.to_csv(OUT / 'leakage_attribution.tsv', sep='\\t', index=False)
    def pick(q, token):
        m = (leak_test.question == q) & leak_test.feature_set.str.contains(token, regex=False)
        return leak_test[m].iloc[0]
    e_held, e_back = pick('easy', 'HELD'), pick('easy', 'PUT BACK')
    h_back = pick('hard', 'PUT BACK')
    LEAK_NOTE = ('On the NEW labels, putting the declared dependencies back moves the easy question from '
                 'AUC %.3f / bound %.4f to AUC %.3f / bound %.4f, and the hard question to AUC %.3f / '
                 'bound %.4f%s. So run 1\\'s 0.906 was carried by BOTH changes: those feature families '
                 'did most of the work, and the 22 rejected dense granule proteins are now negatives that '
                 'the model gets wrong either way.'
                 % (e_held.roc_auc, e_held.best_lower90, e_back.roc_auc, e_back.best_lower90,
                    h_back.roc_auc, h_back.best_lower90,
                    ' - still short of 0.90 on the hard question' if not h_back.reaches_bar else ''))
    print(LEAK_NOTE)
    leak_test.round(4)
""")

nb.md(
    '## 12 · Does it do worse on Einar’s calls than on the internal markers?',
    'This is the diagnosis run 1 failed. If the model separates positives from hyperLOPIT internal '
    'markers but not from curated vacuole residents, it is telling inside from outside the '
    '**parasite**, not inside from outside the **vacuole** — and only the second is the biological '
    'question. Each row below pits one positive provenance against one negative provenance, so the '
    'AUC is comparable across rows. The all-labels model is used, since that is the one whose scores '
    'cover every gene.',
    '`false_positive_rate_at_0_5` on the `gold_negative` and `einar_*` rows is the number to look at: '
    'it is the rate at which the model would claim a vacuole-resident protein reaches the host.')
code("""
    o_easy = oofs[('easy', 'within')][0]
    c_easy = cals[('easy', 'within')]
    pv = labels.provenance
    POSP = ['einar_direct', 'einar_curated', 'host_interaction']
    NEGP = ['einar_curated', 'einar_direct', 'gold_negative', 'hyperlopit_internal']
    rows = []
    for pname in POSP + ['all positives']:
        pmask = (y & (pv == pname)) if pname != 'all positives' else y
        for nname in NEGP + ['einar+gold (the hard negatives)', 'all negatives']:
            if nname == 'all negatives':
                nmask = ~y
            elif nname.startswith('einar+gold'):
                nmask = (~y) & pv.isin(['einar_direct', 'einar_curated', 'gold_negative'])
            else:
                nmask = (~y) & (pv == nname)
            idx = y.index[pmask.fillna(False) | nmask.fillna(False)]
            t = y.loc[idx]
            if t.sum() < 3 or (~t).sum() < 3:
                continue
            s = o_easy.loc[idx]
            rows.append(dict(positives=pname, negatives=nname, n_pos=int(t.sum()), n_neg=int((~t).sum()),
                             roc_auc=roc_auc_score(t, s), average_precision=average_precision_score(t, s),
                             recall_at_0_5=float((s[t] >= 0.5).mean()),
                             false_positive_rate_at_0_5=float((s[~t] >= 0.5).mean()),
                             median_score_pos=float(s[t].median()), median_score_neg=float(s[~t].median())))
    by_prov = pd.DataFrame(rows)
    by_prov.to_csv(OUT / 'by_label_provenance.tsv', sep='\\t', index=False)
    by_prov.round(3)
""")
code("""
    print('per-provenance error rates, all three questions, at the 0.5 cut:')
    rows = []
    for key in QUESTIONS:
        for variant in (['within'] if key == 'easy' else ['within', 'all-labels']):
            o, _m, _u, yq, fq = oofs[(key, variant)]
            p = labels.provenance.reindex(yq.index)
            if key == 'reviewed':
                p = p.fillna('einar (his ruling)')
            for (lab_v, prov_v), grp in yq.groupby([yq, p]):
                s = o.loc[grp.index]
                rows.append(dict(question=key, model=variant, label=bool(lab_v), provenance=prov_v,
                                 n=len(grp), median_score=round(float(s.median()), 3),
                                 error_rate_at_0_5=round(float((s >= 0.5).mean() if not lab_v
                                                               else (s < 0.5).mean()), 3)))
    prov_errors = pd.DataFrame(rows)
    prov_errors.to_csv(OUT / 'by_label_provenance_errors.tsv', sep='\\t', index=False)
    print(prov_errors.to_string(index=False))
""")
nb.md(
    '### The hard question with the two contested genes flipped',
    'GRA9 and NTPase II are negatives in the map and positives in Einar’s review (§1). If the '
    'verdict moves when two of 47 hard negatives change side, the verdict is not robust and the '
    'reader should know.')
code("""
    y_flip = y.copy()
    y_flip.loc[EINAR_FLIPS] = ~y_flip.loc[EINAR_FLIPS]
    g_hard = QUESTIONS['hard']['genes']
    yq_f, fq_f = y_flip.loc[g_hard], fold_all.loc[g_hard]
    o_f, _, _ = out_of_fold(FEATURES, yq_f, fq_f)
    c_f = cross_fitted_calibration(o_f, yq_f, fq_f)
    rows = []
    for t in THRESHOLDS:
        sel = c_f >= t
        n, k = int(sel.sum()), int(yq_f[sel].sum())
        rows.append(dict(threshold=t, calls=n, correct=k, precision=(k / n) if n else np.nan,
                         precision_lower90=clopper_pearson_lower(k, n),
                         reaches_bar=bool(n and clopper_pearson_lower(k, n) >= 0.90)))
    flip = pd.DataFrame(rows)
    FLIP_VERDICT = ('with GRA9 and NTPase II as positives the hard question %s (best lower bound %.4f)'
                    % ('STILL reaches the bar' if flip.reaches_bar.any() else 'still does not reach the bar',
                       flip.precision_lower90.max()))
    print(FLIP_VERDICT)
    flip.round(4)
""")

# --------------------------------------------------------------------------- 13 importance
nb.md(
    '## 13 · Permutation importance, with the family of every feature',
    '`sklearn.inspection.permutation_importance` on the **held-out fold only**: a column is shuffled '
    'and the loss of average precision on genes the model never saw is recorded. Five shuffles per '
    'feature per fold, seed 0, averaged over folds. Run for the hard question (the headline) and for '
    'the all-label model.',
    'With ~20 held-out genes per fold in the hard question the individual importances are small and '
    'their signs flip between folds. Read the top of the table as a direction, not as a ranking to be '
    'trusted feature by feature — and read the **family** column rather than the feature names, '
    'because that is the unit of evidence.')
code("""
    from sklearn.inspection import permutation_importance
    imp_rows = []
    for key, variant in (('hard', 'within'), ('easy', 'within')):
        o, models, used, yq, fq = oofs[(key, variant)]
        per_fold = {}
        for f, m in models.items():
            te = yq.index[fq == f]
            if len(te) < 5 or yq.loc[te].nunique() < 2:
                continue
            r = permutation_importance(m, X.loc[te, used[f]], yq.loc[te], scoring='average_precision',
                                       n_repeats=5, random_state=0, n_jobs=1)
            per_fold[f'fold_{f}'] = pd.Series(r.importances_mean, index=used[f])
        d = pd.DataFrame(per_fold)
        d['importance_mean'] = d.mean(axis=1)
        d['importance_sd'] = d.std(axis=1)
        d['family'] = fam.reindex(d.index)
        d['question'] = key
        d.index.name = 'feature'
        imp_rows.append(d.sort_values('importance_mean', ascending=False))
    importance = pd.concat(imp_rows)
    importance.to_csv(OUT / 'importance.tsv', sep='\\t')
    IMP_HARD = imp_rows[0]
    print('HARD question, top 20 features by permutation importance:')
    print(IMP_HARD.head(20)[['family', 'importance_mean', 'importance_sd']].round(4).to_string())
    print()
    print('summed over families, top 12:')
    print(IMP_HARD.head(40).groupby('family').importance_mean.agg(['size', 'sum'])
          .sort_values('sum', ascending=False).head(12).round(4).to_string())
    IMP_HARD.head(20)[['family', 'importance_mean', 'importance_sd']].round(4)
""")
code("""
    print('EASY question, top 15 - kept so run 1 can be compared (its top two features are now held out):')
    print(imp_rows[1].head(15)[['family', 'importance_mean', 'importance_sd']].round(4).to_string())
    print()
    for c in ('has_signal_peptide', 'n_proximity_partners', 'tdi_coverage'):
        if c in IMP_HARD.index:
            print('  %-24s hard-question rank %d of %d, mean importance %+.4f'
                  % (c, list(IMP_HARD.index).index(c) + 1, len(IMP_HARD), IMP_HARD.loc[c, 'importance_mean']))
""")
nb.md(
    '### Robustness: does the verdict survive removing a family?',
    'Run 1’s single most important lesson was that its answer moved when one feature family was '
    'removed. So the hard question is refitted with each worrying block taken out, and the best lower '
    'bound is reported each time. **A verdict that depends on one family is not a verdict.**')
code("""
    ATT_ADJACENT = sorted(fam.index[fam.isin(['bioid_corpus_membership', 'ipms_corpus_membership',
                                              'ipms_baits', 'pmc_oa', 'pride_proximity'])])
    STRUCT_AVAIL = [c for c in ('has_structure', 'tdi_coverage', 'tdi_source', 'tdi_n_fragments',
                                'tdi_align_identity', 'tdi_align_coverage') if c in FEATURES]
    HOSTADJ = sorted(fam.index[fam.isin(['secretome_partition', 'toxodb_epitopes', 'iedb_bcell',
                                         'crosslink_interactome', 'cyst_wall_interactome',
                                         'spacr_escrt_screen', 'gra12', 'gra17_synthlethal'])])
    sets = {
        'all features': FEATURES,
        'minus foldseek_3di': [c for c in FEATURES if c not in set(FOLDSEEK_3DI)],
        'minus sequence/topology traits': [c for c in FEATURES if c not in set(SEQ_FEATURES)],
        'minus attention-adjacent (BioID/IP-MS corpus membership, PMC, proximity counts)':
            [c for c in FEATURES if c not in set(ATT_ADJACENT)],
        'minus structure availability': [c for c in FEATURES if c not in set(STRUCT_AVAIL)],
        'minus host-facing-adjacent screens': [c for c in FEATURES if c not in set(HOSTADJ)],
        'minus all of the above': [c for c in FEATURES if c not in
                                   set(ATT_ADJACENT) | set(STRUCT_AVAIL) | set(HOSTADJ)],
    }
    rows = []
    for name, cols in sets.items():
        for variant in ('within', 'all-labels'):
            train = QUESTIONS['hard']['y'] if variant == 'within' else y
            tf = QUESTIONS['hard']['fold'] if variant == 'within' else fold_all.loc[y.index]
            o, _, _ = out_of_fold(cols, train, tf, score_on=QUESTIONS['hard']['genes'])
            yq = QUESTIONS['hard']['y']
            c = cross_fitted_calibration(o, yq, QUESTIONS['hard']['fold'])
            best = max(((clopper_pearson_lower(int(yq[c >= t].sum()), int((c >= t).sum())), t,
                         int((c >= t).sum()), int(yq[c >= t].sum()))
                        for t in THRESHOLDS if int((c >= t).sum()) > 0), default=(np.nan,) * 4)
            rows.append(dict(feature_set=name, model=variant, n_features=len(cols),
                             roc_auc=roc_auc_score(yq, o), average_precision=average_precision_score(yq, o),
                             best_lower90=best[0], at_threshold=best[1], calls=best[2], correct=best[3],
                             reaches_bar=bool(best[0] >= 0.90)))
    robust = pd.DataFrame(rows)
    robust.to_csv(OUT / 'robustness_hard.tsv', sep='\\t', index=False)
    ROBUST_NOTE = ('hard question: %d of %d feature sets x training variants reach the bar; best lower bound '
                   'ranges %.3f to %.3f' % (int(robust.reaches_bar.sum()), len(robust),
                                            robust.best_lower90.min(), robust.best_lower90.max()))
    print(ROBUST_NOTE)
    robust.round(4)
""")

# --------------------------------------------------------------------------- 14 predictions
nb.md(
    '## 14 · Predictions for every gene, and what pushed each one up',
    'A calibrated probability for all 8,140 genes from each question’s model: fitted on that '
    'question’s training labels, mapped through the isotonic curve fitted on its pooled '
    'out-of-fold scores. For labelled genes the out-of-fold probability is given too — **that is the '
    'honest one**, and it is what every score in this notebook uses.',
    '**The contributing features are measured, not guessed.** For each gene, each of the top 20 '
    'features of the hard model is set to missing (which the tree handles natively) and the fall in '
    'that gene’s own probability is recorded; the three largest falls are reported. So '
    '`top_features` means "take this measurement away and the call weakens most", per gene, not a '
    'global ranking.',
    '**Nothing here is written into the information map.** Whether any of it may be is decided by the '
    'bounds in §9–§11.')
code("""
    key_model = 'hard'
    raw, final, cols = all_raw_scores[(key_model, 'all-labels')]
    occlude = [c for c in IMP_HARD.head(20).index if c not in CATEGORICAL and c in cols]
    base = raw.values
    deltas = {}
    for c in occlude:
        Xc = X[cols].copy(); Xc[c] = np.nan
        deltas[c] = base - final.predict_proba(Xc)[:, 1]
    deltas = pd.DataFrame(deltas, index=X.index)

    def top_three(row):
        s = row.sort_values(ascending=False)
        return '; '.join(f'{k}={v:+.3f}' for k, v in s.head(3).items() if abs(v) > 1e-4) or 'none above 0.0001'

    predictions = pd.DataFrame(index=X.index)
    predictions.index.name = 'gene_id'
    for (k, variant), c in all_cal.items():
        predictions[f'prob_{k}_{variant}'] = c.round(4)
    for (k, variant), c in cals.items():
        predictions[f'oof_{k}_{variant}'] = c.reindex(predictions.index).round(4)
    predictions['label_map'] = y.reindex(predictions.index)
    predictions['label_einar'] = y_einar.reindex(predictions.index)
    predictions['provenance'] = labels.provenance.reindex(predictions.index)
    predictions['fold'] = fold_all.reindex(predictions.index)
    predictions['in_secretory'] = in_secretory.reindex(predictions.index)
    predictions['has_signal_peptide'] = sp_all.reindex(predictions.index)
    predictions['hyperlopit_marker'] = marker.reindex(predictions.index)
    predictions['top_features'] = deltas.apply(top_three, axis=1)
    predictions.sort_values('prob_hard_all-labels', ascending=False).to_csv(
        OUT / 'predictions.tsv', sep='\\t')
    print('predictions written for', len(predictions), 'genes')
    print()
    print('highest-scoring UNTESTED secretory genes (exploratory, not facts):')
    cand = predictions.loc[[g for g in empty_sec]].sort_values('prob_hard_all-labels', ascending=False)
    print(cand.head(15)[['prob_hard_all-labels', 'has_signal_peptide', 'hyperlopit_marker', 'top_features']].to_string())
""")
code("""
    print('distribution of the hard model probability over untested genes:')
    print(all_cal[('hard', 'all-labels')].loc[empty_genes].describe(percentiles=[.5, .9, .99, .999]).round(4).to_string())
    print()
    print('untested genes at or above each threshold (all / secretory only):')
    print(pd.DataFrame({t: [int((all_cal[('hard', 'all-labels')].loc[empty_genes] >= t).sum()),
                            int((all_cal[('hard', 'all-labels')].loc[empty_sec] >= t).sum())]
                        for t in THRESHOLDS}, index=['all empty', 'empty secretory']).to_string())
""")

# --------------------------------------------------------------------------- 15 limits
nb.md(
    '## 15 · Limitations, and what would actually move the hard question',
    'None of these are softened. The first four are the ones that decide the answer.',
    '1. **The hard set is 105 genes with 58 positives and 47 negatives.** Over five folds that is '
    '~12 positives per fold. The Clopper–Pearson bound cannot reach 0.90 on fewer than ~22 calls even '
    'if every call is correct, so for most of the threshold range the bar is unreachable as '
    '*arithmetic*, before any question of model quality. **This is a data-quantity result before it '
    'is a modelling result.**',
    '2. **Two of the 47 hard negatives are contested** (GRA9, NTPase II: Einar’s review says '
    'host-facing, the curated gold-negative set says vacuole-resident). §12 reruns the question with '
    'them flipped.',
    '3. **The 43 hard negatives are a curation decision about published papers, not an export '
    'assay.** "Shown to stay in the vacuole" usually means immunogold EM or IFA in one strain at one '
    'time point. A protein can be vacuolar in a tachyzoite and exported in a bradyzoite, and nothing '
    'in these labels would see it.',
    '4. **DeepTMHMM misses signal peptides** — GRA16 (TGME49_208830), a secreted effector, is called '
    'GLOBULAR. The secretory-marker arm of the definition recovers some of those genes, but a protein '
    'that is neither predicted nor a hyperLOPIT marker is absent from the hard question entirely. '
    'SignalP or a curated cleavage table would fix it; recorded as an open gap.',
    '5. **The remaining positives are still mostly one assay.** Most accepted claims trace to the PVM '
    'proximity screen (PMID 34749525) read by a human. Einar’s review is what makes them usable, '
    'and it is one reader’s judgement on one screen plus the literature.',
    '6. **Publication bias is reduced, not removed.** `pubmed`, `hosttx_effectors` and the attention '
    'columns are held out, but a famous gene is also the gene with a melting temperature, a palmitome '
    'entry, a phosphoproteome row and a deposited structure. *Missingness itself carries fame*, and no '
    'hold-out removes that. §13 removes the attention-adjacent corpus-membership columns as a check.',
    '7. **The 3Di composition family is a crude shape summary**, not structural comparison, and '
    '637 genes have no 3Di at all while 520 have one copied from a near-identical structure by '
    'alignment. The quality columns travel with it, as instruction 72 requires, but structure '
    '*availability* is itself fame-correlated.',
    '8. **The all-gene probability columns are not independently validated.** Each comes from a model '
    'fitted on all of its question’s labels and an isotonic curve fitted in-sample on pooled '
    'out-of-fold scores. Only the cross-fitted out-of-fold numbers are honest; the all-gene columns '
    'exist to count what *would* be fillable.',
    '9. **The ortholog-copy baseline has coverage 0** under orthogroup-grouped folds, by '
    'construction, and is reported rather than hidden.',
    '10. **One model, fixed settings, no search.** Stages 2–5 of the instruction 72 plan (frozen '
    'sequence/structure blocks, guilt by association, the joint row model, pairs) are not in this '
    'notebook.',
    '11. **Permutation importance is noisy at this size** — ~20 held-out genes per fold in the hard '
    'question. Signs flip between folds. Read families, not features.')
nb.md(
    '### More labels, or a better measurement? — the arithmetic, not a preference',
    'Two different things can block the bar and they need different fixes, so they are separated here '
    'rather than merged into "we need more data".',
    '* **If the model’s true precision at its operating point is already above 0.90**, the only '
    'obstacle is sample size: the Clopper–Pearson bound converges upward to the true precision, so '
    'enough calls will demonstrate it. The table below says how many.',
    '* **If the true precision is below 0.90**, the bound converges to that precision and **no number '
    'of labels will ever clear the bar.** More labels then only measure the failure more precisely. '
    'The fix has to be better evidence or a better model.',
    'The code below works out which case this is.')
code("""
    yq = QUESTIONS['hard']['y']
    rows = []
    for target_prec in (1.00, 0.98, 0.97, 0.95, 0.93, 0.92):
        n = 1
        while n < 20000 and clopper_pearson_lower(int(round(n * target_prec)), n) < 0.90:
            n += 1
        rows.append(dict(true_precision=target_prec,
                         calls_needed=n if n < 20000 else np.nan))
    needed = pd.DataFrame(rows)
    print('calls needed for the lower 90 % bound to reach 0.90, by the model\\'s TRUE precision:')
    print(needed.to_string(index=False))
    print('(at a true precision of 0.90 or below the bound never reaches 0.90, at any number of calls)')
    print()
    hard_best = thr_tables['hard'].loc[thr_tables['hard'].precision_lower90.idxmax()]
    atk_best = at_k[at_k.question == 'hard'].sort_values('precision_lower90').iloc[-1]
    print('hard question, best calibrated threshold: %.2f (%s), %d calls, %d correct, precision %.3f, '
          'lower bound %.4f' % (hard_best.threshold, hard_best.model, hard_best.calls,
                                hard_best.correct, hard_best.precision, hard_best.precision_lower90))
    print('hard question, best point on the raw ranking: top %d, %d correct, precision %.3f, '
          'lower bound %.4f' % (atk_best.top_k, atk_best.correct, atk_best.precision,
                                atk_best.precision_lower90))
    print()
    obs_p = float(max(hard_best.precision, atk_best.precision))
    if obs_p < 0.90:
        WHAT_IS_NEEDED = (
            'WHICH CASE IS THIS? The second one. The best precision this model reaches anywhere - at any '
            'calibrated threshold and at any point on the raw ranking - is %.3f, which is BELOW the 0.90 bar '
            'itself. The lower confidence bound converges upward to the true precision, so with a precision '
            'of %.3f no number of additional labels can ever bring the bound to 0.90: more labels would '
            'simply measure %.3f more precisely. The limiting resource for the hard question is therefore '
            'NOT the number of labels - it is the evidence. Something has to be measured that distinguishes '
            '"in the vacuole" from "in the host cell", because the features available now do not: against '
            'the 43 hard negatives the model has ROC-AUC %.2f, against the hyperLOPIT internal markers it '
            'has %.2f, and that gap is the whole result. Only once a model clears 0.90 precision does the '
            'label count become the binding constraint, and at a true precision of 0.95 that would take '
            'about %d held-out calls, i.e. roughly %d labelled secretory genes at the current call rate '
            'against the %d that exist today.'
            % (obs_p, obs_p, obs_p,
               float(by_prov[by_prov.positives == 'all positives']
                     .set_index('negatives').loc['einar+gold (the hard negatives)', 'roc_auc']),
               float(by_prov[by_prov.positives == 'all positives']
                     .set_index('negatives').loc['hyperlopit_internal', 'roc_auc']),
               int(needed.set_index('true_precision').loc[0.95, 'calls_needed']),
               int(round(int(needed.set_index('true_precision').loc[0.95, 'calls_needed'])
                         / max(hard_best.calls / len(yq), 1e-9))),
               len(yq)))
    else:
        need = int(needed.set_index('true_precision').loc[0.95, 'calls_needed'])
        call_rate = hard_best.calls / len(yq)
        WHAT_IS_NEEDED = (
            'WHICH CASE IS THIS? The first one. The observed precision is %.3f, above the bar, so the block '
            'is sample size: %d calls are needed at a true precision of 0.95, the model calls %.0f %% of '
            'held-out secretory genes, so about %d labelled secretory genes in total against the %d there '
            'are now - roughly %d more reviewed secretory genes.'
            % (obs_p, need, 100 * call_rate, int(round(need / call_rate)), len(yq),
               max(0, int(round(need / call_rate)) - len(yq))))
    print(WHAT_IS_NEEDED)
""")
nb.md(
    '### What would move it, instead of lowering the bar',
    'In the order the numbers above imply, not in order of convenience:',
    '1. **A measurement that separates "in the vacuole" from "in the host cell" directly and '
    'genome-wide.** This is the binding constraint whenever the observed precision sits below the '
    'bar, because no number of labels fixes a model that is wrong 20 % of the time. Everything in '
    'the current feature set that correlates with host-facing also correlates with *being secreted*, '
    'which both classes are. The missing experiment is a **differential** one in the same cells: '
    'selective permeabilisation, or proximity labelling anchored in the host cytosol with a '
    'PV-lumen-anchored control, scored per gene with an FDR. The two existing genome-wide host-side '
    'assays are held out here for good reason — the PVM proximity screen *defines* the positives, '
    'and the host-transcription effector screen is downstream of the conclusion — so neither can '
    'serve as the independent measurement.',
    '2. **More reviewed labels of the same kind as Einar’s 98** — published, secretory, ruled on '
    'by a human, and especially more hard *negatives*. These are what make the bound demonstrable '
    'once a model is good enough, and they are cheap in the sense that no new experiment is needed. '
    'The cell above gives the count for a model at 95 % precision.',
    '3. **SignalP (or curated cleavage sites)** so the secretory set is defined correctly rather than '
    'by a predictor that misses known effectors such as GRA16.',
    'What is **not** on this list: lowering the bar, dropping the lower-bound requirement, filling on '
    'the point estimate, or quoting the easy question. The 90 % bar *with a lower bound* is the whole '
    'point of the fill rule, and the right action when it is not met is to fill nothing.')

# --------------------------------------------------------------------------- 16 README
code("""
    def fmt(df, index=True):
        return '```\\n' + df.to_string(index=index) + '\\n```\\n'

    top10 = IMP_HARD.head(10)[['family', 'importance_mean']].round(4)
    R = []
    R.append('# Host-facing inference v2 on Einar\\'s reviewed labels (2026-11-04)\\n')
    R.append('Built by `scripts/build_host_facing_inference_v2_notebook.py`; the executed notebook is\\n'
             '`notebooks/host_facing_inference_v2_2026_11_04.ipynb`. Rerun of\\n'
             '`results/host_facing_inference_2026_11_04/` after Einar reviewed all 140 published claims.\\n'
             'Run 1\\'s outputs are untouched.\\n')
    R.append('## Why this rerun\\n')
    R.append('Run 1 reported ROC-AUC 0.995 and a lower bound of 0.906, and it was not trustworthy: its\\n'
             'negatives were 477 markers of compartments inside the parasite against 9 curated vacuole\\n'
             'residents, so it separated "inside the parasite" from "secreted" rather than "stays in the\\n'
             'vacuole" from "exported into the host"; against the 9 hard negatives its AUC fell to 0.855 with\\n'
             'a 33 % false-positive rate; and its top two features were a host-transcription effector screen\\n'
             'and a publication count. Einar then rejected 22 dense granule proteins the PVM proximity screen\\n'
             'had called positive, so the map now holds 43 vacuole-resident or curated negatives.\\n')
    R.append('## The three questions and the answer to each\\n')
    R.append(fmt(VERDICTS[['reaches_bar', 'model', 'threshold', 'calls', 'correct', 'precision',
                           'lower90', 'best_lower90', 'fillable_all', 'fillable_secretory']].round(4)))
    R.append('\\n**1 · HARD (the headline).** ' + HARD_ANSWER + '\\n')
    R.append('\\n**2 · HARDEST (Einar-reviewed genes only).** ' + REVIEWED_ANSWER + '\\n')
    R.append('\\n**3 · easy (comparison only, the number run 1 reported).** ' + EASY_ANSWER + '\\n')
    R.append('\\nRobustness of the hard question across feature sets: ' + ROBUST_NOTE + '.\\n')
    R.append('\\n**Where run 1\\'s 0.906 went.** ' + LEAK_NOTE + ' `leakage_attribution.tsv` (a leakage\\n'
             'demonstration, not a result - nothing from it may fill anything).\\n')
    R.append('\\n**Not a calibration artefact.** Ranking the genes by raw out-of-fold score and taking the\\n'
             'top k directly, for every k from 10 upward — which is what *any* threshold rule amounts to —\\n'
             'the bar is cleared by no k on any question. ' + ATK_NOTE + '. `precision_at_k.tsv`.\\n')
    R.append('With the two contested genes flipped (GRA9, NTPase II): ' + FLIP_VERDICT + '.\\n')
    R.append('\\n' + WHAT_IS_NEEDED + '\\n')
    R.append('## Definitions\\n')
    R.append('* **Secretory-pathway set (the hard question):** a predicted DeepTMHMM signal peptide OR a\\n'
             '  hyperLOPIT marker of dense granules / rhoptries 1 / rhoptries 2 / micronemes. %d genes\\n'
             '  genome-wide, %d of them labelled (%d host-facing / %d not), %d untested. %d of the %d\\n'
             '  host-facing positives are outside it (no predicted signal peptide, no secretory marker).\\n'
             % (int(in_secretory.sum()), len(QUESTIONS['hard']['y']), int(QUESTIONS['hard']['y'].sum()),
                int((~QUESTIONS['hard']['y']).sum()), len(empty_sec),
                int((~in_secretory.loc[y.index[y]]).sum()), int(y.sum())))
    R.append('* **Reviewed set (the hardest question):** the 98 genes Einar ruled on, with his decisions as\\n'
             '  labels (66 yes / 32 no).\\n')
    R.append('* **Easy set:** all %d labelled genes (%d / %d), hyperLOPIT internal markers included.\\n'
             % (len(y), int(y.sum()), int((~y).sum())))
    R.append('* **Label provenance** travels with every gene: ' +
             ', '.join('%s (%d)' % (k, v) for k, v in labels.provenance.value_counts().items()) + '.\\n')
    R.append('* **Conflicts:** all 33 genes in `conflicts.tsv` have Einar\\'s review as a source, so all 33\\n'
             '  are resolved and kept. None dropped (run 1 dropped 8).\\n')
    R.append('* **Two genes where the map and his review disagree:** GRA9 (TGME49_251540) and NTPase II\\n'
             '  (TGME49_277270). The map says not host-facing (curated gold negative, 0.90); his review\\n'
             '  accepted the PVM proximity claim (0.90). Questions 1 and 3 use the map column as the\\n'
             '  instruction names it; question 2 uses his decisions; the flip test is in the notebook.\\n')
    R.append('## Method\\n')
    R.append('* Features: **%d** — numeric/boolean `nodes.parquet` columns, the DeepTMHMM sequence and\\n'
             '  topology traits, the 3Di quality columns (`tdi_coverage`, `tdi_align_identity`,\\n'
             '  `tdi_align_coverage`, `tdi_n_fragments`, `tdi_source`) and a new `foldseek_3di` family of\\n'
             '  %d 3Di composition features. Full list with families is printed in the notebook and in\\n'
             '  `importance.tsv`.\\n'
             '* Held out as declared dependencies of `host_facing` (**%d columns**): every column of\\n'
             '  `hosttx_effectors` (%d), `organelle_surface_turboid` (%d), `pubmed` (%d) and\\n'
             '  `pvm_proximity` (%d); `n_host_targets`; the attention columns (%s); and everything matching\\n'
             '  `lopit|compartment|ortholopit|tagm` (%d).\\n'
             '* Folds: 5, grouped by OrthoMCL orthogroup so paralogs never split across folds, seed 0,\\n'
             '  assigned once over all labelled genes and inherited by every question.\\n'
             '* Model: `HistGradientBoostingClassifier`, missing values passed through with no imputation,\\n'
             '  300 rounds at learning rate 0.05, 15 leaves, early stopping OFF, identical for every fold\\n'
             '  and question. No hyperparameter search; nothing tuned on a held-out fold.\\n'
             '* Two training variants per question: `within` (fitted on that question\\'s genes only) and\\n'
             '  `all-labels` (fitted on all %d map labels, scored only on that question\\'s genes). Both are\\n'
             '  out-of-fold; both are reported.\\n'
             '* Calibration: isotonic (`parascape.claims._isotonic`), cross-fitted across folds.\\n'
             '* Bound: exact one-sided Clopper-Pearson, `scipy.stats.beta.ppf(0.10, k, n-k+1)`.\\n'
             % (len(FEATURES), len(FOLDSEEK_3DI), len(HELD),
                len([c for c in held_family if node_family.get(c) == 'hosttx_effectors']),
                len([c for c in held_family if node_family.get(c) == 'organelle_surface_turboid']),
                len([c for c in held_family if node_family.get(c) == 'pubmed']),
                len([c for c in held_family if node_family.get(c) == 'pvm_proximity']),
                ', '.join('`%s`' % c for c in held_attention), len(held_named), len(y)))
    R.append('## Scoreboards (mean over 5 folds)\\n')
    for key in QUESTIONS:
        sub = scoreboard[scoreboard.question == key]
        R.append('### %s\\n' % QUESTIONS[key]['title'])
        R.append(fmt(sub.pivot_table(index='strategy', columns='metric', values='mean')[
            [c for c in ['accuracy', 'precision', 'recall', 'roc_auc', 'average_precision', 'calls', 'coverage']
             if c in sub.metric.unique()]].round(3)))
    R.append('## Calibration\\n')
    R.append(fmt(calibration.groupby(['question', 'model'])[['ece', 'ece_uncalibrated']].first().round(4)))
    R.append('## Top 10 features, hard question, permutation importance on held-out folds\\n')
    R.append(fmt(top10))
    R.append('## Performance by label provenance (does it do worse on the hard negatives?)\\n')
    R.append(fmt(by_prov.round(3), index=False))
    R.append('## Robustness of the hard question\\n')
    R.append(fmt(robust.round(4), index=False))
    R.append('## Files\\n')
    R.append('| File | Contents |\\n|---|---|\\n'
             '| `scoreboard.tsv` | question x strategy x metric, mean plus every fold |\\n'
             '| `thresholds_hard.tsv` | the headline: threshold, calls, correct, precision, lower 90 % bound, genes fillable, reaches 0.90 |\\n'
             '| `thresholds_reviewed.tsv` | the same on the 98 Einar-reviewed genes |\\n'
             '| `thresholds_easy.tsv` | the same on all labels (comparison only) |\\n'
             '| `calibration.tsv` | reliability table and ECE per question and training variant |\\n'
             '| `importance.tsv` | permutation importance per feature with its provenance family, hard and easy questions |\\n'
             '| `by_label_provenance.tsv` | AUC / FPR for each positive provenance against each negative provenance |\\n'
             '| `by_label_provenance_errors.tsv` | per-provenance median score and error rate, every question |\\n'
             '| `robustness_hard.tsv` | the hard-question bound with each worrying feature block removed |\\n'
             '| `precision_at_k.tsv` | precision and lower bound of the top k genes by raw out-of-fold score, every k from 10 up |\\n'
             '| `leakage_attribution.tsv` | the same questions with the declared dependencies put back - a leakage demonstration, not a result |\\n'
             '| `predictions.tsv` | all 8,140 genes: calibrated probability per question model, out-of-fold probability, label, provenance, secretory flag, per-gene occlusion attributions |\\n')
    R.append('\\n## Limitations (not softened)\\n')
    R.append('1. **The hard set is %d genes, %d positives.** ~%d positives per fold. The Clopper-Pearson\\n'
             '   bound cannot reach 0.90 on fewer than ~22 calls even with perfect precision, so over most of\\n'
             '   the threshold range the bar is unreachable as arithmetic, before any question of model\\n'
             '   quality. This is a data-quantity result before it is a modelling result.\\n'
             % (len(QUESTIONS['hard']['y']), int(QUESTIONS['hard']['y'].sum()),
                int(QUESTIONS['hard']['y'].sum()) // 5))
    R.append('2. **Two of the hard negatives are contested** (GRA9, NTPase II), and the contested direction\\n'
             '   is the one that decides the question. Flip test in the notebook.\\n'
             '3. **The 43 hard negatives are a curation decision about published papers, not an export\\n'
             '   assay.** Usually immunogold EM or IFA in one strain at one time point. A protein vacuolar in\\n'
             '   a tachyzoite and exported in a bradyzoite would be invisible to these labels.\\n'
             '4. **DeepTMHMM misses signal peptides** (GRA16 / TGME49_208830 is called GLOBULAR). The\\n'
             '   secretory-marker arm recovers some, but a protein that is neither predicted nor a hyperLOPIT\\n'
             '   marker is absent from the hard question entirely. SignalP would fix it; open gap.\\n'
             '5. **The positives are still mostly one assay** read by one human: the PVM proximity screen\\n'
             '   (PMID 34749525) plus the literature.\\n'
             '6. **Publication bias is reduced, not removed.** The declared dependencies are held out, but a\\n'
             '   famous gene is also the gene with a melting temperature, a palmitome entry and a deposited\\n'
             '   structure. Missingness itself carries fame.\\n'
             '7. **`foldseek_3di` is a crude shape summary**, 637 genes have no 3Di and 520 have one copied\\n'
             '   by alignment; structure availability is fame-correlated.\\n'
             '8. **The all-gene probability columns are not independently validated** (isotonic curve fitted\\n'
             '   in-sample on pooled out-of-fold scores). Only the cross-fitted out-of-fold numbers are honest.\\n'
             '9. **The ortholog-copy baseline has coverage 0** under orthogroup-grouped folds, by construction.\\n'
             '10. **One model, fixed settings.** Stages 2-5 of the instruction 72 plan are not in this notebook.\\n'
             '11. **Permutation importance is noisy at this size** (~20 held-out genes per fold in the hard\\n'
             '    question); signs flip between folds. Read families, not features.\\n')
    (OUT / 'README.md').write_text('\\n'.join(R))
    print((OUT / 'README.md').read_text()[:2500])
""")

nb.write(str(ROOT / 'notebooks/host_facing_inference_v2_2026_11_04.ipynb'))
print('wrote notebook and', ROOT / 'results/host_facing_inference_v2_2026_11_04')
