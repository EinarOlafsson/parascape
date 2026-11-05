"""Instruction 72, target trait 3: binding partners. First scoring run.

Nothing has been scored for the partner trait before this. The lessons from the two host-facing runs
and the domain columns shape every choice here:

* a wrong or easy LABEL set flatters a model far more than a leaky feature does, so the negative set
  is built explicitly and the easy version is reported beside the hard one;
* every feature family is reported against the HARD comparison (cold-node folds, degree-matched
  decoys), never against the flattering one;
* "cannot be scored yet" is a valid answer, and the Clopper-Pearson lower 90 % bound, not the point
  estimate, decides whether anything may be filled.

Positives are the 2,843 StarPath XL-MS crosslink pairs (probability = 1 - local FDR): direct physical
evidence. The 33,646 bait-member pull-down pairs in the information map sit at a provisional 0.5
because no Toxoplasma background frequency background exists; they are scored as a separate, clearly weaker set
and never mixed into the headline positives.

Writes notebooks/partner_inference_2026_11_04.ipynb (executed) and results/partner_inference_2026_11_04/.
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

nb = ExecutedNotebook('Binding-partner inference: cold-node folds and degree-matched decoys')


def code(src: str):
    return nb.code(dedent(src).strip('\n'))


# --------------------------------------------------------------------------- intro
nb.md(
    '# Binding partners — the third target trait, scored for the first time',
    'Instruction 72 (Einar, 2026-11-04). Trait 1 (host-facing) has been scored twice and fills '
    'nothing; trait 3 (function) has a truth set but too few labels per class. **Nothing had been '
    'scored for binding partners at all.** This notebook builds the pair-level question and answers '
    'it, or says why it cannot be answered.',
    '**The fill bar is unchanged:** a predicted pair may be written as a fact only at a calibrated '
    'probability whose held-out precision has a **Clopper-Pearson lower 90 % confidence bound** of '
    '0.90 or better. Not the point estimate.',
    '**Three lessons from the earlier runs are built into the design, not bolted on.**',
    '1. *A wrong or easy label set flatters a model more than a leaky feature does* (the correction '
    'in instruction 72: run 1 of host-facing passed the bar because 22 vacuole-resident proteins had '
    'been counted as positives, worth +0.22, while the near-circular features were worth +0.013). '
    'For pairs there is no verified negative at all, so the negative construction is the single most '
    'consequential decision in this notebook and it is made in the open, three ways.\n'
    '2. *Every feature is reported against the hard comparison* (the `domain_novel` table: odds '
    'ratio 15.1 against all negatives, 0.45 against the negatives that answer the question). Here '
    'the hard comparison is **cold-node folds with degree-matched decoys**; the easy one (random '
    'pairs, known-node folds) is reported beside it so the gap is visible.\n'
    '3. *"Cannot be scored yet" is a result.*',
    '**Positives: direct physical evidence only.** The 2,843 StarPath crosslink pairs '
    '(`parascape/data/crosslink_models.parquet`), each with its own local FDR, so the positive '
    'probability is `1 - local_fdr` rather than a tier default.',
    '**Why the pull-down pairs are not in the headline positives.** '
    '`results/information_map_v0_2026_11_04/pairs.parquet` also holds 33,646 bait-member pairs '
    'parsed from IP-MS and BioID studies, all at a **provisional 0.5**. Without a background '
    'frequency table — a Toxoplasma background frequency — **a sticky contaminant is indistinguishable from a '
    'real partner**: ribosomal proteins, chaperones, tubulin and the proteasome appear in nearly '
    'every pull-down, and nothing in the parsed table says whether a protein was found with this '
    'bait or with all baits. Training on them would teach the model "abundant and sticky", and that '
    'is the pair analogue of the easy contrast. They are scored separately in §12.',
    'A reader who wants only the answer should read §8 (the headline scoreboard), §10 (the bound and '
    'what would fill), §11 (known complexes) and §14 (verdict and limitations).')

# --------------------------------------------------------------------------- 1 data
nb.md(
    '## 1 · The pair universe and the positives',
    'All identifiers are ME49 gene IDs over the 8,140 rows of `nodes.parquet`. The crosslink table '
    'is keyed by StarPath pair id; both endpoints of all 2,843 pairs are in the node table, so '
    'nothing is lost to mapping.')
code("""
    import os
    # One thread per BLAS/OpenMP pool: the tables are small and the machine is shared.
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[v] = '1'
    import json, warnings
    from pathlib import Path
    import numpy as np, pandas as pd
    from scipy.stats import beta
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    ROOT = Path('.').resolve()
    MAP = ROOT / 'results/information_map_v0_2026_11_04'
    OUT = ROOT / 'results/partner_inference_2026_11_04'; OUT.mkdir(parents=True, exist_ok=True)
    SEED = 0

    nodes = pd.read_parquet(ROOT / 'parascape/data/nodes.parquet').set_index('gene_id')
    xl = pd.read_parquet(ROOT / 'parascape/data/crosslink_models.parquet')
    pairs_map = pd.read_parquet(MAP / 'pairs.parquet')
    values = pd.read_parquet(MAP / 'values.parquet')
    seqt = pd.read_parquet(ROOT / 'results/sequence_traits_2026_11_04/sequence_traits.parquet')
    domt = pd.read_parquet(ROOT / 'results/domain_traits_2026_11_04/domain_traits.parquet')
    coex = pd.read_parquet(ROOT / 'results/coexpression_localization_2026_11_04/coexpression_neighbours.parquet')

    GENES = nodes.index
    gi = pd.Series(np.arange(len(GENES)), index=GENES)          # gene -> row number

    def key(a, b):
        a = np.asarray(a, dtype=object); b = np.asarray(b, dtype=object)
        return np.where(a < b, a + '|' + b, b + '|' + a)

    pos = xl[['gene_a', 'gene_b', 'local_fdr', 'n_crosslinks', 'identification_score']].copy()
    pos['probability'] = 1.0 - pos['local_fdr']
    pos['pkey'] = key(pos.gene_a, pos.gene_b)
    in_nodes = pos.gene_a.isin(GENES) & pos.gene_b.isin(GENES)
    print('crosslink pairs: %d, both endpoints in nodes.parquet: %d' % (len(pos), int(in_nodes.sum())))
    pos = pos[in_nodes & (pos.gene_a != pos.gene_b)].drop_duplicates('pkey').reset_index(drop=True)
    xl_nodes = pd.Index(sorted(set(pos.gene_a) | set(pos.gene_b)))
    print('unique undirected positives: %d over %d proteins' % (len(pos), len(xl_nodes)))
    print('positive probability (1 - local FDR): median %.4f, min %.4f, %d pairs above 0.99'
          % (pos.probability.median(), pos.probability.min(), int((pos.probability > 0.99).sum())))
    pos[['gene_a', 'gene_b', 'probability', 'n_crosslinks']].head(8)
""")

nb.md(
    'The pull-down and curated host pairs, for the record and for the exclusion set below. Every '
    'pair that has *any* recorded evidence of interaction — crosslink, pull-down, BioID or a curated '
    'host interaction — is excluded from every decoy set, because a decoy must be a pair for which '
    'nothing has been reported, not merely a pair the crosslinking missed.')
code("""
    pairs_map['pkey'] = key(pairs_map.gene_a.astype(str), pairs_map.gene_b.astype(str))
    print(pairs_map.groupby(['source_tier', 'partner_side']).size().to_string())
    pull = pairs_map[(pairs_map.source_tier == 'observed-screen') & (pairs_map.partner_side == 'parasite')].copy()
    pull = pull[pull.gene_a.isin(GENES) & pull.gene_b.isin(GENES) & (pull.gene_a != pull.gene_b)]
    pull = pull.drop_duplicates('pkey').reset_index(drop=True)
    ANY_EVIDENCE = set(pairs_map.pkey) | set(pos.pkey)
    print('pull-down (bait-member, provisional 0.5) pairs usable: %d over %d proteins'
          % (len(pull), len(set(pull.gene_a) | set(pull.gene_b))))
    print('pairs excluded from all decoy sets (any recorded evidence): %d' % len(ANY_EVIDENCE))
    print('overlap crosslink x pull-down: %d pairs' % len(set(pos.pkey) & set(pull.pkey)))
""")

# --------------------------------------------------------------------------- 2 negatives
nb.md(
    '## 2 · The negatives, which are the hard problem',
    '**There is no experimentally verified "these two proteins do not interact" set for Toxoplasma.** '
    'Nothing in the data says a pair was tested and found not to interact. So every negative here is '
    'a *decoy*: a pair for which no interaction has been reported. That is an assumption, it is '
    'wrong for some unknown fraction of the decoys, and it is the main reason the numbers below are '
    'an upper bound on honesty rather than a measurement of truth.',
    'Three constructions are built and all three are scored. The headline is the hardest defensible '
    'one.',
    '| set | construction | what it removes | what it still allows |\n|---|---|---|---|\n'
    '| `random` *(easy, for comparison)* | both endpoints drawn uniformly from all 8,140 genes | '
    'nothing | abundance, popularity, "was this protein ever seen by mass spectrometry", compartment '
    '— most random pairs differ in everything, so the task is easy |\n'
    '| **`degree_matched` (headline)** | for each positive, one endpoint is kept and the other is '
    'replaced by a protein with a similar number of recorded crosslink partners, drawn from the same '
    'cold-node block | the **popular-protein shortcut** and, because the pool is the crosslink-'
    'detected proteome, the "detectable by XL-MS at all" shortcut (abundance, lysine content, '
    'solubility) | genuine pair-level signal only |\n'
    '| `same_compartment` | both endpoints are hyperLOPIT/TAGM residents of the **same** compartment, '
    'no recorded evidence | the compartment shortcut | popularity |',
    'Why degree matching is the headline: the crosslink graph is extremely uneven — the median '
    'detected protein has 1 partner and the 99th percentile has 28 — so "both proteins are popular" '
    'alone separates positives from random pairs. That is exactly the pair analogue of '
    '"has a signal peptide" in the host-facing run, and the degree-product baseline in §7 measures '
    'how much it is worth in each set.',
    '*A cost of degree matching, stated rather than hidden:* matching on crosslink degree uses the '
    'label source to build the controls. It cannot leak a label into a feature (the degree columns '
    'are held out, §4), but it does mean the decoys are drawn from the 1,630 crosslink-detected '
    'proteins, so the model is never asked about the 6,510 proteins XL-MS never saw. The result '
    'therefore applies to pairs of detectable proteins, which is also the only population where a '
    'prediction could be checked by the same assay.')
code("""
    import collections
    deg = collections.Counter()
    for a, b in zip(pos.gene_a, pos.gene_b):
        deg[a] += 1; deg[b] += 1
    degree = pd.Series(deg, dtype=float).reindex(xl_nodes).fillna(0.0)
    print('crosslink degree over the %d detected proteins: median %.0f, mean %.1f, p90 %.0f, max %.0f'
          % (len(degree), degree.median(), degree.mean(), degree.quantile(0.9), degree.max()))
    pd.DataFrame({'degree': degree.describe()}).T
""")

# --------------------------------------------------------------------------- 3 folds
nb.md(
    '## 3 · Folds: cold-node splits, so a partner is never predicted from a protein seen in training',
    'Instruction 71: *"Pairs (trait 2) need pair-level folds: hold out whole proteins (cold-node '
    'split, `splits.cold_node_pairs`) so a partner is never predicted from the same protein seen in '
    'training."* Two regimes are scored.',
    '* **cold-node (headline)** — whole proteins are held out. Five blocks of proteins, grouped by '
    '**orthogroup** so paralogs cannot sit on both sides, seed 0. For block *f*, the test pairs are '
    'the pairs with **both** endpoints in *f* and the training pairs are those with **neither** '
    'endpoint in *f*. Pairs with exactly one endpoint in the held-out block are dropped for that '
    'fold — that is `splits.cold_node_pairs`\' own rule, and it is why cold-node scoring keeps only '
    'about a fifth of the pairs: the question is genuinely harder and also smaller.\n'
    '* **known-node (easier, reported beside it)** — `splits.known_node_pairs`, which first reserves '
    'a set of training edges covering every protein and then partitions the remaining edges, so '
    'every held-out endpoint has been seen in training with some other partner. This is edge '
    'completion, not generalisation to new proteins.',
    'The manifests are the library ones (`parascape.splits`), not hand-rolled: `make_split`, '
    '`SplitManifest`, `cold_node_pairs`, `known_node_pairs` and `pair_id`. The notebook asserts that '
    'the test set from the library equals the "both endpoints in the held-out block" set and that '
    'pair it drops really does cross a partition.')
code("""
    from parascape import splits as S
    from parascape.splits import Assignment, SplitManifest

    # Blocks are defined over ALL 8,140 proteins, because the random and same-compartment decoy sets
    # reach outside the crosslink-detected proteome and every endpoint must have a block.
    og = nodes['orthogroup'].reindex(GENES).astype(str)
    rng = np.random.default_rng(SEED)
    groups = sorted(og.unique())
    perm = rng.permutation(len(groups))
    block_of_group = {groups[int(j)]: int(i % 5) for i, j in enumerate(perm)}
    node_block = og.map(block_of_group)
    print('cold-node blocks (all proteins per block):', node_block.value_counts().sort_index().tolist())
    print('crosslink-detected proteins per block:',
          node_block.reindex(xl_nodes).value_counts().sort_index().tolist())
    print('orthogroups: %d over %d proteins; no orthogroup crosses a block by construction'
          % (len(groups), len(GENES)))

    # One library SplitManifest per fold: the held-out block is `test`, the other four blocks fill
    # train/tune/calibration. Protected group = orthogroup, so no paralog crosses the boundary.
    ROLE_ROTATION = {0: 'calibration', 1: 'tune', 2: 'train', 3: 'train'}

    def node_split_for(fold):
        rows = []
        for gene, blk in node_block.items():
            role = 'test' if blk == fold else ROLE_ROTATION[(blk - fold) % 5 - 1]
            rows.append(Assignment(gene, og[gene], role))
        return SplitManifest('Tg', 'partner_pairs_cold_node', 'homology', tuple(rows), SEED, 'transductive')

    node_splits = {f: node_split_for(f) for f in range(5)}
    print('manifest identities:', {f: m.identity[:10] for f, m in node_splits.items()})
""")

nb.md(
    'Decoys are drawn **inside the fold structure**: for each positive, the endpoint that is '
    'resampled is replaced by a protein from the *same* cold-node block, so a decoy sits in exactly '
    'the same fold as the positive it mirrors and every fold stays balanced at 1:1. Which endpoint '
    'is resampled is a seeded coin, so the construction is not systematically biased to one side.')
code("""
    logdeg = np.log1p(degree)
    xl_block = node_block.reindex(xl_nodes)
    by_block = {f: xl_block.index[xl_block.values == f] for f in range(5)}

    def degree_matched_decoys(positives, tol=0.25, seed=SEED):
        r = np.random.default_rng(seed)
        rows, taken = [], set(ANY_EVIDENCE)
        for a, b in zip(positives.gene_a.values, positives.gene_b.values):
            keep, swap = (a, b) if r.random() < 0.5 else (b, a)
            blk = int(node_block[swap])
            pool = by_block[blk]
            cand = []
            d = abs(logdeg[pool].values - logdeg[swap])
            width = tol
            for _ in range(6):
                cand = pool[(d <= width)]
                cand = [c for c in cand if c != keep and key([c], [keep])[0] not in taken]
                if cand:
                    break
                width *= 2
            if not cand:                      # nearest degree in the block, ignoring tolerance
                order = np.argsort(d, kind='stable')
                cand = [pool[int(j)] for j in order if pool[int(j)] != keep
                        and key([pool[int(j)]], [keep])[0] not in taken][:1]
            if not cand:
                continue
            pick = cand[int(r.integers(len(cand)))]
            taken.add(key([pick], [keep])[0])
            rows.append(dict(gene_a=keep, gene_b=pick, mirrored_a=a, mirrored_b=b,
                             swapped_degree=float(degree[swap]), decoy_degree=float(degree[pick])))
        out = pd.DataFrame(rows)
        out['pkey'] = key(out.gene_a, out.gene_b)
        return out.drop_duplicates('pkey').reset_index(drop=True)

    dec_degree = degree_matched_decoys(pos)
    print('degree-matched decoys: %d for %d positives' % (len(dec_degree), len(pos)))
    md = pd.DataFrame({'replaced protein degree': dec_degree.swapped_degree.describe(),
                       'decoy protein degree': dec_degree.decoy_degree.describe()})
    print('median |log1p degree difference| between the replaced protein and its decoy: %.3f'
          % np.median(np.abs(np.log1p(dec_degree.swapped_degree) - np.log1p(dec_degree.decoy_degree))))
    md.T
""")
code("""
    def random_decoys(n, seed=SEED + 1):
        r = np.random.default_rng(seed)
        rows, taken = [], set(ANY_EVIDENCE)
        allg = GENES.values
        while len(rows) < n:
            a, b = allg[r.integers(len(allg))], allg[r.integers(len(allg))]
            if a == b:
                continue
            k = key([a], [b])[0]
            if k in taken:
                continue
            taken.add(k); rows.append((a, b))
        out = pd.DataFrame(rows, columns=['gene_a', 'gene_b'])
        out['pkey'] = key(out.gene_a, out.gene_b)
        return out

    comp_obs = values['compartment_hyperlopit'].astype('object')
    comp_pred = nodes['compartment'].astype('object').replace({'unassigned': np.nan})
    comp_any = comp_obs.reindex(GENES).fillna(comp_pred.reindex(GENES))

    def same_compartment_decoys(n, seed=SEED + 2):
        r = np.random.default_rng(seed)
        rows, taken = [], set(ANY_EVIDENCE)
        pools = {c: g.index.values for c, g in comp_any.dropna().to_frame('c').groupby('c')
                 if len(g) >= 4}
        names = sorted(pools)
        weights = np.array([len(pools[c]) for c in names], dtype=float); weights /= weights.sum()
        while len(rows) < n:
            c = names[int(r.choice(len(names), p=weights))]
            pool = pools[c]
            a, b = pool[r.integers(len(pool))], pool[r.integers(len(pool))]
            if a == b:
                continue
            k = key([a], [b])[0]
            if k in taken:
                continue
            taken.add(k); rows.append((a, b, c))
        out = pd.DataFrame(rows, columns=['gene_a', 'gene_b', 'compartment'])
        out['pkey'] = key(out.gene_a, out.gene_b)
        return out

    dec_random = random_decoys(len(pos))
    dec_compart = same_compartment_decoys(len(pos))
    print('random decoys %d; same-compartment decoys %d (compartments used: %d)'
          % (len(dec_random), len(dec_compart), dec_compart.compartment.nunique()))
    print('random-decoy endpoints that were ever seen by XL-MS: %.1f%% of proteins'
          % (100 * pd.Index(np.r_[dec_random.gene_a.values, dec_random.gene_b.values]).isin(xl_nodes).mean()))
    print('degree-matched decoy endpoints seen by XL-MS: 100 % by construction')
""")

# --------------------------------------------------------------------------- 4 features
nb.md(
    '## 4 · Features, one family at a time, and what is held out',
    'Every feature is a property of the **pair**, computed the same way for a positive and for a '
    'decoy, and every one belongs to exactly one family so the family can be held out or permuted '
    'as a block.',
    '| family | features | source |\n|---|---|---|\n'
    '| `coexpression` | the combined percentile of the pair, study support, whether the neighbour '
    'relation is mutual | `results/coexpression_localization_2026_11_04/coexpression_neighbours.parquet` '
    '(top 50 neighbours per gene, so a pair outside both top-50 lists is missing, not zero) |\n'
    '| `fitness` | for each CRISPR fitness column, the **absolute difference** and the **minimum** of '
    'the pair; plus the correlation of the two fitness profiles across columns | the 28 `fit_*` '
    'columns of `nodes.parquet` |\n'
    '| `domain` | number of shared InterPro/Pfam families, Jaccard, min/max domain count, both '
    'domain-free | `results/domain_traits_2026_11_04/domain_traits.parquet` |\n'
    '| `compartment_observed` | same hyperLOPIT marker compartment, same coarse compartment, same '
    'experimental GO cellular component, both assigned | `values.parquet` (observed facts only) |\n'
    '| `compartment_predicted` | the same, from the TAGM/unified `compartment` column | '
    '`nodes.parquet` — **prediction-grade**, kept in a family of its own and reported separately |\n'
    '| `host_facing` | both host-facing, either host-facing, the two labels agree | `values.parquet` |\n'
    '| `sequence` | both have a signal peptide and whether they agree, transmembrane counts, length '
    'difference and minimum, same DeepTMHMM N-terminal side | '
    '`results/sequence_traits_2026_11_04/sequence_traits.parquet` |\n'
    '| `orthogroup` | same orthogroup, both have a Plasmodium ortholog | `nodes.parquet` |\n'
    '| `structure_3di` | cosine similarity of 3Di 3-mer profiles, plus the alignment-quality columns '
    'that instruction 72 requires to travel with any use of 3Di | `sequence_traits.parquet` |',
    '**Held out, because it is derived from the crosslinks themselves — the label source.** '
    '`n_crosslink_partners` and `n_xlink_partners` in `nodes.parquet`; every column of '
    '`crosslink_models.parquet` beyond the label and its FDR (`frac_satisfied`, `min_ca`, '
    '`chai_iptm`, `chai_plddt`, `n_models`, `model_trustworthy`, `identification_score` — these are '
    'the crosslink-derived structure models, and they exist only for positives, so they would be '
    'perfect separators of positives from decoys); and `loc_a` / `loc_b`, the StarPath '
    'localisation annotation, which likewise exists only for crosslinked pairs. Same-compartment '
    'features are therefore taken from the information map and `nodes.parquet`, never from the '
    'crosslink table.',
    '**Also held out of the headline model, for the same reason as the host-facing attention '
    'columns:** `n_ipms_partners`, `n_ipms_studies`, `n_publications`, `lit_tier`, `attention_depth` '
    'and the proteome abundance columns. A pull-down partner count is a partner count from another '
    'assay — the popularity shortcut in a different coat — and a publication count is the '
    'publication-bias shortcut. They are used only as explicit baselines (§7) and in one labelled '
    'feature-set variant.')
code("""
    FIT = [c for c in nodes.columns if c.startswith('fit_') and pd.api.types.is_numeric_dtype(nodes[c])]
    HELD_OUT = {
        'n_crosslink_partners': 'crosslink-derived partner count (label source)',
        'n_xlink_partners': 'crosslink-derived partner count (label source)',
        'identification_score': 'crosslink identification score (exists only for positives)',
        'n_crosslinks': 'number of crosslinks (exists only for positives)',
        'frac_satisfied': 'crosslink-derived structure model', 'min_ca': 'crosslink-derived structure model',
        'chai_iptm': 'crosslink-derived structure model', 'chai_plddt': 'crosslink-derived structure model',
        'n_models': 'crosslink-derived structure model', 'model_trustworthy': 'crosslink-derived structure model',
        'loc_a': 'StarPath localisation annotation, only for crosslinked pairs',
        'loc_b': 'StarPath localisation annotation, only for crosslinked pairs',
        'n_ipms_partners': 'pull-down partner count = the popularity shortcut from another assay',
        'n_ipms_studies': 'pull-down study count = popularity', 'n_publications': 'publication-bias shortcut',
        'lit_tier': 'publication-bias shortcut', 'attention_depth': 'publication-bias shortcut',
        'n_papers_focal': 'publication-bias shortcut', 'n_papers_substantive': 'publication-bias shortcut',
        'n_papers_incidental': 'publication-bias shortcut',
    }
    pd.Series(HELD_OUT, name='reason').to_frame().rename_axis('held-out column')
""")
code("""
    # ---- per-gene lookups, all aligned to GENES -------------------------------------------------
    Z = nodes[FIT].reindex(GENES).astype(float)
    Z = (Z - Z.mean()) / Z.std(ddof=0)                  # z-score per column across all genes
    FITM = Z.values
    DOMCOLS = [c for c in domt.columns if c.startswith('dom_')]
    DOM = domt[DOMCOLS].reindex(GENES).fillna(False).values.astype(bool)
    DOMN = DOM.sum(1)

    def coarse(x):
        if not isinstance(x, str):
            return np.nan
        return x.split(' - ')[0].strip().rstrip('0123456789 ').strip()

    g_comp_obs = comp_obs.reindex(GENES).values
    g_comp_obs_coarse = np.array([coarse(x) for x in g_comp_obs], dtype=object)
    g_comp_pred = comp_pred.reindex(GENES).values
    g_comp_pred_coarse = np.array([coarse(x) for x in g_comp_pred], dtype=object)
    g_go_cc = values['compartment_go_cc'].reindex(GENES).astype('object').values
    g_host = values['host_facing'].reindex(GENES).astype('Float64').astype(float).values
    g_sp = seqt['has_signal_peptide'].reindex(GENES).astype(float).values
    g_ntm = seqt['n_tm'].reindex(GENES).astype(float).values
    g_len = seqt['seq_len'].reindex(GENES).astype(float).values
    g_nside = seqt['nterm_side'].reindex(GENES).astype('object').values
    g_og = nodes['orthogroup'].reindex(GENES).astype('object').values
    g_pf = nodes['has_pf_ortholog'].reindex(GENES).astype(float).values
    g_tdi_id = seqt['tdi_align_identity'].reindex(GENES).astype(float).values
    g_tdi_cov = seqt['tdi_coverage'].reindex(GENES).astype(float).values

    # 3Di 3-mer profiles, L2-normalised, as a sparse matrix: cosine = one row dot product per pair.
    from scipy import sparse
    tdi = seqt['tdi_full'].reindex(GENES).fillna('').astype(str).values
    vocab, rows, cols, data = {}, [], [], []
    for r, s in enumerate(tdi):
        if len(s) < 3:
            continue
        cnt = collections.Counter(s[i:i + 3] for i in range(len(s) - 2) if 'X' not in s[i:i + 3])
        if not cnt:
            continue
        for kmer, c in cnt.items():
            j = vocab.setdefault(kmer, len(vocab))
            rows.append(r); cols.append(j); data.append(float(c))
    K3 = sparse.csr_matrix((data, (rows, cols)), shape=(len(GENES), len(vocab)))
    norm = np.sqrt(K3.multiply(K3).sum(1)).A.ravel(); norm[norm == 0] = 1.0
    K3 = sparse.diags(1.0 / norm) @ K3
    g_tdi_ok = np.asarray(K3.getnnz(axis=1) > 0)
    print('3Di 3-mer profiles for %d of %d genes, %d distinct 3-mers' % (g_tdi_ok.sum(), len(GENES), len(vocab)))

    # co-expression: a lookup on the undirected key, plus how many directions recorded it
    coex = coex[coex.gene.isin(GENES) & coex.partner.isin(GENES)].copy()
    coex['pkey'] = key(coex.gene.values, coex.partner.values)
    cox = coex.groupby('pkey').agg(coexpr_percentile=('combined_percentile', 'max'),
                                   coexpr_study_support=('study_support', 'max'),
                                   coexpr_mutual=('gene', 'size'))
    cox['coexpr_mutual'] = (cox['coexpr_mutual'] >= 2).astype(float)
    print('co-expression neighbour pairs (undirected): %d' % len(cox))
""")
code("""
    FAMILIES = {}

    def register(family, names):
        for n in names:
            FAMILIES[n] = family

    register('coexpression', ['coexpr_percentile', 'coexpr_study_support', 'coexpr_mutual', 'coexpr_present'])
    register('fitness', [f'fitdiff_{c}' for c in FIT] + [f'fitmin_{c}' for c in FIT]
             + ['fit_profile_corr', 'fit_n_shared', 'fit_mean_absdiff', 'fit_max_absdiff'])
    register('domain', ['dom_shared', 'dom_jaccard', 'dom_min_count', 'dom_max_count', 'dom_both_none'])
    register('compartment_observed', ['same_comp_obs', 'same_comp_obs_coarse', 'both_comp_obs', 'same_go_cc'])
    register('compartment_predicted', ['same_comp_pred', 'same_comp_pred_coarse', 'both_comp_pred'])
    register('host_facing', ['both_host_facing', 'any_host_facing', 'host_facing_agree'])
    register('sequence', ['both_signal_peptide', 'signal_peptide_agree', 'len_absdiff', 'len_min',
                          'ntm_min', 'ntm_absdiff', 'both_tm', 'same_nterm_side'])
    register('orthogroup', ['same_orthogroup', 'both_pf_ortholog'])
    register('structure_3di', ['tdi_kmer_cosine', 'tdi_both_measured', 'tdi_min_identity', 'tdi_min_coverage'])
    FEATURES = list(FAMILIES)

    def _eq(a, b):
        ok = pd.notna(a) & pd.notna(b)
        out = np.where(ok, (a == b).astype(float), np.nan)
        return out, ok.astype(float)

    def pair_features(ga, gb, chunk=40000):
        ga, gb = np.asarray(ga, dtype=object), np.asarray(gb, dtype=object)
        out = []
        for start in range(0, len(ga), chunk):
            A, B = ga[start:start + chunk], gb[start:start + chunk]
            ia, ib = gi[A].values, gi[B].values
            d = {}
            fa, fb = FITM[ia], FITM[ib]
            both = np.isfinite(fa) & np.isfinite(fb)
            absd = np.abs(fa - fb)
            for j, c in enumerate(FIT):
                d[f'fitdiff_{c}'] = np.where(both[:, j], absd[:, j], np.nan)
                d[f'fitmin_{c}'] = np.where(both[:, j], np.minimum(fa[:, j], fb[:, j]), np.nan)
            n_shared = both.sum(1)
            d['fit_n_shared'] = n_shared.astype(float)
            with np.errstate(all='ignore'):
                d['fit_mean_absdiff'] = np.where(n_shared > 0, np.nansum(np.where(both, absd, 0), 1) / np.maximum(n_shared, 1), np.nan)
                d['fit_max_absdiff'] = np.where(n_shared > 0, np.nanmax(np.where(both, absd, -np.inf), 1), np.nan)
                xa = np.where(both, fa, np.nan); xb = np.where(both, fb, np.nan)
                ma, mb = np.nanmean(xa, 1, keepdims=True), np.nanmean(xb, 1, keepdims=True)
                ca, cbv = xa - ma, xb - mb
                num = np.nansum(ca * cbv, 1)
                den = np.sqrt(np.nansum(ca ** 2, 1) * np.nansum(cbv ** 2, 1))
                d['fit_profile_corr'] = np.where((n_shared >= 5) & (den > 0), num / np.where(den == 0, 1, den), np.nan)
            da, db = DOM[ia], DOM[ib]
            shared = (da & db).sum(1).astype(float)
            union = (da | db).sum(1).astype(float)
            d['dom_shared'] = shared
            d['dom_jaccard'] = np.where(union > 0, shared / np.where(union == 0, 1, union), np.nan)
            d['dom_min_count'] = np.minimum(DOMN[ia], DOMN[ib]).astype(float)
            d['dom_max_count'] = np.maximum(DOMN[ia], DOMN[ib]).astype(float)
            d['dom_both_none'] = ((DOMN[ia] == 0) & (DOMN[ib] == 0)).astype(float)
            d['same_comp_obs'], d['both_comp_obs'] = _eq(g_comp_obs[ia], g_comp_obs[ib])
            d['same_comp_obs_coarse'], _ = _eq(g_comp_obs_coarse[ia], g_comp_obs_coarse[ib])
            d['same_comp_pred'], d['both_comp_pred'] = _eq(g_comp_pred[ia], g_comp_pred[ib])
            d['same_comp_pred_coarse'], _ = _eq(g_comp_pred_coarse[ia], g_comp_pred_coarse[ib])
            d['same_go_cc'], _ = _eq(g_go_cc[ia], g_go_cc[ib])
            ha, hb = g_host[ia], g_host[ib]
            okh = np.isfinite(ha) & np.isfinite(hb)
            d['both_host_facing'] = np.where(okh, np.fmin(ha, hb), np.nan)
            d['any_host_facing'] = np.fmax(ha, hb)
            d['host_facing_agree'] = np.where(okh, (ha == hb).astype(float), np.nan)
            sa, sb = g_sp[ia], g_sp[ib]
            d['both_signal_peptide'] = sa * sb
            d['signal_peptide_agree'] = (sa == sb).astype(float)
            d['len_absdiff'] = np.abs(g_len[ia] - g_len[ib])
            d['len_min'] = np.minimum(g_len[ia], g_len[ib])
            d['ntm_min'] = np.minimum(g_ntm[ia], g_ntm[ib])
            d['ntm_absdiff'] = np.abs(g_ntm[ia] - g_ntm[ib])
            d['both_tm'] = ((g_ntm[ia] > 0) & (g_ntm[ib] > 0)).astype(float)
            d['same_nterm_side'], _ = _eq(g_nside[ia], g_nside[ib])
            d['same_orthogroup'], _ = _eq(g_og[ia], g_og[ib])
            d['both_pf_ortholog'] = g_pf[ia] * g_pf[ib]
            cos = np.asarray(K3[ia].multiply(K3[ib]).sum(1)).ravel()
            okt = g_tdi_ok[ia] & g_tdi_ok[ib]
            d['tdi_kmer_cosine'] = np.where(okt, cos, np.nan)
            d['tdi_both_measured'] = okt.astype(float)
            d['tdi_min_identity'] = np.fmin(g_tdi_id[ia], g_tdi_id[ib])
            d['tdi_min_coverage'] = np.fmin(g_tdi_cov[ia], g_tdi_cov[ib])
            frame = pd.DataFrame(d, index=key(A, B))
            cj = cox.reindex(frame.index)
            frame['coexpr_percentile'] = cj.coexpr_percentile.values
            frame['coexpr_study_support'] = cj.coexpr_study_support.values
            frame['coexpr_mutual'] = cj.coexpr_mutual.values
            frame['coexpr_present'] = cj.coexpr_percentile.notna().astype(float).values
            out.append(frame[FEATURES].astype(np.float32))
        return pd.concat(out) if out else pd.DataFrame(columns=FEATURES)

    fam_counts = pd.Series(FAMILIES).value_counts().rename('features')
    print('%d pair features in %d families' % (len(FEATURES), fam_counts.size))
    fam_counts.to_frame()
""")
nb.md('Every feature with its family, written out in full, as instruction 72 requires.')
code("""
    feat_table = pd.DataFrame({'feature': FEATURES, 'family': [FAMILIES[f] for f in FEATURES]})
    feat_table.to_csv(OUT / 'features.tsv', sep='\\t', index=False)
    print(feat_table.groupby('family').feature.apply(lambda s: ', '.join(s)).to_string())
    feat_table.head(10)
""")

# --------------------------------------------------------------------------- 5 assemble
nb.md(
    '## 5 · The labelled pair tables',
    'Six labelled sets: three negative constructions x two split regimes. Feature tables are built '
    'for the labelled pairs and their decoys only — the full 8,140 x 8,140 space is 33 million pairs '
    'and is never materialised; the candidate set for filling (§10) is restricted to pairs that '
    'already have some reason to be considered.')
code("""
    DECOYS = {'degree_matched': dec_degree, 'random': dec_random, 'same_compartment': dec_compart}
    LAB = {}
    for name, dec in DECOYS.items():
        frame = pd.concat([
            pos[['gene_a', 'gene_b']].assign(y=True),
            dec[['gene_a', 'gene_b']].assign(y=False)], ignore_index=True)
        frame['pkey'] = key(frame.gene_a, frame.gene_b)
        frame = frame.drop_duplicates('pkey').set_index('pkey')
        LAB[name] = frame
        print('%-18s %d pairs (%d positive, %d decoy)'
              % (name, len(frame), int(frame.y.sum()), int((~frame.y).sum())))

    X = {}
    for name, frame in LAB.items():
        X[name] = pair_features(frame.gene_a.values, frame.gene_b.values).reindex(frame.index)
    print('feature coverage, degree-matched set (fraction of pairs with a value):')
    cov = X['degree_matched'].notna().mean().sort_values()
    print(cov.head(8).round(3).to_string())
    print('...')
    print(cov.tail(6).round(3).to_string())
    cov.describe().round(3)
""")
nb.md(
    '**Cold-node fold assignment.** The library decides which pairs are testable in each fold; the '
    'notebook checks it against the explicit rule.')
code("""
    fold_of = {}
    for name, frame in LAB.items():
        test_fold = pd.Series(np.nan, index=frame.index, dtype=float)
        dropped = {}
        for f in range(5):
            part = S.cold_node_pairs(list(zip(frame.gene_a.values, frame.gene_b.values)), node_splits[f])
            held = {key([a], [b])[0] for a, b, role in part.pairs if role == 'test'}
            expect = set(frame.index[(node_block.reindex(frame.gene_a).values == f)
                                     & (node_block.reindex(frame.gene_b).values == f)])
            assert held == expect, (name, f, len(held), len(expect))
            for a, b in part.excluded_cross_partition:
                assert (node_block[a] == f) != (node_block[b] == f) or node_block[a] != node_block[b]
            test_fold.loc[sorted(held)] = f
            dropped[f] = len(part.excluded_cross_partition)
        fold_of[name] = test_fold
        n = int(test_fold.notna().sum())
        print('%-18s cold-node testable pairs %d of %d (%.0f%%), positives %d'
              % (name, n, len(frame), 100 * n / len(frame),
                 int(frame.y[test_fold.notna()].sum())))
    cold_counts = pd.DataFrame({name: fold_of[name].value_counts().sort_index() for name in LAB})
    cold_counts.index.name = 'fold'
    cold_counts
""")
code("""
    # known-node regime: the library reserves endpoint-covering training edges, then partitions the rest.
    known = {}
    for name, frame in LAB.items():
        part = S.known_node_pairs(list(zip(frame.gene_a.values, frame.gene_b.values)),
                                  organism='Tg', benchmark_id='partner_pairs_known_node', seed=SEED)
        roles = pd.Series({key([a], [b])[0]: role for a, b, role in part.pairs})
        known[name] = roles.reindex(frame.index)
        print('%-18s known-node roles: %s' % (name, known[name].value_counts().to_dict()))
    print('known-node held-out positives, degree-matched set: %d'
          % int(LAB['degree_matched'].y[known['degree_matched'] == 'test'].sum()))
""")

# --------------------------------------------------------------------------- 6 machinery
nb.md(
    '## 6 · Machinery: model, metrics, bound, calibration',
    '`HistGradientBoostingClassifier`, missing values passed straight through with **no imputation** '
    '(a pair outside both top-50 co-expression lists has no percentile; that absence is itself '
    'informative and the tree learns a direction for it). Fixed settings for every fold and every '
    'question — 300 rounds, learning rate 0.05, 15 leaves, 20 samples per leaf, L2 1.0, early '
    'stopping off, seed 0 — so there is **no hyperparameter search and nothing was tuned on a '
    'held-out fold**.',
    'The bound is the exact one-sided **Clopper-Pearson** lower 90 % limit, '
    '`beta.ppf(0.10, k, n-k+1)`. The table printed below is the arithmetic that governs the whole '
    'notebook: below about 22 held-out calls the bound cannot reach 0.90 even if every call is right.')
code("""
    from sklearn.metrics import roc_auc_score, average_precision_score
    from sklearn.ensemble import HistGradientBoostingClassifier
    from parascape.claims import _isotonic, _ece

    def make_model():
        return HistGradientBoostingClassifier(
            max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=20,
            l2_regularization=1.0, early_stopping=False, random_state=0)

    def clopper_pearson_lower(k, n, confidence=0.90):
        if n == 0:
            return np.nan
        if k == 0:
            return 0.0
        return float(beta.ppf(1 - confidence, k, n - k + 1))

    THRESHOLDS = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97, 0.99]

    def usable(Xn, feats, rows):
        '''Drop columns with no observed value in the training rows: a wholly missing column
        carries nothing and sklearn 1.9 cannot bin one.'''
        return [c for c in feats if Xn.loc[rows, c].notna().any()]

    def cold_node_oof(name, features=None, labels=None):
        '''Out-of-fold scores under the cold-node regime: train on pairs with NEITHER endpoint in
        the held-out protein block, score the pairs with BOTH endpoints in it.'''
        frame = LAB[name] if labels is None else labels
        feats = FEATURES if features is None else features
        Xn = X[name]
        y = frame.y
        blk_a, blk_b = node_block.reindex(frame.gene_a).values, node_block.reindex(frame.gene_b).values
        oof = pd.Series(np.nan, index=frame.index, dtype=float)
        models, train_sizes = {}, {}
        for f in range(5):
            te = frame.index[(blk_a == f) & (blk_b == f)]
            tr = frame.index[(blk_a != f) & (blk_b != f)]
            if not len(te) or len(set(y[tr])) < 2:
                continue
            cols = usable(Xn, feats, tr)
            if not cols:
                continue
            m = make_model().fit(Xn.loc[tr, cols], y.loc[tr])
            oof.loc[te] = m.predict_proba(Xn.loc[te, cols])[:, 1]
            models[f], train_sizes[f] = m, cols
        return oof, models, train_sizes

    def known_node_oof(name, features=None):
        '''Known-node regime: five folds over the non-reserved edges, the reserved
        endpoint-covering edges always in training, so every held-out endpoint was seen.'''
        frame, feats, Xn = LAB[name], FEATURES if features is None else features, X[name]
        roles = known[name]
        oof = pd.Series(np.nan, index=frame.index, dtype=float)
        r = np.random.default_rng(SEED)
        movable = frame.index[roles != 'train']
        fold = pd.Series(r.integers(0, 5, len(movable)), index=movable)
        for f in range(5):
            te = fold.index[fold == f]
            tr = frame.index.difference(te)
            if not len(te) or len(set(frame.y[tr])) < 2:
                continue
            cols = usable(Xn, feats, tr)
            if not cols:
                continue
            m = make_model().fit(Xn.loc[tr, cols], frame.y.loc[tr])
            oof.loc[te] = m.predict_proba(Xn.loc[te, cols])[:, 1]
        return oof, fold.reindex(frame.index)

    def cross_fitted_isotonic(scores, y, fold):
        cal = pd.Series(np.nan, index=scores.index, dtype=float)
        idx = scores.dropna().index
        for f in sorted(fold.dropna().unique()):
            tr = [i for i in idx if fold.get(i) != f]
            te = [i for i in idx if fold.get(i) == f]
            if len(tr) < 20 or not te:
                continue
            xs, ys = _isotonic(scores.loc[tr].values, y.loc[tr].astype(float).values)
            cal.loc[te] = np.interp(scores.loc[te].values, xs, ys)
        return cal

    pd.DataFrame([dict(calls=n, lower90_if_every_call_correct=round(clopper_pearson_lower(n, n), 3))
                  for n in (4, 8, 12, 16, 20, 22, 25, 30, 50, 100, 300)]).set_index('calls').T
""")

# --------------------------------------------------------------------------- 7 baselines
nb.md(
    '## 7 · Baselines, in the same table as the model',
    '**random** — a uniform score, seed 0; its precision is the prevalence and that is the floor.\n'
    '**both proteins are popular** — the product of the two crosslink partner counts '
    '(`log1p(deg_a) * log1p(deg_b)`). This is the shortcut degree matching exists to remove, and the '
    'difference between its AUC on the random decoys and on the degree-matched decoys *is* the '
    'measurement of whether the matching worked.\n'
    '**same compartment alone** — the prediction is the same-compartment feature (observed where '
    'available, otherwise the predicted compartment).\n'
    '**shared orthogroup alone** — the prediction is `same_orthogroup`.',
    'Baselines are scored on exactly the pairs the model is scored on, fold by fold.')
code("""
    def baseline_scores(name):
        frame = LAB[name]
        a, b = frame.gene_a.values, frame.gene_b.values
        d_a = np.log1p(degree.reindex(a).fillna(0.0).values)
        d_b = np.log1p(degree.reindex(b).fillna(0.0).values)
        same_comp = X[name]['same_comp_obs'].fillna(X[name]['same_comp_pred']).fillna(0.0).values
        return {
            'baseline: random': pd.Series(np.random.default_rng(SEED).random(len(frame)), index=frame.index),
            'baseline: both popular (degree product)': pd.Series(d_a * d_b, index=frame.index),
            'baseline: same compartment alone': pd.Series(same_comp.astype(float), index=frame.index),
            'baseline: shared orthogroup alone': pd.Series(X[name]['same_orthogroup'].fillna(0.0).values.astype(float),
                                                           index=frame.index),
        }

    BASE = {name: baseline_scores(name) for name in LAB}
    for name in LAB:
        m = fold_of[name].notna()
        y = LAB[name].y[m]
        print(name)
        for bname, s in BASE[name].items():
            print('   %-42s cold-node AUC %.3f' % (bname, roc_auc_score(y, s[m])))
""")
nb.md(
    'The degree-product line is the check on the whole design: on random decoys it should be strong, '
    'on degree-matched decoys it should collapse to chance. Whatever the model achieves above the '
    'degree-matched version of this baseline is pair-level signal rather than popularity.')

# --------------------------------------------------------------------------- 8 scoreboard
nb.md(
    '## 8 · The scoreboard — headline first',
    '**The headline question is cold-node folds with degree-matched decoys.** Everything else in '
    'this table is reported to make the gap visible, exactly as the `domain_novel` table in '
    'instruction 72 reports an odds ratio against both the flattering and the hard comparison.',
    'Calls are made at the calibrated probability 0.5 for the precision/recall columns; the '
    'threshold sweep that the fill bar actually uses is §10.')
code("""
    rows, SCORES = [], {}
    for name in LAB:
        frame = LAB[name]
        for regime in ('cold_node', 'known_node'):
            if regime == 'cold_node':
                s, models, tsz = cold_node_oof(name)
                fold = fold_of[name]
            else:
                s, fold = known_node_oof(name)
                models, tsz = {}, {}
            used_cols = tsz
            m = s.notna() & fold.notna() if regime == 'cold_node' else s.notna()
            y = frame.y[m]
            cal = cross_fitted_isotonic(s[m], y, fold[m])
            SCORES[(name, regime)] = dict(raw=s[m], cal=cal, y=y, fold=fold[m], models=models,
                                          used=used_cols)
            per = []
            for f in sorted(fold[m].dropna().unique()):
                sel = fold[m] == f
                if sel.sum() < 10 or len(set(y[sel])) < 2:
                    continue
                called = (cal[sel].fillna(s[m][sel]) >= 0.5).values
                t = y[sel].values
                tp, fp, fn = int((called & t).sum()), int((called & ~t).sum()), int((~called & t).sum())
                per.append(dict(fold=int(f), n=int(sel.sum()),
                                precision=tp / (tp + fp) if tp + fp else np.nan,
                                recall=tp / (tp + fn) if tp + fn else np.nan,
                                roc_auc=roc_auc_score(t, s[m][sel]),
                                average_precision=average_precision_score(t, s[m][sel])))
            per = pd.DataFrame(per, columns=['fold', 'n', 'precision', 'recall', 'roc_auc',
                                             'average_precision'])
            for metric in ('precision', 'recall', 'roc_auc', 'average_precision'):
                rows.append(dict(negatives=name, regime=regime, strategy='model: all families',
                                 metric=metric, mean=per[metric].mean(),
                                 fold_min=per[metric].min(), fold_max=per[metric].max(),
                                 n_scored=int(m.sum()), n_positive=int(y.sum())))
            for bname, bs in BASE[name].items():
                sb = bs[m]
                per_b = []
                for f in sorted(fold[m].dropna().unique()):
                    sel = fold[m] == f
                    if sel.sum() < 10 or len(set(y[sel])) < 2:
                        continue
                    thr = np.nanmedian(sb[sel]) if bname.startswith('baseline: both') or bname.endswith('random') else 0.5
                    called = (sb[sel] >= max(thr, 1e-9)).values
                    t = y[sel].values
                    tp, fp, fn = int((called & t).sum()), int((called & ~t).sum()), int((~called & t).sum())
                    per_b.append(dict(fold=int(f), precision=tp / (tp + fp) if tp + fp else np.nan,
                                      recall=tp / (tp + fn) if tp + fn else np.nan,
                                      roc_auc=roc_auc_score(t, sb[sel]),
                                      average_precision=average_precision_score(t, sb[sel])))
                per_b = pd.DataFrame(per_b, columns=['fold', 'precision', 'recall', 'roc_auc',
                                                     'average_precision'])
                for metric in ('precision', 'recall', 'roc_auc', 'average_precision'):
                    rows.append(dict(negatives=name, regime=regime, strategy=bname, metric=metric,
                                     mean=per_b[metric].mean(), fold_min=per_b[metric].min(),
                                     fold_max=per_b[metric].max(), n_scored=int(m.sum()),
                                     n_positive=int(y.sum())))
    scoreboard = pd.DataFrame(rows)
    scoreboard.to_csv(OUT / 'scoreboard.tsv', sep='\\t', index=False)
    head = scoreboard[(scoreboard.negatives == 'degree_matched') & (scoreboard.regime == 'cold_node')]
    head.pivot_table(index='strategy', columns='metric', values='mean').round(3)
""")
nb.md('The same strategies on every negative set and both regimes — the gap is the point.')
code("""
    wide = scoreboard[scoreboard.metric == 'roc_auc'].pivot_table(
        index='strategy', columns=['regime', 'negatives'], values='mean').round(3)
    print('ROC-AUC, fold mean:'); print(wide.to_string())
    prec = scoreboard[scoreboard.metric == 'precision'].pivot_table(
        index='strategy', columns=['regime', 'negatives'], values='mean').round(3)
    print(); print('precision at a calibrated 0.5:'); print(prec.to_string())
    wide
""")

# --------------------------------------------------------------------------- 9 feature families
nb.md(
    '## 9 · Each family alone, and each family held out — against the hard comparison only',
    'Instruction 72: *"report it against the hard negatives, never against all genes or all '
    'negatives"*. Every family is therefore scored on the cold-node, degree-matched question. '
    '`alone` fits on that family only; `without` fits on everything else, so the drop is what the '
    'family adds that nothing else supplies.')
code("""
    fam_rows = []
    for fam in sorted(set(FAMILIES.values())):
        cols = [f for f in FEATURES if FAMILIES[f] == fam]
        s_alone, _, _ = cold_node_oof('degree_matched', features=cols)  # family alone
        s_wo, _, _ = cold_node_oof('degree_matched', features=[f for f in FEATURES if FAMILIES[f] != fam])
        m = fold_of['degree_matched'].notna() & s_alone.notna() & s_wo.notna()
        y = LAB['degree_matched'].y[m]
        fam_rows.append(dict(family=fam, n_features=len(cols), n_scored=int(m.sum()),
                             auc_alone=roc_auc_score(y, s_alone[m]) if m.any() else np.nan,
                             ap_alone=average_precision_score(y, s_alone[m]) if m.any() else np.nan,
                             auc_without=roc_auc_score(y, s_wo[m]) if m.any() else np.nan,
                             coverage=float(X['degree_matched'][cols].notna().any(axis=1).mean())))
    full_auc = roc_auc_score(LAB['degree_matched'].y[fold_of['degree_matched'].notna()],
                             SCORES[('degree_matched', 'cold_node')]['raw'])
    fam_table = pd.DataFrame(fam_rows)
    fam_table['auc_drop_when_removed'] = full_auc - fam_table.auc_without
    fam_table = fam_table.sort_values('auc_alone', ascending=False)
    fam_table.to_csv(OUT / 'family_ablation.tsv', sep='\\t', index=False)
    print('full model cold-node degree-matched AUC: %.3f' % full_auc)
    fam_table.round(3)
""")
nb.md(
    '**Permutation importance by family and by feature**, on the held-out fold only: the column (or '
    'the whole family block) is shuffled inside the held-out pairs and the loss of average precision '
    'is recorded. 5 repeats, seed 0.')
code("""
    from sklearn.inspection import permutation_importance
    frame = LAB['degree_matched']
    m = fold_of['degree_matched'].notna()
    models = SCORES[('degree_matched', 'cold_node')]['models']
    used = SCORES[('degree_matched', 'cold_node')]['used']
    imp_rows = []
    for f, mod in models.items():
        te = frame.index[(node_block.reindex(frame.gene_a).values == f) & (node_block.reindex(frame.gene_b).values == f)]
        if len(te) < 20:
            continue
        cols = used[f]
        r = permutation_importance(mod, X['degree_matched'].loc[te, cols], frame.y.loc[te],
                                   scoring='average_precision', n_repeats=5, random_state=0, n_jobs=1)
        for feat, mean, sd in zip(cols, r.importances_mean, r.importances_std):
            imp_rows.append(dict(fold=f, feature=feat, family=FAMILIES[feat], importance=mean, sd=sd))
    importance = (pd.DataFrame(imp_rows).groupby(['feature', 'family'])
                  .agg(importance=('importance', 'mean'), sd=('importance', 'std'),
                       folds=('fold', 'nunique')).reset_index()
                  .sort_values('importance', ascending=False))
    importance.to_csv(OUT / 'importance.tsv', sep='\\t', index=False)
    fam_imp = importance.groupby('family').importance.sum().sort_values(ascending=False)
    print('permutation importance summed by family:'); print(fam_imp.round(4).to_string())
    print(); print('top 15 features:')
    print(importance.head(15).round(4).to_string(index=False))
    importance.head(15)
""")

# --------------------------------------------------------------------------- 10 calibration + thresholds
nb.md(
    '## 10 · Calibration, the threshold sweep, the bound, and how many pairs would fill',
    'Out-of-fold scores are mapped to probabilities by **cross-fitted isotonic regression** (the '
    'calibrator for a fold is fitted only on the other folds), and the expected calibration error is '
    'reported on 10 equal-count bins.')
code("""
    cal_rows = []
    for (name, regime), d in SCORES.items():
        cal, y = d['cal'].dropna(), d['y']
        if len(cal) < 50:
            continue
        yy = y.loc[cal.index].astype(float)
        edges = np.quantile(cal.values, np.linspace(0, 1, 11))
        idx = np.clip(np.searchsorted(edges, cal.values, side='right') - 1, 0, 9)
        for b in range(10):
            sel = idx == b
            if not sel.any():
                continue
            cal_rows.append(dict(negatives=name, regime=regime, bin=b, n=int(sel.sum()),
                                 mean_probability=float(cal.values[sel].mean()),
                                 observed_fraction=float(yy.values[sel].mean())))
        cal_rows.append(dict(negatives=name, regime=regime, bin='ECE', n=len(cal),
                             mean_probability=np.nan, observed_fraction=_ece(cal.values, yy.values)))
    calibration = pd.DataFrame(cal_rows)
    calibration.to_csv(OUT / 'calibration.tsv', sep='\\t', index=False)
    print('expected calibration error:')
    print(calibration[calibration.bin == 'ECE'][['negatives', 'regime', 'n', 'observed_fraction']]
          .rename(columns={'observed_fraction': 'ECE'}).round(4).to_string(index=False))
    calibration[(calibration.negatives == 'degree_matched') & (calibration.regime == 'cold_node')].round(3)
""")
code("""
    thr_rows = []
    for (name, regime), d in SCORES.items():
        cal, y = d['cal'], d['y']
        use = cal.fillna(d['raw'])
        for t in THRESHOLDS:
            sel = use >= t
            n, k = int(sel.sum()), int(y[sel].sum())
            thr_rows.append(dict(negatives=name, regime=regime, threshold=t, calls=n, correct=k,
                                 precision=k / n if n else np.nan,
                                 precision_lower90=clopper_pearson_lower(k, n),
                                 recall=k / int(y.sum()) if int(y.sum()) else np.nan,
                                 reaches_bar=bool(n and clopper_pearson_lower(k, n) >= 0.90)))
    thresholds = pd.DataFrame(thr_rows)
    thresholds.to_csv(OUT / 'thresholds.tsv', sep='\\t', index=False)

    def verdict(name, regime):
        t = thresholds[(thresholds.negatives == name) & (thresholds.regime == regime) & thresholds.reaches_bar]
        if len(t):
            r = t.sort_values('threshold').iloc[0]
            return float(r.threshold), float(r.precision_lower90), int(r.calls)
        t = thresholds[(thresholds.negatives == name) & (thresholds.regime == regime)]
        r = t.loc[t.precision_lower90.idxmax()]
        return np.nan, float(r.precision_lower90), int(r.calls)

    vrows = []
    for (name, regime) in SCORES:
        thr, bound, calls = verdict(name, regime)
        vrows.append(dict(negatives=name, regime=regime, fill_threshold=thr,
                          best_lower90=bound, calls_at_that_point=calls,
                          reaches_0_90=bool(np.isfinite(thr))))
    verdicts = pd.DataFrame(vrows).sort_values(['regime', 'negatives'])
    verdicts.to_csv(OUT / 'verdicts.tsv', sep='\\t', index=False)
    print(verdicts.round(4).to_string(index=False))
    thresholds[(thresholds.negatives == 'degree_matched') & (thresholds.regime == 'cold_node')].round(4)
""")
nb.md(
    '**How many pairs would be filled.** The candidate universe is deliberately not all 33 million '
    'gene pairs; it is the pairs that already have a reason to be considered and are not already '
    'crosslink positives: the 33,646 pull-down pairs plus every co-expression top-50 neighbour pair. '
    'A model fitted on all labelled pairs scores them, the cross-fitted isotonic map from the '
    'headline question converts scores to probabilities, and the count above the fill threshold is '
    'the answer. If no threshold reaches the bar the honest count is **zero**, and the counts at '
    'lower bounds are reported as exploratory only.',
    '**A prior shift that must be stated.** The labelled set is 1:1 positive to decoy by '
    'construction, while the true prevalence of interacting pairs among candidates is far lower, so '
    'a calibrated probability transferred from the labelled set to the candidate set is **not** a '
    'posterior probability for a candidate pair. It is a calibrated *ranking* score at the labelled '
    'prevalence. The fill count below is therefore an upper bound on what could be filled, and it '
    'is reported only because the threshold question has to be asked; nothing is filled unless the '
    'held-out bound licenses it.')
code("""
    cand = pd.Index(sorted((set(pull.pkey) | set(cox.index)) - set(pos.pkey)))
    ca = np.array([k.split('|')[0] for k in cand], dtype=object)
    cbv = np.array([k.split('|')[1] for k in cand], dtype=object)
    ok = pd.Index(ca).isin(GENES) & pd.Index(cbv).isin(GENES)
    cand, ca, cbv = cand[ok], ca[ok], cbv[ok]
    print('candidate pairs scored: %d (%d pull-down, %d co-expression, overlap removed)'
          % (len(cand), len(set(pull.pkey) - set(pos.pkey)), len(set(cox.index) - set(pos.pkey))))

    FULL_COLS = usable(X['degree_matched'], FEATURES, LAB['degree_matched'].index)
    full = make_model().fit(X['degree_matched'][FULL_COLS], LAB['degree_matched'].y)
    Xc = pair_features(ca, cbv)
    cand_raw = pd.Series(full.predict_proba(Xc[FULL_COLS])[:, 1], index=Xc.index)
    d = SCORES[('degree_matched', 'cold_node')]
    xs, ys = _isotonic(d['raw'].values, d['y'].astype(float).values)
    cand_cal = pd.Series(np.interp(cand_raw.values, xs, ys), index=cand_raw.index)
    fill_rows = []
    for t in THRESHOLDS:
        row = thresholds[(thresholds.negatives == 'degree_matched') & (thresholds.regime == 'cold_node')
                         & (thresholds.threshold == t)]
        fill_rows.append(dict(threshold=t, held_out_lower90=float(row.precision_lower90.iloc[0]),
                              reaches_bar=bool(row.reaches_bar.iloc[0]),
                              candidate_pairs_at_or_above=int((cand_cal >= t).sum())))
    fill = pd.DataFrame(fill_rows)
    fill.to_csv(OUT / 'fill_counts.tsv', sep='\\t', index=False)
    print(fill.round(4).to_string(index=False))
    print('pairs that may be written as FACTS under the bar: %d'
          % int(fill[fill.reaches_bar].candidate_pairs_at_or_above.max() if fill.reaches_bar.any() else 0))
    fill
""")

# --------------------------------------------------------------------------- 11 known complexes
nb.md(
    '## 11 · Known complexes: are their pairs ranked highly?',
    'The complexes instruction 71 names for exactly this check: the **MyoA glideosome** (MyoA, MLC1, '
    'ELC1, GAP45, GAP50, GAP40), **AMA1-RON2** and **MIC2-M2AP**. Every within-complex pair is '
    'looked up. A pair that is itself a crosslink positive is reported with its **out-of-fold** '
    'score, so the number is honest; a pair that is not a positive is reported with the full-fit '
    'model score, which is not out of fold for the complex members and is labelled as such. Ranks '
    'are percentiles among all candidate pairs scored in §10.')
code("""
    def resolve(pattern):
        hit = nodes.index[nodes['product'].astype(str).str.contains(pattern, case=False, regex=True, na=False)]
        return hit[0] if len(hit) else None

    MEMBERS = {'MyoA': 'myosin A$', 'MLC1': 'myosin light chain MLC1', 'ELC1': 'essential light chain ELC1',
               'GAP45': 'gliding-associated protein GAP45', 'GAP50': 'gliding-associated protein GAP50',
               'GAP40': 'GAP40 protein', 'AMA1': 'apical membrane antigen AMA1',
               'RON2': 'rhoptry neck protein RON2$', 'MIC2': 'microneme protein MIC2$',
               'M2AP': 'MIC2-associated protein M2AP'}
    resolved = {k: resolve(v) for k, v in MEMBERS.items()}
    print('resolved gene IDs:', {k: v for k, v in resolved.items()})
    COMPLEXES = {'MyoA glideosome': ['MyoA', 'MLC1', 'ELC1', 'GAP45', 'GAP50', 'GAP40'],
                 'AMA1-RON2': ['AMA1', 'RON2'], 'MIC2-M2AP': ['MIC2', 'M2AP']}
    import itertools
    rows = []
    POSKEYS = set(pos.pkey)
    kn = SCORES[('degree_matched', 'known_node')]
    kn_cal = kn['cal'].fillna(kn['raw']).dropna()
    for cname, members in COMPLEXES.items():
        for u, v in itertools.combinations(members, 2):
            a, b = resolved[u], resolved[v]
            if a is None or b is None:
                rows.append(dict(complex=cname, pair=f'{u}-{v}', status='gene id unresolved')); continue
            k = key([a], [b])[0]
            is_pos = k in POSKEYS
            if is_pos and np.isfinite(d['cal'].get(k, np.nan)):
                score, how = float(d['cal'][k]), 'crosslink positive, cold-node OUT-OF-FOLD score'
            elif is_pos and np.isfinite(kn_cal.get(k, np.nan)):
                score, how = float(kn_cal[k]), ('crosslink positive, dropped by the cold-node rule; '
                                                'known-node OUT-OF-FOLD score instead')
            elif is_pos:
                score, how = np.nan, 'crosslink positive, held out by neither regime (always in training)'
            elif k in cand_cal.index:
                score, how = float(cand_cal[k]), 'candidate, full-fit score (not out of fold)'
            else:
                f1 = pair_features([a], [b])
                score, how = float(np.interp(full.predict_proba(f1[FULL_COLS])[:, 1], xs, ys)[0]), \
                    'not a candidate (no pull-down and not a co-expression neighbour), full-fit score'
            rows.append(dict(complex=cname, pair=f'{u}-{v}', gene_a=a, gene_b=b,
                             crosslink_positive=is_pos, probability=score, how=how,
                             percentile_among_candidates=float((cand_cal < score).mean()) if np.isfinite(score) else np.nan,
                             in_pulldown=bool(k in set(pull.pkey))))
    known_complexes = pd.DataFrame(rows)
    known_complexes.to_csv(OUT / 'known_complexes.tsv', sep='\\t', index=False)
    _g = known_complexes[(known_complexes.complex == 'MyoA glideosome')
                         & known_complexes.percentile_among_candidates.notna()]
    glide_scored = int(len(_g))
    glide_high = int((_g.percentile_among_candidates >= 0.90).sum())
    _a = known_complexes[known_complexes.pair == 'AMA1-RON2'].percentile_among_candidates
    ama_pct = float(_a.iloc[0]) if len(_a) and np.isfinite(_a.iloc[0]) else float('nan')
    print('glideosome pairs scored %d, of which %d above the 90th candidate percentile; '
          'AMA1-RON2 at the %.1fth percentile' % (glide_scored, glide_high, 100 * ama_pct))
    print(known_complexes[['complex', 'pair', 'crosslink_positive', 'in_pulldown', 'probability',
                           'percentile_among_candidates', 'how']].round(4).to_string(index=False))
    known_complexes
""")

# --------------------------------------------------------------------------- 12 pull-down set
nb.md(
    '## 12 · The pull-down pairs as a separate, weaker question — and what a background frequency would add',
    'Same machinery, positives = the 33,646 bait-member pairs at a provisional 0.5, degree-matched '
    'decoys on the pull-down degree, cold-node folds. This is **not** a headline result and it is '
    'not comparable to the crosslink question: the positives include every sticky contaminant that '
    'co-purified with any bait. If this set scores *better* than the crosslink set, that is evidence '
    'of exactly that contamination — the model would be learning abundance.')
code("""
    pdeg = collections.Counter()
    for a, b in zip(pull.gene_a, pull.gene_b):
        pdeg[a] += 1; pdeg[b] += 1
    pull_nodes = pd.Index(sorted(set(pull.gene_a) | set(pull.gene_b)))
    pull_degree = pd.Series(pdeg, dtype=float).reindex(pull_nodes).fillna(0.0)
    pog = nodes['orthogroup'].reindex(pull_nodes).astype(str)
    rng2 = np.random.default_rng(SEED)
    pg = sorted(pog.unique()); pperm = rng2.permutation(len(pg))
    pblock_of_group = {pg[int(j)]: int(i % 5) for i, j in enumerate(pperm)}
    pull_block = pog.map(pblock_of_group)

    # one degree-matched decoy per positive, same block, nearest pull-down degree
    plogdeg = np.log1p(pull_degree)
    pb_genes = {f: np.asarray(pull_degree.index[pull_block.values == f], dtype=object) for f in range(5)}
    pb_logdeg = {f: plogdeg.reindex(pb_genes[f]).values for f in range(5)}
    pblock = {g: int(f) for g, f in pull_block.items()}
    plog = {g: float(v) for g, v in plogdeg.items()}
    r = np.random.default_rng(SEED)
    rows, taken = [], set(ANY_EVIDENCE)
    for a, b in zip(pull.gene_a.values, pull.gene_b.values):
        keep, swap = (a, b) if r.random() < 0.5 else (b, a)
        f = pblock[swap]
        genes_f, dd = pb_genes[f], np.abs(pb_logdeg[f] - plog[swap])
        width, pick = 0.25, None
        for _ in range(7):
            sel = np.flatnonzero(dd <= width)
            if len(sel):
                r.shuffle(sel)
                for j in sel[:60]:
                    c = genes_f[j]
                    if c != keep and (c + '|' + keep if c < keep else keep + '|' + c) not in taken:
                        pick = c
                        break
            if pick is not None:
                break
            width *= 2
        if pick is None:
            continue
        taken.add(pick + '|' + keep if pick < keep else keep + '|' + pick)
        rows.append((keep, pick))
    pdec = pd.DataFrame(rows, columns=['gene_a', 'gene_b'])
    pdec['pkey'] = key(pdec.gene_a, pdec.gene_b)
    pdec = pdec.drop_duplicates('pkey')
    pframe = pd.concat([pull[['gene_a', 'gene_b']].assign(y=True),
                        pdec[['gene_a', 'gene_b']].assign(y=False)], ignore_index=True)
    pframe['pkey'] = key(pframe.gene_a, pframe.gene_b)
    pframe = pframe.drop_duplicates('pkey').set_index('pkey')
    LAB['pulldown_degree_matched'] = pframe
    X['pulldown_degree_matched'] = pair_features(pframe.gene_a.values, pframe.gene_b.values).reindex(pframe.index)
    print('pull-down question: %d pairs (%d positive, %d decoy) over %d proteins'
          % (len(pframe), int(pframe.y.sum()), int((~pframe.y).sum()), len(pull_nodes)))
""")
code("""
    blk_a = pull_block.reindex(pframe.gene_a).values
    blk_b = pull_block.reindex(pframe.gene_b).values
    s = pd.Series(np.nan, index=pframe.index, dtype=float)
    pf = pd.Series(np.nan, index=pframe.index, dtype=float)
    for f in range(5):
        te = pframe.index[(blk_a == f) & (blk_b == f)]
        tr = pframe.index[(blk_a != f) & (blk_b != f)]
        if not len(te) or len(set(pframe.y[tr])) < 2:
            continue
        pcols = usable(X['pulldown_degree_matched'], FEATURES, tr)
        mod = make_model().fit(X['pulldown_degree_matched'].loc[tr, pcols], pframe.y.loc[tr])
        s.loc[te] = mod.predict_proba(X['pulldown_degree_matched'].loc[te, pcols])[:, 1]
        pf.loc[te] = f
    mm = s.notna()
    yy = pframe.y[mm]
    cal_p = cross_fitted_isotonic(s[mm], yy, pf[mm])
    use_p = cal_p.fillna(s[mm])
    cands = [(clopper_pearson_lower(int(yy[use_p >= t].sum()), int((use_p >= t).sum())), t)
             for t in THRESHOLDS if int((use_p >= t).sum()) > 0]
    best = max(cands) if cands else (float('nan'), float('nan'))
    print('pull-down question, cold-node degree-matched: %d scored pairs, AUC %.3f, best lower 90%% bound %.3f at threshold %.2f'
          % (int(mm.sum()), roc_auc_score(yy, s[mm]), best[0], best[1]))
    pull_board = pd.DataFrame([dict(negatives='pulldown_degree_matched', regime='cold_node',
                                    strategy='model: all families', metric='roc_auc',
                                    mean=roc_auc_score(yy, s[mm]), fold_min=np.nan, fold_max=np.nan,
                                    n_scored=int(mm.sum()), n_positive=int(yy.sum()))])
    scoreboard = pd.concat([scoreboard, pull_board], ignore_index=True)
    scoreboard.to_csv(OUT / 'scoreboard.tsv', sep='\\t', index=False)
    print()
    print('Compare: crosslink question cold-node degree-matched AUC %.3f on %d pairs.'
          % (full_auc, int(fold_of['degree_matched'].notna().sum())))
""")

# --------------------------------------------------------------------------- 13 predictions
nb.md(
    '## 13 · The top predicted pairs, with their features',
    'Written to `predictions.tsv` with the features that drove them, their calibrated probability, '
    'and an explicit `status` column saying whether the pair may be written as a fact under the bar. '
    'No pair is presented as a fact unless §10 licenses it.')
code("""
    top = cand_cal.sort_values(ascending=False).head(500)
    pf_top = Xc.reindex(top.index)
    pred = pd.DataFrame({'gene_a': [k.split('|')[0] for k in top.index],
                         'gene_b': [k.split('|')[1] for k in top.index],
                         'probability': top.values}, index=top.index)
    pred['product_a'] = nodes['product'].reindex(pred.gene_a).values
    pred['product_b'] = nodes['product'].reindex(pred.gene_b).values
    pred['in_pulldown'] = pred.index.isin(set(pull.pkey))
    pred['coexpr_percentile'] = pf_top['coexpr_percentile'].values
    pred['fit_profile_corr'] = pf_top['fit_profile_corr'].values
    pred['dom_shared'] = pf_top['dom_shared'].values
    pred['same_comp_obs'] = pf_top['same_comp_obs'].values
    pred['same_comp_pred'] = pf_top['same_comp_pred'].values
    pred['same_orthogroup'] = pf_top['same_orthogroup'].values
    pred['tdi_kmer_cosine'] = pf_top['tdi_kmer_cosine'].values
    bar = verdicts[(verdicts.negatives == 'degree_matched') & (verdicts.regime == 'cold_node')].iloc[0]
    bar_open = bool(bar.reaches_0_90) and bool(np.isfinite(bar.fill_threshold))
    pred['status'] = np.where((pred.probability >= bar.fill_threshold) if bar_open else False,
                              'fillable', 'exploratory')
    pred.to_csv(OUT / 'predictions.tsv', sep='\\t', index=True, index_label='pkey')
    print('predictions.tsv: %d rows; fillable %d' % (len(pred), int((pred.status == 'fillable').sum())))
    pred.head(15).round(3)
""")

# --------------------------------------------------------------------------- 14 verdict
nb.md(
    '## 14 · Verdict, what constrains it, and what would change it')
code("""
    hb = verdicts[(verdicts.negatives == 'degree_matched') & (verdicts.regime == 'cold_node')].iloc[0]
    easy_r = verdicts[(verdicts.negatives == 'random') & (verdicts.regime == 'cold_node')].iloc[0]
    easy_k = verdicts[(verdicts.negatives == 'degree_matched') & (verdicts.regime == 'known_node')].iloc[0]
    comp_r = verdicts[(verdicts.negatives == 'same_compartment') & (verdicts.regime == 'cold_node')].iloc[0]
    summary = pd.DataFrame([
        dict(question='HEADLINE cold-node, degree-matched', bound=hb.best_lower90,
             reaches=hb.reaches_0_90, n=int(scoreboard[(scoreboard.negatives == 'degree_matched') &
                                                       (scoreboard.regime == 'cold_node')].n_scored.iloc[0])),
        dict(question='cold-node, same-compartment decoys', bound=comp_r.best_lower90, reaches=comp_r.reaches_0_90,
             n=int(scoreboard[(scoreboard.negatives == 'same_compartment') & (scoreboard.regime == 'cold_node')].n_scored.iloc[0])),
        dict(question='easy: cold-node, random decoys', bound=easy_r.best_lower90, reaches=easy_r.reaches_0_90,
             n=int(scoreboard[(scoreboard.negatives == 'random') & (scoreboard.regime == 'cold_node')].n_scored.iloc[0])),
        dict(question='easy: known-node, degree-matched', bound=easy_k.best_lower90, reaches=easy_k.reaches_0_90,
             n=int(scoreboard[(scoreboard.negatives == 'degree_matched') & (scoreboard.regime == 'known_node')].n_scored.iloc[0])),
    ])
    # how many held-out calls a bound of 0.90 needs at the observed best precision
    bestp = thresholds[(thresholds.negatives == 'degree_matched') & (thresholds.regime == 'cold_node')
                       & (thresholds.calls >= 20)].precision.max()
    need = np.nan
    if np.isfinite(bestp) and bestp > 0.90:
        n = 20
        while n < 20000 and clopper_pearson_lower(int(round(n * bestp)), n) < 0.90:
            n += 1
        need = n
    xp = X['degree_matched']
    ypos = LAB['degree_matched'].y
    sc = xp['same_comp_obs'].fillna(xp['same_comp_pred'])
    samecomp_pos = float(sc[ypos].mean(skipna=True))
    samecomp_dec = float(sc[~ypos].mean(skipna=True))
    print('same compartment: %.1f %% of positives, %.1f %% of degree-matched decoys'
          % (100 * samecomp_pos, 100 * samecomp_dec))
    print(summary.round(4).to_string(index=False))
    print()
    print('best precision at any threshold with >=20 calls (headline question): %.3f' % bestp)
    print('held-out calls a 0.90 bound would need at that precision: %s'
          % ('unreachable at this precision' if not np.isfinite(need) else int(need)))
    summary
""")
nb.md(
    'The written verdict, the constraint and the candidate remedies are in the README below, '
    'generated from these numbers rather than typed.')
code("""
    def fmt(df, index=False):
        return '```\\n' + df.round(4).to_string(index=index) + '\\n```\\n'

    reaches = bool(hb.reaches_0_90)
    fam_imp_top = importance.head(10)[['feature', 'family', 'importance']]
    filled = int(fill[fill.reaches_bar].candidate_pairs_at_or_above.max()) if fill.reaches_bar.any() else 0
    R = []
    R.append('# Binding-partner inference, first scoring run (instruction 72, trait 2)\\n')
    R.append('Notebook: `notebooks/partner_inference_2026_11_04.ipynb`. '
             'Builder: `scripts/build_partner_inference_notebook.py`. Seed 0.\\n')
    R.append('## Question\\n')
    R.append('Positives: the **%d StarPath XL-MS crosslink pairs** over %d proteins, probability '
             '1 - local FDR (median %.4f). Direct physical evidence.\\n'
             % (len(pos), len(xl_nodes), pos.probability.median()))
    R.append('The %d bait-member pull-down pairs in the information map are **not** headline '
             'positives: they sit at a provisional 0.5 because no Toxoplasma background frequency background '
             'exists, and without a background frequency a sticky contaminant is indistinguishable '
             'from a real partner. They are scored separately (§12).\\n' % len(pull))
    R.append('**There is no verified non-interaction set.** Negatives are decoys — pairs with no '
             'recorded evidence — built three ways, and the headline is the hardest: for each '
             'positive, one endpoint kept and the other replaced by a protein of similar crosslink '
             'degree from the same held-out protein block, which removes the popular-protein '
             'shortcut and the "detectable by XL-MS at all" shortcut.\\n')
    R.append('## Scoreboard (lower 90 % Clopper-Pearson bound on held-out precision)\\n')
    R.append(fmt(summary) + '\\n')
    R.append('**Verdict: partner prediction does %sreach the 0.90 bar on the cold-node, '
             'degree-matched question.** Best lower bound **%.3f**%s.\\n'
             % ('' if reaches else 'NOT ', hb.best_lower90,
                ' at threshold %.2f, %d held-out calls' % (hb.fill_threshold, hb.calls_at_that_point)
                if reaches else ''))
    R.append('Pairs that may be written as facts: **%d**.\\n' % filled)
    R.append('## Negative-set and split gap (the pair analogue of the easy contrast)\\n')
    R.append(fmt(scoreboard[(scoreboard.metric == 'roc_auc')]
                 .pivot_table(index='strategy', columns=['regime', 'negatives'], values='mean'),
                 index=True) + '\\n')
    R.append('## Feature families, hard question only\\n')
    R.append(fmt(fam_table) + '\\n')
    R.append('## Top 10 features by permutation importance (held-out folds, hard question)\\n')
    R.append(fmt(fam_imp_top) + '\\n')
    R.append('## Known complexes\\n')
    R.append(fmt(known_complexes[['complex', 'pair', 'crosslink_positive', 'in_pulldown',
                                  'probability', 'percentile_among_candidates']]) + '\\n')
    R.append('## Calibration\\n')
    R.append('Expected calibration error, cross-fitted isotonic:\\n\\n'
             + fmt(calibration[calibration.bin == 'ECE'][['negatives', 'regime', 'n', 'observed_fraction']]
                   .rename(columns={'observed_fraction': 'ECE'})) + '\\n')
    R.append('## What the signal is, and what to check\\n')
    R.append('* **Co-compartment carries most of it.** The strongest family on the hard question is '
             'the predicted compartment (AUC alone %.3f) with the observed marker compartment close '
             'behind (%.3f); permutation importance is dominated by `same_comp_pred_coarse`. '
             'Crosslink positives share a compartment in %.0f %% of pairs against %.0f %% of '
             'degree-matched decoys. Two things follow: the model is largely a co-localisation '
             'predictor, and **whether the pair FDR used by StarPath relied on localisation '
             'agreement must be checked** — if it did, the strongest feature is partly circular. '
             'The crosslink table '
             'does carry `loc_a`/`loc_b` columns, which is why they are held out here.\\n'
             % (float(fam_table.set_index('family').loc['compartment_predicted', 'auc_alone']),
                float(fam_table.set_index('family').loc['compartment_observed', 'auc_alone']),
                100 * samecomp_pos, 100 * samecomp_dec))
    R.append('* **The CRISPR-fitness pair features do not reproduce the StarPath finding.** The '
             'fitness family — the absolute difference and the minimum of each of the %d `fit_*` '
             'columns, plus the profile correlation — reaches AUC %.3f alone and removing it changes '
             'the full model by %+.4f. Instruction 72 records the difference as highly predictive; '
             'on this question, with degree-matched decoys and cold-node folds, it is not. The '
             'mouse-immune-dependence note in 72 found the same thing for single genes: subtracting '
             'two relative scores amplifies noise.\\n'
             % (len(FIT), float(fam_table.set_index('family').loc['fitness', 'auc_alone']),
                -float(fam_table.set_index('family').loc['fitness', 'auc_drop_when_removed'])))
    R.append('* **The known-complex check is mixed and one miss is flat.** %d of the %d glideosome '
             'pairs that could be scored sit above the %.0fth percentile of candidates, but '
             '**AMA1-RON2 is placed at the %.1fth percentile** — the model puts one of the '
             'best-established complexes in Toxoplasma near the bottom, because AMA1 (microneme) and '
             'RON2 (rhoptry neck) are in different compartments and the model is mostly a '
             'co-localisation predictor. That is the clearest single demonstration of what it has '
             'and has not learned.\\n' % (glide_high, glide_scored, 90.0, 100 * ama_pct))
    R.append('## What constrains the answer\\n')
    R.append('* Best precision at any threshold with at least 20 held-out calls: **%.3f**. '
             'A 0.90 lower bound at that precision needs %s held-out calls, against %d scored pairs '
             'in the headline question.\\n'
             % (bestp, 'more than any label count can supply (the precision is below the bar itself)'
                if not np.isfinite(need) else '~%d' % int(need),
                int(hb.calls_at_that_point) if np.isfinite(hb.calls_at_that_point) else 0))
    R.append('* Cold-node scoring keeps only %.0f %% of the labelled pairs (%d of %d): a pair is '
             'testable only when both endpoints land in the same held-out protein block. That is '
             'inherent to the regime, not a bug.\\n'
             % (100 * fold_of['degree_matched'].notna().mean(),
                int(fold_of['degree_matched'].notna().sum()), len(LAB['degree_matched'])))
    R.append('* The decoys are assumed non-interacting. Some fraction of them interact and were '
             'never reported, which depresses measured precision by an unknown amount. Nothing in '
             'the data bounds that fraction.\\n')
    R.append('## What would change it\\n')
    R.append('* **A Toxoplasma background frequency background** (instruction 71, product 4) is the obvious '
             'candidate and the largest single gain available: it would turn %d provisional '
             'bait-member pairs into graded co-purification evidence. The positives would grow from '
             '%d to the specific subset of those %d — and, more important for the scoring, '
             'background frequency-filtered co-purification is the one *independent* assay that could serve as '
             'held-out truth for crosslink-trained predictions, which nothing in this notebook has. '
             'How much it would add cannot be stated as a number before the background exists; what '
             'can be stated is the ceiling: %d pairs are currently unusable as labels for want of '
             'it.\\n' % (len(pull), len(pos), len(pull), len(pull)))
    R.append('* Reciprocal pull-downs or solved complexes would give the first **verified '
             'negatives**; every number here rests on decoys instead.\\n')
    R.append('* AF3 pair confidence on a shortlist is the per-pair evidence the feature table lacks '
             'entirely: every feature here is a property of the two proteins separately, combined. '
             'The crosslink-derived Chai models exist but are held out as label-derived.\\n')
    R.append('## Files\\n\\n| file | contents |\\n|---|---|\\n')
    for f, what in [('scoreboard.tsv', 'every strategy x negative set x split regime, per-fold spread'),
                    ('thresholds.tsv', 'threshold sweep with the Clopper-Pearson lower 90 % bound'),
                    ('verdicts.tsv', 'the lowest threshold reaching the bar, per question'),
                    ('calibration.tsv', 'reliability bins and expected calibration error'),
                    ('importance.tsv', 'permutation importance per feature with its family'),
                    ('family_ablation.tsv', 'each family alone and held out, hard question only'),
                    ('known_complexes.tsv', 'the MyoA glideosome, AMA1-RON2, MIC2-M2AP pair by pair'),
                    ('predictions.tsv', 'top 500 candidate pairs with features and fillable/exploratory status'),
                    ('fill_counts.tsv', 'candidate pairs above each threshold'),
                    ('features.tsv', 'every feature with its family')]:
        R.append('| `%s` | %s |\\n' % (f, what))
    (OUT / 'README.md').write_text(''.join(R))
    print(''.join(R)[:4000])
""")

nb.md(
    '## 15 · Limitations, stated plainly',
    '1. **No verified negatives exist.** Every number rests on decoys. An unknown fraction of the '
    'decoys are real partners that were never reported, so the measured precision is a lower bound '
    'on the truth for the positives and an upper bound on honesty about the negatives — and the two '
    'cannot be separated with this data.\n'
    '2. **Degree matching uses the label source to build the controls.** It cannot leak a label into '
    'a feature (the degree columns are held out) but it restricts the question to pairs of '
    'XL-MS-detectable proteins.\n'
    '3. **Cold-node folds keep only about a fifth of the pairs.** Fewer held-out calls make the '
    'Clopper-Pearson bound harder to clear for a given true precision; the §6 table says exactly how '
    'much harder.\n'
    '4. **Every feature is a pair of single-protein properties combined.** There is no per-pair '
    'measurement in the feature table at all — no AF3 interface confidence, no co-fractionation '
    'profile, no co-purification specificity — because the only per-pair measurements available are '
    'the crosslinks themselves and the structure models built from them, which are the labels.\n'
    '5. **The co-expression table is top-50 per gene**, so most pairs have no co-expression value. '
    'Absence is passed to the tree as missing, but it is missing-not-at-random: it means "not in '
    'either top 50", which is itself weak evidence against.\n'
    '6. **The compartment features mix grades.** `compartment_observed` is marker-grade; '
    '`compartment_predicted` is a TAGM classifier output. They are separate families for that '
    'reason and the ablation shows what each is worth.\n'
    '7. **The 1:1 class balance is a construction, not the prevalence.** Interacting pairs are a '
    'tiny fraction of all pairs, so probabilities calibrated on the labelled set overstate the '
    'posterior for a candidate pair by the prior ratio. Held-out precision is measured at the '
    'labelled balance and means "precision among pairs like these", not "precision over the '
    'proteome".\n'
    '8. **The pull-down question in §12 is not comparable** to the headline and must not be quoted '
    'as a partner-prediction result.')

nb.write(str(ROOT / 'notebooks/partner_inference_2026_11_04.ipynb'))
print('wrote notebooks/partner_inference_2026_11_04.ipynb')
