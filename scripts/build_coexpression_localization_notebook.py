"""Instruction 71, product 1: cross-study co-expression as localization evidence, scored the instruction 70 way.

Writes notebooks/coexpression_localization_2026_11_04.ipynb (executed) and its result tables under
results/coexpression_localization_2026_11_04/. Inputs are the uniform reprocessing only (latest verified app
snapshot, salmon against ME49_host, parasite genes with TPM reclosed after dropping host genes).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from notebook_runner import ExecutedNotebook  # noqa: E402

nb = ExecutedNotebook('Co-expression across studies as localization evidence')
nb.md('# Co-expression across studies as localization evidence',
      'Instruction 71, product 1 (Einar, 2026-11-04). Genes of one compartment are made together in the cell cycle '
      'and respond together to perturbation, so genes that are co-expressed across many independent studies may share '
      'a location. This notebook builds that co-expression from the uniform RNA-seq reprocessing and asks, the '
      'instruction 70 way, whether it recovers held-out hyperLOPIT markers better than baselines.',
      '**What this is not:** a deployment. The truth here is the curation-grade marker set (718 markers, 26 niches), '
      'which mixes prior knowledge with microscopy; scores against it measure agreement with curation, and no '
      'hyperLOPIT-derived value is used as a feature (whole-source hold-out).')

nb.md('## 1 · Inputs',
      'The newest app snapshot with verified runs; one matrix per study (salmon, ME49_host reference, parasite genes '
      'only). Samples with fewer than one million parasite read counts are dropped (mock-infected or host-dominated '
      'libraries carry little parasite signal); studies keep at least six samples after that.')
nb.code('import json, re, hashlib',
        'from pathlib import Path',
        'import numpy as np, pandas as pd',
        "SNAP_ROOT = Path('/nas_mnt/carruthers/Einar/parascape_omics/app_snapshots')",
        "snap = sorted(p for p in SNAP_ROOT.iterdir() if re.fullmatch(r'\\d{8}T\\d{6}Z-\\d+', p.name) and (p / 'rna/snapshot.json').exists())[-1]",
        "meta = json.loads((snap / 'rna/snapshot.json').read_text())",
        "OUT = Path('results/coexpression_localization_2026_11_04'); OUT.mkdir(parents=True, exist_ok=True)",
        "print(snap.name, meta['verified_runs'], 'verified runs;', meta['normalization'])")
nb.code('MIN_READS, MIN_SAMPLES = 1e6, 6',
        'studies = {}',
        "for f in sorted((snap / 'rna/matrices').glob('*/salmon_ME49_host.parquet')):",
        '    m = pd.read_parquet(f)',
        "    counts, tpm = m['count'], m['tpm']",
        '    keep = counts.columns[counts.sum() >= MIN_READS]',
        '    if len(keep) >= MIN_SAMPLES:',
        '        studies[f.parent.name] = np.log2(tpm[keep].astype(float) + 1)',
        "summary = pd.DataFrame({'samples': {s: v.shape[1] for s, v in studies.items()}}).sort_values('samples', ascending=False)",
        "print(len(studies), 'studies,', int(summary.samples.sum()), 'samples')",
        'summary.T')

nb.md('## 2 · Co-expression per study, then combined by rank',
      'Within each study: genes expressed (mean log2 TPM ≥ 1) and variable (top 5,000 by variance), Pearson correlation '
      'across samples. Each study correlation matrix is turned into percentile ranks over its own gene pairs, so a '
      'study with many samples or strong batch structure does not outweigh the others. The combined score of a gene '
      'pair is its mean percentile over the studies that measured both genes; `support` counts those studies.',
      'Studies are the unit of independence: a pair supported by one study is not replicated.')
nb.code('genes = sorted(set().union(*[set(v.index) for v in studies.values()]))',
        'index = {g: i for i, g in enumerate(genes)}; G = len(genes)',
        'total = np.zeros((G, G), dtype=np.float32); support = np.zeros((G, G), dtype=np.uint8)',
        'def study_percentiles(x):',
        '    x = x[(x.mean(axis=1) >= 1)]',
        '    x = x.loc[x.var(axis=1).sort_values(ascending=False).index[:5000]]',
        '    x = x[x.var(axis=1) > 0]',
        '    c = np.corrcoef(x.values)',
        '    iu = np.triu_indices_from(c, 1)',
        '    ranks = pd.Series(c[iu]).rank(pct=True).values.astype(np.float32)',
        '    p = np.zeros_like(c, dtype=np.float32); p[iu] = ranks; p = p + p.T',
        '    return list(x.index), p',
        'def combine(names):',
        '    tot = np.zeros((G, G), dtype=np.float32); sup = np.zeros((G, G), dtype=np.uint8)',
        '    for s in names:',
        '        ids, p = study_percentiles(studies[s]); ix = np.array([index[g] for g in ids])',
        '        tot[np.ix_(ix, ix)] += p; sup[np.ix_(ix, ix)] += 1',
        '        np.fill_diagonal(sup, 0)',
        '    with np.errstate(invalid="ignore", divide="ignore"):',
        '        return np.where(sup > 0, tot / np.maximum(sup, 1), np.nan).astype(np.float32), sup',
        'score, support = combine(list(studies))',
        "print(G, 'genes; pairs supported by >= 3 studies:', int((support >= 3).sum() // 2))")

nb.md('## 3 · Truth: the 718 curation-grade hyperLOPIT markers',
      'From mmc5.xls (Barylyuk 2020, sha256 pinned in instruction 70). Reported at the 26 refined niches and at a coarse '
      'level that merges numbered sub-niches (rhoptries 1/2, ER 1/2, ...), ribosome and proteasome subunits, the '
      'nucleus and mitochondrion sub-compartments and the plasma membrane classes. `unknown` is unlabeled, never a '
      'negative.')
nb.code("MMC5 = Path('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/sources/barylyuk2020_mmc5/mmc5.xls')",
        "assert hashlib.sha256(MMC5.read_bytes()).hexdigest() == 'a58a8778714c42e344099de9429daa3be4e9c70021b0aeab3de65578c165b8f2'",
        "lopit = pd.read_excel(MMC5, 'S4a - TAGM-MAP classification')",
        "markers = lopit[lopit.markers != 'unknown'].set_index('Accession').markers",
        'def coarse(label):',
        "    label = re.sub(r' \\d+$', '', label)",
        "    for key, value in (('ribosome', 'ribosome'), ('proteasome', 'proteasome'), ('nucleus', 'nucleus'),",
        "                       ('mitochondrion', 'mitochondrion'), ('PM -', 'plasma membrane')):",
        '        if key in label: return value',
        '    return label',
        'truth = pd.DataFrame({"fine": markers, "coarse": markers.map(coarse)})',
        'truth = truth[truth.index.isin(index)]',
        "print(len(markers), 'markers;', len(truth), 'measured in the co-expression;', truth.fine.nunique(), 'fine,', truth.coarse.nunique(), 'coarse labels')",
        'truth.coarse.value_counts().to_frame().T')

nb.md('## 4 · Five folds grouped by orthogroup, one frozen seed',
      'Paralogs share an orthogroup and often a location, so a paralog in the training fold would leak the answer. '
      'Whole orthogroups go to one fold (OrthoMCL 6.21; genes without a group are their own group). Seed 0, recorded.')
nb.code("og = pd.read_csv('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/MASTER_parasite_wide_by_orthogroup.csv', usecols=['orthomcl_group', 'tgon_gene_id']).dropna()",
        "group = {g: row.orthomcl_group for row in og.itertuples() for g in str(row.tgon_gene_id).split(';')}",
        "groups = truth.index.map(lambda g: group.get(g, 'self:' + g))",
        'rng = np.random.default_rng(0); unique = np.array(sorted(set(groups))); rng.shuffle(unique)',
        'fold_of_group = {g: i % 5 for i, g in enumerate(unique)}',
        "truth['fold'] = [fold_of_group[g] for g in groups]",
        "truth['group'] = list(groups)",
        "truth.fold.value_counts().sort_index().to_frame().T")

nb.md('## 5 · The strategy and the baselines',
      '**Co-expression kNN:** for a held-out gene, take its 25 strongest partners among training-fold markers, '
      'requiring support from at least 3 studies; each partner votes for its label with weight = combined percentile '
      'minus 0.5 (only above-median pairs vote). Confidence = winning vote share.',
      '**Baselines, same folds:** training-fold majority label; prevalence-weighted random (expected accuracy Σp²); '
      '**random partners** -- the same kNN with partners chosen at random, which isolates the value of the co-expression '
      'itself. Instruction 70 reports improvement over the *best* baseline.')
nb.code('def knn(score, support, level, k=25, min_support=3, random_partners=False, seed=0):',
        '    rng = np.random.default_rng(seed); rows = []',
        '    for fold in range(5):',
        '        train = truth[truth.fold != fold]; test = truth[truth.fold == fold]',
        '        t_idx = np.array([index[g] for g in train.index]); t_lab = train[level].values',
        '        for g, true in zip(test.index, test[level]):',
        '            s = score[index[g], t_idx].copy(); ok = (support[index[g], t_idx] >= min_support) & np.isfinite(s)',
        '            if random_partners: s = np.where(ok, rng.random(len(s)), np.nan)',
        '            s[~ok] = -np.inf',
        '            top = np.argsort(-s)[:k]; top = top[np.isfinite(s[top])]',
        '            w = (score[index[g], t_idx][top] - 0.5).clip(min=0) if not random_partners else np.ones(len(top))',
        '            if len(top) == 0 or w.sum() == 0:',
        '                rows.append((g, fold, true, None, 0.0)); continue',
        '            votes = pd.Series(w).groupby(t_lab[top]).sum()',
        '            rows.append((g, fold, true, votes.idxmax(), votes.max() / votes.sum()))',
        "    return pd.DataFrame(rows, columns=['gene', 'fold', 'true', 'pred', 'confidence'])",
        'def summarize(pred, name):',
        "    called = pred.dropna(subset=['pred'])",
        "    per_fold = called.groupby('fold').apply(lambda d: (d.true == d.pred).mean(), include_groups=False)",
        "    return dict(strategy=name, coverage=len(called) / len(pred), accuracy=(called.true == called.pred).mean(),",
        "                fold_min=per_fold.min(), fold_max=per_fold.max())",
        'def baselines(level):',
        '    out = []',
        '    for fold in range(5):',
        '        train = truth[truth.fold != fold]; test = truth[truth.fold == fold]',
        '        p = train[level].value_counts(normalize=True)',
        "        out.append(dict(fold=fold, majority=(test[level] == p.idxmax()).mean(), prevalence_random=float((p ** 2).sum())))",
        '    return pd.DataFrame(out)')

nb.md('## 6 · Results at both label levels')
nb.code('tables, preds = [], {}',
        "for level in ('coarse', 'fine'):",
        '    co = knn(score, support, level); rnd = knn(score, support, level, random_partners=True)',
        "    preds[level] = co",
        '    b = baselines(level)',
        "    rows = [summarize(co, 'co-expression kNN'), summarize(rnd, 'random partners (baseline)'),",
        "            dict(strategy='training majority (baseline)', coverage=1.0, accuracy=b.majority.mean(), fold_min=b.majority.min(), fold_max=b.majority.max()),",
        "            dict(strategy='prevalence random (baseline)', coverage=1.0, accuracy=b.prevalence_random.mean(), fold_min=b.prevalence_random.min(), fold_max=b.prevalence_random.max())]",
        "    t = pd.DataFrame(rows); t.insert(0, 'level', level); tables.append(t)",
        'board = pd.concat(tables, ignore_index=True)',
        "board.to_csv(OUT / 'scoreboard.tsv', sep='\\t', index=False)",
        'board.round(3)')

nb.md('## 7 · Does the strategy beat the best baseline in each fold?',
      'Instruction 70 initial deployment rule needs: better than the best baseline in at least 4 of 5 folds, and '
      '≥ 80 % precision of confident calls on direct-evidence genes. The second half needs the direct tier (pending, '
      'instruction 70); here the precision of confident calls on the curation markers is reported instead, '
      'clearly labelled.')
nb.code('rule = []',
        "for level in ('coarse', 'fine'):",
        '    co = preds[level]; b = baselines(level); rnd = knn(score, support, level, random_partners=True)',
        '    for fold in range(5):',
        "        c = co[(co.fold == fold) & co.pred.notna()]; r = rnd[(rnd.fold == fold) & rnd.pred.notna()]",
        '        best = max(b.majority[fold], b.prevalence_random[fold], (r.true == r.pred).mean())',
        "        rule.append(dict(level=level, fold=fold, coexpression=(c.true == c.pred).mean(), best_baseline=best))",
        'rule = pd.DataFrame(rule); rule["beats"] = rule.coexpression > rule.best_baseline',
        "rule.groupby('level').beats.sum().to_frame('folds beating best baseline (of 5)')")
nb.code('cal = []',
        "for level in ('coarse', 'fine'):",
        "    c = preds[level].dropna(subset=['pred'])",
        '    for t in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8):',
        '        sel = c[c.confidence >= t]',
        "        cal.append(dict(level=level, confidence_at_least=t, calls=len(sel), precision_vs_curation=(sel.true == sel.pred).mean() if len(sel) else np.nan))",
        'cal = pd.DataFrame(cal); cal.to_csv(OUT / "confidence.tsv", sep="\\t", index=False)',
        'cal.round(3)')

nb.md('## 8 · Which compartments does co-expression place?',
      'Recall per coarse compartment (held-out markers called correctly / held-out markers). Compartments built in a '
      'burst during division (IMC, rhoptries, micronemes, apical) are expected to do best; housekeeping locations with '
      'flat expression (cytosol, ER) worst.')
nb.code("c = preds['coarse']",
        "per = c.assign(hit=c.true == c.pred).groupby('true').agg(markers=('hit', 'size'), recall=('hit', 'mean'))",
        "per.sort_values('recall', ascending=False).round(2).T")

nb.md('## 9 · Replication across independent studies',
      'The studies are split into two halves (alternating by sample count, so both halves get large and small studies); '
      'the network is rebuilt from each half alone and scored on the same folds. If the signal is real it should hold '
      'in both halves, not come from one study.')
nb.code('order = list(summary.index); half_a, half_b = order[0::2], order[1::2]',
        'rep = []',
        "for name, half in (('half A', half_a), ('half B', half_b)):",
        '    s, sup = combine(half)',
        "    rep.append(dict(studies=name, n_studies=len(half), **{k: v for k, v in summarize(knn(s, sup, 'coarse', min_support=2), 'co-expression kNN').items() if k != 'strategy'}))",
        'pd.DataFrame(rep).round(3)')

nb.md('## 10 · Saved for instruction 70',
      'The per-gene neighbour list (top 50 partners with combined percentile and study support) is written as an evidence '
      'block with its provenance (snapshot, studies, rule) so the scoreboard can register it and keep it out of any '
      'feature set it would leak into. No prediction is written into the knowledge matrix from this notebook.')
nb.code('rows = []',
        'for i, g in enumerate(genes):',
        '    s = np.where(support[i] >= 3, score[i], -np.inf); top = np.argsort(-s)[:50]',
        '    rows += [(g, genes[j], float(score[i, j]), int(support[i, j])) for j in top if np.isfinite(s[j])]',
        "neighbours = pd.DataFrame(rows, columns=['gene', 'partner', 'combined_percentile', 'study_support'])",
        "neighbours.to_parquet(OUT / 'coexpression_neighbours.parquet', index=False)",
        "(OUT / 'provenance.json').write_text(json.dumps(dict(snapshot=snap.name, studies=list(studies), min_reads=MIN_READS,",
        "    min_samples=MIN_SAMPLES, rule='per-study Pearson on log2(TPM+1), top 5000 variable expressed genes, percentile ranks, mean over studies',",
        "    provenance_family='uniform_rnaseq_coexpression', excluded_sources=['lopit_tgon (truth only)'],",
        "    truth='Barylyuk 2020 mmc5 markers (curation grade)'), indent=2) + '\\n')",
        "print(len(neighbours), 'neighbour rows for', neighbours.gene.nunique(), 'genes')")

nb.write(str(ROOT / 'notebooks/coexpression_localization_2026_11_04.ipynb'))
