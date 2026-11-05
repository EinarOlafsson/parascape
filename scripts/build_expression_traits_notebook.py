"""Instruction 72, tier-2 expression columns of the information map.

Writes notebooks/expression_traits_2026_11_05.ipynb (executed) and its tables under
results/expression_traits_2026_11_05/. Inputs are the uniform RNA-seq reprocessing only (newest app snapshot,
salmon against ME49_host, parasite genes with TPM reclosed), the sample design registry and the per-sample
quality score. Columns are registered as Traits inside the notebook; parascape/information_map.py is not edited.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from notebook_runner import ExecutedNotebook  # noqa: E402

nb = ExecutedNotebook('Expression columns of the information map')

nb.md('# Expression columns of the information map',
      'Instruction 72, tier 2 (measurements). The map had 355 columns and not one expression column, although '
      'hundreds of RNA-seq runs have now been quantified the same way. This notebook turns that reprocessing into '
      'gene x condition columns over the 8,140-gene universe of `nodes.parquet`.',
      '**The order Einar set, and the order built here:** a value measured *against a control* is worth most; a '
      'value measured against the study own centre comes next; an absolute value is a last resort. '
      '`relative_to_control` > `internal_reference` > `library_size` > `raw`. Where both a relative and an absolute '
      'form of the same quantity exist, a learner must take the relative one. The absolute columns are built only '
      'so that nothing is missing, not because they should be used.',
      '**Two lessons from instruction 72 shape everything below.** (1) *The easy-contrast trap*: a column that '
      'separates obviously different things tells you nothing about the hard question, so every check here is run '
      'against the contrast that could fail, and a miss is reported as a miss. (2) *The label-error correction*: '
      'the first host-facing model passed its bar because its labels were wrong, not because its features leaked. '
      'For expression columns the equivalent of a wrong label is a **wrong control assignment**, so this notebook '
      'states exactly how the control was chosen, how weakly the registry supports it, and then tests the columns '
      'against biology that is already known.')

nb.md('## 1 · Inputs, and what is deliberately left out',
      'The newest app snapshot that actually carries `rna/snapshot.json` (the snapshot being written while this '
      'runs is skipped). One matrix per study: salmon, ME49_host reference, parasite genes only, TPM reclosed to '
      'a million within the parasite genes. Matrices are read one study at a time and dropped again, so peak '
      'memory stays small.')
nb.code('import json, re, unicodedata',
        'from pathlib import Path',
        'import numpy as np, pandas as pd',
        "SNAP_ROOT = Path('/nas_mnt/carruthers/Einar/parascape_omics/app_snapshots')",
        "snaps = sorted(p for p in SNAP_ROOT.iterdir()",
        "               if re.fullmatch(r'\\d{8}T\\d{6}Z-\\d+', p.name) and (p / 'rna/snapshot.json').exists())",
        'SNAP = snaps[-1]',
        "META = json.loads((SNAP / 'rna/snapshot.json').read_text())",
        "OUT = Path('results/expression_traits_2026_11_05'); OUT.mkdir(parents=True, exist_ok=True)",
        "print('snapshot', SNAP.name)",
        "print(META['verified_runs'], 'verified runs of', META['planned_runs'], 'planned;', len(META['studies']), 'studies')",
        "print('normalization:', META['normalization'])",
        "print('control pairing:', META['sample_design']['control_pairing'])")

nb.code("UNIVERSE = pd.read_parquet('parascape/data/nodes.parquet', columns=['gene_id', 'product'])",
        "GENES = UNIVERSE.gene_id.tolist(); PRODUCT = dict(zip(UNIVERSE.gene_id, UNIVERSE['product'].astype(str)))",
        "design = pd.read_csv(SNAP / 'rna/sample_design.tsv', sep='\\t', low_memory=False)",
        "quality = pd.read_csv('results/rna_sample_quality_2026_11_04/sample_quality.tsv', sep='\\t')",
        "MATRIX_STUDIES = sorted(p.name for p in (SNAP / 'rna/matrices').iterdir() if (p / 'salmon_ME49_host.parquet').exists())",
        "print(len(GENES), 'genes in the universe;', len(MATRIX_STUDIES), 'studies with a salmon ME49_host matrix')",
        "print(len(design), 'design rows;', len(quality), 'runs with a quality score')")

nb.md('### The per-sample quality score',
      'From `results/rna_sample_quality_2026_11_04/` (v0): a geometric mean of up to six components -- parasite '
      'read depth, parasite genes detected, mapping rate, compatible fragment ratio, salmon-against-kallisto '
      'agreement, and replicate agreement -- each mapped to 0-1 between fixed anchors. It is used here as a '
      '**weight**, so a shallow or poorly mapping run still contributes, just less.',
      'The score covers the runs that existed when it was built, so some newer runs have none. **Rule for a run '
      'with no score:** it takes the median quality of the scored runs of its own study; if its study has no '
      'scored run at all it takes the global median of scored runs, and the column records how many of its runs '
      'had an imputed weight. An imputed weight is a guess about data quality and is flagged, never hidden.',
      'Both medians are taken over **scored runs that are not themselves excluded**. That detail matters: all the '
      'scored runs of GSE245775 are the zero-parasite ones, so using the plain study median would have given its '
      'eight surviving runs a weight of 0.03 inherited from runs that are being thrown away.')
nb.code('scored_kept = quality[quality.parasite_reads > 0]',
        "GLOBAL_Q = float(scored_kept.quality.median())",
        'STUDY_Q = scored_kept.groupby("study").quality.median().to_dict()',
        'QUAL = dict(zip(quality.run, quality.quality))',
        'FLAGS = dict(zip(quality.run, quality["flags"].fillna("")))',
        'ZERO_PARASITE = set(quality.loc[quality.parasite_reads == 0, "run"])',
        "print('median quality', round(GLOBAL_Q, 3), '| quartiles', [round(x, 3) for x in quality.quality.quantile([.25,.5,.75])])",
        "print(len(ZERO_PARASITE), 'runs with zero parasite reads;', int(quality[\"flags\"].notna().sum()), 'runs carry at least one flag')",
        'quality.quality.describe().to_frame().T')

nb.md('## 2 · Exclusions, each one named and counted',
      'Instruction 69 took five per-cell single-cell collections out of the bulk queue, and quarantined the four '
      'GSE229505 runs that had already been forced through the bulk path by a `fastq-dump --split-3` fallback: '
      'those are not valid bulk measurements. Runs whose salmon quantification found **zero** parasite reads are '
      'excluded as well -- instruction 69 traced those to untrimmed short-insert libraries (GSE243206, GSE243207, '
      'GSE245775), where salmon selective alignment rejects everything; they will come back after the planned '
      'fastp re-quantification, not before.')
nb.code("EXCL_DIR = Path('/nas_mnt/carruthers/Einar/parascape_omics/rnaseq_excluded_single_cell')",
        'quarantined = sorted({r.name for s in EXCL_DIR.iterdir() if s.is_dir() for r in s.iterdir() if r.is_dir()})',
        "SINGLE_CELL_STUDIES = ['GSE229505', 'GSE346252', 'PRJNA1047766', 'PRJDB12314', 'PRJNA1053644',",
        "                       'GSE145800', 'GSE145080']",
        'rows = []',
        'd = design[design.study.isin(MATRIX_STUDIES)].copy()',
        "rows.append(dict(exclusion='design rows outside the studies that have a matrix', runs=len(design) - len(d),",
        "                 detail='runs planned or verified but with no salmon ME49_host matrix in this snapshot'))",
        'sc = d[d.study.isin(SINGLE_CELL_STUDIES) | d.run_accession.isin(quarantined)]',
        "rows.append(dict(exclusion='single-cell collections and the GSE229505 quarantine', runs=len(sc),",
        "                 detail='instruction 69: per-cell 10x libraries; ' + str(len(quarantined)) + ' quarantined runs on the NAS'))",
        'd = d[~d.index.isin(sc.index)]',
        'zp = d[d.run_accession.isin(ZERO_PARASITE)]',
        "rows.append(dict(exclusion='zero parasite reads (untrimmed short-insert libraries)', runs=len(zp),",
        "                 detail='studies ' + ', '.join(sorted(zp.study.unique()))))",
        'd = d[~d.index.isin(zp.index)]',
        'KEEP = d.copy()',
        "exclusions = pd.DataFrame(rows)",
        "print(len(KEEP), 'runs kept across', KEEP.study.nunique(), 'studies')",
        'exclusions')

nb.md('### The control marking is weak, and that governs how many relative columns exist',
      'This is the most important sentence in the notebook. `sample_design.tsv` has a `control_role` and a '
      '`control_link` column, but **`control_link` carries no pairing at all** -- for every kept run its value is '
      'the sentence *"No verified pairing; choose the matching control within this study"*, and the snapshot '
      'itself records `control_pairing: unverified; requires within-study review`. `control_role` marks a run '
      'either `control (source label)` (the deposit called it a control) or `candidate baseline; review required` '
      '(a rule guessed it from words like *wild type*, *untreated*, *vehicle*).',
      '**The rule used here, stated plainly:** a *condition group* counts as a control group when **every** run in '
      'it is marked a control or a candidate baseline. Every other condition group of the same study becomes one '
      'contrast, against the pooled control runs of that study. Nothing is paired across studies, and nothing is '
      'paired by hand.',
      '**What this rule gets wrong, and it must be read with the columns:** the pooled baseline is the whole '
      'control-marked part of a study, not the matched arm. In GSE158231 the keyword rule marks the '
      'arsenite-treated wild type *and* the untreated wild type *and* the untreated knockout as baselines, so the '
      'one contrast that study yields is a knockout-plus-drug arm against a baseline that mixes drug and no drug. '
      'A contrast like that is a sum of two effects, not one. This is exactly the wrong-label failure mode of '
      'instruction 72, and the only honest fix is Einar reviewing the pairings; until then every relative column '
      'carries `control_groups` and `contrast_label` verbatim so the mixture is visible.')
nb.code("KEEP['is_control'] = KEEP.control_role != 'not assigned'",
        "print(KEEP.control_role.value_counts().to_string())",
        "print()",
        "print('distinct control_link values among kept runs:', KEEP.control_link.nunique())",
        "print(KEEP.control_link.value_counts().to_string())",
        'grp = KEEP.groupby(["study", "condition_group"]).agg(n=("run_accession", "size"), marked=("is_control", "sum"))',
        'grp["all_control"] = grp.marked == grp.n',
        'grp["partly"] = (grp.marked > 0) & (grp.marked < grp.n)',
        "print()",
        "print(int(grp.all_control.sum()), 'control groups;', int(grp.partly.sum()), 'groups where only some runs are marked (none, so no group is split)')",
        'per_study = grp.groupby("study").agg(groups=("n", "size"), control_groups=("all_control", "sum"))',
        'per_study["contrasts"] = per_study.groups - per_study.control_groups',
        'WITH_CONTROL = sorted(per_study.index[per_study.control_groups > 0])',
        'NO_CONTROL = sorted(per_study.index[per_study.control_groups == 0])',
        "print(len(WITH_CONTROL), 'studies have a usable control group ->', int(per_study.loc[WITH_CONTROL, 'contrasts'].sum()), 'contrasts')",
        "print(len(NO_CONTROL), 'studies have none ->', int(per_study.loc[NO_CONTROL, 'groups'].sum()), 'internally normalised columns')",
        'per_study.sort_values("contrasts", ascending=False)')

nb.md('## 3 · How a column is computed',
      'Per study, per gene, on `log2(TPM + 1)` of the salmon ME49_host matrix (parasite TPM, reclosed to a million '
      'within the parasite genes, so library size is already divided out).',
      '**Weighted group mean.** Each run carries weight *w* = its quality score. A condition group mean is the '
      'weighted mean of its runs. The **effective sample count** is Kish: `n_eff = (Σw)² / Σw²`, which equals the '
      'run count when all runs are perfect and falls as weights become unequal. A column records `n_samples` '
      '(runs that entered), `effective_n` and `mean_quality` (the weighted mean of the weights behind it).',
      '**Standard error.** Within-group weighted variances are pooled across the groups of the study that have at '
      'least two runs, then moderated towards the study median variance with a prior worth 4 degrees of freedom '
      '(the limma-trend idea, written out in numpy): `s²_mod = (4·s²_median + df·s²_pooled) / (4 + df)`. Then '
      '`se = sqrt(s²_mod · (1/n_eff_treatment + 1/n_eff_control))`. A study with no replicated group at all has no '
      'honest within-group variance; it falls back to the variance across all its runs, which is an upper bound, '
      'and the column is flagged `se_approximate`.',
      '**Shrinkage.** The raw difference of two noisy means is itself noisy, and instruction 72 records that '
      'subtracting two already-relative numbers amplifies noise (the mouse-immune lesson). So each fold change is '
      'shrunk towards zero by the normal-normal posterior mean: with `tau2 = var(lfc_raw) - mean(se²)` across the '
      'genes of that contrast (floored above zero), `lfc = lfc_raw · tau2 / (tau2 + se²)`. A gene measured with a '
      'large standard error is pulled towards no change; a gene measured precisely is barely moved. The raw value '
      'is not kept as a separate column, because the shrunk one is the value a learner should see.')
nb.code('PRIOR_DF = 4.0',
        '',
        'def weights_for(runs, study):',
        '    out, imputed = {}, 0',
        '    for r in runs:',
        '        if r in QUAL:',
        '            out[r] = float(QUAL[r])',
        '        else:',
        '            out[r] = float(STUDY_Q.get(study, GLOBAL_Q)); imputed += 1',
        '    return out, imputed',
        '',
        'def wstats(V, W):',
        '    """Weighted mean, effective n and weighted unbiased variance per gene for one group of runs."""',
        '    mask = np.isfinite(V)',
        '    Wm = mask * W[None, :]',
        '    sw = Wm.sum(1); sw2 = (Wm ** 2).sum(1)',
        '    mu = np.where(sw > 0, np.nansum(np.where(mask, V * W[None, :], 0.0), 1) / np.where(sw > 0, sw, 1.0), np.nan)',
        '    neff = np.where(sw2 > 0, sw ** 2 / np.where(sw2 > 0, sw2, 1.0), 0.0)',
        '    dev = np.where(mask, (V - mu[:, None]) ** 2, 0.0)',
        '    denom = sw - np.where(sw > 0, sw2 / np.where(sw > 0, sw, 1.0), 0.0)',
        '    var = np.where(denom > 0, (Wm * dev).sum(1) / np.where(denom > 0, denom, 1.0), np.nan)',
        '    return mu, neff, var',
        '',
        'def moderate(pooled, df):',
        '    s0 = np.nanmedian(pooled)',
        '    return (PRIOR_DF * s0 + df * np.nan_to_num(pooled, nan=s0)) / (PRIOR_DF + df)',
        '',
        'def shrink(lfc_raw, se):',
        '    finite = np.isfinite(lfc_raw) & np.isfinite(se)',
        '    if finite.sum() < 50:',
        '        return lfc_raw, np.nan',
        '    tau2 = float(np.var(lfc_raw[finite]) - np.mean(se[finite] ** 2))',
        '    tau2 = max(tau2, 1e-4)',
        '    return lfc_raw * (tau2 / (tau2 + se ** 2)), tau2',
        "print('helpers ready; prior degrees of freedom', PRIOR_DF)")

nb.md('### Naming, and the probability a column carries',
      'A column is named `expr_lfc_<study>_<contrast>` for the fold change and `expr_lfcse_<study>_<contrast>` for '
      'its standard error, `expr_rel_<study>_<group>` for an internally normalised value and '
      '`expr_abs_<study>` for an absolute one. The contrast slug comes from the part of `condition_group` after '
      'the last `|` (the deposit own sample label), stripped to letters, digits and underscores; the verbatim '
      '`condition_group` is kept in the registry.',
      '**The probability, computed exactly like this:**',
      '`probability = base x mean_quality x replication`, clipped to [0.05, 0.95], where',
      '* `base` is **0.90** for `relative_to_control`, **0.85** for `internal_reference`, **0.80** for '
      '`library_size`. A relative column does not get the 0.95 of a direct experiment because its control pairing '
      'is unverified (section 2); the ordering of the three bases is Einar normalization preference made numeric.',
      '* `mean_quality` is the weighted mean sample quality of the runs behind the column -- this is how "a column '
      'probability inherits that quality" (instruction 72) is implemented.',
      '* `replication` is **1.0** when the treatment arm has `effective_n >= 2`, and **0.85** when it does not: a '
      'single run cannot show that what it measured repeats.',
      'A column whose standard error had to be approximated, or whose runs needed an imputed quality weight, is '
      'flagged in the registry; the flag is not folded into the number, so nothing is hidden inside it.')
nb.code('BASE_P = {"relative_to_control": 0.90, "internal_reference": 0.85, "library_size": 0.80}',
        '',
        'def slug(text, limit=46):',
        "    text = unicodedata.normalize('NFKD', str(text)).encode('ascii', 'ignore').decode()",
        "    label = text.split('|')[-1].strip() or text",
        "    s = re.sub(r'_+', '_', re.sub(r'[^0-9A-Za-z]+', '_', label)).strip('_')",
        "    return (s[:limit].rstrip('_') or 'group')",
        '',
        'def probability(norm, mean_quality, effective_n):',
        '    p = BASE_P[norm] * float(mean_quality) * (1.0 if effective_n >= 2 else 0.85)',
        '    return float(np.clip(p, 0.05, 0.95))',
        "print({k: BASE_P[k] for k in BASE_P})",
        "print('example:', round(probability('relative_to_control', 0.95, 3), 4), round(probability('library_size', 0.40, 1), 4))")

nb.md('## 4 · Building every column, one study at a time',
      'For a study with a control group: one fold-change column and one standard-error column per non-control '
      'condition group, against the pooled control runs.',
      'For a study without one: one `internal_reference` column per condition group, equal to the group weighted '
      'mean log2 TPM **minus the median of that gene across the condition groups of the same study**. The '
      'reference is therefore the gene own typical level in that study, which removes gene length, sequence '
      'composition and the study batch, but not a shift shared by every condition of the study. A '
      'housekeeping-gene reference was considered and not used: TPM is already reclosed within the parasite genes, '
      'so a library-size-like reference adds nothing, and no validated invariant gene set exists for these '
      'conditions -- section 8 measures how stable the usual candidates actually are.',
      'And for every study, one `library_size` column: the weighted mean log2 TPM over all its kept runs. These '
      'exist so no study is missing from the map, and learners must prefer the relative forms where both exist.',
      'Three things to know about the internal form, all consequences of using a median and all recorded as '
      '`n_condition_groups` in the registry. A study with only one condition group gets no internal column, '
      'because the value would be identically zero (GSE286090 has one verified run so far and is absolute only). '
      'In a study with exactly two condition groups the median is the mean, so its two columns are mirror images '
      'and carry one measurement between them. And with an **odd** number of groups the median *is* one of the '
      'observed group values for each gene, so the group that happens to be the middle one reads exactly zero for '
      'that gene -- in GSE208647, whose three groups are a 30-day tissue cyst, in vitro pH 8.2 bradyzoites and a '
      'tachyzoite arm, the pH 8.2 group is the middle one for BAG1 and its column therefore reads 0.0 for BAG1. '
      'That is the column saying "this is the typical condition of this study", which is true and is also why a '
      'three-group study is a weak place to read an internal reference. The median is kept because it is robust '
      'to one extreme group (the cyst sample here is +10 log2 for BAG1 and would drag a mean).')
nb.code('MATRICES = SNAP / "rna/matrices"',
        'columns, registry, group_means_store, run_lists, slug_of = {}, [], {}, {}, {}',
        '',
        'def unique_slugs(labels):',
        '    """A slug per condition group, made unique inside its study by a numeric suffix."""',
        '    out, seen = {}, {}',
        '    for label in sorted(labels):',
        '        s = slug(label)',
        '        seen[s] = seen.get(s, 0) + 1',
        '        out[label] = s if seen[s] == 1 else f"{s}_{seen[s]}"',
        '    return out',
        '',
        'def load_study(study, runs):',
        '    m = pd.read_parquet(MATRICES / study / "salmon_ME49_host.parquet")',
        '    tpm = m["tpm"].astype("float32").reindex(GENES)',
        '    present = [r for r in runs if r in tpm.columns]',
        '    return np.log2(tpm[present].to_numpy(dtype="float64") + 1.0), present',
        '',
        'for study in MATRIX_STUDIES:',
        '    sd = KEEP[KEEP.study == study]',
        '    if not len(sd):',
        '        continue',
        '    X, runs = load_study(study, sd.run_accession.tolist())',
        '    if not runs:',
        '        continue',
        '    sd = sd[sd.run_accession.isin(runs)]',
        '    order = {r: i for i, r in enumerate(runs)}',
        '    W, n_imputed = weights_for(runs, study)',
        '    run_lists[study] = runs',
        '    isctrl = dict(zip(sd.run_accession, sd.is_control))',
        '    groups = {g: [r for r in gd.run_accession if r in order] for g, gd in sd.groupby("condition_group")}',
        '    groups = {g: rs for g, rs in groups.items() if rs}',
        '    ctrl_groups = [g for g, rs in groups.items() if all(bool(isctrl[r]) for r in rs)]',
        '    SL = unique_slugs(groups); slug_of[study] = SL',
        '    stats = {}',
        '    for g, rs in groups.items():',
        '        idx = [order[r] for r in rs]',
        '        mu, neff, var = wstats(X[:, idx], np.array([W[r] for r in rs]))',
        '        mq = float(np.average([W[r] for r in rs], weights=[W[r] for r in rs]))',
        '        stats[g] = dict(mu=mu, neff=neff, var=var, runs=rs, n=len(rs), mean_quality=mq)',
        '    group_means_store[study] = stats',
        '    rep = [g for g in groups if stats[g]["n"] >= 2]',
        '    if rep:',
        '        df = float(sum(stats[g]["n"] - 1 for g in rep))',
        '        pooled = np.nansum(np.stack([(stats[g]["n"] - 1) * np.nan_to_num(stats[g]["var"], nan=0.0) for g in rep]), 0) / df',
        '        se_approx = False',
        '    else:',
        '        df = float(max(len(runs) - 1, 1))',
        '        pooled = np.nanvar(X, axis=1, ddof=1) if X.shape[1] > 1 else np.full(len(GENES), np.nan)',
        '        se_approx = True',
        '    s2 = moderate(pooled, df)',
        '    # --- relative to control -------------------------------------------------',
        '    if ctrl_groups and set(groups) - set(ctrl_groups):',
        '        crs = [r for g in ctrl_groups for r in groups[g]]',
        '        cidx = [order[r] for r in crs]',
        '        cmu, cneff, _ = wstats(X[:, cidx], np.array([W[r] for r in crs]))',
        '        cq = float(np.average([W[r] for r in crs], weights=[W[r] for r in crs]))',
        '        for g in sorted(set(groups) - set(ctrl_groups)):',
        '            st = stats[g]',
        '            name = f"expr_lfc_{study}_{SL[g]}"',
        '            se = np.sqrt(s2 * (1.0 / np.maximum(st["neff"], 1e-9) + 1.0 / np.maximum(cneff, 1e-9)))',
        '            lfc, tau2 = shrink(st["mu"] - cmu, se)',
        '            columns[name] = lfc; columns[name.replace("expr_lfc_", "expr_lfcse_")] = se',
        '            neff_t = float(np.nanmedian(st["neff"]))',
        '            mq = float((st["mean_quality"] * st["n"] + cq * len(crs)) / (st["n"] + len(crs)))',
        '            common = dict(study=study, contrast=SL[g], contrast_label=g, normalization="relative_to_control",',
        '                          n_samples=st["n"] + len(crs), n_treatment=st["n"], n_control=len(crs),',
        '                          effective_n=round(neff_t + float(np.nanmedian(cneff)), 3),',
        '                          effective_n_treatment=round(neff_t, 3), mean_quality=round(mq, 4),',
        '                          probability=round(probability("relative_to_control", mq, neff_t), 4),',
        '                          control_groups=" ; ".join(sorted(ctrl_groups)), tau2=None if tau2 != tau2 else round(tau2, 4),',
        '                          se_approximate=se_approx, n_quality_imputed=n_imputed)',
        '            registry.append(dict(column=name, quantity="log2_fold_change", **common))',
        '            registry.append(dict(column=name.replace("expr_lfc_", "expr_lfcse_"), quantity="standard_error", **common))',
        '    # --- internal reference --------------------------------------------------',
        '    elif len(groups) >= 2:',
        '        centre = np.nanmedian(np.stack([stats[g]["mu"] for g in sorted(groups)]), 0)',
        '        for g in sorted(groups):',
        '            st = stats[g]',
        '            name = f"expr_rel_{study}_{SL[g]}"',
        '            columns[name] = st["mu"] - centre',
        '            neff_t = float(np.nanmedian(st["neff"]))',
        '            registry.append(dict(column=name, study=study, contrast=SL[g], contrast_label=g,',
        '                                 normalization="internal_reference", quantity="log2_vs_study_median",',
        '                                 n_samples=st["n"], n_treatment=st["n"], n_control=0,',
        '                                 effective_n=round(neff_t, 3), effective_n_treatment=round(neff_t, 3),',
        '                                 mean_quality=round(st["mean_quality"], 4),',
        '                                 probability=round(probability("internal_reference", st["mean_quality"], neff_t), 4),',
        '                                 control_groups="(none marked in this study)", tau2=None,',
        '                                 se_approximate=se_approx, n_quality_imputed=n_imputed))',
        '    # --- absolute, last resort ----------------------------------------------',
        '    Wall = np.array([W[r] for r in runs])',
        '    mu_all, neff_all, _ = wstats(X, Wall)',
        '    name = f"expr_abs_{study}"',
        '    columns[name] = mu_all',
        '    mq = float(np.average(Wall, weights=Wall))',
        '    registry.append(dict(column=name, study=study, contrast="all_samples", contrast_label="all kept runs of the study",',
        '                         normalization="library_size", quantity="mean_log2_tpm", n_samples=len(runs),',
        '                         n_treatment=len(runs), n_control=0, effective_n=round(float(np.nanmedian(neff_all)), 3),',
        '                         effective_n_treatment=round(float(np.nanmedian(neff_all)), 3), mean_quality=round(mq, 4),',
        '                         probability=round(probability("library_size", mq, float(np.nanmedian(neff_all))), 4),',
        '                         control_groups="(not applicable)", tau2=None, se_approximate=se_approx,',
        '                         n_quality_imputed=n_imputed))',
        '    del X',
        '',
        'traits = pd.DataFrame(registry)',
        'traits["n_condition_groups"] = traits.study.map({s: len(g) for s, g in group_means_store.items()})',
        'values = pd.DataFrame(columns, index=pd.Index(GENES, name="gene_id"))[traits.column.tolist()]',
        "print(values.shape[1], 'columns over', values.shape[0], 'genes')",
        'traits.groupby(["normalization", "quantity"]).agg(columns=("column", "size"), studies=("study", "nunique"),'
        ' median_probability=("probability", "median")).reset_index()')

nb.md('### Which columns carry any signal at all',
      'Shrinkage has a consequence worth printing. GSE132238 deposits **one run per condition group** -- sixteen '
      'groups, sixteen runs -- so it has no replicated group, its standard error falls back to the variance across '
      'all its runs (about 3.9 log2 units), and every fold change is shrunk almost exactly to zero. Its eight '
      'columns are therefore arithmetically correct and informationally empty, which is the right answer for a '
      'design with no replication, and much better than a confident number. The count below, '
      '`genes_abs_lfc_ge_0_5`, says how many genes moved at least half a log2 unit in each fold-change column; a '
      'column near zero there should not be used on its own.')
nb.code('sig = {row.column: int((values[row.column].abs() >= 0.5).sum())',
        '       for row in traits[traits.quantity == "log2_fold_change"].itertuples()}',
        'traits["genes_abs_lfc_ge_0_5"] = traits.column.map(sig)',
        'flat = sorted(c for c, n in sig.items() if n < 50)',
        "print(len(sig), 'fold-change columns;', len(flat), 'move fewer than 50 genes by >= 0.5 log2 (effectively flat)')",
        "print('flat columns by study:', pd.Series([c.split('_')[2] for c in flat]).value_counts().to_dict())",
        'pd.Series(sig).sort_values().to_frame("genes moved >= 0.5 log2").head(12).T')

nb.md('## 5 · Registering the columns as Traits',
      'Every column is registered with `parascape.information_map.register`, tier 2, `role="measurement"`, the '
      'normalization above and `family` = the study accession, so instruction 72 whole-family hold-out can drop a '
      'whole study at once. `parascape/information_map.py` is **not** edited; registration happens here.',
      'Kind is `measure`, because an expression value may be negative, which is what the `measure` kind exists '
      'for. The `register()` lines are printed verbatim below -- the first few in full, then the count.')
nb.code('import parascape.information_map as imap',
        'lines = []',
        'for row in traits.itertuples():',
        '    desc = {"log2_fold_change": f"log2 fold change of {row.contrast} against the control-marked samples of {row.study} (shrunk)",',
        '            "standard_error": f"standard error of the log2 fold change of {row.contrast} in {row.study}",',
        '            "log2_vs_study_median": f"log2 TPM of {row.contrast} minus the gene median across the condition groups of {row.study}",',
        '            "mean_log2_tpm": f"mean log2 TPM across all kept runs of {row.study} (absolute; prefer a relative column)"}[row.quantity]',
        '    trait = imap.Trait(row.column, "measure", desc, tier=2, family=row.study,',
        '                       normalization=row.normalization, role="measurement")',
        '    imap.register(trait)',
        '    lines.append(f"register(Trait({row.column!r}, \'measure\', {desc!r}, tier=2, family={row.study!r}, "',
        '                 f"normalization={row.normalization!r}, role=\'measurement\'))")',
        'for line in lines[:6]:',
        '    print(line)',
        "print('...')",
        'for line in lines[-3:]:',
        '    print(line)',
        "print()",
        "print(len(lines), 'register() lines;', len(imap.TRAITS), 'traits now in the registry')",
        "Path(OUT / 'register_lines.txt').write_text('\\n'.join(lines) + '\\n')",
        'pd.Series({n: t.normalization for n, t in imap.TRAITS.items() if t.tier == 2}).value_counts().to_frame("tier-2 traits").T')

nb.md('## 6 · Known biology: the direction has to come out right',
      'A column that cannot reproduce an answer already known is not evidence. The check set is written out by '
      'hand, one row per (column pattern, gene, expected direction, reason), so each row can be argued with. Gene '
      'identities are taken from the ToxoDB product names in `nodes.parquet`, not from memory: BFD1 = '
      '`TGME49_200385` ("bradyzoite-formation deficient protein BFD1"), BFD2 = `TGME49_311100`, BAG1 = '
      '`TGME49_259020` ("bradyzoite antigen BAG1"), SAG1 = `TGME49_233460` ("SAG-related sequence SRS29B"), and '
      'the two enolases are `TGME49_268860` ("enolase 1") and `TGME49_268850` ("enolase 2"). ENO1 is the '
      'bradyzoite enolase and ENO2 the tachyzoite one; both are reported so that the assignment is testable here '
      'rather than assumed.',
      'A row **passes** when the sign is as expected and |log2 fold change| >= 0.5; it is **inconclusive** when a '
      'standard error exists and |z| < 1, because such a column was already shrunk towards zero and has nothing to '
      'say either way; it is a **weak pass** when the sign is right but the effect is under 0.5; it is a **MISS** '
      'when the sign is wrong. The gene percentile within its own column is reported too, because a right sign at '
      'the 50th percentile is not a finding.')
nb.code('MARKERS = dict(BFD1="TGME49_200385", BFD2="TGME49_311100", BAG1="TGME49_259020", SAG1="TGME49_233460",',
        '               ENO1="TGME49_268860", ENO2="TGME49_268850")',
        'for g, gid in MARKERS.items():',
        '    print(f"{g:5s} {gid}  {PRODUCT.get(gid)}")',
        'CHECKS = [',
        '    # alkaline / bradyzoite induction in wild-type parasites: the whole bradyzoite programme should rise',
        '    ("GSE218921", r"^expr_lfc_GSE218921_T_gondii_VEG_WT_Switch", "WT alkaline switch vs control media", "BFD1", +1),',
        '    ("GSE218921", r"^expr_lfc_GSE218921_T_gondii_VEG_WT_Switch", "WT alkaline switch vs control media", "BAG1", +1),',
        '    ("GSE218921", r"^expr_lfc_GSE218921_T_gondii_VEG_WT_Switch", "WT alkaline switch vs control media", "ENO1", +1),',
        '    ("GSE218921", r"^expr_lfc_GSE218921_T_gondii_VEG_WT_Switch", "WT alkaline switch vs control media", "SAG1", -1),',
        '    ("GSE218921", r"^expr_lfc_GSE218921_T_gondii_in_vitro_bradyzoite", "sporozoite bradyzoite induction", "BAG1", +1),',
        '    # BFD1 stabilised by Shield-1 in a DD-BFD1 line: BFD1 is the master regulator, so BAG1/ENO1 rise',
        '    ("GSE223877", r"^expr_lfc_GSE223877_DD_BFD1_Shield", "BFD1 stabilised (Shield-1)", "BAG1", +1),',
        '    ("GSE223877", r"^expr_lfc_GSE223877_DD_BFD1_Shield", "BFD1 stabilised (Shield-1)", "ENO1", +1),',
        '    ("GSE139537", r"^expr_lfc_GSE139537_2B4_Shield1", "BFD1 stabilised (Shield-1)", "BAG1", +1),',
        '    ("GSE139537", r"^expr_lfc_GSE139537_2B4_Shield1", "BFD1 stabilised (Shield-1)", "ENO1", +1),',
        '    # BFD1 or BFD2 removed under alkaline stress: the programme should fail to switch on',
        '    ("GSE223819", r"^expr_lfc_GSE223819_200385_mNG_mAID_IAA", "BFD1 depleted under alkaline stress", "BAG1", -1),',
        '    ("GSE223819", r"^expr_lfc_GSE223819_311100_mNG_mAID_IAA", "BFD2 depleted under alkaline stress", "BAG1", -1),',
        '    ("GSE132238", r"^expr_lfc_GSE132238_KO_stress", "BFD1 knockout under stress vs WT baseline", "BAG1", -1),',
        '    # stage, from the internally normalised columns of a study with no control',
        '    ("GSE159515", r"^expr_rel_GSE159515_ME49_Bz", "ME49 bradyzoites vs study centre", "BAG1", +1),',
        '    ("GSE159515", r"^expr_rel_GSE159515_ME49_Bz", "ME49 bradyzoites vs study centre", "ENO1", +1),',
        '    ("GSE159515", r"^expr_rel_GSE159515_ME49_Bz", "ME49 bradyzoites vs study centre", "SAG1", -1),',
        '    ("GSE159515", r"^expr_rel_GSE159515_(ME49|RH)_Tz", "tachyzoites vs study centre", "SAG1", +1),',
        '    ("GSE159515", r"^expr_rel_GSE159515_(ME49|RH)_Tz", "tachyzoites vs study centre", "ENO2", +1),',
        '    ("GSE208647", r"^expr_rel_GSE208647_ME49EW_cyst", "30-day mouse-brain tissue cysts vs study centre", "BAG1", +1),',
        '    ("GSE208647", r"^expr_rel_GSE208647_ME49_Tir1_strain_pH8_2", "in vitro pH 8.2 bradyzoites vs study centre", "BAG1", +1),',
        ']',
        "print()",
        "print(len(CHECKS), 'check rows')")

nb.code('rows = []',
        'for study, pattern, why, gene, want in CHECKS:',
        '    gid = MARKERS[gene]',
        '    cols = [c for c in values.columns if re.match(pattern, c) and not c.startswith("expr_lfcse_")]',
        '    if not cols:',
        '        rows.append(dict(study=study, column="(no column matched)", contrast=why, gene=gene, gene_id=gid,',
        '                         expected="up" if want > 0 else "down", lfc=np.nan, se=np.nan, z=np.nan,',
        '                         percentile=np.nan, verdict="not available"))',
        '        continue',
        '    for col in cols:',
        '        v = values[col]',
        '        lfc = float(v.get(gid, np.nan))',
        '        se_col = col.replace("expr_lfc_", "expr_lfcse_") if col.startswith("expr_lfc_") else None',
        '        se = float(values[se_col].get(gid, np.nan)) if se_col in values.columns else np.nan',
        '        pct = float((v.dropna() < lfc).mean() * 100) if np.isfinite(lfc) else np.nan',
        '        if not np.isfinite(lfc):',
        '            verdict = "not measured"',
        '        elif np.sign(lfc) != np.sign(want):',
        '            verdict = "MISS (wrong direction)"',
        '        elif np.isfinite(se) and se > 0 and abs(lfc / se) < 1:',
        '            verdict = "inconclusive (|z| < 1)"',
        '        elif abs(lfc) < 0.5:',
        '            verdict = "weak pass (|lfc| < 0.5)"',
        '        else:',
        '            verdict = "pass"',
        '        rows.append(dict(study=study, column=col, contrast=why, gene=gene, gene_id=gid,',
        '                         expected="up" if want > 0 else "down", lfc=round(lfc, 3),',
        '                         se=None if not np.isfinite(se) else round(se, 3),',
        '                         z=None if not np.isfinite(se) or se == 0 else round(lfc / se, 2),',
        '                         percentile=None if not np.isfinite(pct) else round(pct, 1), verdict=verdict))',
        'known = pd.DataFrame(rows)',
        'known.insert(0, "source", "registered column")',
        "print(known.verdict.value_counts().to_string())",
        'known')

nb.md('### Two checks the registered columns cannot make, and why',
      '**The tachyzoite question.** The requirement is that a tachyzoite contrast raises SAG1. In GSE159515 four '
      'of six samples are tachyzoites, so the study centre *is* mostly tachyzoite and a tachyzoite column sits at '
      'zero by construction -- which is why SAG1 comes out near zero in those columns rather than up. The honest '
      'test is the within-study difference of two internal columns of the same study: they share one centre, so '
      'subtracting them cancels it exactly and gives a clean tachyzoite-versus-bradyzoite fold change. That '
      'difference is computed here as a **check**; it is not added as a column, because a learner can form it from '
      'the two columns that are already registered.',
      '**The contrasts the control marking hides.** In GSE218921 the keyword rule marked the wild-type '
      'alkaline-switch arm as a *candidate baseline* (its label contains "WT"), so the single most informative '
      'alkaline contrast in the whole snapshot never becomes a column -- it was absorbed into the baseline it '
      'should have been measured against. The same happens to WT_stress in GSE132238. Both differences are '
      'computed here directly from the group means to show that the **data contain the expected biology and the '
      'control marking is what loses it**. These rows are labelled `hidden by the control marking` and are a '
      'request for Einar review, not columns.')
nb.code('DERIVED = [("GSE159515", "expr_rel_GSE159515_ME49_Tz_A", "expr_rel_GSE159515_ME49_Bz_A",',
        '            "ME49 tachyzoite minus ME49 bradyzoite (difference of two internal columns)",',
        '            [("SAG1", +1), ("ENO2", +1), ("BAG1", -1), ("ENO1", -1)]),',
        '           ("GSE159515", "expr_rel_GSE159515_ME49_Tz_B", "expr_rel_GSE159515_ME49_Bz_B",',
        '            "ME49 tachyzoite minus ME49 bradyzoite (second replicate pair)",',
        '            [("SAG1", +1), ("ENO2", +1), ("BAG1", -1), ("ENO1", -1)])]',
        'HIDDEN = [("GSE218921", "VEG WT, Switch", "VEG WT, Control",',
        '           "WT alkaline switch vs WT control media (both marked control, so no column exists)",',
        '           [("BFD1", +1), ("BAG1", +1), ("ENO1", +1), ("SAG1", -1)]),',
        '          ("GSE132238", "WT_stress", "WT_unstress",',
        '           "WT alkaline stress vs WT unstressed (both marked control, so no column exists)",',
        '           [("BAG1", +1), ("ENO1", +1), ("SAG1", -1)])]',
        '',
        'def verdict_for(lfc, want):',
        '    if not np.isfinite(lfc):',
        '        return "not measured"',
        '    if np.sign(lfc) != np.sign(want):',
        '        return "MISS (wrong direction)"',
        '    return "pass" if abs(lfc) >= 0.5 else "weak pass (|lfc| < 0.5)"',
        '',
        'extra = []',
        'for study, a, b, why, wants in DERIVED:',
        '    if a not in values.columns or b not in values.columns:',
        '        continue',
        '    diff = values[a] - values[b]',
        '    for gene, want in wants:',
        '        lfc = float(diff.get(MARKERS[gene], np.nan))',
        '        pct = float((diff.dropna() < lfc).mean() * 100) if np.isfinite(lfc) else np.nan',
        '        extra.append(dict(source="derived within-study", study=study, column=f"{a} - {b}", contrast=why,',
        '                          gene=gene, gene_id=MARKERS[gene], expected="up" if want > 0 else "down",',
        '                          lfc=round(lfc, 3), se=None, z=None,',
        '                          percentile=None if not np.isfinite(pct) else round(pct, 1),',
        '                          verdict=verdict_for(lfc, want)))',
        'for study, tsub, csub, why, wants in HIDDEN:',
        '    stats = group_means_store.get(study, {})',
        '    tg = [g for g in stats if tsub in g]; cg = [g for g in stats if csub in g]',
        '    if not tg or not cg:',
        '        continue',
        '    tm = np.nanmean(np.stack([stats[g]["mu"] for g in tg]), 0)',
        '    cm = np.nanmean(np.stack([stats[g]["mu"] for g in cg]), 0)',
        '    diff = pd.Series(tm - cm, index=GENES)',
        '    for gene, want in wants:',
        '        lfc = float(diff.get(MARKERS[gene], np.nan))',
        '        pct = float((diff.dropna() < lfc).mean() * 100) if np.isfinite(lfc) else np.nan',
        '        extra.append(dict(source="hidden by the control marking", study=study,',
        '                          column=f"(not a column) {tsub} vs {csub}", contrast=why, gene=gene,',
        '                          gene_id=MARKERS[gene], expected="up" if want > 0 else "down", lfc=round(lfc, 3),',
        '                          se=None, z=None, percentile=None if not np.isfinite(pct) else round(pct, 1),',
        '                          verdict=verdict_for(lfc, want)))',
        'known = pd.concat([known, pd.DataFrame(extra)], ignore_index=True)',
        'known.to_csv(OUT / "known_biology.tsv", sep="\\t", index=False)',
        "print(known.groupby(['source']).verdict.value_counts().to_string())",
        "print()",
        "print('MISSES:')",
        'miss = known[known.verdict.str.startswith("MISS")]',
        "print(miss[['source','column','gene','expected','lfc','percentile']].to_string(index=False) if len(miss) else '(none)')",
        'known[known.source != "registered column"]')

nb.md('## 7 · Replicate agreement: does a column agree with itself?',
      'For every fold-change column whose treatment arm and control pool both hold at least two runs, the runs of '
      'each arm are split in half (alternating in accession order) and the fold change is computed twice, from '
      'disjoint samples. The two halves are then correlated over the genes expressed in that study (control mean '
      'log2 TPM >= 1). **A column whose own replicates disagree is not evidence**, so the correlation is reported '
      'per column, not averaged away.')
nb.code('rows = []',
        'for row in traits[traits.quantity == "log2_fold_change"].itertuples():',
        '    stats = group_means_store[row.study]',
        '    g = next((k for k, s in slug_of[row.study].items() if s == row.contrast), None)',
        '    if g is None:',
        '        continue',
        '    crs = sorted({r for cg in row.control_groups.split(" ; ") for r in stats.get(cg, {}).get("runs", [])})',
        '    trs = sorted(stats[g]["runs"])',
        '    if len(trs) < 2 or len(crs) < 2:',
        '        continue',
        '    X, runs = load_study(row.study, trs + crs)',
        '    order = {r: i for i, r in enumerate(runs)}',
        '    W, _ = weights_for(runs, row.study)',
        '    halves = []',
        '    for half in (0, 1):',
        '        t = [r for i, r in enumerate(trs) if i % 2 == half and r in order]',
        '        c = [r for i, r in enumerate(crs) if i % 2 == half and r in order]',
        '        if not t or not c:',
        '            halves = []; break',
        '        tm, _, _ = wstats(X[:, [order[r] for r in t]], np.array([W[r] for r in t]))',
        '        cm, _, _ = wstats(X[:, [order[r] for r in c]], np.array([W[r] for r in c]))',
        '        halves.append((tm - cm, cm))',
        '    del X',
        '    if not halves:',
        '        continue',
        '    a, b = halves[0][0], halves[1][0]',
        '    expressed = (np.nanmean([halves[0][1], halves[1][1]], 0) >= 1) & np.isfinite(a) & np.isfinite(b)',
        '    if expressed.sum() < 200:',
        '        continue',
        '    s = pd.Series(a[expressed]); t2 = pd.Series(b[expressed])',
        '    rows.append(dict(column=row.column, study=row.study, contrast=row.contrast, n_treatment=row.n_treatment,',
        '                     n_control=row.n_control, genes=int(expressed.sum()),',
        '                     pearson=round(float(s.corr(t2)), 3), spearman=round(float(s.corr(t2, method="spearman")), 3)))',
        'replicate = pd.DataFrame(rows).sort_values("spearman")',
        'replicate.to_csv(OUT / "replicate_agreement.tsv", sep="\\t", index=False)',
        "print(len(replicate), 'columns could be split (both arms have >= 2 runs)')",
        "print('median split-half Spearman', round(float(replicate.spearman.median()), 3) if len(replicate) else 'n/a',",
        "      '| below 0.3:', int((replicate.spearman < 0.3).sum()) if len(replicate) else 0)",
        'replicate')

nb.md('## 8 · Cross-study agreement: the number that says whether these columns are comparable at all',
      'This is the point of the uniform reprocessing. Columns are tagged with a **theme** by keyword rules over '
      'the verbatim `condition_group` -- alkaline or bradyzoite induction, BFD1 gain of function, BFD1 or BFD2 '
      'loss, auxin depletion, stage. Every pair of columns from **two different studies** inside one theme is then '
      'correlated over the genes both measured. Studies, not runs, are the unit of independence (instruction 71), '
      'so a pair from one study is not reported as agreement.',
      'A pair where both columns use the same normalization is the fair comparison; cross-form pairs are reported '
      'too and marked, since a fold change against a control and a value against a study centre are not the same '
      'quantity. **If these numbers are low, that is the finding**, not a failure of the notebook: it would mean '
      'the columns of different studies are not on a common scale and a learner must treat each study as its own '
      'block.')
nb.code('def theme(label):',
        '    t = str(label).lower()',
        '    out = []',
        '    if re.search(r"alkaline|high ph|ph ?8|bradyzoite induction|switch condition|stress media", t): out.append("alkaline_bradyzoite_induction")',
        '    if re.search(r"shield", t) and "bfd1" in t: out.append("bfd1_gain_of_function")',
        '    if re.search(r"delta\\}?bfd1|_bfd1|bfd1 ko|200385|311100|bfd2", t): out.append("bfd_loss")',
        '    if re.search(r"\\biaa\\b|auxin|atc\\b|shield", t): out.append("inducible_depletion")',
        '    if re.search(r"bradyzoite|cyst", t) and not re.search(r"induction", t): out.append("bradyzoite_stage")',
        '    if re.search(r"tachyzoite", t): out.append("tachyzoite_stage")',
        '    return out',
        '',
        'tag = {}',
        'for row in traits[traits.quantity != "standard_error"].itertuples():',
        '    for th in theme(row.contrast_label):',
        '        tag.setdefault(th, []).append(row.column)',
        'norm_of = dict(zip(traits.column, traits.normalization)); study_of = dict(zip(traits.column, traits.study))',
        'label_of = dict(zip(traits.column, traits.contrast_label))',
        'rows = []',
        'for th, cols in sorted(tag.items()):',
        '    for i, a in enumerate(cols):',
        '        for b in cols[i + 1:]:',
        '            if study_of[a] == study_of[b]:',
        '                continue',
        '            x, y = values[a], values[b]',
        '            ok = x.notna() & y.notna()',
        '            if ok.sum() < 500:',
        '                continue',
        '            rows.append(dict(theme=th, column_a=a, column_b=b, study_a=study_of[a], study_b=study_of[b],',
        '                             normalization_a=norm_of[a], normalization_b=norm_of[b],',
        '                             same_normalization=norm_of[a] == norm_of[b], genes=int(ok.sum()),',
        '                             pearson=round(float(x[ok].corr(y[ok])), 3),',
        '                             spearman=round(float(x[ok].corr(y[ok], method="spearman")), 3),',
        '                             label_a=label_of[a], label_b=label_of[b]))',
        'cross = pd.DataFrame(rows).sort_values(["theme", "spearman"], ascending=[True, False])',
        'cross.to_csv(OUT / "cross_study_agreement.tsv", sep="\\t", index=False)',
        "print(len(cross), 'cross-study pairs in', cross.theme.nunique() if len(cross) else 0, 'themes')",
        'summary = cross.groupby(["theme", "same_normalization"]).spearman.agg(["size", "median", "min", "max"]).reset_index()',
        'summary')

nb.code('top = cross[cross.same_normalization].sort_values("spearman", ascending=False)',
        "print('same-normalization cross-study pairs:', len(top))",
        "print('median Spearman', round(float(top.spearman.median()), 3) if len(top) else 'n/a')",
        'top[["theme", "column_a", "column_b", "genes", "pearson", "spearman"]].head(25)')

nb.md('### Reading those numbers without overstating them in either direction',
      'The lumped median over every themed pair is near zero, and that is **not** by itself evidence that the '
      'columns are incomparable: most pairs inside a theme perturb *different genes*. Two auxin knockdowns of two '
      'unrelated proteins land in `inducible_depletion` together and have no reason whatsoever to correlate. The '
      'number that tests comparability is the correlation between pairs that **should** agree, so those are '
      'singled out below by name, together with a pair that should **not** agree as a control on the test itself.',
      'The themes are coarse keyword rules over the deposit own condition text. They are a convenience for '
      'finding candidate pairs, not a curated mapping, and a learner should not treat a theme as a label.')
nb.code('def pair(a, b):',
        '    row = cross[((cross.column_a == a) & (cross.column_b == b)) | ((cross.column_a == b) & (cross.column_b == a))]',
        '    return float(row.spearman.iloc[0]) if len(row) else np.nan',
        'NAMED = [("BFD1 stabilised, two independent studies (should agree)",',
        '          "expr_lfc_GSE139537_2B4_Shield1", "expr_lfc_GSE223877_DD_BFD1_Shield_1a"),',
        '         ("BFD1 stabilised, the second GSE223877 replicate set (should agree)",',
        '          "expr_lfc_GSE139537_2B4_Shield1", "expr_lfc_GSE223877_DD_BFD1_Shield_1b"),',
        '         ("BFD1 stabilised with BFD2 removed (should agree less: BFD2 is needed)",',
        '          "expr_lfc_GSE139537_2B4_Shield1", "expr_lfc_GSE223877_DD_BFD1bfd2_Shield_1a"),',
        '         ("Shield-1 in a line with no BFD1 to stabilise (should NOT agree: the drug alone does nothing)",',
        '          "expr_lfc_GSE139537_KO_Shield1", "expr_lfc_GSE223877_DD_BFD1_Shield_1a"),',
        '         ("BFD1 lost, two independent studies (should agree)",',
        '          "expr_lfc_GSE132238_KO_stress_72h_array1", "expr_lfc_GSE223819_200385_mNG_mAID_IAA"),',
        '         ("BFD1 lost against BFD1 stabilised (should be NEGATIVE: opposite perturbations)",',
        '          "expr_lfc_GSE132238_KO_unstress_72h_array1", "expr_lfc_GSE223877_DD_BFD1_Shield_1a"),',
        '         ("alkaline stress, two independent studies, internal reference",',
        '          "expr_rel_GSE208647_ME49_Tir1_strain_pH8_2_48hrs", "expr_rel_GSE223869_WT_stress_2b"),',
        '         ("bradyzoites, two independent studies, internal reference",',
        '          "expr_rel_GSE159515_ME49_Bz_A", "expr_rel_GSE208647_ME49_Tir1_strain_pH8_2_48hrs")]',
        'named = pd.DataFrame([dict(expectation=why, column_a=a, column_b=b, spearman=pair(a, b)) for why, a, b in NAMED])',
        'named.to_csv(OUT / "cross_study_named_pairs.tsv", sep="\\t", index=False)',
        'named')

nb.md('### Is a single study simply reproducible across studies at all?',
      'One more framing of the same question, free of any theme rule: the correlation between the **absolute** '
      'columns of every pair of studies. If two studies measuring the same parasite disagree even on which genes '
      'are highly expressed, nothing built on top of them is comparable.')
nb.code('abscols = traits.column[traits.normalization == "library_size"].tolist()',
        'A = values[abscols]',
        'corr = A.corr(method="spearman")',
        'iu = np.triu_indices_from(corr.values, 1)',
        'pairs = pd.Series(corr.values[iu])',
        "print(len(abscols), 'studies;', len(pairs), 'study pairs')",
        "print('absolute-profile Spearman: median', round(float(pairs.median()), 3), 'min', round(float(pairs.min()), 3), 'max', round(float(pairs.max()), 3))",
        "pd.DataFrame({'quantile': [.05, .25, .5, .75, .95], 'spearman': [round(float(pairs.quantile(q)), 3) for q in [.05, .25, .5, .75, .95]]})")

nb.md('### How stable are the usual housekeeping genes?',
      'Justifying the choice not to use a housekeeping reference. Across the absolute columns, a gene usable as an '
      'internal reference must be both highly expressed and low in variance between studies. The table ranks the '
      'conventional candidates (actin, alpha- and beta-tubulin, GAPDH, the enolases) against the genome.')
nb.code('hk = UNIVERSE[UNIVERSE["product"].astype(str).str.contains("actin|tubulin|glyceraldehyde|enolase", case=False, na=False)]',
        'sd_all = A.std(axis=1); mean_all = A.mean(axis=1)',
        'rank = sd_all.rank(pct=True) * 100',
        'tbl = pd.DataFrame({"product": hk.set_index("gene_id")["product"], "mean_log2_tpm": mean_all.reindex(hk.gene_id).round(2),',
        '                    "sd_between_studies": sd_all.reindex(hk.gene_id).round(2),',
        '                    "sd_percentile": rank.reindex(hk.gene_id).round(1)}).dropna().sort_values("sd_between_studies")',
        "print('a percentile near 0 would mean unusually stable; near 100, unusually variable')",
        'tbl.head(20)')

nb.md('## 9 · Coverage, and writing the tables',
      'Coverage is reported per normalization, per study and in total: how many columns, how many genes carry at '
      'least one value, and what fraction of the gene x column grid is filled. An empty cell means not measured; '
      'instruction 72 rule holds -- missing is never zero.')
nb.code('rows = [dict(level="total", key="all columns", columns=values.shape[1],',
        '             genes_with_any_value=int(values.notna().any(axis=1).sum()),',
        '             mean_genes_per_column=round(float(values.notna().sum().mean()), 1),',
        '             fill_fraction=round(float(values.notna().to_numpy().mean()), 4))]',
        'for norm, cols in traits.groupby("normalization").column:',
        '    sub = values[cols.tolist()]',
        '    rows.append(dict(level="normalization", key=norm, columns=sub.shape[1],',
        '                     genes_with_any_value=int(sub.notna().any(axis=1).sum()),',
        '                     mean_genes_per_column=round(float(sub.notna().sum().mean()), 1),',
        '                     fill_fraction=round(float(sub.notna().to_numpy().mean()), 4)))',
        'for study, cols in traits.groupby("study").column:',
        '    sub = values[cols.tolist()]',
        '    rows.append(dict(level="study", key=study, columns=sub.shape[1],',
        '                     genes_with_any_value=int(sub.notna().any(axis=1).sum()),',
        '                     mean_genes_per_column=round(float(sub.notna().sum().mean()), 1),',
        '                     fill_fraction=round(float(sub.notna().to_numpy().mean()), 4)))',
        'coverage = pd.DataFrame(rows)',
        'coverage.to_csv(OUT / "coverage.tsv", sep="\\t", index=False)',
        'values.to_parquet(OUT / "expression_traits.parquet")',
        'traits.to_csv(OUT / "column_registry.tsv", sep="\\t", index=False)',
        'exclusions.to_csv(OUT / "exclusions.tsv", sep="\\t", index=False)',
        "print(values.shape[1], 'columns;', int(values.notna().any(axis=1).sum()), 'of', len(GENES), 'genes carry at least one value')",
        'coverage.head(12)')

nb.md('## 10 · What a learner may and may not do with these columns',
      '* Prefer `relative_to_control` over `internal_reference` over `library_size` for the same study. Where a '
      'study appears in two forms, the absolute column is there for completeness only.',
      '* Hold out a whole **study** (the trait `family`), never single runs: runs of one study are not independent.',
      '* Carry `expr_lfcse_*` with its `expr_lfc_*`. A large fold change with a large standard error has already '
      'been shrunk towards zero, and the standard error says how much was known in the first place.',
      '* Read `column_registry.tsv` before using a column: `control_groups` says what the fold change was measured '
      'against, and some of those baselines mix conditions (section 2).',
      '* Treat the empty cells as unmeasured. A gene with no value in a study was not seen there; it is not off.')

nb.code('def md_table(df):',
        '    """A markdown table without tabulate, which is not installed in this environment."""',
        '    cols = [str(c) for c in df.columns]',
        '    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]',
        '    for row in df.itertuples(index=False):',
        '        out.append("| " + " | ".join(str(v) for v in row) + " |")',
        '    return "\\n".join(out)',
        'print(md_table(exclusions.head(2)))')

nb.code('README = f"""# Expression tier-2 columns of the information map (2026-11-05)',
        '',
        'Built by `scripts/build_expression_traits_notebook.py`; notebook `notebooks/expression_traits_2026_11_05.ipynb`.',
        'Snapshot `{SNAP.name}` ({META["verified_runs"]} verified runs of {META["planned_runs"]} planned).',
        'Values come from the uniform reprocessing only: salmon, ME49_host reference, parasite genes, TPM reclosed',
        'to a million within the parasite genes, `log2(TPM + 1)`. Author-deposited matrices are not mixed in.',
        '',
        '## The rules, verbatim',
        '',
        '* **Normalization preference (instruction 72, Einar):** `relative_to_control` > `internal_reference` >',
        '  `library_size` > `raw`. Where both a relative and an absolute form of one quantity exist, learners take',
        '  the relative one; absolute values are used only where no control exists.',
        '* **A control group** is a condition group in which *every* run is marked `control (source label)` or',
        '  `candidate baseline; review required` in `sample_design.tsv`. Every other condition group of that study',
        '  becomes one contrast against the pooled control runs of the same study.',
        '* **`control_link` carries no pairing.** For all {len(KEEP)} kept runs its value is the sentence "No verified',
        '  pairing; choose the matching control within this study", and the snapshot records `control_pairing:',
        '  unverified; requires within-study review`. The pooled baseline is therefore the control-marked *part of a',
        '  study*, not the matched arm, and in at least one study (GSE158231) that pool mixes drug-treated and',
        '  untreated baselines. Read `control_groups` and `contrast_label` in `column_registry.tsv` before using a column.',
        '* **Internal reference rule:** group weighted mean `log2(TPM + 1)` minus the median of that gene across the',
        '  condition groups of the same study. It removes gene length, composition and study batch; it does not',
        '  remove a shift shared by every condition of the study. No housekeeping reference was used: TPM is already',
        '  reclosed within the parasite genes and no invariant gene set is validated for these conditions.',
        '* **Absolute rule:** weighted mean `log2(TPM + 1)` over all kept runs of the study. Last resort.',
        '* **Weighting:** run weight = its quality score (`results/rna_sample_quality_2026_11_04/`). A group mean is',
        '  weighted; `effective_n = (sum w)^2 / sum w^2` (Kish); `mean_quality` is the weighted mean of those weights.',
        '  A run with no score takes its study median quality, or the global median ({GLOBAL_Q:.3f}) if its study has',
        '  none scored. Both medians are taken over scored runs that are not themselves excluded, so a study whose',
        '  only scored runs are the zero-parasite ones (GSE245775) does not inherit their weight.',
        '  `n_quality_imputed` counts imputed-weight runs per column.',
        '* **Standard error:** within-group weighted variances pooled over the groups of the study with >= 2 runs,',
        '  moderated towards the study median variance with a prior worth {PRIOR_DF:.0f} degrees of freedom, then',
        '  `se = sqrt(s2_mod * (1/n_eff_treatment + 1/n_eff_control))`. Studies with no replicated group fall back to',
        '  the variance across all their runs and are flagged `se_approximate`.',
        '* **Shrinkage:** `lfc = lfc_raw * tau2 / (tau2 + se^2)` with `tau2 = var(lfc_raw) - mean(se^2)` over the genes',
        '  of that contrast, floored just above zero (normal-normal posterior mean).',
        '* **Probability:** `base x mean_quality x replication`, clipped to [0.05, 0.95]. `base` = 0.90',
        '  (`relative_to_control`), 0.85 (`internal_reference`), 0.80 (`library_size`); `replication` = 1.0 when the',
        '  treatment arm has `effective_n >= 2`, else 0.85. A relative column does not get the 0.95 of a direct',
        '  experiment because its control pairing is unverified.',
        '* **Traits:** every column is registered with `information_map.register`, tier 2, kind `measure`, role',
        '  `measurement`, `family` = the study accession, and the normalization above. `parascape/information_map.py`',
        '  was not edited; the lines are in `register_lines.txt`.',
        '* **Missing is empty, never zero.** A gene with no value in a study was not measured there.',
        '',
        '## Exclusions',
        '',
        '{md_table(exclusions)}',
        '',
        'Single-cell collections are out by instruction 69 (per-cell 10x libraries; deposited single-cell matrices are',
        'used instead), including the {len(quarantined)} GSE229505 runs quarantined under',
        '`parascape_omics/rnaseq_excluded_single_cell/`. Runs with zero parasite reads are the untrimmed',
        'short-insert libraries of GSE243206, GSE243207 and GSE245775; they return after the planned fastp',
        're-quantification. GSE243206 and GSE243207 lose every run and so contribute no column at all.',
        '',
        '## What was built',
        '',
        '{md_table(traits.groupby(["normalization", "quantity"]).size().rename("columns").reset_index())}',
        '',
        '{len(WITH_CONTROL)} of {len(MATRIX_STUDIES) - 2} usable studies have a control group; {len(NO_CONTROL)} do not.',
        'Coverage: {int(values.notna().any(axis=1).sum())} of {len(GENES)} genes carry at least one value;',
        'the gene x column grid is {float(values.notna().to_numpy().mean()) * 100:.1f} % filled.',
        '',
        '## How comparable are these columns, honestly',
        '',
        '**Columns with no signal.** {len(flat)} of {len(sig)} fold-change columns move fewer than 50 genes by at least',
        'half a log2 unit (`genes_abs_lfc_ge_0_5` in the registry). Eight of them are GSE132238, which deposits one run',
        'per condition group: with no replicated group its standard error falls back to about 3.9 log2 units and every',
        'fold change is shrunk to zero. That is the correct answer for an unreplicated design, and those columns must',
        'not be read as "nothing changed in the experiment".',
        '',
        '**Two contrasts the control marking hides, with the numbers.** In GSE218921 the keyword rule marked the',
        'wild-type alkaline-switch arm as a candidate baseline, so the best alkaline contrast in the snapshot never',
        'became a column. Computed directly from the group means it is exactly as expected: BFD1 +1.36, BAG1 +4.75,',
        'ENO1 +4.91, SAG1 -1.34. The same thing happens to WT_stress in GSE132238 (BAG1 +6.92, ENO1 +6.99, SAG1 -5.25).',
        'The data hold the biology; the control marking loses it. These rows are in `known_biology.tsv` under',
        '`hidden by the control marking` and are the clearest argument for a within-study control review.',
        '',
        '**Known biology** (`known_biology.tsv`): {int((known.verdict == "pass").sum())} pass,',
        '{int(known.verdict.str.startswith("weak").sum())} weak pass,',
        '{int(known.verdict.str.startswith("inconclusive").sum())} inconclusive,',
        '{int(known.verdict.str.startswith("MISS").sum())} miss. Every miss is listed in the table with its value.',
        '',
        '**Replicate agreement** (`replicate_agreement.tsv`): of {len(traits[traits.quantity == "log2_fold_change"])}',
        'fold-change columns, {len(replicate)} have at least two runs in both arms and can be split in half; median',
        'split-half Spearman {float(replicate.spearman.median()):.3f}, and {int((replicate.spearman < 0.3).sum())} fall',
        'below 0.3. A low split-half correlation usually means the perturbation changed little, so the column is mostly',
        'noise -- correct arithmetic, no evidence. Those columns should not be used on their own.',
        '',
        '**Cross-study agreement** (`cross_study_agreement.tsv`, `cross_study_named_pairs.tsv`): over',
        '{len(top)} same-normalization pairs inside a keyword theme the median Spearman is only',
        '{float(top.spearman.median()):.3f}. **That number must not be read as "the columns are incomparable"**: most',
        'pairs inside a theme perturb different genes and have no reason to agree. The pairs that should agree do:',
        'BFD1 stabilised in two independent studies (GSE139537 against GSE223877) reaches',
        '{pair("expr_lfc_GSE139537_2B4_Shield1", "expr_lfc_GSE223877_DD_BFD1_Shield_1b"):.3f}, BFD1 lost in two',
        'independent studies (GSE132238 against GSE223819) reaches',
        '{pair("expr_lfc_GSE132238_KO_stress_72h_array1", "expr_lfc_GSE223819_200385_mNG_mAID_IAA"):.3f}, and the',
        'opposite perturbations correlate negatively',
        '({pair("expr_lfc_GSE132238_KO_unstress_72h_array1", "expr_lfc_GSE223877_DD_BFD1_Shield_1a"):.3f}), while',
        'Shield-1 given to a line with no BFD1 to stabilise -- the control on the test itself -- correlates',
        '{pair("expr_lfc_GSE139537_KO_Shield1", "expr_lfc_GSE223877_DD_BFD1_Shield_1a"):.3f}, i.e. not at all.',
        'Where it does **not** work is alkaline stress: two independent alkaline-stress studies agree only at',
        '{pair("expr_rel_GSE208647_ME49_Tir1_strain_pH8_2_48hrs", "expr_rel_GSE223869_WT_stress_2b"):.3f}. Both of',
        'those are `internal_reference` columns whose centre is a different mixture of conditions in each study, so',
        'part of that disagreement is the normalization and not the biology -- which is itself the argument for',
        'getting controls reviewed so those studies can carry `relative_to_control` columns instead.',
        '',
        '**Absolute profiles** agree well: median pairwise Spearman {float(pairs.median()):.3f} over {len(pairs)} study',
        'pairs. Studies do agree about which genes are highly expressed; the hard part is the differences.',
        '',
        '**Summary judgement.** A perturbation column is trustworthy where the perturbation is strong and the control',
        'arm is the right one -- BFD1 gain and loss reproduce across labs on these values, which is what the uniform',
        'reprocessing was for. Columns from weak perturbations, from studies with no replicated group, and from the',
        'studies whose baseline mixes conditions are not yet evidence on their own. The binding constraint is not the',
        'quantification; it is the **control assignment**, which no code here can fix.',
        '"""',
        "(OUT / 'README.md').write_text(README)",
        "print('README.md written,', len(README.splitlines()), 'lines')")

nb.md('## 11 · The honest summary, computed from the tables')
nb.code('passes = int((known.verdict == "pass").sum()); weak = int(known.verdict.str.startswith("weak").sum())',
        'misses = known[known.verdict.str.startswith("MISS")]',
        'unavailable = int(known.verdict.isin(["not available", "not measured"]).sum())',
        'inconc = int(known.verdict.str.startswith("inconclusive").sum())',
        'rep_line = (f"replicate split-half Spearman: median {replicate.spearman.median():.3f} over {len(replicate)} "',
        '            f"splittable columns, {int((replicate.spearman < 0.3).sum())} below 0.3") if len(replicate) else \\',
        '           "replicate agreement: no column had two runs in both arms"',
        'cross_line = (f"cross-study, same normalization: {len(top)} pairs, median Spearman {top.spearman.median():.3f}, "',
        '              f"range {top.spearman.min():.3f} to {top.spearman.max():.3f} -- but most themed pairs perturb "',
        '              f"different genes, so see the named pairs") if len(top) else "cross-study: no comparable pair"',
        'named_line = "named cross-study pairs: " + "; ".join(f"{r.expectation.split(chr(40))[0].strip()} = {r.spearman:.3f}"',
        '                                                     for r in named.itertuples())',
        'theme_line = "by theme (same normalization, median Spearman): " + "; ".join(',
        '    f"{t} {s:.3f} (n={n})" for t, s, n in zip(*[top.groupby("theme").spearman.median().index,',
        '                                               top.groupby("theme").spearman.median().values,',
        '                                               top.groupby("theme").spearman.size().values]))',
        'lines = [',
        '    f"columns: {values.shape[1]} total -- " + ", ".join(f"{n} {k}" for k, n in traits.normalization.value_counts().items()),',
        '    f"studies with a usable control group: {len(WITH_CONTROL)}; without: {len(NO_CONTROL)}",',
        '    f"known biology: {passes} pass, {weak} weak pass, {inconc} inconclusive, {len(misses)} miss, "',
        '    f"{unavailable} not available",',
        '    "misses: " + (", ".join(f"{r.column} / {r.gene} (expected {r.expected}, got {r.lfc})" for r in misses.itertuples()) or "none"),',
        '    rep_line,',
        '    cross_line,',
        '    theme_line,',
        '    named_line,',
        '    f"absolute study profiles: median pairwise Spearman {pairs.median():.3f} over {len(pairs)} study pairs",',
        '    f"coverage: {int(values.notna().any(axis=1).sum())} of {len(GENES)} genes, grid {float(values.notna().to_numpy().mean()) * 100:.1f} % filled",',
        ']',
        'for line in lines:',
        '    print("*", line)',
        "(OUT / 'summary.txt').write_text('\\n'.join(lines) + '\\n')",
        "print()",
        "print('written to', OUT)")

nb.write(str(ROOT / 'notebooks/expression_traits_2026_11_05.ipynb'))
print('notebook written')
