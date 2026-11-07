"""Condition fingerprints (instruction 75), first two axes, as an executed notebook.

A fingerprint is built from a study that samples a condition axis in distinct samples, then every other
bulk RNA-seq run is scored against it -- annotating datasets with metadata their authors never published.
Scores are within-sample ranks (never absolute values, so a fingerprint learned in one study does not just
score that study's batch), the output is a mixture rather than a label, and every axis is validated on a
held-out study that states its own condition before it ships. An annotation is INFERRED: it carries
`filled:<fingerprint>@<version>` and the study it depends on, and is never evidence for a gene trait.

    python scripts/build_condition_fingerprints_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_06"
VERSION = "v1"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"condition_fingerprints_{DATE}")
NB = os.path.join(REPO, "notebooks", f"condition_fingerprints_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    nb = ExecutedNotebook("Condition fingerprints: cell cycle and bradyzoite conversion",
                          namespace={"REPO": REPO, "OUT": OUT, "VERSION": VERSION})
    nb.md(
        "**What this is.** Two condition fingerprints, each built on one study and used to annotate every "
        "bulk RNA-seq run of the uniform reprocessing with something its authors did not publish: where "
        "in the cell cycle the culture sits, and how far it has converted to bradyzoites (instruction 75).",
        "**How it scores.** Within each sample, genes are turned into percentile ranks; a fingerprint "
        "score is the mean rank of its genes. Ranks make a fingerprint usable across studies, platforms "
        "and normalizations, because nothing study-specific survives the ranking except the order.",
        "**How it is checked.** Each axis is tested on studies it was not built on and that state their "
        "own condition. An axis that fails is published as exploratory or not at all.",
        "**What an annotation is.** Inferred, never observed: `filled:<fingerprint>@" + VERSION + "`, "
        "with the study it depends on. It must never become evidence for a gene, or the model learns "
        f"from its own guess. Built {_dt.date.today().isoformat()}.")

    nb.code(
        "import json, os, re",
        "from pathlib import Path",
        "import numpy as np, pandas as pd",
        "from scipy.stats import rankdata, mannwhitneyu",
        "nodes = pd.read_parquet(os.path.join(REPO, 'parascape', 'data', 'nodes.parquet'))",
        "GENES = nodes.gene_id.astype(str).tolist()",
        "BASE = Path('/nas_mnt/carruthers/Einar/parascape_omics')",
        "snap = sorted(p for p in (BASE / 'app_snapshots').iterdir()",
        "              if re.fullmatch(r'\\d{8}T\\d{6}Z-\\d+', p.name) and (p / 'rna/snapshot.json').exists())[-1]",
        "design = pd.read_csv(snap / 'rna/sample_design.tsv', sep='\\t', dtype=str)",
        "design = design[design.processing_status == 'verified']",
        "SINGLE_CELL = {'GSE229505', 'GSE346252', 'PRJNA1047766', 'PRJDB12314', 'PRJNA1053644', 'GSE145800', 'GSE145080'}",
        "design = design[~design.study.isin(SINGLE_CELL)]",
        "print(snap.name, len(design), 'verified bulk runs in', design.study.nunique(), 'studies')")

    nb.md("## 1 · Every run, as within-sample ranks",
          "Salmon TPM on the ME49 + host reference, parasite genes only. A gene not detected in a run "
          "(TPM 0) ties at the bottom. The rank is a percentile in [0, 1].")
    nb.code(
        "def percentile_ranks(X):",
        "    \"\"\"Columns of X -> within-column percentile ranks (ties averaged), NaN kept as NaN.\"\"\"",
        "    out = np.full(X.shape, np.nan)",
        "    for j in range(X.shape[1]):",
        "        col = X[:, j]; ok = ~np.isnan(col)",
        "        if ok.sum() > 1:",
        "            out[ok, j] = (rankdata(col[ok]) - 1) / (ok.sum() - 1)",
        "    return out",
        "frames = []",
        "for study, sd in design.groupby('study'):",
        "    path = snap / 'rna/matrices' / study / 'salmon_ME49_host.parquet'",
        "    if not path.exists():",
        "        continue",
        "    tpm = pd.read_parquet(path)['tpm'].reindex(GENES)",
        "    runs = [r for r in sd.run_accession if r in tpm.columns]",
        "    if runs:",
        "        frames.append(pd.DataFrame(percentile_ranks(tpm[runs].to_numpy(dtype=float)), index=GENES, columns=runs))",
        "R = pd.concat(frames, axis=1)",
        "# A run deposited under two studies is ranked once, under the first.",
        "R = R.loc[:, ~R.columns.duplicated()]",
        "# A run with no value for any gene (GSE245775: verified, but zero counts in the snapshot) cannot be scored.",
        "empty_runs = R.columns[R.isna().all()]",
        "R = R.drop(columns=empty_runs)",
        "print(len(empty_runs), 'runs with no values dropped:', sorted(design.set_index('run_accession').loc[empty_runs].study.unique()))",
        "meta = design.drop_duplicates('run_accession').set_index('run_accession').loc[R.columns, ['study', 'condition_group']]",
        "print(R.shape[1], 'runs ranked over', R.shape[0], 'genes')")

    nb.md("## 2 · Cell cycle: built on single cells, tested on a synchronised bulk series",
          "**Built on** Xue et al. 2020 (PMID 32065584): single parasites placed along the cycle, each of "
          "873 genes assigned the phase it peaks in (G1a, G1b, S, M, C). **Tested on** the synchronised "
          "tachyzoite series in this project (`cellcycle19092_*`: hourly from 1 to 12 h, in duplicate, "
          "plus asynchronous and blocked cultures) -- a different study, a different platform, and a "
          "condition it states: the hour after release.",
          "Each sample gets a score per phase (mean rank of that phase genes, minus the mean over all "
          "phase genes, so the phases compete), and a **phase angle**: the five phases placed around a "
          "circle in cycle order, weighted by their positive scores. A culture moving through the cycle "
          "should have an angle that advances hour by hour.")
    nb.code(
        "PHASES = ['G1a', 'G1b', 'S', 'M', 'C']",
        "phase_of = nodes.set_index('gene_id').cellcycle_phase.dropna()",
        "phase_of = phase_of[phase_of.isin(PHASES)]",
        "print(phase_of.value_counts().reindex(PHASES).to_dict())",
        "ANGLE = {p: 2 * np.pi * i / len(PHASES) for i, p in enumerate(PHASES)}",
        "",
        "def phase_scores(ranks):",
        "    \"\"\"ranks: genes x samples percentile ranks -> samples x phases, centred across phases.\"\"\"",
        "    s = pd.DataFrame({p: ranks.reindex(phase_of.index[phase_of == p]).mean() for p in PHASES})",
        "    return s.sub(s.mean(axis=1), axis=0)",
        "",
        "def phase_mixture(scores):",
        "    \"\"\"The positive part of each phase score, closed to one: a relative mixture, not a cell count.\"\"\"",
        "    pos = scores.clip(lower=0)",
        "    return pos.div(pos.sum(axis=1).replace(0, np.nan), axis=0)",
        "",
        "def phase_angle(scores):",
        "    mix = phase_mixture(scores).fillna(0)",
        "    z = sum(mix[p] * np.exp(1j * ANGLE[p]) for p in PHASES)",
        "    return pd.Series(np.angle(z) % (2 * np.pi), index=scores.index), pd.Series(np.abs(z), index=scores.index)")
    nb.code(
        "cc_cols = [c for c in nodes.columns if c.startswith('cellcycle19092_')]",
        "series = nodes.set_index('gene_id')[cc_cols].astype(float)",
        "series_ranks = pd.DataFrame(percentile_ranks(series.to_numpy()), index=series.index, columns=cc_cols)",
        "cc = phase_scores(series_ranks)",
        "angle, strength = phase_angle(cc)",
        "cc['angle_deg'] = np.degrees(angle).round(1); cc['strength'] = strength.round(3)",
        "cc['condition'] = [c.split('_')[1] for c in cc.index]",
        "cc.round(3)")

    nb.md("### Does the inferred phase advance with time?",
          "Hours 1 to 12, the two replicates averaged. For each step from one hour to the next, the "
          "angle change is taken the short way round the circle; a cycle that advances gives mostly "
          "positive steps. The null shuffles the hour labels 10,000 times.")
    nb.code(
        "hours = [f'{h}h' for h in range(1, 13)]",
        "hourly = cc[cc.condition.isin(hours)].groupby('condition')[PHASES].mean().reindex(hours)",
        "a, _ = phase_angle(hourly)",
        "def forward_share(angles):",
        "    d = np.angle(np.exp(1j * np.diff(angles)))",
        "    return float((d > 0).mean()), float(d.sum() / (2 * np.pi))",
        "obs_forward, obs_turns = forward_share(a.to_numpy())",
        "rng = np.random.default_rng(0)",
        "null = np.array([forward_share(rng.permutation(a.to_numpy()))[0] for _ in range(10000)])",
        "p_cc = float((null >= obs_forward).mean())",
        "dominant = hourly.idxmax(axis=1)",
        "cc_validation = pd.DataFrame([dict(axis='cell cycle', built_on='Xue 2020 single-cell phases',",
        "    tested_on='synchronised series cellcycle19092 (hours 1-12)', measure='share of hourly steps that advance the phase',",
        "    value=round(obs_forward, 3), chance=round(float(null.mean()), 3), p=p_cc,",
        "    detail=f'{obs_turns:.2f} turns of the cycle over 12 h; dominant phase by hour: ' + ' '.join(dominant))])",
        "print(dict(zip(hours, dominant))); print('turns over 12 h', round(obs_turns, 2))",
        "cc_validation")

    nb.md("## 3 · Bradyzoite conversion: built on one study, tested on every other that states it",
          "**Built on** GSE223869, which has wild type (BFD1-Ty) and a line without the second "
          "converter, BFD2, each in standard and in alkaline stress medium. Two fingerprints come from it:",
          "* **Alkaline response**: the 150 genes that rise most and the 150 that fall most in wild type "
          "under stress. It reads everything stress does, including what does not need conversion.",
          "* **Conversion**: genes that move in wild type under stress AND fail to move as far without "
          "BFD2 (the smaller of the two changes, per gene). It should read the state of the parasite, "
          "not the medium.",
          "A score is the mean rank of the risers minus the mean rank of the fallers, **converted to a "
          "fraction** between the build study own wild-type standard (0 %) and stressed (100 %) means, "
          "clipped to [0, 1]: \"68 % bradyzoite\" means 68 % of the way from that tachyzoite culture "
          "to that converted one.",
          "**Tested on** every other study whose condition text states bradyzoite-inducing or "
          "tachyzoite conditions, and on two held-out studies that **remove the converter under stress** "
          "(BFD1 knocked out in GSE132238, BFD1 depleted by auxin in GSE223819): there the alkaline "
          "response should stay high and the conversion score should fall.")
    nb.code(
        "BUILD = 'GSE223869'",
        "bm = meta[meta.study == BUILD]",
        "bfd2_gone = bm.condition_group.str.contains('bfd2', case=False)",
        "stress = bm.condition_group.str.contains('alkaline', case=False)",
        "tpm_b = pd.read_parquet(snap / 'rna/matrices' / BUILD / 'salmon_ME49_host.parquet')['tpm'].reindex(GENES)",
        "L = np.log2(tpm_b + 1)",
        "wt_up, wt_down = bm.index[~bfd2_gone & stress].tolist(), bm.index[~bfd2_gone & ~stress].tolist()",
        "ko_up = bm.index[bfd2_gone & stress].tolist()",
        "expressed = L[wt_up + wt_down].mean(axis=1) > 1",
        "lfc_wt = (L[wt_up].mean(axis=1) - L[wt_down].mean(axis=1))[expressed]",
        "lfc_dep = (L[wt_up].mean(axis=1) - L[ko_up].mean(axis=1))[expressed]",
        "conv_up = np.minimum(lfc_wt, lfc_dep); conv_down = np.maximum(lfc_wt, lfc_dep)",
        "SETS = {'response': (lfc_wt.nlargest(150).index, lfc_wt.nsmallest(150).index),",
        "        'conversion': (conv_up.nlargest(150).index, conv_down.nsmallest(150).index)}",
        "MARKERS = {'TGME49_259020': 'BAG1', 'TGME49_268860': 'ENO1', 'TGME49_200385': 'BFD1', 'TGME49_233460': 'SAG1'}",
        "score, fraction = {}, {}",
        "for name, (up, down) in SETS.items():",
        "    sc = R.reindex(up).mean() - R.reindex(down).mean()",
        "    lo, hi = sc[wt_down].mean(), sc[wt_up].mean()",
        "    score[name], fraction[name] = sc, ((sc - lo) / (hi - lo)).clip(0, 1)",
        "    print(name, 'markers among risers', [MARKERS[g] for g in MARKERS if g in up],",
        "          'fallers', [MARKERS[g] for g in MARKERS if g in down],",
        "          '; overlap with response risers', len(set(up) & set(SETS['response'][0])))",
        "print(len(wt_up), 'stressed and', len(wt_down), 'standard wild-type runs;', len(ko_up), 'stressed without BFD2')")
    nb.code(
        "BRADY = r'bradyzoite|alkaline|high ph|ph ?8|cyst|switch|stress media|\\bbz\\b|(?<!un)stress_\\d'",
        "TACHY = r'tachyzoite|ph ?7|standard culture|control media|unstress|\\btz\\b'",
        "text = meta.condition_group.fillna('').str.lower()",
        "is_b = text.str.contains(BRADY, regex=True); is_t = text.str.contains(TACHY, regex=True)",
        "stated = pd.Series(np.where(is_b & ~is_t, 'bradyzoite', np.where(is_t & ~is_b, 'tachyzoite', '')), index=meta.index)",
        "# The converter-removal runs, named one by one: a keyword rule caught the wrong lines here once.",
        "CONVERTER = [('GSE132238', 'BFD1 knocked out', r'\\| WT_stress', r'\\| KO_stress'),",
        "             ('GSE223819', 'BFD1 depleted by auxin', r'\\| 200385-mNG-mAID_Veh', r'\\| 200385-mNG-mAID_\\+IAA')]",
        "for study, _, _, removed in CONVERTER:",
        "    stated[(meta.study == study) & meta.condition_group.str.contains(removed, regex=True)] = 'converter removed'",
        "held = meta.study != BUILD",
        "def auc(a, b):",
        "    a, b = a.dropna(), b.dropna()",
        "    return float(mannwhitneyu(a, b).statistic / (len(a) * len(b))) if len(a) and len(b) else np.nan",
        "rows = []",
        "for study in sorted(meta[held & stated.isin(['bradyzoite', 'tachyzoite'])].study.unique()):",
        "    s = stated[meta.study == study]",
        "    b, t = s.index[s == 'bradyzoite'], s.index[s == 'tachyzoite']",
        "    if len(b) and len(t):",
        "        row = dict(study=study, bradyzoite=len(b), tachyzoite=len(t))",
        "        for name in SETS:",
        "            row[f'auc_{name}'] = round(auc(score[name][b], score[name][t]), 3)",
        "            row[f'fraction_b_{name}'] = round(float(fraction[name][b].median()), 2)",
        "            row[f'fraction_t_{name}'] = round(float(fraction[name][t].median()), 2)",
        "        rows.append(row)",
        "per_study = pd.DataFrame(rows)",
        "hb, ht = stated.index[held & (stated == 'bradyzoite')], stated.index[held & (stated == 'tachyzoite')]",
        "pooled = {name: auc(score[name][hb], score[name][ht]) for name in SETS}",
        "print({k: round(v, 3) for k, v in pooled.items()}, 'pooled AUC over', len(hb), 'bradyzoite and', len(ht), 'tachyzoite held-out runs')",
        "per_study")

    nb.md("### Removing the converter under stress",
          "Same medium, same study, the converter present or removed. The alkaline response should "
          "barely move; the conversion score should fall. `auc` is how often a converter-present run "
          "outscores a converter-removed one (1 = always).")
    nb.code(
        "conv_rows = []",
        "for study, what, present, removed in CONVERTER:",
        "    m = meta[meta.study == study]",
        "    pr = m.index[m.condition_group.str.contains(present, regex=True)]",
        "    rm = m.index[m.condition_group.str.contains(removed, regex=True)]",
        "    for name in SETS:",
        "        conv_rows.append(dict(study=study, removed=what, fingerprint=name, present_runs=len(pr), removed_runs=len(rm),",
        "                              fraction_present=round(float(fraction[name][pr].median()), 2),",
        "                              fraction_removed=round(float(fraction[name][rm].median()), 2),",
        "                              auc=round(auc(score[name][pr], score[name][rm]), 3)))",
        "converter = pd.DataFrame(conv_rows)",
        "converter.to_csv(os.path.join(OUT, 'converter_removal.tsv'), sep='\\t', index=False)",
        "converter")

    nb.md("### A metadata swap the fingerprint exposes",
          "In GSE223819 two tagged lines carry a genotype field and a sample title that name different "
          "genes (208020 against 311100, which is BFD2). The line titled 311100 converts poorly even on "
          "vehicle -- what a tag that weakens BFD2 would do -- and the line titled 208020 converts "
          "fully. The fingerprint therefore sides with the **titles**. Shown as evidence, not applied: "
          "the registry keeps the deposit text.")
    nb.code(
        "m = meta[meta.study == 'GSE223819'].copy()",
        "m['genotype_field'] = m.condition_group.str.extract(r'genotype: ME49/TIR1/?([0-9]*)', expand=False).replace('', 'parental')",
        "m['sample_title'] = m.condition_group.str.extract(r'\\| ([0-9a-z]+)', expand=False)",
        "m['treatment'] = np.where(m.condition_group.str.contains(r'IAA_|\\+IAA', regex=True), 'auxin', 'vehicle')",
        "m['conversion'] = fraction['conversion'].reindex(m.index).values",
        "swap = m.groupby(['genotype_field', 'sample_title', 'treatment']).conversion.median().round(2).unstack()",
        "swap.to_csv(os.path.join(OUT, 'gse223819_label_evidence.tsv'), sep='\\t')",
        "swap")

    nb.md("## 4 · The annotations, and what ships",
          "An axis ships when its held-out test clears chance clearly: the cell-cycle phase advances "
          "more often than shuffled hours do (p < 0.01); a bradyzoite fingerprint separates stated "
          "conditions in held-out studies (pooled AUC at least 0.8 and most testable studies above "
          "0.7). The conversion fingerprint must also fall when the converter is removed (AUC at least "
          "0.8 in both converter studies), or it is no more than the alkaline response under another "
          "name. Anything that fails is written `exploratory:` and should not be used.")
    nb.code(
        "cc_ok = p_cc < 0.01",
        "ok, rows = {}, [cc_validation]",
        "for name in SETS:",
        "    col = per_study[f'auc_{name}']",
        "    conv_auc = converter[converter.fingerprint == name].auc",
        "    ok[name] = pooled[name] >= 0.8 and (col > 0.7).mean() > 0.5",
        "    if name == 'conversion':",
        "        ok[name] = ok[name] and bool((conv_auc >= 0.8).all())",
        "    rows.append(pd.DataFrame([dict(axis=f'bradyzoite {name}', built_on=f'{BUILD} alkaline vs standard' + (' and without BFD2' if name == 'conversion' else ''),",
        "        tested_on=f'{len(per_study)} other studies that state the condition', measure='pooled AUC, stated bradyzoite vs tachyzoite runs',",
        "        value=round(pooled[name], 3), chance=0.5, p=np.nan,",
        "        detail=f'{int((col > 0.7).sum())} of {len(col)} studies above 0.7; converter removed vs present AUC ' + ', '.join(f'{v:.2f}' for v in conv_auc))]))",
        "validation = pd.concat(rows, ignore_index=True)",
        "validation['ships'] = [bool(cc_ok)] + [bool(ok[n]) for n in SETS]",
        "validation.to_csv(os.path.join(OUT, 'validation.tsv'), sep='\\t', index=False)",
        "per_study.to_csv(os.path.join(OUT, 'bradyzoite_per_study.tsv'), sep='\\t', index=False)",
        "tier = lambda good, name: ('filled' if good else 'exploratory') + f':{name}@{VERSION}'",
        "ccs = phase_scores(R); ang, stren = phase_angle(ccs); mix = phase_mixture(ccs)",
        "ann = pd.DataFrame({'run': R.columns, 'study': meta.study.values, 'condition_group': meta.condition_group.values})",
        "for p in PHASES:",
        "    ann[f'cellcycle_{p}'] = mix[p].round(3).values",
        "ann['cellcycle_angle_deg'] = np.degrees(ang).round(1).values",
        "ann['cellcycle_source'] = tier(cc_ok, 'cellcycle_xue2020')",
        "ann['cellcycle_depends_on'] = 'PMID 32065584'",
        "for name in SETS:",
        "    ann[f'bradyzoite_{name}_fraction'] = fraction[name].reindex(R.columns).round(3).values",
        "    ann[f'bradyzoite_{name}_source'] = tier(ok[name], f'bradyzoite_{name}_gse223869')",
        "ann['bradyzoite_stated'] = stated.reindex(R.columns).values",
        "ann['bradyzoite_depends_on'] = BUILD",
        "ann.to_csv(os.path.join(OUT, 'sample_annotations.tsv'), sep='\\t', index=False)",
        "fp = [pd.DataFrame({'fingerprint': 'cellcycle_xue2020', 'gene': phase_of.index, 'set': phase_of.values})]",
        "for name, (up, down) in SETS.items():",
        "    fp += [pd.DataFrame({'fingerprint': f'bradyzoite_{name}_gse223869', 'gene': up, 'set': 'rises'}),",
        "           pd.DataFrame({'fingerprint': f'bradyzoite_{name}_gse223869', 'gene': down, 'set': 'falls'})]",
        "pd.concat(fp).to_csv(os.path.join(OUT, 'fingerprints.tsv'), sep='\\t', index=False)",
        "validation")

    nb.md("### What the annotations say about datasets nobody labelled",
          "The studies that state neither condition, ranked by their median inferred conversion -- the "
          "new metadata this item exists to produce. Read with the validation above.")
    nb.code(
        "unl = ann[(ann.bradyzoite_stated == '') & (ann.study != BUILD)]",
        "unlabelled = unl.groupby('study').agg(runs=('run', 'size'),",
        "    conversion=('bradyzoite_conversion_fraction', 'median'), response=('bradyzoite_response_fraction', 'median'),",
        "    cellcycle_angle=('cellcycle_angle_deg', 'median')).sort_values('conversion', ascending=False)",
        "unlabelled.to_csv(os.path.join(OUT, 'unlabelled_studies.tsv'), sep='\\t')",
        "unlabelled.head(15)")
    nb.code(
        "def md(frame):",
        "    cols = [str(c) for c in frame.columns]",
        "    lines = ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join('---' for _ in cols) + '|']",
        "    return '\\n'.join(lines + ['| ' + ' | '.join(str(v) for v in r) + ' |' for r in frame.itertuples(index=False)])",
        "readme = [f'# Condition fingerprints {VERSION} (instruction 75)', '',",
        "    'Notebook: `notebooks/condition_fingerprints_" + DATE + ".ipynb` (executed). Builder: '",
        "    '`scripts/build_condition_fingerprints_notebook.py`. Snapshot ' + snap.name + f', {R.shape[1]} bulk runs scored.', '',",
        "    '## Validation on held-out studies', '', md(validation[['axis', 'built_on', 'tested_on', 'value', 'chance', 'p', 'detail', 'ships']]), '',",
        "    '## Removing the converter under stress', '', md(converter), '',",
        "    'Both bradyzoite fingerprints fall when BFD1 is removed; the conversion-specific one falls only slightly '",
        "    'further, so most of the alkaline response already depends on BFD1.', '',",
        "    '## GSE223819: genotype field and sample title disagree', '', md(swap.reset_index()), '',",
        "    'The line titled 311100 (BFD2) converts poorly even on vehicle; the fingerprint sides with the titles. '",
        "    'Evidence only; the registry keeps the deposit text.', '',",
        "    'Files: `sample_annotations.tsv` (every run: phase mixture and angle, both bradyzoite fractions, the stated '",
        "    'condition, `filled:`/`exploratory:` source and `depends_on`), `fingerprints.tsv`, `validation.tsv`, '",
        "    '`bradyzoite_per_study.tsv`, `converter_removal.tsv`, `unlabelled_studies.tsv`. Annotations are inferred '",
        "    'and are never evidence for a gene trait.']",
        "open(os.path.join(OUT, 'README.md'), 'w').write('\\n'.join(readme) + '\\n')",
        "print('README written')")
    nb.write(NB)


if __name__ == "__main__":
    main()
