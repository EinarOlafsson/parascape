"""Condition fingerprints (instruction 75), the strain-type axis: a positive control, as an executed notebook.

Every bulk RNA-seq run whose metadata states its parasite strain is labelled with the strain type (I: RH, GT1;
II: ME49, Prugniaud/Pru; III: VEG, CTG). A fingerprint learned from the other studies -- the genes whose
within-sample ranks separate the types in training, and each type centroid -- must recover the stated type of a
held-out study. Leave-one-study-out: no run of a held-out study is used for gene selection or centroids.
Strain is already known per run, so this axis is a positive control for the fingerprint method rather than new
metadata; annotations for runs that do not state a strain are written as INFERRED (`filled:strain_type@v1`).

Stated limit: reads are mapped to the ME49 (type II) reference, so part of the signal is mapping of polymorphic
genes, which is a real consequence of strain origin but not purely expression.

    python scripts/build_strain_fingerprint_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
VERSION = "v1"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"strain_fingerprint_{DATE}")
NB = os.path.join(REPO, "notebooks", f"strain_fingerprint_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    nb = ExecutedNotebook("Condition fingerprints: strain type (positive control)",
                          namespace={"REPO": REPO, "OUT": OUT, "VERSION": VERSION})
    nb.md("**The axis (instruction 75).** Parasite strain type, read from expression. Strain is stated in the "
          "metadata of most runs, so this is a positive control: if a fingerprint cannot recover a stated, strong "
          "property of a held-out study, the method is not ready for the axes nobody states.",
          "**Validation.** Leave one study out: gene selection and type centroids use every other study; the "
          "held-out study runs are then assigned the nearest centroid. Compared with always guessing the "
          "commonest training type.",
          "**Stated limit.** Reads are mapped to the ME49 (type II) reference; polymorphic genes map less well "
          "from type I and III parasites, so part of the signal is mapping, which follows strain origin.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, re",
            "from pathlib import Path",
            "import numpy as np, pandas as pd",
            "from scipy.stats import rankdata",
            "nodes = pd.read_parquet(os.path.join(REPO, 'parascape', 'data', 'nodes.parquet'))",
            "GENES = nodes.gene_id.astype(str).tolist()",
            "BASE = Path('/nas_mnt/carruthers/Einar/parascape_omics')",
            "snap = sorted(p for p in (BASE / 'app_snapshots').iterdir()",
            "              if re.fullmatch(r'\\d{8}T\\d{6}Z-\\d+', p.name) and (p / 'rna/snapshot.json').exists())[-1]",
            "design = pd.read_csv(snap / 'rna/sample_design.tsv', sep='\\t', dtype=str)",
            "design = design[design.processing_status == 'verified']",
            "SINGLE_CELL = {'GSE229505', 'GSE346252', 'PRJNA1047766', 'PRJDB12314', 'PRJNA1053644', 'GSE145800', 'GSE145080'}",
            "design = design[~design.study.isin(SINGLE_CELL)].drop_duplicates('run_accession')",
            "print(snap.name, len(design), 'verified bulk runs in', design.study.nunique(), 'studies')")
    nb.md("## 1 · The stated strain type of each run",
          "From the run attributes whose key names the parasite strain (keys naming a mouse or host strain are "
          "skipped). Type I: RH, GT1, 'type 1'. Type II: ME49, Prugniaud/Pru/PTG, 'type 2'. Type III: VEG, CTG, "
          "CEP, 'type 3'. Anything else (field isolates, unnamed clones) stays unstated.")
    nb.code("TYPES = [('I', r'\\b(RH|GT-?1|type ?(1|I))\\b'), ('III', r'\\b(VEG|CTG|CEP|type ?(3|III))\\b'),",
            "         ('II', r'\\b(ME-?49|PRU|PRUGNIAUD|PTG|type ?(2|II))\\b')]",
            "def stated_type(attributes):",
            "    try:",
            "        j = json.loads(attributes) if isinstance(attributes, str) else {}",
            "    except Exception:",
            "        return None",
            "    if not isinstance(j, dict):",
            "        return None",
            "    for k, v in j.items():",
            "        key = k.lower()",
            "        if 'strain' not in key or any(w in key for w in ('mice', 'mouse', 'host', 'cell')):",
            "            continue",
            "        text = str(v)",
            "        hits = {t for t, rx in TYPES if re.search(rx, text, flags=re.IGNORECASE)}",
            "        if len(hits) == 1:",
            "            return hits.pop()",
            "    return None",
            "design['strain_type'] = design.attributes_json.map(stated_type)",
            "print(design.strain_type.value_counts(dropna=False).to_dict())",
            "per_study = design.dropna(subset=['strain_type']).groupby('study').strain_type.agg(lambda s: ','.join(sorted(set(s))))",
            "print(len(per_study), 'studies state a type;', int(per_study.str.contains(',').sum()), 'state more than one')")
    nb.md("## 2 · Every run as within-sample ranks", "Salmon TPM on ME49 + host, parasite genes, percentile ranks.")
    nb.code("def percentile_ranks(X):",
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
            "R = R.loc[:, ~R.columns.duplicated()]",
            "R = R.drop(columns=R.columns[R.isna().all()])",
            "meta = design.set_index('run_accession').loc[R.columns, ['study', 'strain_type']]",
            "print(R.shape[1], 'runs ranked;', int(meta.strain_type.notna().sum()), 'with a stated type in',",
            "      meta.dropna().study.nunique(), 'studies')")
    nb.md("## 3 · Leave one study out",
          "For each study that states a type: in the training studies only, keep the 300 genes whose ranks differ "
          "most between types (F statistic on training runs), form one centroid per type, assign each held-out "
          "run to the nearest centroid (correlation). Chance = the commonest training type.")
    nb.code("labelled = meta.dropna(subset=['strain_type'])",
            "X = R[labelled.index].fillna(0.5)",
            "def select_genes(Xtr, ytr, n=300):",
            "    groups = [Xtr.loc[:, ytr == t] for t in sorted(set(ytr))]",
            "    if len(groups) < 2:",
            "        return Xtr.index[:0]",
            "    grand = Xtr.mean(axis=1)",
            "    between = sum(g.shape[1] * (g.mean(axis=1) - grand) ** 2 for g in groups) / (len(groups) - 1)",
            "    within = sum(((g.sub(g.mean(axis=1), axis=0)) ** 2).sum(axis=1) for g in groups) / max(Xtr.shape[1] - len(groups), 1)",
            "    f = between / within.replace(0, np.nan)",
            "    return f.dropna().nlargest(n).index",
            "def predict(Xtr, ytr, Xte):",
            "    genes = select_genes(Xtr, ytr)",
            "    cents = {t: Xtr.loc[genes, ytr == t].mean(axis=1) for t in sorted(set(ytr))}",
            "    C = pd.DataFrame(cents)",
            "    corr = pd.DataFrame({t: Xte.loc[genes].corrwith(C[t]) for t in C}, index=Xte.columns)",
            "    return corr.idxmax(axis=1), corr",
            "rows = []",
            "for study in sorted(labelled.study.unique()):",
            "    te = labelled.index[labelled.study == study]",
            "    tr = labelled.index[labelled.study != study]",
            "    ytr = labelled.loc[tr, 'strain_type']",
            "    if ytr.nunique() < 2:",
            "        continue",
            "    pred, _ = predict(X[tr], ytr, X[te])",
            "    truth = labelled.loc[te, 'strain_type']",
            "    chance = ytr.value_counts().index[0]",
            "    for r in te:",
            "        rows.append(dict(study=study, run=r, truth=truth[r], predicted=pred[r], chance_guess=chance))",
            "loso = pd.DataFrame(rows)",
            "loso['right'] = loso.truth == loso.predicted",
            "loso['chance_right'] = loso.truth == loso.chance_guess",
            "loso.to_csv(os.path.join(OUT, 'leave_one_study_out.tsv'), sep='\\t', index=False)",
            "by_study = loso.groupby('study').agg(type=('truth', lambda s: ','.join(sorted(set(s)))), runs=('run', 'size'),",
            "                                      right=('right', 'mean'), chance=('chance_right', 'mean')).round(2)",
            "summary = dict(runs=len(loso), studies=loso.study.nunique(), accuracy=round(float(loso.right.mean()), 3),",
            "               chance=round(float(loso.chance_right.mean()), 3),",
            "               study_majority_right=round(float((by_study.right > 0.5).mean()), 3))",
            "print(summary)",
            "pd.crosstab(loso.truth, loso.predicted)")
    nb.code("by_study")
    nb.md("### Per type, and what cannot be validated",
          "Recall per stated type. A type stated by a single study cannot be validated by leaving studies out: "
          "with that study held out, the type is absent from training and cannot be predicted at all. Such a "
          "type is reported as unvalidated, and inferred calls of it are flagged.")
    nb.code("recall = loso.groupby('truth').right.mean().round(3)",
            "studies_per_type = labelled.groupby('strain_type').study.nunique()",
            "validated = [t for t in recall.index if studies_per_type.get(t, 0) >= 2]",
            "summary.update(recall=recall.to_dict(), studies_per_type=studies_per_type.to_dict(), validated_types=validated,",
            "               accuracy_validated_types=round(float(loso[loso.truth.isin(validated)].right.mean()), 3))",
            "summary")
    nb.md("## 4 · Annotations for runs that do not state a strain",
          "Fitted on every labelled run, written as inferred with the correlation margin; never evidence for a "
          "gene trait.")
    nb.code("unl = meta.index[meta.strain_type.isna()]",
            "pred, corr = predict(X, labelled.strain_type, R[unl].fillna(0.5)) if len(unl) else (pd.Series(dtype=str), pd.DataFrame())",
            "ann = pd.DataFrame({'run_accession': unl, 'study': meta.loc[unl, 'study'].values,",
            "                    'strain_type_inferred': pred.reindex(unl).values,",
            "                    'margin': (corr.max(axis=1) - corr.apply(lambda r: r.nlargest(2).iloc[-1], axis=1)).reindex(unl).round(3).values,",
            "                    'provenance': f'filled:strain_type@{VERSION}', 'depends_on': 'all labelled studies'})",
            "ann['type_validated'] = ann.strain_type_inferred.isin(validated)",
            "ann.to_csv(os.path.join(OUT, 'annotations.tsv'), sep='\\t', index=False)",
            "json.dump(dict(summary=summary, snapshot=snap.name, version=VERSION, annotated_runs=len(ann),",
            "               stated=design.strain_type.value_counts().to_dict()), open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "ann.groupby(['study', 'strain_type_inferred']).size().unstack(fill_value=0)")
    nb.write(NB)


if __name__ == "__main__":
    main()
