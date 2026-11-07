"""Co-expression -> localization through the harness (instruction 71 product 1; 70.08), as an executed notebook.

The first product 1 build (`build_coexpression_localization_notebook.py`, 2026-11-04, 441 runs) scored
co-expression nearest-neighbour votes against the 26-niche markers (coarse accuracy 0.466 against 0.132 for the
best baseline). This one is its complement on the 735-run snapshot: label-free co-expression components, scored
per coarse niche through the harness with a frozen partition each, side by side with all allowed evidence.

Genes of one compartment are made together (rhoptry, microneme and IMC waves), so co-expression across
the uniformly reprocessed RNA-seq should place proteins without hyperLOPIT. Built per study (within-study
gene-gene Spearman on log TPM, studies with at least 6 bulk runs), combined across studies by mean rank so
no single large study dominates, then reduced to its leading components -- **label-free features**: a
feature built from neighbours' labels would carry frozen-test labels into training. Scored on the
localization niche questions of `build_localization_questions_notebook.py` with the co-expression
components as the only evidence, a fresh frozen partition per question (seed 2).

    python scripts/build_coexpression_components_localization_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_06"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"coexpression_components_localization_{DATE}")
NB = os.path.join(REPO, "notebooks", f"coexpression_components_localization_{DATE}.ipynb")
MARKERS = "/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/sources/barylyuk2020_mmc5/mmc5.xls"


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Co-expression places proteins: localization from RNA alone",
                          namespace={"REPO": REPO, "OUT": OUT, "MARKERS": MARKERS})
    nb.md(
        "**The question.** Can co-expression across the uniformly reprocessed RNA-seq, and nothing "
        "else, tell which compartment a protein is in? (instruction 71, product 1; instruction 70.08).",
        "**How.** Within each study with at least 6 bulk runs, the Spearman correlation of every pair "
        "of parasite genes on log2(TPM + 1); each study matrix turned into ranks and averaged, so a "
        "pair is close when it is close in many studies, not in one big one; the consensus reduced to "
        "its 20 leading components. Those 20 numbers per gene are the only evidence the model sees.",
        "**Why components and not neighbours' labels.** A feature such as \"the commonest niche among "
        "my co-expression neighbours\" uses other genes' labels, including the frozen-test genes', and "
        "would leak them into training. The components use no label at all.",
        f"Built {_dt.date.today().isoformat()}.")

    nb.code(
        "import os, re",
        "from pathlib import Path",
        "import numpy as np, pandas as pd",
        "from scipy.stats import rankdata",
        "from parascape import harness as H",
        "from parascape.question_registry import Question",
        "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet')).set_index('gene_id')",
        "GENES = nodes.index.astype(str).tolist()",
        "BASE = Path('/nas_mnt/carruthers/Einar/parascape_omics')",
        "snap = sorted(p for p in (BASE / 'app_snapshots').iterdir()",
        "              if re.fullmatch(r'\\d{8}T\\d{6}Z-\\d+', p.name) and (p / 'rna/snapshot.json').exists()",
        "              and (p / 'rna/sample_design.tsv').exists())[-1]   # a snapshot still being written is skipped",
        "design = pd.read_csv(snap / 'rna/sample_design.tsv', sep='\\t', dtype=str)",
        "design = design[design.processing_status == 'verified'].drop_duplicates('run_accession')",
        "SINGLE_CELL = {'GSE229505', 'GSE346252', 'PRJNA1047766', 'PRJDB12314', 'PRJNA1053644', 'GSE145800', 'GSE145080'}",
        "design = design[~design.study.isin(SINGLE_CELL)]",
        "print(snap.name, len(design), 'verified bulk runs in', design.study.nunique(), 'studies')")

    nb.md("## 1 · Co-expression per study, combined by mean rank",
          "Genes expressed (TPM >= 1) in at least half the runs of a study enter that study; a pair "
          "missing from a study is left out of its average rather than counted as uncorrelated. The "
          "per-study matrix is ranked within the study, scaled to [0, 1], and the ranks are averaged.")
    nb.code(
        "n = len(GENES)",
        "total = np.zeros((n, n), dtype=np.float32)",
        "count = np.zeros((n, n), dtype=np.uint8)",
        "used = []",
        "for study, sd in design.groupby('study'):",
        "    path = snap / 'rna/matrices' / study / 'salmon_ME49_host.parquet'",
        "    if not path.exists():",
        "        continue",
        "    tpm = pd.read_parquet(path)['tpm'].reindex(GENES)",
        "    runs = [r for r in sd.run_accession if r in tpm.columns and tpm[r].notna().any() and tpm[r].sum() > 0]",
        "    if len(runs) < 6:",
        "        continue",
        "    X = np.log2(tpm[runs].to_numpy(dtype=float) + 1)",
        "    keep = (np.nan_to_num(X) >= 1).mean(axis=1) >= 0.5",
        "    if keep.sum() < 1000:",
        "        continue",
        "    R = np.apply_along_axis(rankdata, 1, X[keep])",
        "    R = (R - R.mean(axis=1, keepdims=True))",
        "    R /= np.clip(np.linalg.norm(R, axis=1, keepdims=True), 1e-9, None)",
        "    C = (R @ R.T).astype(np.float32)",
        "    ranks = (rankdata(C.ravel()).reshape(C.shape) / C.size).astype(np.float32)",
        "    idx = np.flatnonzero(keep)",
        "    total[np.ix_(idx, idx)] += ranks",
        "    count[np.ix_(idx, idx)] += 1",
        "    used.append(dict(study=study, runs=len(runs), genes=int(keep.sum())))",
        "    del C, ranks, R",
        "used = pd.DataFrame(used)",
        "consensus = np.where(count > 0, total / np.maximum(count, 1), np.nan).astype(np.float32)",
        "del total",
        "print(len(used), 'studies combined;', int((count.diagonal() > 0).sum()), 'genes in at least one;',",
        "      'median studies per gene', int(np.median(count.diagonal()[count.diagonal() > 0])))",
        "used.to_csv(os.path.join(OUT, 'studies_used.tsv'), sep='\\t', index=False)",
        "used.sort_values('runs', ascending=False).head(10)")

    nb.md("## 2 · The leading components: 20 label-free numbers per gene",
          "Genes seen in at least 3 studies; missing pairs filled with the neutral rank 0.5; the "
          "consensus centred and its 20 leading eigenvectors taken. A gene seen in fewer studies has "
          "no components (missing, never zero).")
    nb.code(
        "from scipy.sparse.linalg import eigsh",
        "seen = count.diagonal() >= 3",
        "S = np.nan_to_num(consensus[np.ix_(seen, seen)], nan=0.5) - 0.5",
        "vals, vecs = eigsh(S.astype(np.float64), k=20, which='LA')",
        "order = np.argsort(vals)[::-1]",
        "comps = pd.DataFrame(np.nan, index=GENES, columns=[f'coexp_pc{i + 1}' for i in range(20)])",
        "comps.loc[np.array(GENES)[seen], :] = vecs[:, order] * np.sqrt(np.clip(vals[order], 0, None))",
        "comps.to_parquet(os.path.join(OUT, 'coexpression_components.parquet'))",
        "print(int(seen.sum()), 'genes with components; leading eigenvalues', np.round(vals[order][:5], 1))")

    nb.md("## 3 · Localization from co-expression alone",
          "The 13 coarse niches with at least 25 hyperLOPIT markers, every other marker a curated "
          "\"no\"; evidence = the 20 components only; boosted trees; nested cross-validation; a fresh "
          "frozen partition per question (seed 2; seeds 0 and 1 hold the all-evidence and genome-wide "
          "runs of the same questions).")
    nb.code(
        "m = pd.read_excel(MARKERS, 'S4a - TAGM-MAP classification')",
        "m = m[m.markers.astype(str) != 'unknown'].set_index('Accession')",
        "COARSE = {'rhoptries 1': 'rhoptries', 'rhoptries 2': 'rhoptries', 'apical 1': 'apical', 'apical 2': 'apical',",
        "          'ER 1': 'ER', 'ER 2': 'ER', 'PM - peripheral 1': 'plasma membrane', 'PM - peripheral 2': 'plasma membrane',",
        "          'PM - integral': 'plasma membrane', 'nucleus - chromatin': 'nucleus', 'nucleus - non-chromatin': 'nucleus',",
        "          'nucleolus': 'nucleus', 'mitochondrion - membranes': 'mitochondrion', 'mitochondrion - soluble': 'mitochondrion',",
        "          '40S ribosome': 'ribosome', '60S ribosome': 'ribosome', '19S proteasome': 'proteasome', '20S proteasome': 'proteasome'}",
        "niche = m.markers.astype(str).map(lambda v: COARSE.get(v, v))",
        "niche = niche[niche.index.isin(nodes.index)]",
        "NICHES = [k for k, v in niche.value_counts().items() if v >= 25]",
        "features = comps.loc[niche.index]",
        "families = pd.Series('coexpression_uniform_rnaseq', index=features.columns)",
        "groups = nodes['orthogroup'].reindex(niche.index)",
        "groups = groups.where(groups.notna(), pd.Series({g: 'self:' + g for g in niche.index}))",
        "rows = []",
        "for n_ in NICHES:",
        "    q = Question(id='localization_' + n_.replace(' ', '_'), text=f'Is this protein in the {n_}?',",
        "                 answer_type='boolean', question_class='B', label_source='barylyuk2020 mmc5 markers',",
        "                 label_family='lopit_tgon', label_tier='observed-curated', label_category='place',",
        "                 evidence_forbidden=('lopit_tgon',))",
        "    y = (niche == n_)",
        "    r = H.run_question(q, labels=y, features=features, families=families, groups=groups,",
        "                       strategies=[H.Strategy.of('boosted trees, uniform co-expression components', 'hist_gbt', max_iter=150)],",
        "                       seed=2, importance_repeats=1, log_path=os.path.join(REPO, H.HOLDOUT_LOG))",
        "    cv, ft = r.cross_validated, r.frozen_test or {}",
        "    rows.append(dict(niche=n_, markers=int(y.sum()), with_components=int(features.loc[y[y].index].notna().all(axis=1).sum()),",
        "                     coexp_cv_auc=round(cv['roc_auc'], 3), coexp_test_auc=round(ft.get('roc_auc', np.nan), 3),",
        "                     coexp_test_ap=round(ft.get('average_precision', np.nan), 3), may_fill=H.may_fill(r)))",
        "scores = pd.DataFrame(rows)",
        "prior = os.path.join(REPO, 'results/localization_questions_2026_11_06/scores.tsv')",
        "if os.path.exists(prior):",
        "    allev = pd.read_csv(prior, sep='\\t')[['niche', 'cv_auc', 'test_auc']].rename(",
        "        columns={'cv_auc': 'all_evidence_cv_auc', 'test_auc': 'all_evidence_test_auc'})",
        "    scores = scores.merge(allev, on='niche', how='left')",
        "scores = scores.sort_values('coexp_cv_auc', ascending=False)",
        "scores.to_csv(os.path.join(OUT, 'scores.tsv'), sep='\\t', index=False)",
        "scores")

    nb.md("## Reading it",
          "Co-expression alone against all allowed evidence, niche by niche, on the cross-validated "
          "estimate (the steadier one; frozen partitions hold 4-13 positives). A niche co-expression "
          "places well is one whose proteins are made together; one it cannot place is made on no "
          "shared schedule across these studies.")
    nb.code(
        "summary = pd.DataFrame([dict(niches=len(scores), studies=len(used),",
        "    median_coexp_cv_auc=float(scores.coexp_cv_auc.median()),",
        "    median_all_evidence_cv_auc=float(scores.all_evidence_cv_auc.median()) if 'all_evidence_cv_auc' in scores else np.nan,",
        "    coexp_above_0_8=int((scores.coexp_cv_auc >= 0.8).sum()), may_fill=int(scores.may_fill.sum()))])",
        "summary.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\\t', index=False)",
        "summary")
    nb.write(NB)


if __name__ == "__main__":
    main()
