"""Localization as harness questions (instruction 70, action 70.04), as an executed notebook.

The harness scores yes/no labels, so localization is asked one coarse niche at a time: "is this
protein in the dense granules?", with every other curation-grade hyperLOPIT marker as a "no". Each
question is built here, not registered in the application, so that nothing the application shows
changes before Einar has seen the numbers. Every column of category `place` is held out (the
hyperLOPIT experiment and everything transferred from it), so the question is whether the other
evidence places a protein.

    python scripts/build_localization_questions_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_06"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"localization_questions_{DATE}")
NB = os.path.join(REPO, "notebooks", f"localization_questions_{DATE}.ipynb")
MARKERS = "/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/sources/barylyuk2020_mmc5/mmc5.xls"


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Localization, one niche at a time, through the harness",
                          namespace={"REPO": REPO, "OUT": OUT, "MARKERS": MARKERS})
    nb.md(
        "**What this is.** The first localization scores on a frozen partition (instruction 70). "
        "Localization is the target Einar chose first; the harness scores yes/no questions, so each "
        "coarse niche is one question.",
        "**The labels.** The 718 final hyperLOPIT markers of Barylyuk et al. 2020 (mmc5, sha256 "
        "`a58a8778...`), merged from 26 niches to coarse ones (rhoptries 1 and 2 are rhoptries). For "
        "a niche, its markers are \"yes\" and every other marker is \"no\" -- a curated negative, not "
        "an unlabelled gene. **Caveat, stated first:** the paper chose 656 of the markers from prior "
        "knowledge, inferred function and profile consistency, so marker status is curation, not "
        "independent experiment, and evidence that encodes prior knowledge may recover it.",
        "**What is held out.** Every column of the `place` category, whatever its family -- the "
        "hyperLOPIT experiment, its posteriors, and transfers built from it.",
        f"Built {_dt.date.today().isoformat()}.")

    nb.code(
        "import os",
        "import numpy as np, pandas as pd",
        "from parascape import harness as H, datasets as ds, organisms as O",
        "from parascape.question_registry import Question",
        "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet')).set_index('gene_id')",
        "m = pd.read_excel(MARKERS, 'S4a - TAGM-MAP classification')",
        "m = m[m.markers.astype(str) != 'unknown'].set_index('Accession')",
        "COARSE = {'rhoptries 1': 'rhoptries', 'rhoptries 2': 'rhoptries', 'apical 1': 'apical', 'apical 2': 'apical',",
        "          'ER 1': 'ER', 'ER 2': 'ER', 'PM - peripheral 1': 'plasma membrane', 'PM - peripheral 2': 'plasma membrane',",
        "          'PM - integral': 'plasma membrane', 'nucleus - chromatin': 'nucleus', 'nucleus - non-chromatin': 'nucleus',",
        "          'nucleolus': 'nucleus', 'mitochondrion - membranes': 'mitochondrion', 'mitochondrion - soluble': 'mitochondrion',",
        "          '40S ribosome': 'ribosome', '60S ribosome': 'ribosome', '19S proteasome': 'proteasome', '20S proteasome': 'proteasome'}",
        "niche = m.markers.astype(str).map(lambda v: COARSE.get(v, v))",
        "niche = niche[niche.index.isin(nodes.index)]",
        "counts = niche.value_counts()",
        "print(len(m), 'markers;', len(niche), 'in the node table;', counts.size, 'coarse niches')",
        "counts.to_frame('markers')")

    nb.md("## The questions, and the evidence each may use",
          "Coarse niches with at least 25 markers. Each question holds out the `place` category.")
    nb.code(
        "def family_of(column):",
        "    found = ds.provenance(column, O.TOXOPLASMA)",
        "    return found.key if found is not None else 'unregistered'",
        "numeric = [c for c in nodes.columns if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]",
        "families = pd.Series({c: family_of(c) for c in numeric}, name='family')",
        "full = nodes.loc[niche.index, numeric].astype('float64')",
        "groups = nodes['orthogroup'].reindex(niche.index)",
        "groups = groups.where(groups.notna(), pd.Series({g: 'self:' + g for g in niche.index}))",
        "NICHES = [n for n, k in counts.items() if k >= 25]",
        "def question(n):",
        "    return Question(id='localization_' + n.replace(' ', '_'), text=f'Is this protein in the {n}?',",
        "                    answer_type='boolean', question_class='B', label_source='barylyuk2020 mmc5 markers',",
        "                    label_family='lopit_tgon', label_tier='observed-curated', label_category='place',",
        "                    evidence_forbidden=('lopit_tgon',))",
        "features = H.allowed_features(question(NICHES[0]), full, families)",
        "print(len(NICHES), 'questions;', features.shape[1], 'allowed columns of', full.shape[1],",
        "      '(', full.shape[1] - features.shape[1], 'held out as place )')")

    nb.md("## The runs",
          "Boosted trees over every allowed column, nested 3 x 3 cross-validation inside train + tune, "
          "calibration, then the frozen partition once per question (seed 0, logged).")
    nb.code(
        "rows = []",
        "for n in NICHES:",
        "    q = question(n)",
        "    y = (niche == n)",
        "    r = H.run_question(q, labels=y, features=features, families=families, groups=groups,",
        "                       strategies=[H.Strategy.of('boosted trees', 'hist_gbt', max_iter=150)], seed=0,",
        "                       importance_repeats=2, log_path=os.path.join(REPO, H.HOLDOUT_LOG))",
        "    cv, ft = r.cross_validated, r.frozen_test or {}",
        "    base = r.baselines[r.baselines.partition == 'frozen_test']",
        "    top = r.importance.head(3).family.tolist() if len(r.importance) else []",
        "    rows.append(dict(niche=n, markers=int(y.sum()), cv_auc=round(cv['roc_auc'], 3), cv_ap=round(cv['average_precision'], 3),",
        "                     test_n=ft.get('n'), test_positives=ft.get('positives'), test_auc=round(ft.get('roc_auc', np.nan), 3),",
        "                     test_ap=round(ft.get('average_precision', np.nan), 3), test_calls=ft.get('calls'),",
        "                     test_precision=ft.get('precision'), test_lower90=round(ft.get('precision_lower90', np.nan), 3),",
        "                     best_baseline_lower90=round(float(base.precision_lower90.max()), 3) if len(base) else np.nan,",
        "                     may_fill=H.may_fill(r), top_families=', '.join(top), flags=';'.join(r.flags)))",
        "scores = pd.DataFrame(rows).sort_values('test_auc', ascending=False)",
        "scores.to_csv(os.path.join(OUT, 'scores.tsv'), sep='\\t', index=False)",
        "scores")

    nb.md("## Sensitivity: genome-wide evidence only",
          "Some evidence families were measured on a chosen subset of genes -- `gra12` is a CRISPR "
          "screen of 236 genes, a library of secretory candidates. Whether a gene was measured at all "
          "then encodes the choice someone made, which is prior knowledge of location by another door: "
          "the trees can read the missingness, not the measurement. The same questions are run again "
          "with only the families measured on at least half the genes, on a fresh seed (a different "
          "strategy may not look at the first partition). A niche that holds up here is placed by "
          "measurement; one that collapses was placed by who chose to measure it.")
    nb.code(
        "coverage = nodes[numeric].notna().mean()",
        "family_cover = coverage.groupby(families.reindex(numeric).values).max()",
        "WIDE = sorted(f for f, c in family_cover.items() if c >= 0.5 and f in set(families.reindex(features.columns)))",
        "print(len(WIDE), 'genome-wide families of', families.reindex(features.columns).nunique())",
        "rows = []",
        "for n in NICHES:",
        "    q = question(n)",
        "    y = (niche == n)",
        "    r = H.run_question(q, labels=y, features=features, families=families, groups=groups,",
        "                       strategies=[H.Strategy.of('boosted trees, genome-wide families', 'hist_gbt', families=WIDE, max_iter=150)],",
        "                       seed=1, importance_repeats=2, log_path=os.path.join(REPO, H.HOLDOUT_LOG))",
        "    cv, ft = r.cross_validated, r.frozen_test or {}",
        "    rows.append(dict(niche=n, wide_cv_auc=round(cv['roc_auc'], 3), wide_test_auc=round(ft.get('roc_auc', np.nan), 3),",
        "                     wide_top_families=', '.join(r.importance.head(3).family.tolist()) if len(r.importance) else ''))",
        "wide = scores[['niche', 'cv_auc', 'test_auc']].merge(pd.DataFrame(rows), on='niche')",
        "wide['cv_change'] = (wide.wide_cv_auc - wide.cv_auc).round(3)",
        "wide.to_csv(os.path.join(OUT, 'genome_wide_sensitivity.tsv'), sep='\\t', index=False)",
        "wide")

    nb.md("## Reading it",
          "AUC says whether the evidence ranks a niche's proteins above the other markers; the "
          "frozen-test lower bound says whether any calls are good enough to write into the map. With "
          "a handful of positives per frozen partition the bound is wide by construction -- more "
          "markers, not a better model, is what tightens it. The top families show what carries each "
          "niche; a family that restates prior knowledge (a curated annotation) is the caveat above, "
          "showing.")
    nb.code(
        "summary = pd.DataFrame([dict(questions=len(scores), median_test_auc=float(scores.test_auc.median()),",
        "    above_0_8=int((scores.test_auc >= 0.8).sum()), may_fill=int(scores.may_fill.sum()),",
        "    median_genome_wide_cv_auc=float(wide.wide_cv_auc.median()), median_cv_change=float(wide.cv_change.median()))])",
        "summary.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\\t', index=False)",
        "summary")
    nb.write(NB)


if __name__ == "__main__":
    main()
