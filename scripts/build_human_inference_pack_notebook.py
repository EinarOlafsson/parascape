"""Human inference pilot (instruction 64.28): local benchmarks on human truth, from the human gene pack.

Pilot targets: five Human Protein Atlas subcellular locations (Nucleoplasm, Plasma membrane, Mitochondria, Golgi
apparatus, Endoplasmic reticulum), truth only where HPA grades the location Enhanced, Supported or Approved.
Evidence: expression only (GTEx medians, HPA RNA specificity/distribution/single-cell specificity, and the shipped
infection responses of the gene's one reviewed protein). Location-like columns (secretome location, protein
class, which includes predicted membrane/secreted classes) are excluded as the label's own category. Each question
runs through the harness: nested CV, a frozen test touched once and logged. Then cached/fresh parity, a cold-cache
install of the pack, and coverage and gaps.

    python scripts/build_human_inference_pack_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"human_inference_pilot_{DATE}")
NB = os.path.join(REPO, "notebooks", f"human_inference_pilot_{DATE}.ipynb")
LOCATIONS = ["Nucleoplasm", "Plasma membrane", "Mitochondria", "Golgi apparatus", "Endoplasmic reticulum"]


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Human inference pilot", namespace={"REPO": REPO, "OUT": OUT, "LOCATIONS": LOCATIONS})
    nb.md("**The pilot (instruction 64.28).** Human performance measured on human truth, never borrowed from the "
          "parasites. Truth: Human Protein Atlas subcellular location where HPA grades it Enhanced, Supported or "
          "Approved (Uncertain and unannotated genes are left out). Evidence: expression only -- GTEx fibroblast "
          "and hepatocyte medians, HPA RNA tissue specificity, distribution and single-cell specificity, and the "
          "infection responses ParaScape ships for the gene's one reviewed protein. Secretome location and protein "
          "class carry location information (predicted membrane/secreted classes) and are excluded.",
          "**Frozen scope (declared before running):** five locations, one strategy (boosted trees, harness "
          "defaults), seed 0, each gene its own group (no human homology groups in the pack: paralogs may sit on "
          "both sides of a split -- a stated gap).",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, shutil, numpy as np, pandas as pd",
            "from parascape import harness as H, packs, paths, space_catalogue as SC",
            "from parascape.question_registry import Question",
            "genes = pd.read_parquet(os.path.join(REPO, 'results/human_gene_space_2026_11_07/hs_genes.parquet')).set_index('gene_id')",
            "links = pd.read_parquet(os.path.join(REPO, 'results/human_gene_space_2026_11_07/hs_protein_gene_links.parquet'))",
            "prot = pd.read_parquet(os.path.join(REPO, 'parascape/data/hs_host_proteins.parquet')).set_index('host_id')",
            "graded = genes[genes.HPA_subcellular_reliability.isin(['Enhanced', 'Supported', 'Approved'])]",
            "print(len(genes), 'genes;', len(graded), 'with a graded location')")
    nb.md("## 1 · Evidence table")
    nb.code("feat = pd.DataFrame(index=genes.index)",
            "for c in ('gtex_fibroblast_median_tpm', 'gtex_lcm_hepatocyte_median_tpm'):",
            "    feat['log_' + c] = np.log1p(genes[c])",
            "for c in ('HPA_RNA_tissue_specificity', 'HPA_RNA_tissue_distribution', 'HPA_RNA_singlecell_specificity'):",
            "    feat = feat.join(pd.get_dummies(genes[c], prefix=c).astype(float))",
            "one = links[links.link == 'one gene'].drop_duplicates('host_id')",
            "numeric = [c for c in prot.columns if pd.api.types.is_numeric_dtype(prot[c])]",
            "response = prot.loc[one.host_id, numeric]",
            "response.index = pd.Index(one.gene_id.astype(str).tolist())",
            "response = response[~response.index.duplicated()]",
            "feat = feat.join(response.add_prefix('infection_')).astype(float)  # nullable NA -> NaN for the models",
            "families = pd.Series({c: ('infection response' if c.startswith('infection_') else 'expression') for c in feat.columns})",
            "print(feat.shape[1], 'evidence columns:', families.value_counts().to_dict())")
    nb.md("## 2 · Five locations through the harness",
          "Each question: nested cross-validation on train + tune, calibration on its own partition, and a frozen "
          "test touched once (logged in `results/harness_holdout_log.jsonl`).")
    nb.code("rows, results, questions = [], {}, {}",
            "loc_sets = graded.HPA_subcellular_location.fillna('').str.split(r'[;,]').map(lambda v: {x.strip() for x in v})",
            "for loc in LOCATIONS:",
            "    labels = loc_sets.map(lambda s: loc in s).astype(bool)",
            "    entities = labels.index",
            "    q = Question(id='hs_location_' + loc.lower().replace(' ', '_'), text=f'Is this human protein in the {loc}?',",
            "                 answer_type='boolean', question_class='B', label_source='Human Protein Atlas subcellular location (graded)',",
            "                 label_family='hpa_location', evidence_forbidden=('hpa_location',), label_tier='observed-curated', label_category='place', status='scored')",
            "    questions[loc] = q",
            "    r = H.run_question(q, labels=labels, features=feat.loc[entities], families=families,",
            "                       groups=pd.Series(entities, index=entities), organism='Hs',",
            "                       strategies=[H.Strategy.of('boosted trees, expression only', 'hist_gbt')], seed=0)",
            "    results[loc] = r",
            "    ft = r.frozen_test or {}",
            "    base = r.baselines[(r.baselines.partition == 'frozen_test')] if len(r.baselines) else r.baselines",
            "    rows.append(dict(location=loc, genes=len(labels), positives=int(labels.sum()), cv_auc=round(r.cross_validated['roc_auc'], 3),",
            "                     frozen_auc=round(ft.get('roc_auc', float('nan')), 3), frozen_ap=round(ft.get('average_precision', float('nan')), 3),",
            "                     prevalence=round(float(labels.mean()), 3), precision_lower90=round(ft.get('precision_lower90', float('nan')), 3),",
            "                     may_fill=r.may_fill, flags=', '.join(r.flags)))",
            "    r.rows.assign(location=loc).to_parquet(os.path.join(OUT, f\"rows_{loc.replace(' ', '_')}.parquet\"), index=False)",
            "scores = pd.DataFrame(rows); scores.to_csv(os.path.join(OUT, 'scores.tsv'), sep='\\t', index=False)",
            "scores")
    nb.md("## 3 · Cached and fresh agree",
          "One question re-run from scratch with the identical strategy (the harness flags it a repeat and logs it): "
          "the per-gene rows must be identical.")
    nb.code("loc = LOCATIONS[2]",
            "labels = loc_sets.map(lambda s: loc in s).astype(bool)",
            "again = H.run_question(questions[loc], labels=labels, features=feat.loc[labels.index], families=families,",
            "                       groups=pd.Series(labels.index, index=labels.index), organism='Hs',",
            "                       strategies=[H.Strategy.of('boosted trees, expression only', 'hist_gbt')], seed=0)",
            "a, b = results[loc].rows.reset_index(drop=True), again.rows.reset_index(drop=True)",
            "parity = dict(location=loc, identical_rows=bool(a.equals(b)), flags=list(again.flags))",
            "parity")
    nb.md("## 4 · Cold cache: the pack installs and opens on its own")
    nb.code("state = os.path.join(OUT, '_cold'); shutil.rmtree(state, ignore_errors=True)",
            "os.environ[paths.ENV_STATE] = state",
            "try:",
            "    before = {e.code: e.state for e in SC.catalogue()}",
            "    where = SC.install('Hs')",
            "    after = {e.code: SC.describe(e) for e in SC.catalogue() if e.code == 'Hs'}",
            "    back = pd.read_parquet(where / 'hs_gene_nodes.parquet')",
            "    cold = dict(state_before=before['Hs'], after=after['Hs'], genes=len(back), same_as_built=bool(back.gene_id.equals(genes.reset_index().gene_id)))",
            "finally:",
            "    os.environ.pop(paths.ENV_STATE); shutil.rmtree(state, ignore_errors=True)",
            "cold")
    nb.md("## 5 · Coverage and gaps, stated")
    nb.code("cover = dict(genes=len(genes), graded_location=len(graded), with_infection_response=int(feat.filter(like='infection_').notna().any(axis=1).sum()),",
            "             questions=len(LOCATIONS), fill_allowed=int(scores.may_fill.sum()),",
            "             gaps=['no human homology groups: paralogs can sit on both sides of a split',",
            "                   'only expression evidence: no human interaction, structure or domain columns in the pack yet',",
            "                   'HPA locations are cell-line immunofluorescence; tissue context differs',",
            "                   'catalogue strategies other than the harness models need a human context adapter (not built)'])",
            "json.dump(dict(scores=rows, parity=parity, cold_cache=cold, coverage=cover), open(os.path.join(OUT, 'summary.json'), 'w'), indent=1, default=str)",
            "cover")
    nb.write(NB)


if __name__ == "__main__":
    main()
