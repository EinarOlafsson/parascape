"""Mouse inference pilot (instruction 64.29): a local benchmark on mouse truth, from the mouse gene pack.

Pilot target: is a mouse protein detected on the surface of bone-marrow-derived macrophages (BMDM surface
proteomics, shipped in `mm_host_proteins.parquet`; 1,296 proteins tested, 150 detected). Evidence: mouse
expression only -- BMDM, brain and skeletal-muscle TPM from the shipped table and the Bgee summaries of the
protein's gene from the mouse pack (one reviewed link). Excluded: the surface intensity (the label's own source)
and every location column, including the locations transferred from human orthologs (the label's own category,
and inferred). The proteins are the units: the truth is measured on proteins. Then cached/fresh parity, a
cold-cache install of the pack, and coverage and gaps.

    python scripts/build_mouse_inference_pack_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"mouse_inference_pilot_{DATE}")
NB = os.path.join(REPO, "notebooks", f"mouse_inference_pilot_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Mouse inference pilot", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**The pilot (instruction 64.29).** Mouse performance measured on mouse truth, never borrowed from human "
          "or parasite results. Truth: detection on the BMDM surface by surface proteomics (proteins tested; "
          "untested proteins are unknown, not negative). Evidence: mouse expression only. Every location column is "
          "excluded -- including the locations transferred from human orthologs, which are both the label's own "
          "category and inferred.",
          "**Frozen scope (declared before running):** one question, boosted trees at harness defaults, seed 0, each "
          "protein its own group (no mouse homology groups -- a stated gap).",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, shutil, numpy as np, pandas as pd",
            "from parascape import harness as H, paths, space_catalogue as SC",
            "from parascape.question_registry import Question",
            "prot = pd.read_parquet(os.path.join(REPO, 'parascape/data/mm_host_proteins.parquet')).set_index('host_id')",
            "import zipfile, io",
            "pack = os.path.join(REPO, 'results/mouse_gene_space_2026_11_07/parascape_Mm_2026.11.07.zip')",
            "with zipfile.ZipFile(pack) as z:",
            "    genes = pd.read_parquet(io.BytesIO(z.read('mm_gene_nodes.parquet'))).set_index('gene_id')",
            "    links = pd.read_parquet(io.BytesIO(z.read('mm_protein_gene_links.parquet')))",
            "tested = prot[prot.bmdm_surface_detected.notna()]",
            "labels = tested.bmdm_surface_detected.astype(bool)",
            "print(len(tested), 'proteins tested,', int(labels.sum()), 'detected on the BMDM surface')")
    nb.md("## 1 · Evidence table (mouse expression only)")
    nb.code("feat = np.log1p(prot[['bmdm_tpm', 'brain_tpm', 'skeletal_muscle_tpm']].astype(float)).add_prefix('log_')",
            "one = links[links.link == 'one gene'].drop_duplicates('host_id').set_index('host_id').gene_id",
            "bgee = genes[['bgee_top_score', 'bgee_tissues_present']].reindex(one.values)",
            "bgee.index = one.index",
            "top = pd.get_dummies(genes.bgee_top_tissue.reindex(one.values), prefix='bgee_top').astype(float)",
            "top.index = one.index",
            "top = top.loc[:, top.sum() >= 20]",
            "feat = feat.join(bgee).join(top)",
            "families = pd.Series({c: 'expression' for c in feat.columns})",
            "assert not any('surface' in c or 'location' in c for c in feat.columns)",
            "print(feat.shape[1], 'evidence columns; proteins with Bgee through a gene link:', int(bgee.reindex(tested.index).notna().any(axis=1).sum()))")
    nb.md("## 2 · The question through the harness")
    nb.code("q = Question(id='mm_bmdm_surface', text='Is this mouse protein on the macrophage surface?', answer_type='boolean',",
            "             question_class='B', label_source='mm_host_proteins.parquet:bmdm_surface_detected', label_family='bmdm_surface', evidence_forbidden=('bmdm_surface',),",
            "             label_tier='observed-screen', label_category='place', status='scored', min_labels_for_frozen_test=150)",
            "strategy = [H.Strategy.of('boosted trees, mouse expression only', 'hist_gbt')]",
            "groups = pd.Series(tested.index, index=tested.index)",
            "r = H.run_question(q, labels=labels, features=feat.loc[tested.index], families=families, groups=groups, organism='Mm',",
            "                   strategies=strategy, seed=0)",
            "r.rows.to_parquet(os.path.join(OUT, 'rows.parquet'), index=False)",
            "summary = r.summary(); summary.to_csv(os.path.join(OUT, 'scores.tsv'), sep='\\t', index=False)",
            "print('flags:', r.flags, '; may fill:', r.may_fill)",
            "summary[['partition', 'n', 'positives', 'roc_auc', 'average_precision', 'calls', 'precision', 'precision_lower90']]")
    nb.code("r.baselines[['name', 'partition', 'roc_auc', 'precision']] if 'name' in r.baselines else r.baselines")
    nb.md("## 3 · Cached and fresh agree", "The same question re-run from scratch (flagged a repeat and logged).")
    nb.code("again = H.run_question(q, labels=labels, features=feat.loc[tested.index], families=families, groups=groups,",
            "                       organism='Mm', strategies=strategy, seed=0)",
            "parity = dict(identical_rows=bool(r.rows.reset_index(drop=True).equals(again.rows.reset_index(drop=True))), flags=list(again.flags))",
            "parity")
    nb.md("## 4 · Cold cache: the mouse pack installs and opens on its own")
    nb.code("state = os.path.join(OUT, '_cold'); shutil.rmtree(state, ignore_errors=True)",
            "os.environ[paths.ENV_STATE] = state",
            "try:",
            "    before = {e.code: e.state for e in SC.catalogue()}['Mm']",
            "    where = SC.install('Mm')",
            "    after = [SC.describe(e) for e in SC.catalogue() if e.code == 'Mm'][0]",
            "    back = pd.read_parquet(where / 'mm_gene_nodes.parquet')",
            "    cold = dict(state_before=before, after=after, genes=len(back), same_as_pack=bool(back.gene_id.equals(genes.reset_index().gene_id)))",
            "finally:",
            "    os.environ.pop(paths.ENV_STATE); shutil.rmtree(state, ignore_errors=True)",
            "cold")
    nb.md("## 5 · Coverage and gaps, stated")
    nb.code("cover = dict(proteins=len(prot), tested=len(tested), detected=int(labels.sum()),",
            "             gaps=['one observed mouse label in the pilot; the transferred locations are not truth',",
            "                   'no mouse homology groups: paralogs can sit on both sides of a split',",
            "                   'BMDM surface proteomics measures one cell type in culture; other tissues untested',",
            "                   'catalogue strategies other than the harness models need a mouse context adapter (not built)'])",
            "json.dump(dict(scores=summary.to_dict('records'), flags=list(r.flags), may_fill=r.may_fill, parity=parity, cold_cache=cold, coverage=cover),",
            "          open(os.path.join(OUT, 'summary.json'), 'w'), indent=1, default=str)",
            "cover")
    nb.write(NB)


if __name__ == "__main__":
    main()
