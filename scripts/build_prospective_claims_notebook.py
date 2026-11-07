"""Frozen prospective claims (instruction 63): predictions fixed today, to be graded by data that does not exist yet.

The glucose-withdrawal CRISPR screen measured 7,393 genes; about 750 it did not. A model fitted on every
measured gene ranks the unmeasured ones, and the ranking is frozen -- written, hashed, dated -- before any
new screen can be read. When a later genome-wide screen covers those genes, the frozen file is graded
against it unchanged: the only test that no tuning, however careful, can have seen.

None of these is a fill: the harness found no strategy that clears the 0.90 bar on this question. They
are a ranking with probabilities, offered as a prediction to check.

    python scripts/build_prospective_claims_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_06"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"prospective_claims_{DATE}")
NB = os.path.join(REPO, "notebooks", f"prospective_claims_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Frozen prospective claims: genes the glucose screen did not measure",
                          namespace={"REPO": REPO, "OUT": OUT})
    nb.md(
        "**What this is.** A ranking of the genes the glucose-withdrawal screen did not measure, by how "
        "likely each is to be needed for growth without glucose -- written down, hashed and dated now, "
        "so a future screen can grade it without anyone having been able to adjust it (instruction 63).",
        "**How it was made.** The reference strategy of the harness (boosted trees over every allowed "
        "evidence family; every pooled fitness screen held out, as the question requires) fitted on all "
        "7,393 measured genes, then applied to the unmeasured ones. No frozen test partition is touched: "
        "this is the deployment fit, and the strategy's held-out record is in "
        "`results/borrowed_strategies_2026_11_06/` (frozen-test AUC 0.824).",
        "**What it is not.** A fill. No strategy clears the 0.90 bar on this question, so these are "
        "predictions to check, never facts in the map.",
        f"Frozen {_dt.date.today().isoformat()}.")
    nb.code(
        "import hashlib, json, os, time",
        "import numpy as np, pandas as pd",
        "from parascape import harness as H, question_registry as Q, datasets as ds, organisms as O",
        "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet')).set_index('gene_id')",
        "question = Q.get('required_without_glucose')",
        "measured = nodes['fit_no_glucose'].notna()",
        "labels = (nodes['fit_no_glucose'] <= -2.0)[measured]",
        "def family_of(column):",
        "    found = ds.provenance(column, O.TOXOPLASMA)",
        "    return found.key if found is not None else 'unregistered'",
        "numeric = [c for c in nodes.columns if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]",
        "families = pd.Series({c: family_of(c) for c in numeric})",
        "features = H.allowed_features(question, nodes[numeric].astype('float64'), families)",
        "print(int(measured.sum()), 'measured genes;', int((~measured).sum()), 'unmeasured;',",
        "      features.shape[1], 'allowed columns')")
    nb.code(
        "strategy = H.Strategy.of('boosted trees (reference)', 'hist_gbt', max_iter=150)",
        "usable = [c for c in features.columns if features.loc[labels.index, c].notna().any()]",
        "model = strategy.build().fit(features.loc[labels.index, usable], labels.astype(int))",
        "unmeasured = nodes.index[~measured]",
        "p = model.predict_proba(features.loc[unmeasured, usable])[:, 1]",
        "claims = pd.DataFrame({'gene_id': unmeasured, 'question': question.id,",
        "                       'probability_needed': np.round(p, 4)}).sort_values('probability_needed', ascending=False)",
        "claims['rank'] = np.arange(1, len(claims) + 1)",
        "claims['product'] = nodes.loc[claims.gene_id, 'product'].values if 'product' in nodes else ''",
        "claims.head(15)")
    nb.md("## The freeze",
          "The claims file, its SHA-256, the strategy and evidence hashes, the version and the date go "
          "into `manifest.json`. Grading later means: recompute the hash, confirm it matches, then score "
          "the ranking against the new screen -- precision in the top 50, 100 and 200, and AUC.")
    nb.code(
        "from parascape import __version__",
        "path = os.path.join(OUT, 'claims.tsv')",
        "claims.to_csv(path, sep='\\t', index=False)",
        "digest = hashlib.sha256(open(path, 'rb').read()).hexdigest()",
        "manifest = dict(question=question.id, label=question.label_source, frozen_at=time.strftime('%Y-%m-%dT%H:%M:%S'),",
        "                parascape_version=__version__, strategy=strategy.name, strategy_hash=strategy.identity,",
        "                evidence_hash=H.evidence_identity(question, features), n_training=int(len(labels)),",
        "                n_training_positive=int(labels.sum()), n_claims=int(len(claims)),",
        "                claims_sha256=digest, held_out_record='results/borrowed_strategies_2026_11_06/scores.tsv',",
        "                grade_with='precision in the top 50/100/200 and AUC against the next genome-wide '",
        "                           'glucose-withdrawal (or glycolysis-null) screen covering these genes')",
        "json.dump(manifest, open(os.path.join(OUT, 'manifest.json'), 'w'), indent=1)",
        "manifest")
    nb.write(NB)


if __name__ == "__main__":
    main()
