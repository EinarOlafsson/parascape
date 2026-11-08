"""Frozen prospective localization claims (instruction 63), as an executed notebook.

The companion of the glucose freeze, for the question Einar cares about most: where is a protein? The unified
hyperLOPIT assignment covers 3,827 of the 8,140 genes; the other 4,313 have none. A model per niche, fitted on
the assigned genes with every place column held out, ranks the unassigned ones, and the ranking is frozen --
written, hashed, dated -- before any new localisation data can be read.

None of these is a fill: no localization strategy clears the 0.90 bar (`question_registry`: a frozen partition
holds 4-13 positives). They are predictions to grade, unchanged, against the next independent localisation
source -- a new hyperLOPIT, a tagging campaign, or curated per-protein papers.

    python scripts/build_prospective_localization_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"prospective_localization_{DATE}")
NB = os.path.join(REPO, "notebooks", f"prospective_localization_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Frozen prospective localization claims: the genes hyperLOPIT did not place",
                          namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**What this is.** For every gene with no hyperLOPIT assignment, the niche the data say it is most "
          "likely in, with a probability -- written down, hashed and dated now, so an independent localisation "
          "source can grade it later without anyone having been able to adjust it (instruction 63).",
          "**How it was made.** One model per niche (boosted trees over every allowed evidence family), fitted "
          "on the assigned genes with **every column of the place category held out**, so no localisation "
          "measurement can restate the answer. The deployment fit; no frozen partition is touched.",
          "**What it is not.** A fill. No localization question clears the 0.90 bar, so these are predictions "
          "to check, never facts in the map.",
          f"Frozen {_dt.date.today().isoformat()}.")
    nb.code("import hashlib, json, os, time",
            "import numpy as np, pandas as pd",
            "from parascape import harness as H, datasets as ds, organisms as O",
            "from parascape.question_registry import Question",
            "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet')).set_index('gene_id')",
            "assigned = nodes['lopit_unified'].notna()",
            "niches = nodes.loc[assigned, 'lopit_unified'].value_counts()",
            "def family_of(column):",
            "    found = ds.provenance(column, O.TOXOPLASMA)",
            "    return found.key if found is not None else 'unregistered'",
            "numeric = [c for c in nodes.columns if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]",
            "families = pd.Series({c: family_of(c) for c in numeric})",
            "print(int(assigned.sum()), 'assigned genes;', int((~assigned).sum()), 'unassigned;', len(niches), 'niches')",
            "niches.to_frame('assigned genes')")
    nb.md("## 1 · One model per niche, place held out",
          "A niche with fewer than 50 assigned genes is not modelled and is named. The held-out record of this "
          "question family is in `question_registry` (frozen-test AUC 0.846-0.933 per niche, with every place "
          "column held out); this notebook is the deployment fit on all assigned genes.")
    nb.code("MIN = 50",
            "unassigned = nodes.index[~assigned]",
            "strategy = H.Strategy.of('boosted trees (reference)', 'hist_gbt', max_iter=150)",
            "columns, skipped, scores = None, [], {}",
            "for niche, n in niches.items():",
            "    if n < MIN:",
            "        skipped.append(f'{niche} ({n})')",
            "        continue",
            "    q = Question(id=f'prospective_localization_{niche}', text=f'Is this protein in the {niche}?',",
            "                 answer_type='boolean', question_class='B', label_source='nodes.parquet:lopit_unified',",
            "                 label_family=family_of('lopit_unified'), evidence_forbidden=(family_of('lopit_unified'),),",
            "                 label_tier='observed-curated', label_category='place', status='scored')",
            "    features = H.allowed_features(q, nodes[numeric].astype('float64'), families)",
            "    labels = (nodes.loc[assigned, 'lopit_unified'] == niche)",
            "    usable = [c for c in features.columns if features.loc[labels.index, c].notna().any()]",
            "    model = strategy.build().fit(features.loc[labels.index, usable], labels.astype(int))",
            "    scores[niche] = model.predict_proba(features.loc[unassigned, usable])[:, 1]",
            "    columns = usable",
            "print(len(scores), 'niches modelled;', features.shape[1], 'allowed columns of', len(numeric))",
            "print('too few assigned genes to model:', ', '.join(skipped) if skipped else 'none')",
            "P = pd.DataFrame(scores, index=unassigned)",
            "P.round(3).head()")
    nb.md("## 2 · The claim for each unplaced gene",
          "The niche with the highest probability, and how far it is ahead of the next one -- the margin is the "
          "part to read: a gene whose best two niches are a hair apart is not placed by this evidence.")
    nb.code("order = np.argsort(-P.to_numpy(), axis=1)",
            "best = P.columns.to_numpy()[order[:, 0]]",
            "top = np.take_along_axis(P.to_numpy(), order[:, :2], axis=1)",
            "claims = pd.DataFrame({'gene_id': unassigned.astype(str), 'niche': best,",
            "                       'probability': np.round(top[:, 0], 4), 'margin': np.round(top[:, 0] - top[:, 1], 4),",
            "                       'second_niche': P.columns.to_numpy()[order[:, 1]]})",
            "claims['product'] = nodes.loc[claims.gene_id, 'product'].values if 'product' in nodes else ''",
            "claims = claims.sort_values('probability', ascending=False).reset_index(drop=True)",
            "claims['rank'] = np.arange(1, len(claims) + 1)",
            "print(len(claims), 'claims;', int((claims.margin >= 0.2).sum()), 'with a clear best niche (margin >= 0.2)')",
            "claims.niche.value_counts().to_frame('claims')")
    nb.md("## 3 · The freeze",
          "The claims file, its SHA-256, the strategy and evidence hashes, the version and the date go into "
          "`manifest.json`. Grading later means: recompute the hash, confirm it matches, then score the claims "
          "against the new source -- accuracy on the genes it places, against the commonest niche as chance, and "
          "accuracy among the claims with a margin of at least 0.2.")
    nb.code("from parascape import __version__",
            "path = os.path.join(OUT, 'claims.tsv')",
            "claims.to_csv(path, sep='\\t', index=False)",
            "digest = hashlib.sha256(open(path, 'rb').read()).hexdigest()",
            "manifest = dict(question='prospective_localization (one model per niche)',",
            "                label='nodes.parquet:lopit_unified', frozen_at=time.strftime('%Y-%m-%dT%H:%M:%S'),",
            "                parascape_version=__version__, strategy=strategy.name, strategy_hash=strategy.identity,",
            "                niches=sorted(scores), niches_not_modelled=skipped, n_training=int(assigned.sum()),",
            "                n_claims=int(len(claims)), n_clear_margin=int((claims.margin >= 0.2).sum()),",
            "                evidence_columns=len(columns), claims_sha256=digest,",
            "                held_out_record='parascape/question_registry.py: LOCALIZATION_NICHES (frozen-test AUC per niche)',",
            "                grade_with='accuracy against the commonest-niche baseline on the genes an independent '",
            "                           'localisation source places (new hyperLOPIT, tagging campaign or curated '",
            "                           'per-protein papers), overall and among claims with margin >= 0.2')",
            "json.dump(manifest, open(os.path.join(OUT, 'manifest.json'), 'w'), indent=1)",
            "manifest")
    nb.write(NB)


if __name__ == "__main__":
    main()
