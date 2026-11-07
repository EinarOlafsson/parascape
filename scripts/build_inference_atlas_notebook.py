"""The parasite inference atlases (instructions 64.26 Toxoplasma, 64.27 Plasmodium), as an executed notebook.

For each organism: the manifest (every catalogue method with its role, held-out rows for known genes, deployment
claims for unknown genes, its card, or a specific status), the canonical class hold-outs already cached in the
track record, and a reproducibility check -- randomly sampled cached records rebuilt fresh with the same frozen
configuration and compared call by call.

    python scripts/build_inference_atlas_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"inference_atlas_{DATE}")
NB = os.path.join(REPO, "notebooks", f"inference_atlas_{DATE}.ipynb")
SAMPLES = 4


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("The parasite inference atlases", namespace={"REPO": REPO, "OUT": OUT, "SAMPLES": SAMPLES})
    nb.md("**What this records (instructions 64.26, 64.27).** For Toxoplasma and Plasmodium, what each of the 39 "
          "catalogue methods has precomputed: held-out per-gene rows for genes with a known answer (each hidden "
          "with its orthogroup), deployment claims for genes without one, its calibration card -- or, where there "
          "is none, the specific reason and how it is run on demand. Then a check that the cache is real: a "
          f"random sample of {SAMPLES} cached (method, label) records per organism is rebuilt from scratch with the "
          "same frozen configuration and compared call by call.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, json, shutil, numpy as np, pandas as pd",
            "from parascape import atlas as AT, track_record as T, precompute_track_record as P",
            "manifests = {org: AT.manifest(org) for org in ('Tg', 'Pf')}",
            "for org, m in manifests.items():",
            "    m.to_csv(os.path.join(OUT, f'manifest_{org}.tsv'), sep='\\t', index=False)",
            "pd.DataFrame({org: AT.coverage(m) for org, m in manifests.items()})")
    nb.md("## 1 · Every method, per organism", "Role, what is stored for known and for unknown genes, the card.")
    nb.code("show = ['strategy', 'role', 'known_genes', 'unknown_genes', 'card']",
            "manifests['Tg'][show]")
    nb.code("manifests['Pf'][show]")
    nb.md("## 2 · Canonical class queries already cached",
          "The track record also holds every class of every label hidden as a whole set (mode `set`): the canonical "
          "'what would it say about this class' query, stored once and read by query identity.")
    nb.code("L = T.shipped()",
            "sets = L[L['mode'].astype(str) == 'set']",
            "canon = sets.groupby([sets.organism.astype(str), sets.target.astype(str)], observed=True).agg(",
            "    methods=('strategy', 'nunique'), set_queries=('set_name', 'nunique'), rows=('gene_id', 'size'))",
            "canon.to_csv(os.path.join(OUT, 'canonical_set_queries.tsv'), sep='\\t')",
            "canon")
    nb.md("## 3 · Do cached records match a fresh run?",
          f"{SAMPLES} (method, label) pairs per organism, drawn at random (seed 7), rebuilt with "
          "`precompute_track_record.build_partition` (the production computation) and compared with the shipped "
          "rows: same genes, same call (an abstention on both sides counts as the same), same accuracy. The shipped "
          "record dates from 2026-10-17; any code change since that alters calls shows up here as a mismatch, and "
          "is reported, not hidden.")
    nb.code("rng = np.random.default_rng(7)",
            "together = L[L['mode'].astype(str) == 'together']",
            "checks = []",
            "work = os.path.join(OUT, '_fresh'); shutil.rmtree(work, ignore_errors=True)",
            "for org in ('Tg', 'Pf'):",
            "    pairs = together[together.organism.astype(str) == org][['strategy', 'target']].astype(str).drop_duplicates()",
            "    pick = pairs.iloc[rng.choice(len(pairs), size=min(SAMPLES, len(pairs)), replace=False)]",
            "    for key, target in pick.itertuples(index=False):",
            "        part = P.plan(org, [target], [key])",
            "        P.run(part, os.path.join(work, org))",
            "        fresh = P.assemble(os.path.join(work, org), part)",
            "        fresh = fresh[fresh['mode'].astype(str) == 'together'] if len(fresh) else fresh",
            "        old = together[(together.organism.astype(str) == org) & (together.strategy.astype(str) == key) & (together.target.astype(str) == target)]",
            "        old = old[old.setting_key.astype(str) == old.setting_key.astype(str).iloc[0]]",
            "        j = old.merge(fresh, on='gene_id', suffixes=('_cached', '_fresh')) if len(fresh) else old.iloc[:0]",
            "        same = (j.prediction_cached.astype(str) == j.prediction_fresh.astype(str)) | (j.prediction_cached.isna() & j.prediction_fresh.isna())",
            "        checks.append(dict(organism=org, method=key, label=target, genes_cached=old.gene_id.nunique(),",
            "                           genes_fresh=fresh.gene_id.nunique() if len(fresh) else 0, same_call=round(float(same.mean()), 4) if len(j) else float('nan'),",
            "                           accuracy_cached=round(float(old.correct.fillna(False).mean()), 4),",
            "                           accuracy_fresh=round(float(fresh.correct.fillna(False).mean()), 4) if len(fresh) else float('nan')))",
            "checks = pd.DataFrame(checks); checks.to_csv(os.path.join(OUT, 'reproducibility.tsv'), sep='\\t', index=False)",
            "shutil.rmtree(work, ignore_errors=True)",
            "checks")
    nb.write(NB)


if __name__ == "__main__":
    main()
