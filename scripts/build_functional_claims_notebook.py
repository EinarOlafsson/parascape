"""Precomputed functional claims (instruction 68.02), as an executed notebook.

For each parasite and each EC major class (1 oxidoreductase ... 6 ligase) with at least 50 annotated genes: is an
annotated gene in that class? Graded through the harness (orthogroup-grouped nested cross-validation, the frozen
test touched once), then deployed on the genes with NO EC annotation, which become ranked, calibrated functional
claims. Evidence that could restate how an EC number is assigned is held out: every function and domain column, and
the EC, InterPro, Pfam and EC-by-orthology families. A claim is a FILL only where its class may fill (frozen-test
lower bound >= 0.90) and the gene clears the calibrated threshold; every other claim is exploratory, ranked, never a
fact, never a training label.

    python scripts/build_functional_claims_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"functional_claims_{DATE}")
NB = os.path.join(REPO, "notebooks", f"functional_claims_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Precomputed functional claims", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**What this precomputes (instruction 68.02).** For genes with no EC annotation in Toxoplasma and "
          "Plasmodium: the probability that each EC major class applies, from evidence that cannot restate how EC "
          "numbers are assigned, each with the held-out grade of its class.",
          "**Held out as evidence:** every column of the function and domain categories, and the EC, InterPro, Pfam "
          "and EC-by-orthology families -- EC assignments are commonly made from domains and orthologs, so those "
          "inputs would partly restate the answer.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, re, numpy as np, pandas as pd",
            "from parascape import datasets as ds, harness as H, information_categories as IC, organisms as O",
            "from parascape.question_registry import Question",
            "EC = re.compile(r'\\b([1-7])\\.(?:\\d+|-)\\.(?:\\d+|-)\\.(?:n?\\d+|-)')",
            "NAMES = {'1': 'oxidoreductase', '2': 'transferase', '3': 'hydrolase', '4': 'lyase', '5': 'isomerase', '6': 'ligase', '7': 'translocase'}",
            "def majors(v):",
            "    return sorted({m.group(1) for m in EC.finditer(str(v))}) if isinstance(v, str) else []")
    nb.md("## 1 · Per organism and class: grade, then deploy")
    nb.code("claims, grades = [], []",
            "for org, table in ((O.TOXOPLASMA, 'nodes.parquet'), (O.FALCIPARUM, 'pf_nodes.parquet')):",
            "    nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data', table)).set_index('gene_id')",
            "    numeric = [c for c in nodes.columns if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]",
            "    fam = lambda c, org=org: (ds.provenance(c, org).key if ds.provenance(c, org) else 'unregistered')",
            "    families = pd.Series({c: fam(c) for c in numeric})",
            "    def blocked(c):",
            "        try:",
            "            cats = set(IC.categories_of(c, org))",
            "        except Exception:",
            "            cats = set()",
            "        f = families.get(c, '')",
            "        return bool(cats & {'function', 'domain'}) or any(k in f.lower() for k in ('ec', 'interpro', 'pfam', 'domain'))",
            "    keep = [c for c in numeric if not blocked(c)]",
            "    groups = nodes['orthogroup'].astype(object) if 'orthogroup' in nodes else pd.Series(nodes.index, index=nodes.index)",
            "    groups = groups.where(groups.notna(), pd.Series({g: 'self:' + g for g in nodes.index}))",
            "    sets = nodes['ec_number'].map(majors)",
            "    annotated = sets[sets.map(len) > 0]",
            "    unannotated = nodes.index[sets.map(len) == 0]",
            "    ec_family = fam('ec_number') if 'ec_number' in families else 'ec'",
            "    print(org, len(annotated), 'annotated;', len(unannotated), 'without EC;', len(keep), 'evidence columns of', len(numeric))",
            "    for k in sorted({m for s in annotated for m in s}):",
            "        labels = annotated.map(lambda s, k=k: k in s)",
            "        if labels.sum() < 50 or (~labels).sum() < 50:",
            "            continue",
            "        q = Question(id=f'function_ec{k}_{org}', text=f'Is this {org} enzyme a {NAMES[k]} (EC {k})?', answer_type='boolean',",
            "                     question_class='B', label_source=f'{table}:ec_number', label_family=ec_family,",
            "                     evidence_forbidden=(ec_family,), label_tier='observed-curated', label_category='function', status='scored')",
            "        feats = nodes.loc[labels.index, keep].astype('float64')",
            "        feats = H.allowed_features(q, feats, families)",
            "        r = H.run_question(q, labels=labels.astype(bool), features=feats, families=families, groups=groups.reindex(labels.index),",
            "                           strategies=[H.Strategy.of('boosted trees', 'hist_gbt', max_iter=150)], seed=0, organism=org,",
            "                           importance_repeats=1, log_path=H.HOLDOUT_LOG)",
            "        ft = r.frozen_test or {}",
            "        grades.append(dict(organism=org, ec_class=k, name=NAMES[k], annotated=len(labels), positives=int(labels.sum()),",
            "                           frozen_auc=ft.get('roc_auc'), average_precision=ft.get('average_precision'), share=float(labels.mean()),",
            "                           precision_lower90=ft.get('precision_lower90'), may_fill=bool(r.may_fill),",
            "                           top_evidence=', '.join(r.importance.head(3)['family'].tolist()) if len(r.importance) else ''))",
            "        dep = H.deploy(r, nodes.loc[unannotated, feats.columns].astype('float64'))",
            "        claims.append(pd.DataFrame({'organism': org, 'gene_id': unannotated.astype(str), 'ec_class': k, 'name': NAMES[k],",
            "                                    'probability': dep.probability.round(4).values, 'calibrated': dep.calibrated.round(4).values,",
            "                                    'called': dep.called.values, 'fills': dep.fills.values,",
            "                                    'class_frozen_auc': round(ft.get('roc_auc', np.nan), 3)}))",
            "G = pd.DataFrame(grades); C = pd.concat(claims, ignore_index=True)",
            "G.to_csv(os.path.join(OUT, 'grades.tsv'), sep='\\t', index=False)",
            "G.round(3)")
    nb.md("## 2 · The claims",
          "One row per unannotated gene and class. `called`: yes, at the calibrated threshold that made the class "
          "held-out precision. `fills`: called AND the class may fill. Ranked best first within organism.")
    nb.code("C = C.sort_values(['organism', 'calibrated'], ascending=[True, False]).reset_index(drop=True)",
            "C.to_parquet(os.path.join(OUT, 'functional_claims.parquet'), index=False)",
            "C.to_parquet(os.path.join(REPO, 'parascape/data/functional_claims.parquet'), index=False)",
            "summary = C.groupby(['organism', 'ec_class', 'name']).agg(genes=('gene_id', 'size'), called=('called', 'sum'), fills=('fills', 'sum')).reset_index()",
            "summary.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\\t', index=False)",
            "json.dump(dict(claims=len(C), called=int(C.called.sum()), fills=int(C.fills.sum()), classes=len(G)),",
            "          open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "summary")
    nb.write(NB)


if __name__ == "__main__":
    main()
