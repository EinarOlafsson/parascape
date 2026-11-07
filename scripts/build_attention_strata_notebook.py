"""Do predictions hold for genes nobody has studied? (instruction 64.16), as an executed notebook.

The harness reference strategy on `required_without_glucose`, with its frozen test read again within strata
of study effort (publications per gene). One fresh partition (seed 20), logged; nothing tuned on it.

    python scripts/build_attention_strata_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"attention_strata_{DATE}")
NB = os.path.join(REPO, "notebooks", f"attention_strata_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Do predictions hold for genes nobody has studied?", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**The question.** A model can look good on average by being right about well-studied genes. Here the "
          "frozen test of the glucose question is read again within three strata of study effort -- no "
          "publication, 1-4, 5 or more -- with the threshold and calibration chosen on train + tune applied "
          "unchanged (instruction 64.16). Publication counts are attention and are held out as evidence; they "
          "are used here only to split the test genes.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, numpy as np, pandas as pd",
            "from parascape import harness as H, question_registry as Q, datasets as ds, organisms as O",
            "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet')).set_index('gene_id')",
            "q = Q.get('required_without_glucose')",
            "y = (nodes['fit_no_glucose'] <= -2.0)[nodes['fit_no_glucose'].notna()]",
            "def family_of(c):",
            "    f = ds.provenance(c, O.TOXOPLASMA); return f.key if f is not None else 'unregistered'",
            "numeric = [c for c in nodes.columns if pd.api.types.is_numeric_dtype(nodes[c]) or pd.api.types.is_bool_dtype(nodes[c])]",
            "families = pd.Series({c: family_of(c) for c in numeric})",
            "X = H.allowed_features(q, nodes.loc[y.index, numeric].astype('float64'), families)",
            "groups = nodes['orthogroup'].reindex(y.index)",
            "groups = groups.where(groups.notna(), pd.Series({g: 'self:' + g for g in y.index}))",
            "pubs = nodes.loc[y.index, 'n_publications']",
            "strata = pd.cut(pubs, [-1, 0, 4, np.inf], labels=['no publication', '1-4 publications', '5 or more'])",
            "strata.value_counts()")
    nb.code("r = H.run_question(q, labels=y, features=X, families=families, groups=groups,",
            "                   strategies=[H.Strategy.of('boosted trees (reference)', 'hist_gbt', max_iter=150)],",
            "                   seed=20, importance_repeats=1, strata=strata.astype(str),",
            "                   log_path=os.path.join(REPO, H.HOLDOUT_LOG))",
            "table = pd.DataFrame(r.frozen_test['by_stratum']).T",
            "table.loc['all'] = {k: r.frozen_test[k] for k in table.columns if k in r.frozen_test}",
            "table.to_csv(os.path.join(OUT, 'by_stratum.tsv'), sep='\\t')",
            "table")
    nb.write(NB)


if __name__ == "__main__":
    main()
