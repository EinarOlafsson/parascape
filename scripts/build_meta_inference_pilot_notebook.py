"""Meta-inference pilot on one label (instruction 64.34), as an executed notebook.

Does combining what 14 methods said about a gene beat the best single method and the plain plurality vote?
Base calls are the hold-out record (each gene's call made with its own answer hidden); the meta-model is fitted
on outer-training orthogroups and tested on untouched outer orthogroups. Stated limit: the base calls come from the
record's own folds, so a base call for an outer-training gene could come from a model that saw an outer-test
gene's label -- a strictly nested version recomputes the base methods inside each outer fold.

    python scripts/build_meta_inference_pilot_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"meta_inference_pilot_{DATE}")
NB = os.path.join(REPO, "notebooks", f"meta_inference_pilot_{DATE}.ipynb")
TARGET = "dtm_class"


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Meta-inference pilot: do 14 methods together beat the best one?",
                          namespace={"REPO": REPO, "OUT": OUT, "TARGET": TARGET})
    nb.md(f"**The pilot (instruction 64.34).** Label `{TARGET}` (Toxoplasma). Each gene has a held-out call from up "
          "to 14 methods. A meta-model (multinomial logistic regression on the one-hot calls) is fitted on 4/5 of "
          "the orthogroups and tested on the other 1/5, five times; compared on the same test genes with the best "
          "single method (chosen on the training orthogroups), the plurality vote, and the commonest class.",
          "**Stated limit.** The base calls come from the hold-out record's own folds; a strictly nested version "
          "would recompute every base method inside each outer fold. This pilot says whether that is worth doing.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, numpy as np, pandas as pd",
            "from sklearn.linear_model import LogisticRegression",
            "from sklearn.model_selection import GroupKFold",
            "from scipy.stats import binomtest",
            "from parascape import track_record as T",
            "L = T.shipped('Tg')",
            "L = L[(L['target'].astype(str) == TARGET) & (L['mode'].astype(str) == 'together')].copy()",
            "first = L.groupby(L['strategy'].astype(str), observed=True)['setting_key'].first()",
            "L = L[L['setting_key'].astype(str).values == L['strategy'].astype(str).map(first).astype(str).values]",
            "calls = L.assign(call=np.where(L['abstained'].astype(bool), 'ABSTAIN', L['prediction'].astype(str))) \\",
            "         .pivot_table(index='gene_id', columns='strategy', values='call', aggfunc='first', observed=True).fillna('ABSTAIN')",
            "truth = L.groupby('gene_id', observed=True)['truth'].first().astype(str).reindex(calls.index)",
            "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet')).set_index('gene_id')",
            "groups = nodes['orthogroup'].reindex(calls.index).fillna(pd.Series(calls.index, index=calls.index))",
            "X = pd.get_dummies(calls.astype(str)).astype(float)",
            "print(calls.shape[0], 'genes,', calls.shape[1], 'methods,', truth.nunique(), 'classes')")
    nb.code("rows, paired = [], []",
            "for fold, (tr, te) in enumerate(GroupKFold(5).split(X, truth, groups)):",
            "    ytr, yte = truth.iloc[tr], truth.iloc[te]",
            "    meta = LogisticRegression(max_iter=3000, C=1.0).fit(X.iloc[tr], ytr)",
            "    p_meta = meta.predict(X.iloc[te])",
            "    acc_tr = {m: float((calls.iloc[tr][m] == ytr).mean()) for m in calls.columns}",
            "    best = max(acc_tr, key=acc_tr.get)",
            "    p_best = calls.iloc[te][best].to_numpy()",
            "    plural = calls.iloc[te].apply(lambda r: r[r != 'ABSTAIN'].value_counts().index[0] if (r != 'ABSTAIN').any() else 'ABSTAIN', axis=1).to_numpy()",
            "    common = np.full(len(te), ytr.value_counts().index[0])",
            "    y = yte.to_numpy()",
            "    for name, p in (('meta-model', p_meta), (f'best single ({best})', p_best), ('plurality vote', plural), ('commonest class', common)):",
            "        rows.append(dict(fold=fold, method=name.split(' (')[0], chosen=best if name.startswith('best') else '', accuracy=float((p == y).mean()), n=len(y)))",
            "    paired.append(dict(fold=fold, meta_only=int(((p_meta == y) & (p_best != y)).sum()), best_only=int(((p_best == y) & (p_meta != y)).sum())))",
            "res = pd.DataFrame(rows)",
            "summary = res.groupby('method').agg(accuracy=('accuracy', 'mean'), low=('accuracy', 'min'), high=('accuracy', 'max')).round(3).sort_values('accuracy', ascending=False)",
            "pp = pd.DataFrame(paired).sum()",
            "test = binomtest(int(pp.meta_only), int(pp.meta_only + pp.best_only), 0.5)",
            "summary.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\\t')",
            "res.to_csv(os.path.join(OUT, 'folds.tsv'), sep='\\t', index=False)",
            "print('meta right where the best single is wrong:', int(pp.meta_only), '; the reverse:', int(pp.best_only), '; McNemar p', round(test.pvalue, 4))",
            "summary")
    nb.write(NB)


if __name__ == "__main__":
    main()
