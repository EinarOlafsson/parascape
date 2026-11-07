"""Strictly nested meta-inference on one label (instruction 64.34 admission), as an executed notebook.

For each outer fold of whole orthogroups, every base method is asked about the outer-test genes with their labels
hidden, and about the outer-training genes through inner folds that never see an outer-test label. The meta-model
is fitted on the inner calls of the training genes only and scored on the outer-test genes: no base call it learns
from, and none it is tested on, came from a model that saw an outer-test label.

    python scripts/build_meta_inference_nested_notebook.py              # the notebook (runs workers as needed)
    python scripts/build_meta_inference_nested_notebook.py --strategy K # one base method, all folds (worker)
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARGET = sys.argv[sys.argv.index("--target") + 1] if "--target" in sys.argv else "dtm_class"
_SUFFIX = "" if TARGET == "dtm_class" else f"_{TARGET}"
OUT = os.path.join(REPO, "results", f"meta_inference_nested_{DATE}{_SUFFIX}")
NB = os.path.join(REPO, "notebooks", f"meta_inference_nested_{DATE}{_SUFFIX}.ipynb")
OUTER, INNER, SEED = 5, 3, 100


def worker(key: str) -> None:
    """All outer and inner calls of one base method, written as one parquet file."""
    import numpy as np
    import pandas as pd
    sys.path.insert(0, REPO)
    from parascape import strategies as S, track_record as T
    ctx = S.Context.shipped("Tg")
    strategy = S.get(key)
    truth = ctx.truth(TARGET)
    kept = truth.where(truth.map(truth.value_counts()) >= S.MIN_CLASS)
    labelled = np.flatnonzero(kept.notna().to_numpy())
    outer = T.folds_by_group(ctx, labelled, OUTER, SEED)
    chosen = T._permitted(ctx, dict(strategy.settings(ctx, target=TARGET)), TARGET)
    rows = []

    def ask(visible, query, fold, role):
        pred, _support = T._predict(ctx, strategy, visible, query, chosen, TARGET)
        if pred is None:
            return False
        p = pd.Series(pred).reset_index(drop=True)
        for g in query:
            call = p.iloc[g]
            rows.append(dict(outer_fold=fold, role=role, gene=int(g), gene_id=str(ctx.gene_ids[g]),
                             truth=str(kept.iloc[g]), call=call if isinstance(call, str) else "ABSTAIN"))
        return True

    for k in range(OUTER):
        test = labelled[outer == k]
        train = labelled[outer != k]
        visible_outer = kept.copy()
        visible_outer.iloc[test] = np.nan
        if not ask(visible_outer, test, k, "outer_test"):
            break
        inner = T.folds_by_group(ctx, train, INNER, SEED + 1 + k)
        for j in range(INNER):
            hidden = train[inner == j]
            visible = visible_outer.copy()
            visible.iloc[hidden] = np.nan
            ask(visible, hidden, k, "inner_train")
    pd.DataFrame(rows).to_parquet(os.path.join(OUT, f"calls_{key}.parquet"))


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Strictly nested meta-inference: admitted or not",
                          namespace={"REPO": REPO, "OUT": OUT, "TARGET": TARGET})
    nb.md(f"**The admission test (instruction 64.34).** Label `{TARGET}`. {OUTER} outer folds of whole orthogroups; "
          f"inside each, the base methods' calls on the training genes come from {INNER} inner folds that never "
          "see an outer-test label, and their calls on the test genes from models fitted with the test labels "
          "hidden. The meta-model learns from the inner calls only and is scored once on the outer-test genes, "
          "against the best single method chosen on the training genes, the plurality vote and the commonest class.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, subprocess, sys",
            "from concurrent.futures import ThreadPoolExecutor",
            "import numpy as np, pandas as pd",
            "from sklearn.linear_model import LogisticRegression",
            "from scipy.stats import binomtest",
            "from parascape import track_record as T",
            "KEYS = T.supported('Tg')",
            "SCRIPT = os.path.join(REPO, 'scripts', 'build_meta_inference_nested_notebook.py')",
            "def run(key):",
            "    out = os.path.join(OUT, f'calls_{key}.parquet')",
            "    if not os.path.exists(out):",
            "        subprocess.run([sys.executable, SCRIPT, '--strategy', key, '--target', TARGET], check=False,",
            "                       env={**os.environ, 'OMP_NUM_THREADS': '2'})",
            "    return pd.read_parquet(out).assign(strategy=key) if os.path.exists(out) else None",
            "with ThreadPoolExecutor(max_workers=7) as pool:",
            "    parts = [p for p in pool.map(run, KEYS) if p is not None and len(p)]",
            "calls = pd.concat(parts, ignore_index=True)",
            "print(calls.strategy.nunique(), 'base methods;', calls.gene_id.nunique(), 'genes')")
    nb.code("rows, paired = [], dict(meta_only=0, best_only=0)",
            "for k in sorted(calls.outer_fold.unique()):",
            "    c = calls[calls.outer_fold == k]",
            "    tr = c[c.role == 'inner_train'].pivot_table(index='gene_id', columns='strategy', values='call', aggfunc='first').fillna('ABSTAIN')",
            "    te = c[c.role == 'outer_test'].pivot_table(index='gene_id', columns='strategy', values='call', aggfunc='first').fillna('ABSTAIN')",
            "    cols = sorted(set(tr.columns) & set(te.columns)); tr, te = tr[cols], te[cols]",
            "    ytr = c[c.role == 'inner_train'].groupby('gene_id').truth.first().reindex(tr.index)",
            "    yte = c[c.role == 'outer_test'].groupby('gene_id').truth.first().reindex(te.index)",
            "    Xtr = pd.get_dummies(tr.astype(str)); Xte = pd.get_dummies(te.astype(str)).reindex(columns=Xtr.columns, fill_value=0)",
            "    meta = LogisticRegression(max_iter=3000).fit(Xtr.astype(float), ytr)",
            "    pm = meta.predict(Xte.astype(float))",
            "    best = max(cols, key=lambda m: float((tr[m] == ytr).mean()))",
            "    pb = te[best].to_numpy()",
            "    pv = te.apply(lambda r: r[r != 'ABSTAIN'].value_counts().index[0] if (r != 'ABSTAIN').any() else 'ABSTAIN', axis=1).to_numpy()",
            "    y = yte.to_numpy()",
            "    for name, p in (('meta-model', pm), ('best single', pb), ('plurality vote', pv), ('commonest class', np.full(len(y), ytr.value_counts().index[0]))):",
            "        rows.append(dict(fold=int(k), method=name, chosen=best if name == 'best single' else '', accuracy=float((p == y).mean()), n=len(y)))",
            "    paired['meta_only'] += int(((pm == y) & (pb != y)).sum()); paired['best_only'] += int(((pb == y) & (pm != y)).sum())",
            "res = pd.DataFrame(rows)",
            "summary = res.groupby('method').agg(accuracy=('accuracy', 'mean'), low=('accuracy', 'min'), high=('accuracy', 'max')).round(3).sort_values('accuracy', ascending=False)",
            "p = binomtest(paired['meta_only'], paired['meta_only'] + paired['best_only'], 0.5).pvalue",
            "admitted = bool(summary.loc['meta-model', 'accuracy'] > summary.loc['best single', 'accuracy'] and p < 0.05)",
            "res.to_csv(os.path.join(OUT, 'folds.tsv'), sep='\\t', index=False); summary.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\\t')",
            "print('paired:', paired, '; McNemar p', round(p, 5), '; admitted:', admitted)",
            "summary")
    nb.write(NB)


if __name__ == "__main__":
    if "--strategy" in sys.argv:
        os.makedirs(OUT, exist_ok=True)
        worker(sys.argv[sys.argv.index("--strategy") + 1])
    else:
        main()
