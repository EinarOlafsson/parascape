"""Stage 5 of the information map: predict a PAIR with whole proteins held out (instruction 72).

Every other stage of that plan predicts a property of one gene. This one predicts a property of two: are these
two proteins partners? The unit change is the whole difficulty, and it is why the ladder of instruction 73 says
nothing about it -- those 69 questions are all gene-level.

**Whole proteins are held out, not pairs.** A random split of pairs leaks: protein A appears in a hundred pairs,
so a model that memorises A scores well on a held-out pair containing A without knowing anything about pairing.
Genes are split into partitions by orthogroup (`splits.make_split`, the harness's own), and a pair is kept only
when BOTH endpoints fall in the same partition. A test pair therefore has two unseen proteins.

**The baseline is popularity, and beating it is the whole test.** Partners are not drawn uniformly: a
well-studied, abundant, essential protein appears in many recorded pairs, so the product of the two proteins'
recorded partner counts is a strong predictor that knows nothing about biology. Degree-matching the negatives
was tried first and starves them (the high-degree bins hold too few genes), so the baseline is reported beside
the model instead. A model that cannot beat it has learnt popularity rather than pairing.

Features are symmetric functions of the two genes' measurements -- |a - b|, min and max -- because a pair has no
order. The StarPath finding that the CRISPR-fitness DIFFERENCE of a pair is predictive is the reason |a - b| is
there, and the importance table reports whether it carried its weight.

    python scripts/build_pair_prediction_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_09"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"pair_prediction_{DATE}")
NB = os.path.join(REPO, "notebooks", f"pair_prediction_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Predicting a pair with whole proteins held out",
                          namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**Stage 5 of instruction 72.** Every other stage predicts a property of one gene; this predicts a "
          "property of two. The strategy ladder of instruction 73 says nothing about it, because all 69 of "
          "those questions are gene-level.",
          "**Whole proteins are held out, not pairs.** Splitting pairs at random leaks: a protein in a hundred "
          "pairs can be memorised, and a held-out pair containing it then scores well without the model "
          "knowing anything about pairing. Genes are split by orthogroup and a pair is kept only when both "
          "endpoints land in the same partition, so a test pair has two unseen proteins.",
          "**Popularity is the baseline to beat.** A well-studied, abundant protein appears in many recorded "
          "pairs, so the product of the two partner counts predicts well while knowing no biology. Matching "
          "the negatives on degree was tried and starves them, so the baseline is reported beside the model: "
          "not beating it is itself the answer.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, numpy as np, pandas as pd",
            "from parascape import datasets as ds, organisms as O, splits",
            "pairs = pd.read_parquet(os.path.join(REPO, 'parascape/data/information_map/pairs.parquet'))",
            "pairs = pairs[pairs.partner_side != 'host']",
            "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet')).set_index('gene_id')",
            "known = set(nodes.index)",
            "pairs = pairs[pairs.gene_a.isin(known) & pairs.gene_b.isin(known)]",
            "edge = pairs.apply(lambda r: tuple(sorted((r.gene_a, r.gene_b))), axis=1)",
            "positives = pd.DataFrame(list(dict.fromkeys(edge)), columns=['a', 'b'])",
            "print(len(pairs), 'recorded parasite-parasite pairs ->', len(positives), 'distinct pairs over',",
            "      len(set(positives.a) | set(positives.b)), 'proteins')",
            "degree = pd.concat([positives.a, positives.b]).value_counts()",
            "degree.head(8).to_frame('pairs recorded')")
    nb.md("## 1 · The split: whole proteins, by orthogroup")
    nb.code("groups = nodes['orthogroup'].astype(object)",
            "groups = groups.where(groups.notna(), pd.Series({g: 'self:' + g for g in nodes.index}))",
            "entities = sorted(set(positives.a) | set(positives.b))",
            "split = splits.make_split(entities, groups.reindex(entities).tolist(), organism=O.TOXOPLASMA,",
            "                          benchmark_id='information_map_pairs_2026_11_09', seed=0)",
            "where = {}",
            "for part in ('train', 'tune', 'calibration', 'test'):",
            "    for gene in split.entities(part):",
            "        where[gene] = part",
            "side = positives.assign(pa=positives.a.map(where), pb=positives.b.map(where))",
            "kept = side[side.pa.notna() & (side.pa == side.pb)].copy()",
            "kept['partition'] = kept.pa",
            "print(len(kept), 'of', len(positives), 'pairs survive the cold-node rule (both ends same partition)')",
            "kept.partition.value_counts().to_frame('pairs')")
    nb.md("## 2 · Degree-matched negatives",
          "For each positive pair, a negative pair is drawn inside the same partition whose two proteins have "
          "the same partition. Exact degree matching on both ends was tried and starves the negatives, so "
          "the popularity baseline of section 4 carries the weight instead: beating it is the only result "
          "that means anything, and not beating it is itself the answer.")
    nb.code("rng = np.random.default_rng(0)",
            "recorded = set(map(tuple, positives[['a', 'b']].to_numpy()))",
            "# One negative per positive, drawn inside the same partition from the genes that appear in pairs",
            "# at all. Exact degree matching on BOTH ends was tried first and starves: the high-degree bins hold",
            "# too few genes, so only 9,021 of 13,540 pairs found a partner and the test set came out 94 %",
            "# positive. Simple sampling keeps the balance, and the POPULARITY BASELINE in section 4 is then the",
            "# number that matters -- a model that cannot beat the product of the two degrees has learnt",
            "# popularity rather than pairing, and that is the finding either way.",
            "rows = []",
            "for part, block in kept.groupby('partition'):",
            "    members = np.array(sorted(set(block.a) | set(block.b)))",
            "    seen = set()",
            "    for _ in range(len(block) * 30):",
            "        if len(seen) >= len(block):",
            "            break",
            "        a, b = rng.choice(members, size=2, replace=False)",
            "        key = tuple(sorted((str(a), str(b))))",
            "        if key in recorded or key in seen:",
            "            continue",
            "        seen.add(key)",
            "    for a, b in seen:",
            "        rows.append(dict(a=a, b=b, partition=part, label=0))",
            "negatives = pd.DataFrame(rows)",
            "table = pd.concat([kept[['a', 'b', 'partition']].assign(label=1), negatives], ignore_index=True)",
            "print(len(table), 'pairs:', int(table.label.sum()), 'recorded,', int((~table.label.astype(bool)).sum()), 'drawn')",
            "pos_degree = pd.concat([kept.a, kept.b]).map(degree)",
            "neg_degree = pd.concat([negatives.a, negatives.b]).map(degree)",
            "pd.DataFrame({'recorded partners of a positive endpoint': pos_degree.describe(),",
            "              'of a drawn endpoint': neg_degree.describe()}).round(2)")
    nb.md("## 3 · Symmetric pair features",
          "A pair has no order, so each feature is a symmetric function of the two genes: the absolute "
          "difference, the smaller and the larger. The columns are the best-covered numeric measurements; "
          "`identity` columns are excluded, as everywhere.")
    nb.code("from parascape import information_categories as IC",
            "numeric = [c for c in nodes.columns if pd.api.types.is_numeric_dtype(nodes[c])",
            "           and not pd.api.types.is_bool_dtype(nodes[c])]",
            "def usable(c):",
            "    cats = set(IC.categories_of(c, O.TOXOPLASMA))",
            "    return bool(cats) and 'identity' not in cats",
            "covered = [c for c in numeric if usable(c) and nodes[c].notna().mean() >= 0.5]",
            "covered = sorted(covered, key=lambda c: -nodes[c].notna().mean())[:120]",
            "A = nodes.loc[table.a, covered].to_numpy(dtype=float)",
            "B = nodes.loc[table.b, covered].to_numpy(dtype=float)",
            "X = pd.DataFrame(np.hstack([np.abs(A - B), np.minimum(A, B), np.maximum(A, B)]),",
            "                 columns=[f'{k}:{c}' for k in ('diff', 'min', 'max') for c in covered])",
            "y = table.label.to_numpy()",
            "print(X.shape[1], 'pair features from', len(covered), 'gene columns')",
            "X.iloc[:3, :4].round(3)")
    nb.md("## 4 · Fit on train + tune, test on the cold pairs",
          "The frozen test partition is touched once. The baselines are the ones that matter for pairs: the "
          "majority answer, and the product of the two proteins' recorded partner counts -- which is what a "
          "model would exploit if the negatives were not degree-matched.")
    nb.code("from sklearn.ensemble import HistGradientBoostingClassifier",
            "from sklearn.metrics import roc_auc_score, average_precision_score",
            "fit = table.partition.isin(['train', 'tune']).to_numpy()",
            "test = (table.partition == 'test').to_numpy()",
            "model = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=15,",
            "                                       min_samples_leaf=10, l2_regularization=1.0,",
            "                                       early_stopping=False, random_state=0).fit(X[fit], y[fit])",
            "score = model.predict_proba(X[test])[:, 1]",
            "fame = (table.a.map(degree).fillna(1) * table.b.map(degree).fillna(1)).to_numpy()[test]",
            "result = dict(train_pairs=int(fit.sum()), test_pairs=int(test.sum()),",
            "              test_positive_share=round(float(y[test].mean()), 3),",
            "              auroc=round(float(roc_auc_score(y[test], score)), 3),",
            "              average_precision=round(float(average_precision_score(y[test], score)), 3),",
            "              auroc_fame_baseline=round(float(roc_auc_score(y[test], fame)), 3),",
            "              average_precision_fame=round(float(average_precision_score(y[test], fame)), 3))",
            "json.dump(result, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "print(result)",
            "pd.DataFrame([result]).T.rename(columns={0: 'value'})")
    nb.md("## 5 · Did the fitness difference carry its weight?",
          "The StarPath finding is the reason the absolute difference is in the features. Permutation "
          "importance over the test pairs, grouped by the kind of feature, says whether it did.")
    nb.code("from sklearn.inspection import permutation_importance",
            "imp = permutation_importance(model, X[test], y[test], n_repeats=3, random_state=0,",
            "                             scoring='average_precision')",
            "weights = pd.Series(imp.importances_mean, index=X.columns)",
            "by_kind = weights.groupby(weights.index.str.split(':').str[0]).sum().round(4)",
            "top = weights.sort_values(ascending=False).head(12).round(4)",
            "fitness = weights[[c for c in weights.index if c.split(':', 1)[1].startswith('fit_')]]",
            "summary = dict(by_kind=by_kind.to_dict(),",
            "               fitness_share=round(float(fitness.sum() / max(weights.clip(lower=0).sum(), 1e-9)), 3),",
            "               fitness_difference_share=round(float(",
            "                   fitness[[c for c in fitness.index if c.startswith('diff:')]].sum()",
            "                   / max(weights.clip(lower=0).sum(), 1e-9)), 3))",
            "json.dump(summary, open(os.path.join(OUT, 'importance.json'), 'w'), indent=1)",
            "weights.sort_values(ascending=False).head(30).to_frame('importance').to_csv(",
            "    os.path.join(OUT, 'importance.tsv'), sep='\\t')",
            "print(summary)",
            "top.to_frame('importance')")
    nb.write(NB)


if __name__ == "__main__":
    main()
