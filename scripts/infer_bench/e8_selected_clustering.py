"""E8 of the Infer benchmark: does clustering on only the informative features recover a label better?

Question (Einar, 2026-11-25): with E1 done, choose for one kind of inference the features that
carry information about it, cluster on those, and see whether more is gained than from clustering
on everything.

Selecting features with the label and then grading clusters against the same label is circular,
so the selection never sees the genes it is graded on:

1. **split**: E7's untouched 20% of gene FAMILIES (seed 2026) is the grading set; the other 80%
   (development) is the only place labels are read for selection and for naming clusters;
2. **selection on development genes only**: every feature is scored as E1 scores it (kNN k = 15
   cosine, family 5-fold, real against shuffled, `SEEDS` seeds), on development labels; a feature
   is kept with BH q < 0.05 and a gain above the noise floor. The extra views (sequence models,
   structure), which E1 scores only as whole views, are kept or dropped whole by the same test;
3. **three embeddings**, each prepared as figure 6 (median fill, z-score, at most 30 components
   per view, views side by side):
   - `all views` - everything left after the label exclusions (E5);
   - `selected` - only what step 2 kept;
   - `random` - as many node features and extra views as `selected`, drawn at random (`RANDOM_DRAWS`
     draws): the control for having fewer columns at all;
4. **clusters** (no label used): Louvain on the cosine kNN graph (E5 choice), resolution 0.5 and 1;
5. **grading, on the untouched genes**:
   - recovery: clusters enriched for a class among DEVELOPMENT genes (hypergeometric, BH q < 0.05);
     the share of untouched labelled genes that fall in a cluster enriched for their own class,
     averaged over classes; an untouched gene without enough data to be placed counts as missed;
   - the same with development labels shuffled (`SHUFFLES` times): the chance level;
   - prediction: kNN from development genes to untouched ones (accuracy against chance; for a
     binary label the area under the ROC curve for the rarer class);
   - adjusted mutual information between clusters and the label on the untouched genes.

Writes `selection.parquet`, `results.parquet` and `summary.json`.

    python scripts/infer_bench/e8_selected_clustering.py --out DIR [--labels "a;b"] [--seeds 5]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from sklearn.metrics import adjusted_mutual_info_score, roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e3_preparation as E3  # noqa: E402
import e4_similarity as E4  # noqa: E402
import e5_clustering as E5  # noqa: E402
import e7_inference as E7  # noqa: E402

IM = H.IM
SEEDS = 5
RANDOM_DRAWS = 3
SHUFFLES = 5
RESOLUTIONS = (0.5, 1.0)
DEFAULT_LABELS = ("host-facing (papers)", "effect on host", "hyperLOPIT niche", "compartment (papers)",
                  "enzymatic function (EC class)", "cell-cycle phase",
                  "interaction module (STRING, measured in Toxoplasma)",
                  "interaction module (STRING, copied by orthology)",
                  "complex (CORUM human, copied by orthology)")


def select(nodes, views, extra, y, kind, fam, drop, seeds):
    """Score every node feature and extra view on the given (development-only) labels."""
    rows = []
    for vname, cols in views.items():
        for c in cols:
            if c in drop:
                continue
            X = nodes[[c]].to_numpy(float)
            if (~np.isnan(X[:, 0])).sum() < 100:
                continue
            rows.append({"kind_": "feature", "view": vname, "feature": c,
                         **E1.score_row(E1.one_column(X), y, kind, fam, seeds)})
    for vname, X in extra.items():
        cover = ~np.isnan(X).any(axis=1)
        if cover.sum() < 100:
            continue
        rows.append({"kind_": "extra view", "view": vname, "feature": f"view:{vname}",
                     **E1.score_row(IM.embed(X, cover), y, kind, fam, seeds)})
    t = pd.DataFrame(rows)
    t["q"] = np.nan
    for _k, part in t.groupby("kind_"):
        t.loc[part.index, "q"] = E1.bh(part["p"].to_numpy())
    t["kept"] = (t["q"] < 0.05) & (t["effect"] > t["noise_floor"])
    return t


def embed(nodes, views, extra, features: dict, extras: list) -> tuple:
    """Prepared blocks side by side for the chosen node features (per view) and extra views."""
    blocks = []
    for vname, cols in features.items():
        if not cols:
            continue
        X = nodes[cols].to_numpy(float)
        cover = (~np.isnan(X)).mean(axis=1) >= 0.5
        if cover.sum() >= 100:
            blocks.append(E3.prepare(X, cover))
    for vname in extras:
        X = extra[vname]
        cover = ~np.isnan(X).any(axis=1)
        if cover.sum() >= 100:
            blocks.append(E3.prepare(X, cover))
    if not blocks:
        return None, None
    E = np.hstack(blocks)
    covered = (~np.isnan(E)).mean(axis=1) >= 0.5
    E = np.nan_to_num(E)
    return E[covered], np.flatnonzero(covered)


def recovery(clusters: np.ndarray, y: np.ndarray, dev: np.ndarray, fin: np.ndarray) -> float:
    """Class-averaged share of untouched genes in a cluster enriched (from development labels) for their class."""
    tr = dev & pd.notna(y) & (clusters >= 0)
    N = int(tr.sum())
    counts = pd.Series(y[tr]).value_counts()
    keys, ps = [], []
    for c in np.unique(clusters[tr]):
        in_c = tr & (clusters == c)
        n = int(in_c.sum())
        for cls, K_ in counts.items():
            x = int((y[in_c] == cls).sum())
            if x:
                keys.append((c, cls)); ps.append(hypergeom.sf(x - 1, N, K_, n))
    q = E1.bh(np.array(ps)) if ps else np.array([])
    enriched = {k for k, qq in zip(keys, q) if qq < 0.05}
    per_class: dict = {}
    for i in np.flatnonzero(fin & pd.notna(y)):
        per_class.setdefault(y[i], []).append(clusters[i] >= 0 and (clusters[i], y[i]) in enriched)
    return float(np.mean([np.mean(v) for v in per_class.values()])) if per_class else np.nan


def shuffled_dev(y, dev, seed):
    """Development labels permuted among themselves; untouched labels kept (they are the truth)."""
    yp = y.copy()
    m = dev & pd.notna(y)
    yp[m] = np.random.default_rng(seed).permutation(yp[m])
    return yp


def prediction(E, rows, y, kind, dev, fin, n):
    """kNN from development genes to untouched genes: (accuracy, chance, AUROC of the rarer class, placed share)."""
    full = np.full((n, E.shape[1]), np.nan)
    full[rows] = E
    ok = ~np.isnan(full).any(axis=1)
    tr, te = dev & pd.notna(y) & ok, fin & pd.notna(y) & ok
    counts = pd.Series(y[dev & pd.notna(y)]).value_counts()
    chance = float(counts.iloc[0] / counts.sum())
    placed = float(te.sum() / max(1, (fin & pd.notna(y)).sum()))
    if tr.sum() <= IM.K or not te.any():
        return np.nan, chance, np.nan, placed
    classes = np.unique(y[tr])
    p = E4.knn_predict(full, y, "categorical", tr, te, classes, "cosine")
    acc = float((classes[p.argmax(axis=1)] == y[te]).mean())
    auroc = np.nan
    if kind == "binary":
        rare = counts.index[-1]
        truth = (y[te] == rare).astype(int)
        if 0 < truth.sum() < len(truth) and rare in classes:
            auroc = float(roc_auc_score(truth, p[:, list(classes).index(rare)]))
    return acc, chance, auroc, placed


def grade(E, rows, y, kind, dev, fin, n, name, lname, extra_info, log):
    """Cluster one embedding and grade it on the untouched genes."""
    out = []
    acc, chance, auroc, placed = prediction(E, rows, y, kind, dev, fin, n)
    for res in RESOLUTIONS:
        c0 = E5.louvain(E, res, 0)
        full = np.full(n, -1)
        full[rows] = c0
        rec = recovery(full, y, dev, fin)
        null = np.nanmean([recovery(full, shuffled_dev(y, dev, 200 + s), dev, fin) for s in range(SHUFFLES)])
        lab = fin & pd.notna(y) & (full >= 0)
        ami = adjusted_mutual_info_score(y[lab].astype(str), full[lab]) if lab.sum() > 10 else np.nan
        out.append({"label": lname, "embedding": name, "resolution": res, "columns": int(E.shape[1]),
                    "genes_clustered": int(len(rows)), "clusters": int(len(set(c0))),
                    "untouched_placed": placed, "recovery": rec, "shuffled_recovery": null,
                    "recovery_gain": rec - null, "ami": ami, "knn_accuracy": acc, "chance": chance,
                    "knn_lift": acc - chance if np.isfinite(acc) else np.nan, "auroc_rarer_class": auroc,
                    **extra_info})
        log(f"  {name} res {res}: recovery {rec:.3f} (shuffled {null:.3f}), AMI {ami:.3f}, kNN {acc:.3f} vs {chance:.3f}")
    return out


def main(argv=None) -> int:
    """Select on development genes, cluster three ways, grade on the untouched genes."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--labels", default=";".join(DEFAULT_LABELS))
    ap.add_argument("--seeds", type=int, default=SEEDS)
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    wanted = [x for x in args.labels.split(";") if x]
    missing = set(wanted) - set(L)
    if missing:
        raise SystemExit(f"unknown labels: {sorted(missing)}; known: {list(L)}")
    from parascape import label_catalogue
    label_catalogue.write({k: L[k] for k in wanted}, args.out,
                          {k: IM.exclusions(L[k][2], views) for k in wanted})
    fam = H.families(nodes)
    fin_all = E7.final_split(fam)
    n = len(nodes)
    sel_rows, res_rows = [], []
    for lname in wanted:
        y, kind, tcols = L[lname]
        if kind == "continuous":
            log(f"{lname}: continuous, skipped (clusters need classes)")
            continue
        drop = set(IM.exclusions(tcols, views))
        dev, fin = ~fin_all, fin_all
        y_dev = np.where(dev, y, None).astype(object)
        S = select(nodes, views, extra, y_dev, kind, fam, drop, args.seeds)
        S.insert(0, "label", lname)
        sel_rows.append(S)
        kept = S[S["kept"]]
        feats = {v: g["feature"].tolist() for v, g in kept[kept["kind_"] == "feature"].groupby("view")}
        extras = kept.loc[kept["kind_"] == "extra view", "view"].tolist()
        n_feat = sum(len(v) for v in feats.values())
        log(f"{lname}: kept {n_feat} of {int((S['kind_'] == 'feature').sum())} features and "
            f"{len(extras)} of {int((S['kind_'] == 'extra view').sum())} extra views (development genes only)")
        info = {"selected_features": n_feat, "selected_extra_views": len(extras)}
        allf = {v: [c for c in cols if c not in drop] for v, cols in views.items()}
        allf = {v: c for v, c in allf.items() if len(c) >= 2}
        E, rows = embed(nodes, views, extra, allf, list(extra))
        res_rows += grade(E, rows, y, kind, dev, fin, n, "all views", lname, info, log)
        def save():
            pd.concat(sel_rows).to_parquet(os.path.join(args.out, "selection.parquet"), index=False)
            pd.DataFrame(res_rows).to_parquet(os.path.join(args.out, "results.parquet"), index=False)
        if n_feat + len(extras) == 0:
            log(f"{lname}: nothing selected")
            save()
            continue
        E, rows = embed(nodes, views, extra, feats, extras)
        res_rows += grade(E, rows, y, kind, dev, fin, n, "selected", lname, info, log)
        pool = S.loc[S["kind_"] == "feature", ["view", "feature"]].to_numpy()
        pool_x = S.loc[S["kind_"] == "extra view", "view"].to_numpy()
        for d in range(RANDOM_DRAWS):
            rng = np.random.default_rng(300 + d)
            pick = pool[rng.choice(len(pool), min(n_feat, len(pool)), replace=False)] if n_feat else []
            rf: dict = {}
            for v, c in pick:
                rf.setdefault(v, []).append(c)
            rx = list(rng.choice(pool_x, min(len(extras), len(pool_x)), replace=False)) if len(extras) else []
            E, rows = embed(nodes, views, extra, rf, rx)
            if E is None:
                continue
            res_rows += grade(E, rows, y, kind, dev, fin, n, f"random {d + 1}", lname, info, log)
        save()
    R = pd.DataFrame(res_rows)
    R["embedding_group"] = R["embedding"].str.replace(r"random \d+", "random", regex=True)
    agg = (R.groupby(["label", "embedding_group", "resolution"])
           [["recovery", "shuffled_recovery", "recovery_gain", "ami", "knn_accuracy", "chance",
             "auroc_rarer_class", "untouched_placed", "columns"]].mean().round(3).reset_index())
    json.dump({"labels": wanted, "seeds": args.seeds, "final_seed": E7.FINAL_SEED,
               "table": agg.to_dict("records")},
              open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
