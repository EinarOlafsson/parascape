"""E9 of the Infer benchmark: choose, from the data, how to extract information and how to cluster.

Einar, 2026-11-25: "make data based decisions on how to extract information and how to cluster. you
should test several strategies and you should have positive controls."

Every choice is made on DEVELOPMENT genes (80% of gene families) and graded on the untouched 20%
(E7 split, seed 2026), so a strategy never sees the genes it is graded on.

**Stage A - how to extract information (which features to cluster on).** Feature rankings, each
computed on development labels only:

- `knn` - E1 single-feature neighbour score against shuffled labels (read from the E8 selection
  table when given, else computed);
- `mutual information` - sklearn `mutual_info_classif` of each median-filled feature;
- `kruskal` - Kruskal-Wallis test of the feature across classes (BH q);
- `forest` - random-forest impurity importance over all features.

Each ranking gives sets of its top 10, 30 and 100 features, plus its own significance set (`knn`:
q < 0.05 and above the noise floor; `kruskal`: q < 0.05). Extra views (sequence models, structure)
are added whole when their E1 view score passes on development genes, the same for every
strategy, so strategies differ only in node features. Controls: `all views`; `random` sets of 10,
30 and 100 features (3 draws each). Graded with Louvain (resolution 0.5, the E5 choice): recovery
gain over shuffled development labels, and kNN lift over chance (AUROC for binary labels).

**Stage B - how to cluster.** On the stage A winner for the label and on `all views`: Louvain at
resolution 0.5, 1 and 2; HDBSCAN (minimum cluster 10); k-means, Gaussian mixture and Ward at
k = 2 x classes. Graded the same way, plus stability (adjusted Rand index between seeds 0 and 1).

**Positive controls.**

- `planted modules at s`: 20 groups of 15 random fully-measured genes; the 20 planted features (10
  from the transcript view, 10 from the fitness view) get one random direction per group, of length
  s x sqrt(20) column SDs (s = 0.5, 1, 2, 4; the first build used 0.25-1, where no strategy recovered a group), as E0 plants. The label is the group (other genes
  unlabelled). A strategy should find the planted features (recall and precision of the selected
  set are recorded) and recover the groups.
- known complexes: per-class recovery is kept for every class; the ribosome and proteasome classes of
  the hyperLOPIT niche and the CORUM complex labels are the biological controls.

**Negative controls**: shuffled development labels (the chance level of recovery), random feature
sets (the effect of fewer columns alone).

Writes `stage_a.parquet`, `stage_b.parquet`, `per_class.parquet`, `planted_selection.parquet`,
`rankings.parquet`, `label_catalogue.tsv` and `summary.json`; the figure is
`scripts/infer_bench/fig_e9.py`.

    python scripts/infer_bench/e9_strategies.py --out DIR [--e8 DIR] [--labels "a;b"] [--seeds 5]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import kruskal
from sklearn.cluster import HDBSCAN, AgglomerativeClustering, KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e5_clustering as E5  # noqa: E402
import e7_inference as E7  # noqa: E402
import e8_selected_clustering as E8  # noqa: E402

IM = H.IM
TOPK = (10, 30, 100)
RANDOM_DRAWS = 3
PLANT_STRENGTHS = (0.5, 1.0, 2.0, 4.0)
PLANT_GROUPS, PLANT_SIZE, PLANT_COLS = 20, 15, 10
PLANT_VIEWS = ("transcript", "fitness")
PLANT_SEED = 777
CONTROL_CLASSES = ("40S ribosome", "60S ribosome", "19S proteasome", "20S proteasome",
                   "40S ribosomal subunit, cytoplasmic", "60S ribosomal subunit, cytoplasmic",
                   "39S ribosomal subunit, mitochondrial", "PA700 complex", "PA28-20S proteasome")


# ------------------------------------------------------------------ rankings (development only)
def filled(nodes, cols):
    """Columns median-filled (for the rank-only methods)."""
    X = nodes[cols].to_numpy(float)
    med = np.nanmedian(X, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    holes = np.isnan(X)
    X[holes] = np.take(med, np.where(holes)[1])
    return X


def rankings(nodes, views, extra, y_dev, kind, fam, drop, seeds, knn_table=None):
    """Every feature scored by the four methods on development labels; extra views by kNN."""
    feats, seen = [], set()
    for v, cols in views.items():                           # a column in two views is ranked once
        for c in cols:
            if c not in drop and c not in seen and (~np.isnan(nodes[c].to_numpy(float))).sum() >= 100:
                feats.append((v, c)); seen.add(c)
    if knn_table is None:
        knn_table = E8.select(nodes, views, extra, y_dev, kind, fam, drop, seeds)
    k = knn_table.drop_duplicates("feature").set_index("feature")
    cols = [c for _v, c in feats]
    lab = pd.notna(y_dev)
    X = filled(nodes, cols)[lab]
    yl = np.asarray(y_dev[lab]).astype(str)
    mi = mutual_info_classif(X, yl, random_state=0)
    kr_p = []
    for j in range(X.shape[1]):
        groups = [X[yl == c, j] for c in np.unique(yl)]
        try:
            kr_p.append(kruskal(*groups).pvalue)
        except ValueError:                                  # a constant column
            kr_p.append(1.0)
    forest = RandomForestClassifier(n_estimators=300, min_samples_leaf=5, n_jobs=2, random_state=0)
    fi = forest.fit(X, yl).feature_importances_
    t = pd.DataFrame({"view": [v for v, _c in feats], "feature": cols, "mutual information": mi,
                      "kruskal_p": kr_p, "forest": fi})
    t["kruskal_q"] = E1.bh(t["kruskal_p"].to_numpy())
    t["kruskal"] = -np.log10(np.clip(t["kruskal_p"], 1e-300, 1))
    t["knn"] = t["feature"].map(k["effect"])
    t["knn_kept"] = t["feature"].map(k["kept"]).fillna(False).astype(bool)
    xv = knn_table[knn_table["kind_"] == "extra view"]
    extras = xv.loc[xv["kept"], "view"].tolist()
    return t, extras


def strategy_sets(t: pd.DataFrame, rng_seed: int = 300) -> dict:
    """Named feature sets: top-k per method, each method own significance set, and random sets."""
    sets = {}
    for m in ("knn", "mutual information", "kruskal", "forest"):
        order = t.sort_values(m, ascending=False)["feature"].tolist()
        for k in TOPK:
            sets[f"{m} top {k}"] = order[:k]
    sets["knn significant"] = t.loc[t["knn_kept"], "feature"].tolist()
    sets["kruskal significant"] = t.loc[t["kruskal_q"] < 0.05, "feature"].tolist()
    for k in TOPK:
        for d in range(RANDOM_DRAWS):
            rng = np.random.default_rng(rng_seed + 10 * k + d)
            sets[f"random {k} #{d + 1}"] = list(rng.choice(t["feature"].to_numpy(), min(k, len(t)), replace=False))
    return sets


def by_view(t: pd.DataFrame, cols: list) -> dict:
    """{view: [columns]} for a feature list."""
    v = dict(zip(t["feature"], t["view"]))
    out: dict = {}
    for c in cols:
        out.setdefault(v[c], []).append(c)
    return out


# ------------------------------------------------------------------ grading
def recovery_by_class(clusters, y, dev, fin) -> dict:
    """Per class: share of untouched genes in a cluster enriched (development labels) for that class."""
    from scipy.stats import hypergeom
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
    per: dict = {}
    for i in np.flatnonzero(fin & pd.notna(y)):
        per.setdefault(y[i], []).append(clusters[i] >= 0 and (clusters[i], y[i]) in enriched)
    return {c: float(np.mean(v)) for c, v in per.items()}


def run_cluster(E, method, setting, seed):
    """Cluster labels for every row of E."""
    if method == "louvain":
        return E5.louvain(E, setting, seed)
    if method == "hdbscan":
        return HDBSCAN(min_cluster_size=int(setting), copy=True).fit_predict(E)
    if method == "kmeans":
        return KMeans(n_clusters=int(setting), n_init=4, random_state=seed).fit_predict(E)
    if method == "gmm":
        return GaussianMixture(n_components=int(setting), covariance_type="diag", random_state=seed).fit(E).predict(E)
    if method == "ward":
        return AgglomerativeClustering(n_clusters=int(setting), linkage="ward").fit_predict(E)
    raise ValueError(method)


def grade(E, rows, y, kind, dev, fin, n, method="louvain", setting=0.5, stability=False):
    """Recovery, shuffled recovery, per-class recovery and kNN lift for one embedding and clustering."""
    c0 = run_cluster(E, method, setting, 0)
    full = np.full(n, -1)
    full[rows] = c0
    per = recovery_by_class(full, y, dev, fin)
    rec = float(np.mean(list(per.values()))) if per else np.nan
    null = []
    for s in range(E8.SHUFFLES):
        pn = recovery_by_class(full, E8.shuffled_dev(y, dev, 200 + s), dev, fin)
        null.append(np.mean(list(pn.values())) if pn else np.nan)
    acc, chance, auroc, placed = E8.prediction(E, rows, y, kind, dev, fin, n)
    out = {"recovery": rec, "shuffled_recovery": float(np.nanmean(null)), "clusters": int(len(set(c0[c0 >= 0]))),
           "knn_accuracy": acc, "chance": chance, "auroc_rarer_class": auroc, "untouched_placed": placed,
           "columns": int(E.shape[1])}
    out["recovery_gain"] = out["recovery"] - out["shuffled_recovery"]
    out["knn_lift"] = acc - chance if np.isfinite(acc) else np.nan
    if stability:
        out["stability_ari"] = float(adjusted_rand_score(c0, run_cluster(E, method, setting, 1)))
    return out, per


# ------------------------------------------------------------------ positive control: planting
def planted(nodes, views, strength):
    """A copy of the nodes with PLANT_GROUPS modules planted, the label, and the planted columns."""
    rng = np.random.default_rng(PLANT_SEED)
    cols = []
    for v in PLANT_VIEWS:
        pool = [c for c in views[v] if (~np.isnan(nodes[c].to_numpy(float))).mean() > 0.9]
        cols += list(rng.choice(pool, PLANT_COLS, replace=False))
    X = nodes[cols].to_numpy(float)
    full = np.flatnonzero(~np.isnan(X).any(axis=1))
    members = rng.choice(full, PLANT_GROUPS * PLANT_SIZE, replace=False).reshape(PLANT_GROUPS, PLANT_SIZE)
    y = np.full(len(nodes), None, dtype=object)
    Xp = X.copy()
    scale = strength * np.sqrt(len(cols))
    for g, idx in enumerate(members):
        Xp = H.plant(Xp, idx, scale, PLANT_SEED + 1 + g)
        y[idx] = f"module {g + 1:02d}"
    out = nodes.copy()
    out[cols] = Xp
    return out, y, cols


# ------------------------------------------------------------------ main
def main(argv=None) -> int:
    """Stage A and stage B for every label and for the planted controls."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--e8", default="", help="E8 output folder: reuse its kNN selection table")
    ap.add_argument("--labels", default=";".join(E8.DEFAULT_LABELS))
    ap.add_argument("--seeds", type=int, default=E8.SEEDS)
    ap.add_argument("--no-planted", action="store_true")
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
    label_catalogue.write({k: L[k] for k in wanted}, args.out, {k: IM.exclusions(L[k][2], views) for k in wanted})
    e8sel = None
    # Reused only from a FINISHED E8 run (its summary is written last); otherwise computed here.
    if args.e8 and os.path.exists(os.path.join(args.e8, "summary.json")):
        e8sel = pd.read_parquet(os.path.join(args.e8, "selection.parquet"))
        if int(json.load(open(os.path.join(args.e8, "summary.json")))["seeds"]) != args.seeds:
            raise SystemExit("E8 selection was made with other seeds; rerun or pass --seeds to match")
    fam = H.families(nodes)
    fin = E7.final_split(fam)
    dev = ~fin
    n = len(nodes)
    jobs = [(lname, nodes, L[lname][0], L[lname][1], set(IM.exclusions(L[lname][2], views)), None)
            for lname in wanted if L[lname][1] != "continuous"]
    if not args.no_planted:
        for s in PLANT_STRENGTHS:
            pn, py, pcols = planted(nodes, views, s)
            jobs.append((f"planted modules at {s}", pn, py, "categorical", set(), pcols))
    A, B, PC, PS, R = [], [], [], [], []

    def save():
        for name, rows in (("stage_a", A), ("stage_b", B), ("per_class", PC), ("planted_selection", PS)):
            pd.DataFrame(rows).to_parquet(os.path.join(args.out, f"{name}.parquet"), index=False)
        if R:
            pd.concat(R).to_parquet(os.path.join(args.out, "rankings.parquet"), index=False)

    for lname, nd, y, kind, drop, truth in jobs:
        y_dev = np.where(dev, y, None).astype(object)
        knn_t = None
        if e8sel is not None and lname in set(e8sel["label"]):
            knn_t = e8sel[e8sel["label"] == lname].drop(columns=["label"])
        t, extras = rankings(nd, views, extra, y_dev, kind, fam, drop, args.seeds, knn_t)
        R.append(t.assign(label=lname))
        sets = strategy_sets(t)
        log(f"{lname}: {len(t)} features ranked, {len(extras)} extra views pass, "
            f"knn significant {len(sets['knn significant'])}")
        if truth is not None:
            for sname, cols in sets.items():
                hit = len(set(cols) & set(truth))
                PS.append({"label": lname, "strategy": sname, "selected": len(cols), "planted_found": hit,
                           "recall": hit / len(truth), "precision": hit / len(cols) if cols else np.nan})
        allf = {v: [c for c in cols if c not in drop] for v, cols in views.items()}
        allf = {v: c for v, c in allf.items() if len(c) >= 2}
        embeds = {"all views": E8.embed(nd, views, extra, allf, list(extra))}
        for sname, cols in sets.items():
            if cols:
                embeds[sname] = E8.embed(nd, views, extra, by_view(t, cols), extras)
        for sname, (E, rows) in embeds.items():
            if E is None:
                continue
            g, per = grade(E, rows, y, kind, dev, fin, n)
            A.append({"label": lname, "strategy": sname, "family": sname.split(" #")[0], **g})
            PC += [{"label": lname, "stage": "A", "strategy": sname, "method": "louvain 0.5", "class": c,
                    "recovery": r} for c, r in per.items()]
        a = pd.DataFrame([r for r in A if r["label"] == lname])
        real = a[~a["strategy"].str.startswith("random") & (a["strategy"] != "all views")]
        best = real.sort_values(["recovery_gain", "knn_lift"], ascending=False)["strategy"].iloc[0]
        log(f"{lname}: stage A winner {best} (gain {real['recovery_gain'].max():.3f}; "
            f"all views {a.loc[a.strategy == 'all views', 'recovery_gain'].iloc[0]:.3f})")
        k_cls = 2 * len(pd.unique(pd.Series(y[pd.notna(y)])))
        clusterings = ([("louvain", r) for r in (0.5, 1, 2)] + [("hdbscan", 10)]
                       + [(m, k_cls) for m in ("kmeans", "gmm", "ward")])
        for sname in ("all views", best):
            E, rows = embeds[sname]
            for method, setting in clusterings:
                try:
                    g, per = grade(E, rows, y, kind, dev, fin, n, method, setting, stability=True)
                except Exception as error:  # noqa: BLE001 - one failing method must not stop the run
                    log(f"  {method} {setting}: {error}")
                    continue
                B.append({"label": lname, "features": "all views" if sname == "all views" else "stage A winner",
                          "strategy": sname, "method": method, "setting": setting, **g})
                PC += [{"label": lname, "stage": "B", "strategy": sname, "method": f"{method} {setting}",
                        "class": c, "recovery": r} for c, r in per.items()]
        save()
        log(f"{lname}: stage B done")
    save()
    a, b = pd.DataFrame(A), pd.DataFrame(B)
    a["family"] = a["family"].str.replace(r"random (\d+)", r"random \1", regex=True)
    real_labels = [x for x in a["label"].unique() if not x.startswith("planted")]
    ar = a[a["label"].isin(real_labels)].copy()
    ar["rank"] = ar.groupby("label")["recovery_gain"].rank(ascending=False)
    summary = {"labels": wanted, "seeds": args.seeds, "final_seed": E7.FINAL_SEED,
               "stage_a_mean_rank": ar.groupby("family")["rank"].mean().sort_values().round(2).to_dict(),
               "stage_a_mean_gain": ar.groupby("family")["recovery_gain"].mean().sort_values(ascending=False).round(3).to_dict(),
               "stage_b_mean_gain": b[b["label"].isin(real_labels)].groupby(["features", "method", "setting"])
               ["recovery_gain"].mean().round(3).reset_index().to_dict("records"),
               "planted": a[a["label"].str.startswith("planted")].groupby(["label", "family"])["recovery_gain"]
               .mean().round(3).reset_index().to_dict("records")}
    json.dump(summary, open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
