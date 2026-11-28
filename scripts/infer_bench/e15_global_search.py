"""Infer E15: the best single clustering strategy for ParaScape, read out on known biology (Einar 2026-11-28).

Unsupervised strategies only (one strategy for every question): per-view weights (0 = view dropped, up to 2),
per-view PCA (off/5/10/30), PCA whitening, distance (cosine / euclidean / correlation), kNN 5-50, and the
clustering (Louvain resolution 0.2-10 with 1 or 5 consensus seeds, k-means, Ward, HDBSCAN). Readout: the mean over
the nine E8/E11 labels of cross-validated class recovery F1 (e11_optimize.recovery_f1: clusters enriched for a
class among training genes call it; F1 on validation genes), with the largest cluster <= 10 % of genes.

Stages (all on development genes, E7 final split; untouched genes graded once at the end):
  --worker i/n : random configurations on fold 1 (written to trials_w<i>.parquet);
  --finalize   : best quarter on folds 2-3; best 12 refined locally (3 rounds x 20 perturbations of weights and
                 resolution); stability (5 column redraws) of the best 5; the choice maximises mean F1 x stability;
                 graded on untouched genes for every label, plus CORUM-complex and STRING-pair excess (E10 scores),
                 against the E8 reference and the E10 default.
  --shuffled   : every label permuted inside development genes before the search (negative control).

    python scripts/infer_bench/e15_global_search.py --out DIR --worker 0/8 --trials 80
    python scripts/infer_bench/e15_global_search.py --out DIR --finalize
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e7_inference as E7  # noqa: E402
import e10_ops_ledger as E10  # noqa: E402
import e11_optimize as E11  # noqa: E402

LABELS = ("host-facing (papers)", "effect on host", "hyperLOPIT niche", "compartment (papers)",
          "enzymatic function (EC class)", "cell-cycle phase", "interaction module (STRING, measured in Toxoplasma)",
          "interaction module (STRING, copied by orthology)", "complex (CORUM human, copied by orthology)",
          "KEGG pathway (copied by orthology)", "hyperLOPIT niche (copied by orthology)")


def setup(shuffled: bool):
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    fam = H.families(nodes)
    fin = E7.final_split(fam)
    dev = ~fin
    labels = {}
    for i, name in enumerate(LABELS):
        y, kind, _t = L[name]
        y = np.asarray(y, dtype=object)
        if shuffled:
            m = dev & pd.notna(y)
            y = y.copy()
            y[m] = np.random.default_rng(4242 + i).permutation(y[m])
        labels[name] = (y, kind)
    # One strategy for every label, so every column that is circular for ANY label is left out (the union of the
    # per-label exclusions E8/E11 apply): the place view, the cell-cycle columns, ... A strategy may not win a label
    # by weighting that label's own measurement.
    drop = set()
    for name in LABELS:
        drop |= set(H.IM.exclusions(L[name][2], views))
    views = {v: [c for c in cols if c not in drop] for v, cols in views.items()}
    views = {v: cols for v, cols in views.items() if cols}
    D = E11.Data(nodes, views, extra, drop)
    inner = H.family_folds(fam, 3, 7)
    folds = [(dev & (inner != f), dev & (inner == f)) for f in range(3)]
    view_names = list(views) + list(extra)
    print(f"{len(drop)} circular columns left out; views: {', '.join(views)}", flush=True)
    return nodes, views, extra, labels, D, dev, fin, folds, view_names


def sample(rng, view_names) -> dict:
    w = {v: float(rng.choice([0.0, 0.5, 1.0, 1.0, 1.5, 2.0])) for v in view_names}
    if all(x == 0 for x in w.values()):
        w[view_names[0]] = 1.0
    return {"weights": w, "pca": int(rng.choice([0, 5, 10, 30])), "whiten": bool(rng.random() < 0.5),
            "metric": str(rng.choice(["cosine", "euclidean", "correlation"])), "knn": int(rng.choice([5, 10, 15, 30, 50])),
            "method": str(rng.choice(["louvain"] * 4 + ["kmeans", "ward", "hdbscan"])),
            "resolution": float(np.exp(rng.uniform(np.log(0.2), np.log(10.0)))), "consensus": int(rng.choice([1, 1, 5])),
            "k": int(np.exp(rng.uniform(np.log(20), np.log(500)))), "min_size": int(rng.choice([5, 10, 20, 40]))}


def perturb(cfg, rng, view_names) -> dict:
    c = json.loads(json.dumps(cfg))
    for v in rng.choice(view_names, 2, replace=False):
        c["weights"][v] = float(np.clip(c["weights"][v] + rng.choice([-0.5, 0.5]), 0, 2))
    c["resolution"] = float(np.clip(c["resolution"] * np.exp(rng.normal(0, 0.3)), 0.1, 15))
    c["k"] = int(np.clip(c["k"] * np.exp(rng.normal(0, 0.3)), 10, 800))
    if rng.random() < 0.2:
        c["knn"] = int(rng.choice([5, 10, 15, 30, 50]))
    return c


def embed(D, nodes, views, extra, cfg, col_rng=None):
    blocks = []
    for v, cols in views.items():
        w = cfg["weights"].get(v, 0.0)
        if w <= 0 or len(cols) < 1:
            continue
        cc = list(col_rng.choice(cols, max(1, int(0.8 * len(cols))), replace=False)) if col_rng is not None else cols
        b = D.block(v, cc, cfg["pca"] if len(cc) > 1 else 0)
        blocks.append(b * w / np.sqrt(max(b.shape[1], 1)))
    for v in extra:
        w = cfg["weights"].get(v, 0.0)
        if w > 0:
            b = D.block(v, [], cfg["pca"], extra=True)
            blocks.append(b * w / np.sqrt(max(b.shape[1], 1)))
    if len(D._blocks) > 400:
        D._blocks.clear()
    E = np.hstack(blocks)
    rows = np.flatnonzero((~np.isnan(E)).mean(axis=1) >= 0.5)
    E = np.nan_to_num(E[rows]).astype(np.float64)
    if cfg["whiten"] and E.shape[1] > 2:
        E = PCA(n_components=0.99, whiten=True, random_state=0).fit_transform(E)
    return E, rows


def partition(D, nodes, views, extra, cfg, seed=0, col_rng=None):
    E, rows = embed(D, nodes, views, extra, cfg, col_rng)
    c = E11.cluster(E, cfg, seed)
    full = np.full(D.n, -1)
    full[rows] = c
    return full


def score(full, labels, train, val):
    f1 = {n: E11.recovery_f1(full, y, kind, train, val)[0] for n, (y, kind) in labels.items()}
    on = full[full >= 0]
    largest = float(np.bincount(on).max() / len(on)) if len(on) else 1.0
    return float(np.nanmean(list(f1.values()))), f1, largest


def evaluate(ctx, cfg, folds_idx, seed=0):
    nodes, views, extra, labels, D, dev, fin, folds, vn = ctx
    full = partition(D, nodes, views, extra, cfg, seed)
    out = []
    for f in folds_idx:
        m, per, largest = score(full, labels, *folds[f])
        out.append((m, per, largest))
    return out, full


def worker(args, ctx):
    i, n = (int(x) for x in args.worker.split("/"))
    rng = np.random.default_rng(1000 + i)
    path = os.path.join(args.out, f"trials_w{i}.parquet")
    rows = pd.read_parquet(path).to_dict("records") if os.path.exists(path) else []
    for t in range(len(rows), args.trials):
        cfg = sample(rng, ctx[-1])
        t0 = time.time()
        try:
            (m, per, largest), = evaluate(ctx, cfg, [0])[0]
        except Exception as e:  # noqa: BLE001
            m, per, largest = np.nan, {}, 1.0
        rows.append({"worker": i, "trial": t, "cfg": json.dumps(cfg), "f1_fold0": m, "largest": largest,
                     "seconds": round(time.time() - t0, 1), **{f"f1:{k}": v for k, v in per.items()}})
        pd.DataFrame(rows).to_parquet(path, index=False)
        print(f"w{i} t{t}: {m:.4f} largest {largest:.2f} ({time.time() - t0:.0f}s)", flush=True)


def finalize(args, ctx):
    nodes, views, extra, labels, D, dev, fin, folds, vn = ctx
    T = pd.concat([pd.read_parquet(f) for f in glob.glob(os.path.join(args.out, "trials_w*.parquet"))])
    T = T[(T["largest"] <= 0.10) & T["f1_fold0"].notna()].sort_values("f1_fold0", ascending=False)
    keep = T.head(max(len(T) // 4, 12))
    S = []
    for cfg_s in keep["cfg"]:
        cfg = json.loads(cfg_s)
        res, _ = evaluate(ctx, cfg, [0, 1, 2])
        S.append({"cfg": cfg_s, "cv": float(np.mean([r[0] for r in res])), "largest": max(r[2] for r in res)})
        print("stage2", round(S[-1]["cv"], 4), flush=True)
    S = pd.DataFrame(S).sort_values("cv", ascending=False)
    S.to_parquet(os.path.join(args.out, "stage2.parquet"), index=False)
    pool = [json.loads(c) for c in S["cfg"].head(12)]
    best = [(float(s), json.loads(c)) for s, c in zip(S["cv"].head(12), S["cfg"].head(12))]
    rng = np.random.default_rng(77)
    for r in range(3):                                      # local refinement around the current best
        cand = []
        for s0, cfg in sorted(best, key=lambda x: -x[0])[:4]:
            for _ in range(5):
                c = perturb(cfg, rng, vn)
                res, _ = evaluate(ctx, c, [0, 1, 2])
                if max(x[2] for x in res) <= 0.10:
                    cand.append((float(np.mean([x[0] for x in res])), c))
        best = sorted(best + cand, key=lambda x: -x[0])[:12]
        print(f"refine round {r + 1}: best cv {best[0][0]:.4f}", flush=True)
    import cluster_quality as Q
    finals = []
    for s0, cfg in best[:5]:
        ref = partition(D, nodes, views, extra, cfg, 0)
        reds = [partition(D, nodes, views, extra, cfg, 11 + b, np.random.default_rng(900 + b)) for b in range(5)]
        per, _ = Q.stability(ref, reds)
        ok = ref >= 0
        sizes = np.bincount(ref[ok], minlength=len(per))
        st = float(np.average(per[:len(sizes)], weights=sizes))
        finals.append({"cv": s0, "stability": st, "cfg": cfg, "score": s0 * st})
        print(f"final candidate cv {s0:.4f} stability {st:.2f}", flush=True)
    choice = max(finals, key=lambda d: d["score"])
    refs = {"chosen": choice["cfg"],
            "E8 reference": {**choice["cfg"], "weights": {v: 1.0 for v in vn}, "pca": 30, "whiten": False,
                             "metric": "cosine", "knn": 15, "method": "louvain", "resolution": 0.5, "consensus": 1},
            "E10 default": {**choice["cfg"], "weights": {v: 1.0 for v in vn}, "pca": 0, "whiten": True,
                            "metric": "cosine", "knn": 15, "method": "louvain", "resolution": 4.0, "consensus": 5}}
    true_labels = setup(False)[3] if args.shuffled else labels
    ref = E10.reference(nodes, print)
    out = {}
    for name, cfg in refs.items():
        full = partition(D, nodes, views, extra, cfg, 0)
        m, per, largest = score(full, true_labels, dev, fin)
        ex = E10.score(full, ref)
        test = [r for r in ex if r["split"] == "test"][0]
        out[name] = {"mean_f1_untouched": m, "per_label": per, "largest": largest,
                     "groups": int(len(np.unique(full[full >= 0]))),
                     "complexes_excess": test["complexes_excess"], "string_f1_excess": test["string_f1_excess"]}
        print(name, round(m, 4), {k: round(v, 3) for k, v in per.items()}, flush=True)
    json.dump({"shuffled": args.shuffled, "trials": int(len(T)), "finalists": finals, "results": out},
              open(os.path.join(args.out, "final.json"), "w"), indent=1, default=str)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--worker", default="")
    ap.add_argument("--trials", type=int, default=80)
    ap.add_argument("--finalize", action="store_true")
    ap.add_argument("--shuffled", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    ctx = setup(args.shuffled)
    if args.finalize:
        finalize(args, ctx)
    else:
        worker(args, ctx)
    return 0


if __name__ == "__main__":
    sys.exit(main())
