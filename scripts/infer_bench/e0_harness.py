"""E0 of the Infer benchmark: the harness, checked before any experiment uses it.

Design: parascape-paper `infer_benchmark/DESIGN.md` (Einar, 2026-11-24). E0 asks whether the harness
can be trusted, in five parts, each written to its own table:

1. **figure 6 reproduced** (`fig6_reproduced.tsv`): the figure 6 protocol (random 2-fold, seeds 0 and
   1, five shuffled-label nulls, kNN k = 15 on cosine) re-run through this harness and compared cell
   by cell with the published `panel_a/stats.csv`. The harness passes when every cell agrees.
2. **family-held-out folds** (`folds.tsv`): the same cells scored with 5 folds that keep each
   orthogroup whole, against 5 random folds. The drop is how much a paralog that shares both data
   and label lifts a score.
3. **the null** (`null.tsv`): the spread of the shuffled-label baseline over 5 seeds against the
   seed-to-seed spread of the real score, so the "has information" threshold of E1 is set from
   measured noise.
4. **planted complexes** (`planted.tsv`): the 5 complexes of 8 genes per view and strength from
   `maxproj_inputs/planted.csv`, planted by the manifest recipe; recovery is the share of each
   member's 15 nearest neighbours (in its own view) that are in its complex, against chance.
5. **leak probe** (`leaks.tsv`): every single column scored alone against every label (random 2-fold,
   seed 0). A column or a view of two columns or fewer that scores near-perfectly (AUC >= 0.95 or
   |rho| >= 0.9) is flagged for review before E1.

The scoring code is `scripts/build_information_matrix.py` itself (embed, knn_predict, score, labels,
exclusions), imported, not copied, so figure 6 and the benchmark cannot drift apart.

    python scripts/infer_bench/e0_harness.py [--parts fig6,folds,null,planted,leaks] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import date

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, ROOT)
import build_information_matrix as IM  # noqa: E402

INPUTS = "/media/carruthers/mnt2/parascape/maxproj_inputs"
FIG6 = "/media/carruthers/mnt3/claude/parascape-paper/main_figures/figure_6/panel_a/stats.csv"
PARTS = ("fig6", "folds", "null", "planted", "leaks")
N_FOLDS = 5
LEAK_AUC, LEAK_RHO = 0.95, 0.9


# ------------------------------------------------------------------------------ folds
def families(nodes: pd.DataFrame) -> np.ndarray:
    """Each gene's orthogroup; a gene without one is its own family."""
    og = nodes["orthogroup"].astype(str) if "orthogroup" in nodes else pd.Series("", index=nodes.index)
    bad = og.isin(["", "nan", "None"])
    return np.where(bad, "self:" + nodes["gene_id"].astype(str), og).astype(object)


def family_folds(fam: np.ndarray, n: int, seed: int) -> np.ndarray:
    """A fold per gene such that every family is in exactly one fold (families shuffled by seed)."""
    uniq = pd.unique(fam)
    order = np.random.default_rng(seed).permutation(len(uniq))
    fold_of = {f: i % n for i, f in zip(range(len(uniq)), uniq[order])}
    return np.array([fold_of[f] for f in fam])


def random_folds(n_genes: int, n: int, seed: int) -> np.ndarray:
    """A fold per gene, at random."""
    return np.random.default_rng(seed).permutation(n_genes) % n


def cv_score(E: np.ndarray, y: np.ndarray, kind: str, folds: np.ndarray) -> tuple:
    """Held-out score of the label from neighbours in E, with the given fold of every gene."""
    ok = ~np.isnan(E).any(axis=1) & pd.notna(y)
    idx = np.flatnonzero(ok)
    if len(idx) < 40:
        return np.nan, len(idx)
    yv = y[idx]
    classes = np.unique(yv) if kind != "continuous" else None
    preds = np.full((len(idx),) if kind == "continuous" else (len(idx), len(classes)), np.nan)
    f = folds[idx]
    for k in np.unique(f):
        tr, te = f != k, f == k
        if tr.sum() < IM.K + 1 or te.sum() == 0:
            continue
        m_tr = np.zeros(len(E), bool); m_tr[idx[tr]] = True
        m_te = np.zeros(len(E), bool); m_te[idx[te]] = True
        preds[te] = IM.knn_predict(E, y, kind, m_tr, m_te, classes)
    done = ~np.isnan(preds if preds.ndim == 1 else preds[:, 0])
    return IM.score(yv[done], preds[done], kind), len(idx)


def shuffled(y: np.ndarray, seed: int) -> np.ndarray:
    """The label permuted among the genes that carry one, so every class keeps its size."""
    yp = y.copy()
    ok = pd.notna(yp)
    yp[ok] = np.random.default_rng(seed).permutation(yp[ok])
    return yp


# ------------------------------------------------------------------------------ inputs
def load():
    """Nodes, views, extra views, labels and the per-label view matrices, as figure 6 builds them."""
    nodes = pd.read_parquet(os.path.join(IM.CACHE, "nodes.parquet"))
    imap = pd.read_parquet(os.path.join(IM.CACHE, "information_map", "values.parquet"))
    views = IM.build_views(nodes)
    extra = IM.extra_views(nodes)
    L = IM.labels(nodes, imap)
    return nodes, views, extra, L


def matrices(nodes, views, extra, drop) -> dict:
    """The embedded view matrices, without the columns in `drop` (figure 6 matrix_for)."""
    mats = {}
    for name, cols in views.items():
        keep = [c for c in cols if c not in drop]
        if len(keep) < 2:
            continue
        X = nodes[keep].to_numpy(float)
        cover = (~np.isnan(X)).mean(axis=1) >= 0.5
        if cover.sum() < 100:
            continue
        mats[name] = IM.embed(X, cover)
    for name, X in extra.items():
        cover = ~np.isnan(X).any(axis=1)
        if cover.sum() >= 100:
            mats[name] = IM.embed(X, cover)
    return mats


# ------------------------------------------------------------------------------ parts
def part_fig6(nodes, views, extra, L, log) -> pd.DataFrame:
    """Figure 6 panel A, re-run cell by cell and set against the published numbers."""
    rows = []
    for lname, (y, kind, tcols) in L.items():
        mats = matrices(nodes, views, extra, IM.exclusions(tcols, views))
        for vname, E in mats.items():
            s0, n0 = IM.cell(E, y, kind, 0)
            s1, _ = IM.cell(E, y, kind, 1)
            nulls = [IM.cell(E, shuffled(y, 100 + b), kind, 0)[0] for b in range(IM.N_NULL)]
            rows.append({"label": lname, "view": vname, "genes": n0, "score": np.nanmean([s0, s1]),
                         "null_mean": float(np.nanmean(nulls))})
        log(f"fig6: {lname}")
    got = pd.DataFrame(rows)
    pub = pd.read_csv(FIG6)
    m = got.merge(pub[["label", "view", "genes", "score", "null_mean"]], on=["label", "view"],
                  how="outer", suffixes=("", "_published"), indicator=True)
    m["score_diff"] = (m["score"] - m["score_published"]).abs()
    m["null_diff"] = (m["null_mean"] - m["null_mean_published"]).abs()
    return m


def part_folds(nodes, views, extra, L, log) -> pd.DataFrame:
    """Random 5-fold against family-held-out 5-fold, seed 0, every label x view."""
    fam = families(nodes)
    rows = []
    for lname, (y, kind, tcols) in L.items():
        mats = matrices(nodes, views, extra, IM.exclusions(tcols, views))
        for vname, E in mats.items():
            r, n = cv_score(E, y, kind, random_folds(len(E), N_FOLDS, 0))
            f, _ = cv_score(E, y, kind, family_folds(fam, N_FOLDS, 0))
            rows.append({"label": lname, "view": vname, "genes": n, "random_folds": r,
                         "family_folds": f, "drop": r - f})
        log(f"folds: {lname}")
    return pd.DataFrame(rows)


def part_null(nodes, views, extra, L, log, seeds: int = 5) -> pd.DataFrame:
    """Seed-to-seed spread of the real score and of the shuffled baseline, family folds."""
    fam = families(nodes)
    rows = []
    for lname, (y, kind, tcols) in L.items():
        mats = matrices(nodes, views, extra, IM.exclusions(tcols, views))
        for vname, E in mats.items():
            real = [cv_score(E, y, kind, family_folds(fam, N_FOLDS, s))[0] for s in range(seeds)]
            null = [cv_score(E, shuffled(y, 100 + s), kind, family_folds(fam, N_FOLDS, s))[0]
                    for s in range(seeds)]
            rows.append({"label": lname, "view": vname, "real_mean": np.nanmean(real),
                         "real_sd": np.nanstd(real, ddof=1), "null_mean": np.nanmean(null),
                         "null_sd": np.nanstd(null, ddof=1),
                         "gain": np.nanmean(real) - np.nanmean(null)})
        log(f"null: {lname}")
    return pd.DataFrame(rows)


def plant(X: np.ndarray, members: np.ndarray, strength: float, seed: int) -> np.ndarray:
    """The manifest recipe: one unit direction per complex, strength x column SD, members only."""
    out = X.copy()
    d = np.random.default_rng(seed).normal(size=X.shape[1])
    d /= np.linalg.norm(d)
    sd = np.nanstd(X, axis=0)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    out[members] = X[members] + strength * sd * d
    return out


def part_planted(log, planted_file: str = "planted.csv") -> pd.DataFrame:
    """Recovery of each planted complex in its own view, embedded as figure 6 embeds a view.

    `planted_v2.csv` holds the same complexes with strength x sqrt(columns of the view), so every
    view is planted at the same detectability (E0 finding, 2026-11-25)."""
    genes = pd.read_parquet(os.path.join(INPUTS, "genes.parquet"))
    vdefs = {v["name"]: v["columns"] for v in json.load(open(os.path.join(INPUTS, "views.json")))}
    manifest = json.load(open(os.path.join(INPUTS, "manifest.json")))
    seed0 = int(manifest["planted"]["seed"])
    planted = pd.read_csv(os.path.join(INPUTS, planted_file))
    order = {c: i for i, c in enumerate(pd.unique(planted["complex"]))}
    pos = {g: i for i, g in enumerate(genes["gene_id"].astype(str))}
    rows = []
    nominal = "strength_nominal" if "strength_nominal" in planted else "strength"
    for (view, _s), part in planted.groupby(["view", nominal], sort=False):
        strength = float(part["strength"].iloc[0])
        X = genes[vdefs[view]].to_numpy(float)
        for cname, members in part.groupby("complex", sort=False):
            idx = np.array([pos[g] for g in members["gene"]])
            Xp = plant(X, idx, float(strength), seed0 + order[cname])
            cover = (~np.isnan(Xp)).mean(axis=1) >= 0.5 if Xp.shape[1] > 1 else ~np.isnan(Xp[:, 0])
            E = IM.embed(Xp, cover)
            ok = ~np.isnan(E).any(axis=1)
            cov = np.flatnonzero(ok)
            m = idx[ok[idx]]
            if len(m) < 2:
                continue
            from sklearn.neighbors import NearestNeighbors
            nn = NearestNeighbors(n_neighbors=IM.K + 1, metric="cosine").fit(E[cov])
            _, ind = nn.kneighbors(E[m])
            neigh = cov[ind[:, 1:]]
            share = np.isin(neigh, m).mean()
            chance = (len(m) - 1) / (len(cov) - 1)
            rows.append({"view": view, "strength": _s, "applied_strength": strength, "complex": cname,
                         "members": len(m),
                         "neighbour_share": share, "chance": chance, "fold_over_chance": share / chance})
        log(f"planted: {view} x{_s}")
    return pd.DataFrame(rows)


def part_leaks(nodes, views, extra, L, log) -> pd.DataFrame:
    """Every single column against every label; near-perfect single columns are flagged."""
    rows = []
    for lname, (y, kind, tcols) in L.items():
        drop = IM.exclusions(tcols, views)
        for vname, cols in views.items():
            for c in cols:
                if c in drop:
                    continue
                X = nodes[[c]].to_numpy(float)
                cover = ~np.isnan(X[:, 0])
                if cover.sum() < 100:
                    continue
                E = IM.embed(X, cover)
                E = np.c_[E, np.ones(len(E))]           # cosine on one column needs a second axis
                E[~cover] = np.nan
                s, n = cv_score(E, y, kind, random_folds(len(E), 2, 0))
                rows.append({"label": lname, "view": vname, "column": c, "genes": n, "score": s})
        log(f"leaks: {lname}")
    t = pd.DataFrame(rows)
    kinds = {k: v[1] for k, v in L.items()}
    t["kind"] = t["label"].map(kinds)
    t["flag"] = np.where(t["kind"] == "continuous", t["score"].abs() >= LEAK_RHO, t["score"] >= LEAK_AUC)
    small = [v for v, cols in views.items() if len(cols) <= 2]
    t["small_view"] = t["view"].isin(small)
    return t.sort_values("score", ascending=False)


def main(argv=None) -> int:
    """Run the chosen parts and write their tables and a summary."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--parts", default=",".join(PARTS))
    ap.add_argument("--planted", default="planted.csv", help="planted.csv or planted_v2.csv")
    ap.add_argument("--labels", default="", help="only these labels, separated by ';' (a quick check)")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", f"infer_bench_e0_{date.today():%Y_%m_%d}"))
    args = ap.parse_args(argv)
    parts = [p for p in args.parts.split(",") if p]
    unknown = set(parts) - set(PARTS)
    if unknown:
        raise SystemExit(f"unknown parts: {sorted(unknown)}; choose from {PARTS}")
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    summary = {"parts": parts, "cache": IM.CACHE, "k": IM.K, "n_folds": N_FOLDS}
    if set(parts) - {"planted"}:
        nodes, views, extra, L = load()
        wanted = [x for x in args.labels.split(";") if x]
        missing = set(wanted) - set(L)
        if missing:
            raise SystemExit(f"unknown labels: {sorted(missing)}; known: {list(L)}")
        L = {k: v for k, v in L.items() if not wanted or k in wanted}
    if "fig6" in parts:
        t = part_fig6(nodes, views, extra, L, log)
        t.to_csv(os.path.join(args.out, "fig6_reproduced.tsv"), sep="\t", index=False)
        both = t[t["_merge"] == "both"]
        summary["fig6"] = {"cells": int(len(both)), "only_here": int((t["_merge"] == "left_only").sum()),
                           "only_published": int((t["_merge"] == "right_only").sum()),
                           "max_score_diff": float(both["score_diff"].max()),
                           "cells_off_by_more_than_0.001": int((both["score_diff"] > 1e-3).sum())}
    if "folds" in parts:
        t = part_folds(nodes, views, extra, L, log)
        t.to_csv(os.path.join(args.out, "folds.tsv"), sep="\t", index=False)
        summary["folds"] = {"cells": int(len(t)), "median_drop": float(t["drop"].median()),
                            "largest_drops": t.nlargest(5, "drop")[["label", "view", "drop"]].to_dict("records")}
    if "null" in parts:
        t = part_null(nodes, views, extra, L, log)
        t.to_csv(os.path.join(args.out, "null.tsv"), sep="\t", index=False)
        summary["null"] = {"median_real_sd": float(t["real_sd"].median()),
                           "median_null_sd": float(t["null_sd"].median())}
    if "planted" in parts:
        t = part_planted(log, args.planted)
        t.to_csv(os.path.join(args.out, "planted.tsv"), sep="\t", index=False)
        summary["planted"] = (t.groupby("strength")["fold_over_chance"].median().round(1).to_dict())
    if "leaks" in parts:
        t = part_leaks(nodes, views, extra, L, log)
        t.to_csv(os.path.join(args.out, "leaks.tsv"), sep="\t", index=False)
        summary["leaks"] = {"columns_scored": int(len(t)), "flagged": int(t["flag"].sum()),
                            "flagged_rows": t[t["flag"]].head(20)[["label", "column", "score"]].to_dict("records")}
    json.dump(summary, open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
