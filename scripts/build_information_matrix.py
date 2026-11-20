"""Information class x label predictive matrix (Einar, 2026-11-20; new figure 6).

Each information class (view) is a block of gene features; each label is something we know for some
genes. A cell of the matrix asks: predicting the label of held-out genes from their nearest
neighbours in this view, how well does it do, above a shuffled-label null?

Rules (spaCR paper session advice, adopted):
  * step 0: two seeds of 2-fold cross-fitting per cell; |difference| is reported as the noise
    floor (minimal detectable effect), and cells are only called when the gain exceeds it;
  * a planted known positive per view (top 15 % of the view's first principal component, a label
    that view must recover) and one random label that every view must score at chance;
  * no leakage: columns from the dataset that defines the label are removed from every view, the
    localisation view is removed when the target is localisation, and the cell-cycle time course
    when the target is cell-cycle phase;
  * kNN on cosine distance in each view (k = 15), the same notion of "neighbour" as the gene graph.

Also writes a label x label association table over shared genes (Einar: does a sparse label track a
well-covered one?), pair-level partner recovery per view, and greedy view combinations per label.

Writes results/information_matrix_<date>/.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, chi2_contingency, kruskal
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import NearestNeighbors

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from parascape import datasets as DS                    # noqa: E402
from parascape import information_categories as IC      # noqa: E402

CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
STRUCT = "/media/carruthers/mnt2/parascape/structure_seqs"
K = 15
N_NULL = 5
MIN_COVER = 0.2
VIEW_CATS = ["transcript", "translation", "chromatin", "protein_abundance", "modification",
             "fitness", "place", "interaction", "sequence", "structure_computed", "conservation",
             "domain", "host_response", "phenotype"]
SECRETED_NICHES = {"dense granules", "rhoptries 1", "rhoptries 2", "micronemes", "PM - integral",
                   "PM - peripheral 1", "PM - peripheral 2"}


def dataset_of(col):
    try:
        d = DS.provenance(col)
        return d.key if d is not None else None
    except Exception:
        return None


def build_views(nodes):
    views = {}
    for cat in VIEW_CATS:
        cols = [c for c in nodes.columns if cat in IC.categories_of(c)
                and pd.api.types.is_numeric_dtype(nodes[c])
                and nodes[c].notna().mean() >= MIN_COVER and nodes[c].nunique() > 2]
        if len(cols) >= 2:
            views[cat] = cols
    return views


def extra_views(nodes):
    """Views that are not node columns: ESM-2 embedding, secondary structure and 3Di composition."""
    out = {}
    ids = nodes["gene_id"].astype(str)
    esm = os.path.join(ROOT, "parascape", "data", "esm_features.parquet")
    if os.path.exists(esm):
        e = pd.read_parquet(esm)
        if "gene_id" in e.columns:
            e = e.set_index("gene_id")
        out["esm"] = e.reindex(ids).to_numpy(float)
    ssp = os.path.join(STRUCT, "secondary_structure.parquet")
    if os.path.exists(ssp):
        s = pd.read_parquet(ssp).drop_duplicates("gene_id").set_index("gene_id")
        def comp(seq):
            n = max(len(seq), 1)
            f = [seq.count(c) / n for c in "HEC"]
            di = ["HH", "HE", "HC", "EH", "EE", "EC", "CH", "CE", "CC"]
            pairs = [seq[i:i + 2] for i in range(len(seq) - 1)]
            m = max(len(pairs), 1)
            return f + [pairs.count(d) / m for d in di]
        ss = np.array([comp(x) for x in s["ss"]])
        out["secondary_structure"] = pd.DataFrame(ss, index=s.index).reindex(ids).to_numpy(float)
    fa = os.path.join(STRUCT, "3di.fasta")
    if os.path.exists(fa):
        seqs, name = {}, None
        for ln in open(fa):
            ln = ln.strip()
            if ln.startswith(">"):
                name = ln[1:].split()[0]
            elif name:
                seqs[name] = seqs.get(name, "") + ln
        sp = pd.read_parquet(ssp) if os.path.exists(ssp) else None
        acc2gene = {}
        if sp is not None:
            acc2gene = dict(zip(sp["uniprot"], sp["gene_id"]))
        alpha = "ACDEFGHIKLMNPQRSTVWY"
        rows = {}
        for k_, v in seqs.items():
            acc = k_.split("-")[1] if k_.startswith("AF-") else k_
            g = acc2gene.get(acc)
            if g:
                n = max(len(v), 1)
                rows[g] = [v.count(c) / n for c in alpha]
        if rows:
            out["3di"] = pd.DataFrame.from_dict(rows, orient="index").reindex(ids).to_numpy(float)
    return out


def embed(X, cover):
    """Median-impute, standardise, PCA to <= 30 dims; rows outside `cover` are NaN."""
    Xc = X[cover].copy()
    med = np.nanmedian(Xc, axis=0)
    idx = np.where(np.isnan(Xc))
    Xc[idx] = np.take(med, idx[1])
    sd = Xc.std(axis=0); sd[sd == 0] = 1
    Xc = (Xc - Xc.mean(axis=0)) / sd
    k = min(30, Xc.shape[1], Xc.shape[0] - 1)
    if Xc.shape[1] > k:
        Xc = PCA(n_components=k, random_state=0).fit_transform(Xc)
    out = np.full((X.shape[0], Xc.shape[1]), np.nan)
    out[cover] = Xc
    return out


def score(y_true, pred, kind):
    if kind == "continuous":
        r = spearmanr(y_true, pred).statistic
        return float(r) if np.isfinite(r) else np.nan
    classes = np.unique(y_true)
    if len(classes) < 2:
        return np.nan
    aucs = []
    for i, c in enumerate(classes):
        yt = (y_true == c).astype(int)
        if 0 < yt.sum() < len(yt):
            aucs.append(roc_auc_score(yt, pred[:, i]))
        if len(classes) == 2:
            break
    return float(np.mean(aucs))


def knn_predict(E, y, kind, train, test, classes):
    nn = NearestNeighbors(n_neighbors=min(K, train.sum() - 1), metric="cosine")
    nn.fit(E[train])
    _, ind = nn.kneighbors(E[test])
    ytr = y[train]
    if kind == "continuous":
        return ytr[ind].mean(axis=1)
    return np.stack([(ytr[ind] == c).mean(axis=1) for c in classes], axis=1)


def cell(E, y, kind, seed):
    ok = ~np.isnan(E).any(axis=1) & pd.notna(y)
    idx = np.flatnonzero(ok)
    if len(idx) < 40:
        return np.nan, len(idx)
    yv = y[idx]
    classes = np.unique(yv) if kind != "continuous" else None
    rng = np.random.default_rng(seed)
    fold = rng.permutation(len(idx)) % 2
    preds = np.empty((len(idx),) if kind == "continuous" else (len(idx), len(classes)))
    for f in (0, 1):
        tr, te = fold != f, fold == f
        m_tr = np.zeros(len(E), bool); m_tr[idx[tr]] = True
        m_te = np.zeros(len(E), bool); m_te[idx[te]] = True
        preds[te] = knn_predict(E, y, kind, m_tr, m_te, classes)
    return score(yv, preds, kind), len(idx)


def labels(nodes, imap):
    ids = nodes["gene_id"].astype(str)
    v = imap.reindex(ids)
    L = {}
    hf = v["host_facing"].astype(str) == "True"
    niche = v["compartment_hyperlopit"].astype(str)
    # Positive-unlabelled: host-facing genes against all other genes. Defining negatives as
    # internally localised hyperLOPIT markers (first version) made the label a copy of the niche
    # (Cramer V 0.99 with hyperLOPIT), so any association with localisation was circular.
    y = np.where(hf.to_numpy(), "host-facing", "other").astype(object)
    L["effect on host"] = (y, "binary", ["pvm_proximity", "hosttx_effectors"])
    ec = nodes["ec_number"].astype(str).str.extract(r"^(\d)")[0]
    L["enzymatic function (EC class)"] = (ec.where(ec.notna()).to_numpy(object), "categorical",
                                          ["ec_number"])
    L["hyperLOPIT niche"] = (niche.where(niche != "nan").to_numpy(object), "categorical", ["__place__"])
    L["phenotype score, fibroblasts"] = (nodes["fit_invitro_hff"].to_numpy(float), "continuous",
                                         ["fit_invitro_hff"])
    # Wang 2020 reports the opposite sign to the other screens (registry note); flipped so that
    # negative = loss costs growth in every phenotype score.
    L["phenotype score, IFN-gamma macrophages"] = (-nodes["fit_ifng"].to_numpy(float), "continuous",
                                                   ["fit_ifng"])
    L["phenotype score, in vivo (brain)"] = (nodes["fit_invivo_brain"].to_numpy(float),
                                             "continuous", ["fit_invivo_brain"])
    L["cell-cycle phase"] = (nodes["cellcycle_phase"].astype(object).where(
        nodes["cellcycle_phase"].notna()).to_numpy(object), "categorical", ["__cellcycle__"])
    L["nutrient dependence (glucose)"] = (nodes["fit_glucose_dependence"].to_numpy(float),
                                          "continuous", ["fit_glucose_dependence"])
    img = nodes["screen_any_phenotype"].astype(object)
    L["imaging phenotype"] = (img.where(img.notna()).astype(str).to_numpy(object), "binary",
                              ["screen_any_phenotype"])
    return L


def exclusions(target_cols, views):
    """Columns to drop from every view for this target: same dataset, or whole classes."""
    drop = set()
    for t in target_cols:
        if t == "__place__":
            drop |= set(views.get("place", []))
            drop |= {c for vv in views.values() for c in vv if c.startswith("lopit")}
        elif t == "__cellcycle__":
            drop |= {c for vv in views.values() for c in vv if c.startswith("cellcycle")}
        else:
            key = dataset_of(t)
            for vv in views.values():
                for c in vv:
                    if c == t or (key and dataset_of(c) == key):
                        drop.add(c)
    return drop


def main():
    out = os.path.join(ROOT, "results", f"information_matrix_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    nodes = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    imap = pd.read_parquet(os.path.join(CACHE, "information_map", "values.parquet"))
    views = build_views(nodes)
    extra = extra_views(nodes)
    print("views:", {k: len(v) for k, v in views.items()}, "extra:", list(extra))
    L = labels(nodes, imap)

    def matrix_for(drop):
        mats = {}
        for name, cols in views.items():
            keep = [c for c in cols if c not in drop]
            if len(keep) < 2:
                continue
            X = nodes[keep].to_numpy(float)
            cover = (~np.isnan(X)).mean(axis=1) >= 0.5
            if cover.sum() < 100:
                continue
            mats[name] = embed(X, cover)
        for name, X in extra.items():
            cover = ~np.isnan(X).any(axis=1)
            if cover.sum() >= 100:
                mats[name] = embed(X, cover)
        return mats

    base = matrix_for(set())
    # planted known positive per view and one random label
    rng = np.random.default_rng(7)
    for name, E in base.items():
        ok = ~np.isnan(E).any(axis=1)
        pc1 = np.full(len(E), np.nan); pc1[ok] = E[ok, 0]
        thr = np.nanquantile(pc1, 0.85)
        y = np.where(np.isnan(pc1), None, np.where(pc1 >= thr, "planted", "rest")).astype(object)
        L[f"PLANTED: {name}"] = (y, "binary", [])
    L["RANDOM (must be chance)"] = (rng.choice(["a", "b"], len(nodes)).astype(object), "binary", [])

    rows = []
    for lname, (y, kind, tcols) in L.items():
        mats = matrix_for(exclusions(tcols, views)) if tcols else base
        for vname, E in mats.items():
            s0, n0 = cell(E, y, kind, 0)
            s1, _ = cell(E, y, kind, 1)
            nulls = []
            for b in range(N_NULL):
                yp = y.copy()
                ok = pd.notna(yp)
                yp[ok] = np.random.default_rng(100 + b).permutation(yp[ok])
                nulls.append(cell(E, yp, kind, 0)[0])
            nm, ns = float(np.nanmean(nulls)), float(np.nanstd(nulls, ddof=1))
            sc = np.nanmean([s0, s1])
            rows.append({"label": lname, "kind": kind, "view": vname, "genes": n0,
                         "score": sc, "seed0": s0, "seed1": s1, "replicate_diff": abs(s0 - s1),
                         "null_mean": nm, "null_sd": ns,
                         "gain": sc - nm, "z": (sc - nm) / ns if ns > 0 else np.nan,
                         "called": (sc - nm) > max(2 * abs(s0 - s1), 2 * ns)})
        print(f"  {lname}: done")
    cells = pd.DataFrame(rows)
    cells.to_csv(os.path.join(out, "cells.tsv"), sep="\t", index=False)

    # greedy view combinations per real label (seed 0)
    comb = []
    for lname, (y, kind, tcols) in L.items():
        if lname.startswith(("PLANTED", "RANDOM")):
            continue
        mats = matrix_for(exclusions(tcols, views)) if tcols else base
        chosen, best = [], -np.inf
        for step in range(4):
            cand = []
            for v in mats:
                if v in chosen:
                    continue
                E = np.hstack([mats[x] for x in chosen + [v]])
                cand.append((cell(E, y, kind, 0)[0], v))
            cand = [c for c in cand if np.isfinite(c[0])]
            if not cand:
                break
            s, v = max(cand)
            if s <= best + 0.005:
                break
            chosen.append(v); best = s
            comb.append({"label": lname, "step": step + 1, "views": "+".join(chosen), "score": s})
    pd.DataFrame(comb).to_csv(os.path.join(out, "combinations.tsv"), sep="\t", index=False)

    # label x label association over shared genes
    assoc = []
    real = {k: v for k, v in L.items() if not k.startswith(("PLANTED", "RANDOM"))}
    for a, (ya, ka, _) in real.items():
        for b, (yb, kb, _) in real.items():
            if a >= b:
                continue
            ok = pd.notna(ya) & pd.notna(yb)
            n = int(ok.sum())
            if n < 30:
                assoc.append({"label_a": a, "label_b": b, "shared_genes": n, "measure": "",
                              "value": np.nan}); continue
            A, B = ya[ok], yb[ok]
            if ka == "continuous" and kb == "continuous":
                val, meas = spearmanr(A.astype(float), B.astype(float)).statistic, "Spearman rho"
            elif ka != "continuous" and kb != "continuous":
                t = pd.crosstab(pd.Series(A), pd.Series(B))
                chi = chi2_contingency(t)[0]
                val = np.sqrt(chi / (n * (min(t.shape) - 1))) if min(t.shape) > 1 else np.nan
                meas = "Cramer V"
            else:
                cat, con = (A, B) if ka != "continuous" else (B, A)
                groups = [con[cat == c].astype(float) for c in np.unique(cat)]
                groups = [g for g in groups if len(g) > 2]
                h = kruskal(*groups).statistic if len(groups) > 1 else np.nan
                val, meas = (h - len(groups) + 1) / (n - len(groups)), "epsilon squared"
            assoc.append({"label_a": a, "label_b": b, "shared_genes": n, "measure": meas,
                          "value": float(val)})
    pd.DataFrame(assoc).to_csv(os.path.join(out, "label_association.tsv"), sep="\t", index=False)

    # binding partners, pair level: cosine similarity in each view, positives vs random pairs
    pairs = pd.read_parquet(os.path.join(CACHE, "information_map", "pairs.parquet"))
    pairs = pairs[(pairs.source_tier == "observed-direct") & (pairs.partner_side != "host")]
    pos_ = {g: i for i, g in enumerate(nodes["gene_id"].astype(str))}
    P = np.array([(pos_[a], pos_[b]) for a, b in zip(pairs.gene_a, pairs.gene_b)
                  if a in pos_ and b in pos_ and a != b])
    prow = []
    rng = np.random.default_rng(3)
    for vname, E in base.items():
        ok = ~np.isnan(E).any(axis=1)
        pp = P[ok[P[:, 0]] & ok[P[:, 1]]]
        if len(pp) < 50:
            continue
        cov = np.flatnonzero(ok)
        neg = np.stack([rng.choice(cov, len(pp)), rng.choice(cov, len(pp))], axis=1)
        unit = E / np.linalg.norm(np.nan_to_num(E), axis=1, keepdims=True)
        s_pos = np.nansum(unit[pp[:, 0]] * unit[pp[:, 1]], axis=1)
        s_neg = np.nansum(unit[neg[:, 0]] * unit[neg[:, 1]], axis=1)
        auc = roc_auc_score(np.r_[np.ones(len(pp)), np.zeros(len(neg))], np.r_[s_pos, s_neg])
        prow.append({"view": vname, "pairs": len(pp), "auc": auc})
    pd.DataFrame(prow).to_csv(os.path.join(out, "partners.tsv"), sep="\t", index=False)
    json.dump({"views": {k: len(v) for k, v in views.items()}, "extra_views": list(extra),
               "labels": list(L), "k": K, "n_null": N_NULL},
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print("wrote", out)


if __name__ == "__main__":
    main()
