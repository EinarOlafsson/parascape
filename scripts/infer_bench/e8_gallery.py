"""E8 opened up: for every label x embedding cell of figure E8 panel A, the map, the clusters and their quality.

For each E8 label, one page with three rows: all views (reference), the E1-selected features, and random draw 1
of the same size, each rebuilt exactly as `e8_selected_clustering.py` builds it. Each row has:
  1. UMAP of the embedding coloured by the label's classes (5 largest classes in colour, the rest dark grey,
     unlabelled genes light grey);
  2-3. the same UMAP coloured by Louvain cluster at resolution 0.5 and 1.0: a cluster enriched for a class
     (hypergeometric, Benjamini-Hochberg q < 0.05 over clusters x classes) takes that class's colour, the other
     clusters alternate light greys;
  4-5. per cluster, at resolution 0.5 and 1.0: share of the cluster's labelled genes in its top class (%),
     stability (mean best Jaccard over B redraws of the columns within views, new seed each) and tightness
     (mean cosine silhouette). Clusters are ordered by the top-class share.
Binary labels count only the rarer (positive) class.

Writes page_<label>.png/.svg, clusters.csv (every label x embedding x resolution x cluster), genes.parquet
(per-gene cluster, silhouette and confidence), and one legend .md; skips labels whose page exists (resume).

    python scripts/infer_bench/e8_gallery.py E8_DIR OUT_DIR [--redraws 10]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
import cluster_quality as Q  # noqa: E402
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e5_clustering as E5  # noqa: E402
import e8_selected_clustering as E8  # noqa: E402
from parascape import figstyle as F  # noqa: E402

IM = H.IM
RES = (0.5, 1.0)
MM = 1 / 25.4
GREYS = ("#D7D7D7", "#BEBEBE")
OTHER = F.DARK_GREY
CLASS_COLOURS = (F.ACCENT, F.TEAL, F.NAVY, F.MAGENTA, F.AMBER)


def safe(name: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in name).strip("_")


def embeddings(nodes, views, extra, S, drop):
    """{name: (features by view, extra views)} for all views, selected, random 1, as E8 builds them."""
    allf = {v: [c for c in cols if c not in drop] for v, cols in views.items()}
    allf = {v: c for v, c in allf.items() if len(c) >= 2}
    kept = S[S["kept"]]
    feats = {v: g["feature"].tolist() for v, g in kept[kept["kind_"] == "feature"].groupby("view")}
    extras = kept.loc[kept["kind_"] == "extra view", "view"].tolist()
    n_feat = sum(len(v) for v in feats.values())
    out = {"all views (ref.)": (allf, list(extra))}
    if n_feat + len(extras):
        out["selected"] = (feats, extras)
        pool = S.loc[S["kind_"] == "feature", ["view", "feature"]].to_numpy()
        pool_x = S.loc[S["kind_"] == "extra view", "view"].to_numpy()
        rng = np.random.default_rng(300)                       # E8 random draw 1
        pick = pool[rng.choice(len(pool), min(n_feat, len(pool)), replace=False)] if n_feat else []
        rf: dict = {}
        for v, c in pick:
            rf.setdefault(v, []).append(c)
        rx = list(rng.choice(pool_x, min(len(extras), len(pool_x)), replace=False)) if len(extras) else []
        out["random, same size"] = (rf, rx)
    return out


def embed_full(nodes, views, extra, feats, extras, n, rng=None):
    """The E8 embedding (rows outside coverage dropped); with rng, columns redrawn within every block."""
    ex = extra
    if rng is not None:
        feats = {v: list(rng.choice(c, len(c))) for v, c in feats.items()}
        ex = {k: X[:, rng.integers(0, X.shape[1], X.shape[1])] for k, X in extra.items()}
    return E8.embed(nodes, views, ex, feats, extras)


def enrichment(c: np.ndarray, y: np.ndarray, kind: str) -> pd.DataFrame:
    """Per cluster: the most enriched class, its share of the cluster's labelled genes, p and q."""
    lab = pd.notna(y) & (c >= 0)
    counts = pd.Series(y[lab]).value_counts()
    classes = [counts.index[-1]] if kind == "binary" else list(counts.index)
    N = int(lab.sum())
    rows = []
    for k in np.unique(c[c >= 0]):
        m = lab & (c == k)
        nl = int(m.sum())
        best = None
        for cls in classes:
            x = int((y[m] == cls).sum())
            p = hypergeom.sf(x - 1, N, int(counts[cls]), nl) if x else 1.0
            if best is None or p < best[2]:
                best = (cls, x, p)
        rows.append({"cluster": int(k), "genes": int((c == k).sum()), "labelled": nl, "top_class": best[0],
                     "top_count": best[1], "top_share": best[1] / nl if nl else np.nan, "p": best[2]})
    t = pd.DataFrame(rows)
    t["q"] = E1.bh(t["p"].to_numpy() * (len(classes)))  # Bonferroni over classes inside BH over clusters
    t["q"] = t["q"].clip(upper=1.0)
    return t


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("e8")
    ap.add_argument("out")
    ap.add_argument("--redraws", type=int, default=10)
    ap.add_argument("--labels", default="")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    import umap
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    sel = pd.read_parquet(os.path.join(args.e8, "selection.parquet"))
    labels = [x for x in pd.unique(sel["label"]) if not args.labels or x in args.labels.split(";")]
    n = len(nodes)
    ctab = os.path.join(args.out, "clusters.csv")
    gtab = os.path.join(args.out, "genes.parquet")
    C_all = pd.read_csv(ctab).to_dict("records") if os.path.exists(ctab) else []
    G_all = [pd.read_parquet(gtab)] if os.path.exists(gtab) else []
    F.use()
    for lname in labels:
        page = os.path.join(args.out, f"page_{safe(lname)}.png")
        if os.path.exists(page):
            log(f"{lname}: page exists, skipped")
            continue
        y, kind, tcols = L[lname]
        y = np.asarray(y, dtype=object)
        drop = set(IM.exclusions(tcols, views))
        embs = embeddings(nodes, views, extra, sel[sel["label"] == lname], drop)
        counts = pd.Series(y[pd.notna(y)]).value_counts()
        if kind == "binary":
            shown = [counts.index[-1]] + ([counts.index[0]] if counts.iloc[0] < 3000 else [])
            colour = {shown[0]: F.ACCENT, **({shown[1]: F.NAVY} if len(shown) > 1 else {})}
        else:
            shown = list(counts.index[:5])
            colour = dict(zip(shown, CLASS_COLOURS))
        rows_out = []
        for ename, (feats, extras) in embs.items():
            E, rows = embed_full(nodes, views, extra, feats, extras, n)
            U = umap.UMAP(n_neighbors=15, min_dist=0.3, metric="cosine", random_state=0).fit_transform(E)
            redraw_E = [embed_full(nodes, views, extra, feats, extras, n, np.random.default_rng(900 + b))
                        for b in range(args.redraws)]
            per_res = {}
            for res in RES:
                c = np.full(n, -1)
                c[rows] = E5.louvain(E, res, 0)
                red = []
                for b, (Eb, rb) in enumerate(redraw_E):
                    lab = np.full(n, -1)
                    lab[rb] = E5.louvain(Eb, res, 1 + b)
                    red.append(lab)
                stab, conf = Q.stability(c, red)
                sil_c, sil_g = Q.tightness(E, c[rows])
                en = enrichment(c, y, kind)
                en["stability"] = stab[en["cluster"].to_numpy()]
                en["silhouette"] = sil_c[en["cluster"].to_numpy()]
                en["mean_gene_confidence"] = [float(np.mean(conf[c == k])) for k in en["cluster"]]
                en = en.sort_values("top_share", ascending=False).reset_index(drop=True)
                en["rank"] = np.arange(1, len(en) + 1)
                for r in en.to_dict("records"):
                    C_all.append({"label": lname, "embedding": ename, "resolution": res, **r})
                g = pd.DataFrame({"label": lname, "embedding": ename, "resolution": res,
                                  "gene_id": nodes["gene_id"].to_numpy()[rows], "cluster": c[rows],
                                  "silhouette": sil_g, "confidence": conf[rows]})
                G_all.append(g)
                per_res[res] = (c[rows], en)
                log(f"{lname} | {ename} | res {res}: {len(en)} clusters, stability median "
                    f"{np.nanmedian(en['stability']):.2f}, silhouette median {np.nanmedian(en['silhouette']):.2f}")
            rows_out.append((ename, U, y[rows], per_res))
        draw_page(lname, kind, rows_out, shown, colour, page)
        pd.DataFrame(C_all).to_csv(ctab, index=False)
        pd.concat(G_all).to_parquet(gtab, index=False)
        log(f"{lname}: page written")
    write_legend(args.out, args.redraws)
    log(f"wrote {args.out}")
    return 0


def draw_page(lname, kind, rows_out, shown, colour, path):
    cell = F.SLOT_IN * 25.4
    sq, gap, left, top = 46.0, 9.0, 16.0, 12.0
    kmax = {res: max(len(r[3][res][1]) for r in rows_out) for res in RES}
    bw = {res: cell * kmax[res] for res in RES}
    mini, mgap = 13.0, 4.0
    rowh = max(sq, 3 * mini + 2 * mgap) + 22
    W = left + 3 * (sq + gap) + 34 + bw[0.5] + 18 + bw[1.0] + 6
    Hh = top + 6 + len(rows_out) * rowh + 4
    fig = plt.figure(figsize=(W * MM, Hh * MM))

    def ax_at(x, yt, w, h):
        return fig.add_axes([x / W, 1 - (yt + h) / Hh, w / W, h / Hh])

    for ri, (ename, U, yr, per_res) in enumerate(rows_out):
        yt = top + 6 + ri * rowh
        if ri == 0:
            for j, head in enumerate(("label classes", "clusters, resolution 0.5", "clusters, resolution 1")):
                fig.text((left + j * (sq + gap) + sq / 2) / W, 1 - (yt - 7) / Hh, head, ha="center", va="bottom",
                         fontsize=F.FONT_PT)
        fig.text(2 / W, 1 - (yt - 1) / Hh, ename, fontsize=F.FONT_PT, weight="bold", va="bottom")
        ax = ax_at(left, yt, sq, sq)
        ax.scatter(U[:, 0], U[:, 1], s=0.8, c=GREYS[0], alpha=F.ALPHA, lw=0, rasterized=True)
        lab = pd.notna(yr)
        rest = lab & ~np.isin(yr, shown)
        ax.scatter(U[rest, 0], U[rest, 1], s=2, c=OTHER, alpha=F.ALPHA, lw=0, rasterized=True)
        for cls in reversed(shown):
            m = yr == cls
            ax.scatter(U[m, 0], U[m, 1], s=3, c=colour[cls], alpha=F.ALPHA, lw=0, rasterized=True)
        for j, res in enumerate(RES):
            c, en = per_res[res]
            ax = ax_at(left + (j + 1) * (sq + gap), yt, sq, sq)
            col_of = {}
            for i, r in en.iterrows():
                col_of[r["cluster"]] = (colour.get(r["top_class"], OTHER) if r["q"] < 0.05 else GREYS[i % 2])
            cols = np.array([col_of.get(k, GREYS[0]) for k in c])
            ax.scatter(U[:, 0], U[:, 1], s=0.8, c=cols, alpha=F.ALPHA, lw=0, rasterized=True)
        for ax in fig.axes[-3:]:
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_xlabel("UMAP 1", labelpad=2); ax.set_ylabel("UMAP 2", labelpad=2)
        x0 = left + 3 * (sq + gap) + 34
        for res in RES:
            c, en = per_res[res]
            xs = np.arange(len(en))
            bar_cols = [colour.get(r["top_class"], OTHER) if r["q"] < 0.05 else F.CONTROL for _, r in en.iterrows()]
            for mi, (colname, ylab, ylim, ref) in enumerate((
                    ("top_share", "Top class (%)", (0, 105), None),
                    ("stability", "Stability (Jaccard)", (0, 1.05), 0.75),
                    ("silhouette", "Silhouette", (-0.3, 0.6), 0.0))):
                ax = ax_at(x0, yt + mi * (mini + mgap), bw[res], mini)
                v = en[colname].to_numpy(float) * (100 if colname == "top_share" else 1)
                ax.bar(xs, np.nan_to_num(v), width=F.BAR_FRACTION,
                       color=bar_cols if mi == 0 else F.DARK_GREY, alpha=F.ALPHA)
                if ref is not None:
                    ax.axhline(ref, color=F.INK, lw=0.6, ls=":", alpha=0.6)
                ax.set_xlim(-0.6, kmax[res] - 0.4)
                ax.set_ylim(*ylim)
                ax.set_xticks(xs)
                ax.set_xticklabels([str(i + 1) for i in xs] if mi == 2 else [])
                if res == RES[0]:
                    ax.set_ylabel(ylab, labelpad=2, rotation=0, ha="right", va="center")
                if mi == 2:
                    ax.set_xlabel(f"Cluster (resolution {res:g})", labelpad=2)
            x0 += bw[res] + 18
    from matplotlib.lines import Line2D
    mk = lambda col, lab_: Line2D([], [], ls="", marker="o", ms=3, color=col, alpha=F.ALPHA, label=lab_)
    hs = [mk(colour[c_], str(c_)) for c_ in shown]
    if kind != "binary":
        hs.append(mk(OTHER, "other classes"))
    hs += [mk(GREYS[1], "unlabelled / cluster not enriched"),
           Line2D([], [], color=F.INK, lw=0.6, ls=":", label="0.75 stability / 0 silhouette")]
    fig.legend(handles=hs, loc="lower center", ncol=min(len(hs), 4), frameon=False, fontsize=F.FONT_PT,
               bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 1 - 2 / Hh, lname, ha="center", va="top", fontsize=F.FONT_PT, weight="bold")
    fig.text(0.995, 0.995, "DEVELOPMENT FIGURE - made by a script (e8_gallery.py)", ha="right", va="top",
             color=F.DARK_GREY, fontsize=F.NOTE_PT)
    fig.savefig(path, dpi=300)
    fig.savefig(path[:-4] + ".svg")
    plt.close(fig)


def write_legend(out, redraws):
    open(os.path.join(out, "legend.md"), "w").write(
        "**E8 gallery: the map, clusters and cluster quality behind every cell of figure E8 panel A.** "
        "PRELIMINARY (development figures). One page per label; rows: all views (reference), the features E1 "
        "selected on development genes, and random draw 1 of the same size, each rebuilt as E8 builds it "
        "(per view median fill, z-score, PCA to at most 30 components; blocks side by side). "
        "Column 1: UMAP (cosine, 15 neighbours, min_dist 0.3, seed 0) coloured by the label's five largest "
        "classes (binary labels: the positive class, and the negative class when it is a measured set); dark "
        "grey = other classes; light grey = unlabelled. Columns 2-3: the same map coloured by Louvain cluster "
        "(cosine kNN 15, seed 0) at resolution 0.5 and 1.0; a cluster enriched for a class (one-sided "
        "hypergeometric over the labelled genes, Bonferroni over classes, Benjamini-Hochberg over clusters, "
        "q < 0.05) takes that class's colour, the others alternate light greys. Bars, per cluster ordered by "
        "top-class share: top class (%) = share of the cluster's labelled genes in its most enriched class "
        "(coloured if enriched, grey if not); stability = mean over "
        f"{redraws} redraws (columns resampled with replacement within each view, new seed) of the best Jaccard "
        "overlap with any redraw cluster (Hennig 2007: > 0.75 stable, < 0.5 dissolved; dotted line 0.75); "
        "silhouette = mean cosine silhouette of the cluster's genes (dotted line 0). Per-cluster numbers: "
        "clusters.csv; per-gene cluster, silhouette and confidence (mean Jaccard of the gene's cluster with "
        "where it lands in each redraw): genes.parquet.\n")


if __name__ == "__main__":
    sys.exit(main())
