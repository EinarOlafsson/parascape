"""Quality of a clustering: how stable, how tight, and how certain each gene's call is.

- stability (per cluster): Hennig's clusterboot idea. Re-cluster B times, each time with the columns of every
  view redrawn with replacement and a new seed; for each original cluster take the best Jaccard overlap with any
  cluster of the redraw; average over the B redraws. Guide (Hennig 2007): > 0.75 stable, 0.6-0.75 patterns,
  0.5-0.6 weak, < 0.5 dissolved.
- tightness (per cluster): mean silhouette of its genes (cosine distance), -1 to 1.
- confidence (per gene): the mean over redraws of the Jaccard between the gene's original cluster and the
  cluster it lands in; 1 = always with the same companions.
"""
from __future__ import annotations

import numpy as np
from sklearn.metrics import silhouette_samples


def jaccard_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Jaccard overlap between every cluster of a (rows) and of b (columns); -1 labels ignored."""
    ok = (a >= 0) & (b >= 0)
    ka, kb = a.max() + 1, b.max() + 1
    inter = np.zeros((ka, kb))
    np.add.at(inter, (a[ok], b[ok]), 1)
    na = np.bincount(a[ok], minlength=ka)[:, None]
    nb = np.bincount(b[ok], minlength=kb)[None, :]
    union = na + nb - inter
    return np.where(union > 0, inter / np.maximum(union, 1), 0.0)


def stability(ref: np.ndarray, redraws: list) -> tuple:
    """(per-cluster mean best Jaccard, per-gene confidence) of `ref` over the redraw labelings."""
    k = ref.max() + 1
    per_cluster = np.zeros((len(redraws), k))
    per_gene = np.zeros((len(redraws), len(ref)))
    for r, lab in enumerate(redraws):
        J = jaccard_matrix(ref, lab)
        per_cluster[r] = J.max(axis=1)
        ok = (ref >= 0) & (lab >= 0)
        per_gene[r, ok] = J[ref[ok], lab[ok]]
    return per_cluster.mean(axis=0), per_gene.mean(axis=0)


def tightness(E: np.ndarray, labels: np.ndarray, max_genes: int = 8000, seed: int = 0) -> tuple:
    """(per-cluster mean silhouette, per-gene silhouette), cosine distance; NaN for clusters of one gene."""
    n = len(labels)
    idx = np.arange(n)
    if n > max_genes:
        idx = np.sort(np.random.default_rng(seed).choice(n, max_genes, replace=False))
    s = np.full(n, np.nan)
    if len(np.unique(labels[idx])) > 1:
        s[idx] = silhouette_samples(E[idx], labels[idx], metric="cosine")
    k = labels.max() + 1
    per = np.array([np.nanmean(s[labels == j]) if (labels == j).sum() > 1 else np.nan for j in range(k)])
    return per, s
