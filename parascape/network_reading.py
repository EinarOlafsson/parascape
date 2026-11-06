"""Reading the network: which links count, how far a gene is from another, and how far to trust
what its links suggest (the Network module, instruction 87).

The links are `star_edges` frames (`a`, `b` positions in the gene table, `source` a measured
layer such as co-fitness or crosslink MS, `strength` the link percentile within its layer).
Which layers count and how strong a link must be are the user's choice -- "edges you define" --
and every number here is computed over exactly that choice, so a suggestion and its hit rate
describe the same network.

Pure functions over arrays and frames; no window.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import sparse

from . import map_reading as MR

#: Link steps drawn sharply around the gene in focus; beyond, the network fades.
SHARP_HOPS = 2


def chosen(edges: pd.DataFrame, layers, min_strength: float = 0.0) -> pd.DataFrame:
    """The links of the chosen layers at or above a strength (0-1 percentile within a layer)."""
    if edges is None or not len(edges):
        return pd.DataFrame(columns=["a", "b", "source", "strength"])
    keep = edges["source"].isin(set(layers)) & (edges["strength"].astype(float) >= min_strength)
    return edges.loc[keep, ["a", "b", "source", "strength"]]


def adjacency(links: pd.DataFrame, n: int) -> sparse.csr_matrix:
    """Symmetric gene-by-gene matrix of the strongest link between each pair."""
    if links is None or not len(links):
        return sparse.csr_matrix((n, n))
    a = links["a"].to_numpy(int)
    b = links["b"].to_numpy(int)
    w = links["strength"].to_numpy(float) + 1e-6
    m = sparse.coo_matrix((np.r_[w, w], (np.r_[a, b], np.r_[b, a])), shape=(n, n)).tocsr()
    m.sum_duplicates()
    # A pair joined by several layers keeps one entry; its weight is the strongest of them.
    m.data = np.minimum(m.data, 1.0 + 1e-6)
    return m


def hops(adj: sparse.csr_matrix, sources, max_hops: int = 6) -> np.ndarray:
    """Link steps from the nearest source to every gene (`inf` beyond `max_hops` or unreachable)."""
    n = adj.shape[0]
    dist = np.full(n, np.inf)
    frontier = np.unique(np.asarray([int(s) for s in sources if 0 <= int(s) < n], dtype=int))
    if not len(frontier):
        return dist
    dist[frontier] = 0
    seen = np.zeros(n, bool)
    seen[frontier] = True
    indptr, indices = adj.indptr, adj.indices
    for step in range(1, max_hops + 1):
        nxt = np.unique(np.concatenate([indices[indptr[i]:indptr[i + 1]] for i in frontier])
                        if len(frontier) else np.zeros(0, int))
        nxt = nxt[~seen[nxt]]
        if not len(nxt):
            break
        dist[nxt] = step
        seen[nxt] = True
        frontier = nxt
    return dist


def links_of(links: pd.DataFrame, gene: int) -> pd.DataFrame:
    """One gene's links: partner position, layer and strength, strongest first."""
    if links is None or not len(links):
        return pd.DataFrame(columns=["partner", "source", "strength"])
    a, b = links["a"].to_numpy(int), links["b"].to_numpy(int)
    mine = (a == gene) | (b == gene)
    rows = links[mine]
    partner = np.where(rows["a"].to_numpy(int) == gene, rows["b"].to_numpy(int),
                       rows["a"].to_numpy(int))
    out = pd.DataFrame({"partner": partner, "source": rows["source"].to_numpy(),
                        "strength": rows["strength"].to_numpy(float)})
    return out.sort_values("strength", ascending=False).reset_index(drop=True)


def agreement(adj: sparse.csr_matrix, labels) -> MR.Agreement | None:
    """For every gene of known label with linked genes of known label: do they name it?

    The network counterpart of `map_reading.agreement` -- direct links instead of map
    neighbours, the gene never its own evidence, against the commonest label's share.
    """
    has, values = MR._labelled(labels)
    idx = np.flatnonzero(has)
    if len(idx) < 3:
        return None
    indptr, indices = adj.indptr, adj.indices
    calls, truth = [], []
    for i in idx:
        nbrs = indices[indptr[i]:indptr[i + 1]]
        nbrs = nbrs[has[nbrs] & (nbrs != i)]
        if not len(nbrs):
            continue
        calls.append(Counter(values[nbrs]).most_common(1)[0][0])
        truth.append(values[i])
    if not calls:
        return None
    calls, truth = np.asarray(calls, dtype=object), np.asarray(truth, dtype=object)
    right = calls == truth
    precision, called = {}, {}
    for label in set(calls):
        mask = calls == label
        called[label] = int(mask.sum())
        precision[label] = float(right[mask].mean())
    commonest = Counter(values[idx]).most_common(1)[0][1] / len(idx)
    return MR.Agreement(float(right.mean()), float(commonest), int(len(calls)), precision, called)


def spectral_layout(adj: sparse.csr_matrix, seed: int = 0, dims: int = 10) -> np.ndarray:
    """A 2D arrangement BY THE LINKS, the way UMAP itself starts: a `dims`-dimensional spectral
    embedding of the largest connected part of the network, brought down to 2D by UMAP, so
    densely linked groups sit together. Two eigenvectors alone put this network on a line.
    Genes outside that part are placed in a ring round the outside -- seen as unlinked rather
    than piled in the middle."""
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse.linalg import eigsh
    n = adj.shape[0]
    xy = np.zeros((n, 2))
    rng = np.random.default_rng(seed)
    _k, comp = connected_components(adj, directed=False)
    sizes = np.bincount(comp)
    main = comp == int(np.argmax(sizes)) if len(sizes) else np.zeros(n, bool)
    if main.sum() > dims + 2:
        sub = adj[main][:, main]
        d = np.asarray(sub.sum(axis=1)).ravel()
        dinv = sparse.diags(1.0 / np.sqrt(d))
        lap = sparse.identity(sub.shape[0]) - dinv @ sub @ dinv
        try:
            _vals, vecs = eigsh(lap, k=dims + 1, which="SM", tol=1e-3, maxiter=4000,
                                v0=rng.normal(size=sub.shape[0]))
            emb = vecs[:, 1:]
        except Exception:
            emb = rng.normal(size=(sub.shape[0], dims))
        try:
            import umap
            emb2 = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1,
                             random_state=seed).fit_transform(emb)
        except Exception:
            emb2 = emb[:, :2]
        emb2 = emb2 - np.median(emb2, axis=0)
        scale = np.percentile(np.abs(emb2), 95, axis=0)
        xy[main] = emb2 / np.where(scale > 0, scale, 1.0)
    else:
        main = np.zeros(n, bool)
    out = np.flatnonzero(~main)
    if len(out):
        angle = np.linspace(0, 2 * np.pi, len(out), endpoint=False)
        xy[out] = np.column_stack([1.6 * np.cos(angle), 1.6 * np.sin(angle)])
    return xy


def map_layout(xyz) -> np.ndarray:
    """A 2D arrangement AS ON THE MAP: the map coordinates projected onto their two widest axes."""
    xyz = np.asarray(xyz, dtype=float)
    centred = xyz - np.median(xyz, axis=0)
    _u, _s, vt = np.linalg.svd(centred - centred.mean(axis=0), full_matrices=False)
    xy = centred @ vt[:2].T
    scale = np.percentile(np.abs(xy), 95, axis=0)
    return xy / np.where(scale > 0, scale, 1.0)


@dataclass
class Linkage:
    """A gene's links in words: per layer, and its strongest partners."""

    gene: int
    per_layer: list = field(default_factory=list)    # (layer, n)
    partners: pd.DataFrame | None = None


def linkage(links: pd.DataFrame, gene: int) -> Linkage:
    """A gene's links counted per layer, with its partners strongest first."""
    mine = links_of(links, gene)
    per_layer = Counter(mine["source"].astype(str)).most_common()
    return Linkage(gene, per_layer, mine)
