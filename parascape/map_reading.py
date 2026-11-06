"""Reading a map: what a gene's neighbours are, and how far to trust them.

The Map module's two jobs, in numbers. Synthesis -- "where does my gene sit, and what
are the genes around it known to be?" -- is `neighbours`. Inference -- "a gene among
dense granule proteins is probably one" -- is only worth saying with how often that
kind of call is right on this same map, which is `agreement`: for every gene whose
label is known, do its neighbours' labels name it? Reported against chance (always
guessing the commonest label), because a 60 % agreement means nothing when one label
covers 55 % of the genes.

Pure functions over arrays; no window, so every number on the screen is testable.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

#: How many neighbours a gene is read by.
K = 10
#: Values that mean "no label", never counted as a class.
UNLABELLED = {"", "unassigned", "nan", "none", "unknown"}


def _labelled(labels) -> np.ndarray:
    values = np.asarray([str(v).strip() for v in labels], dtype=object)
    return np.array([v.lower() not in UNLABELLED for v in values]), values


def _index(xyz):
    from sklearn.neighbors import NearestNeighbors
    return NearestNeighbors().fit(np.asarray(xyz, dtype=float))


@dataclass
class Neighbourhood:
    """The genes nearest a selection, and the labels among them."""

    members: list                         # indices of the selected genes found on the map
    neighbours: list                      # indices of their nearest genes, nearest first
    counts: list = field(default_factory=list)   # [(label, n)], commonest first, labelled only
    unlabelled: int = 0

    @property
    def top(self):
        """(label, n, of) for the commonest label, or None when no neighbour is labelled."""
        if not self.counts:
            return None
        label, n = self.counts[0]
        return label, n, sum(c for _, c in self.counts)


def neighbours(xyz, labels, selected, placed=None, k: int = K) -> Neighbourhood:
    """The `k` nearest genes to each selected gene, pooled, and their labels.

    For one gene, its own `k` neighbours. For a list, the union of every member's
    neighbours that are not themselves in the list -- what surrounds the list -- with the
    labels counted once per gene.

    :param xyz: (n, d) coordinates of every gene.
    :param labels: n labels ("" or "unassigned" when unknown).
    :param selected: indices of the selected genes.
    :param placed: optional (n,) bool, the genes this map has a position for.
    """
    xyz = np.asarray(xyz, dtype=float)
    n = len(xyz)
    placed = np.ones(n, bool) if placed is None else np.asarray(placed, bool)
    members = [int(i) for i in selected if 0 <= int(i) < n and placed[int(i)]]
    if not members:
        return Neighbourhood([], [])
    pool = np.flatnonzero(placed)
    model = _index(xyz[pool])
    want = min(len(pool), k + len(members) + 1)
    _, found = model.kneighbors(xyz[members], n_neighbors=want)
    member_set = set(members)
    order, seen = [], set()
    for row in found:
        taken = 0
        for j in pool[row]:
            j = int(j)
            if j in member_set:
                continue
            if j not in seen:
                seen.add(j)
                order.append(j)
            taken += 1
            if taken >= k:
                break
    has, values = _labelled(labels)
    counts = Counter(values[j] for j in order if has[j])
    return Neighbourhood(members, order, counts.most_common(),
                         unlabelled=sum(1 for j in order if not has[j]))


@dataclass
class Agreement:
    """How often a gene's neighbours name its own label, on this map."""

    accuracy: float                       # majority of k neighbours == own label
    chance: float                         # share of the commonest label
    n: int                                # labelled genes it was measured on
    precision: dict = field(default_factory=dict)   # label -> right when neighbours say it
    called: dict = field(default_factory=dict)      # label -> how often neighbours say it

    def for_label(self, label: str):
        """(precision, times called) for calls of `label`, or None when never called."""
        if label not in self.precision:
            return None
        return self.precision[label], self.called[label]


def agreement(xyz, labels, placed=None, k: int = K) -> Agreement | None:
    """For every labelled gene, ask its `k` labelled neighbours for a majority label.

    Leave-one-out: a gene is never its own neighbour. None when fewer than `k + 1`
    labelled genes are on the map.
    """
    xyz = np.asarray(xyz, dtype=float)
    n = len(xyz)
    placed = np.ones(n, bool) if placed is None else np.asarray(placed, bool)
    has, values = _labelled(labels)
    idx = np.flatnonzero(has & placed)
    if len(idx) <= k:
        return None
    model = _index(xyz[idx])
    _, found = model.kneighbors(xyz[idx], n_neighbors=k + 1)
    truth = values[idx]
    calls = []
    for row in found:
        votes = Counter(truth[j] for j in row[1:])
        calls.append(votes.most_common(1)[0][0])
    calls = np.asarray(calls, dtype=object)
    right = calls == truth
    precision, called = {}, {}
    for label in set(calls):
        mask = calls == label
        called[label] = int(mask.sum())
        precision[label] = float(right[mask].mean())
    commonest = Counter(truth).most_common(1)[0][1] / len(truth)
    return Agreement(float(right.mean()), float(commonest), int(len(idx)), precision, called)


def label_agreement(xyz, labels, k: int = K) -> tuple:
    """(accuracy, chance) of `agreement` -- the score a map gets for one label."""
    result = agreement(xyz, labels, k=k)
    return (float("nan"), float("nan")) if result is None else (result.accuracy, result.chance)


# ---------------------------------------------------------------------- groups (clusters)
@dataclass
class GroupRow:
    """One cluster read against a label: its best-matching value, and how well it matches."""

    group: int
    size: int
    labelled: int             # members with a known label
    label: str                # the commonest known label among them ("" when none)
    hits: int                 # members with that label
    precision: float          # hits / labelled: how pure the group is for that label
    recall: float             # hits / all genes with that label: how much of the label it holds
    f1: float
    unknown: int              # members whose label is unknown: the candidates


@dataclass
class LabelRow:
    """One label value and the group that captures it best."""

    label: str
    genes: int                # genes with this label
    group: int                # best group by F1, -1 when no group holds any
    precision: float
    recall: float
    f1: float


@dataclass
class Learned:
    """Could the label be learned from the groups? Leave-one-out within each group."""

    accuracy: float           # a gene's group-mates name its label (majority, itself excluded)
    chance: float             # the commonest label's share
    called: int               # labelled genes that got a call (in a group with other labelled members)
    labelled: int             # labelled genes in all
    candidates: int           # unknown genes in groups whose precision clears CANDIDATE_PRECISION


#: A group must be at least this pure for its unknown members to be counted as candidates.
CANDIDATE_PRECISION = 0.6
#: ...and hold at least this many genes of known label, so purity is not one gene.
CANDIDATE_MIN_LABELLED = 5


def _f1(p: float, r: float) -> float:
    return 0.0 if p + r == 0 else 2 * p * r / (p + r)


def group_rows(groups, labels, noise: int = -1) -> list:
    """Every group against the label, best match first (by F1)."""
    groups = np.asarray(groups)
    has, values = _labelled(labels)
    totals = Counter(values[has])
    rows = []
    for g in sorted(set(groups[groups != noise].tolist())):
        inside = groups == g
        known = values[inside & has]
        unknown = int((inside & ~has).sum())
        if not len(known):
            rows.append(GroupRow(int(g), int(inside.sum()), 0, "", 0, float("nan"),
                                 float("nan"), 0.0, unknown))
            continue
        label, hits = Counter(known).most_common(1)[0]
        p, r = hits / len(known), hits / totals[label]
        rows.append(GroupRow(int(g), int(inside.sum()), int(len(known)), str(label), int(hits),
                             float(p), float(r), float(_f1(p, r)), unknown))
    rows.sort(key=lambda row: (-row.f1, -row.size))
    return rows


def label_rows(groups, labels, noise: int = -1) -> list:
    """Every label value and the group that captures it best (by F1), best first."""
    groups = np.asarray(groups)
    has, values = _labelled(labels)
    out = []
    for label, n in Counter(values[has]).items():
        mine = has & (values == label)
        best = (-1, 0.0, 0.0, 0.0)
        for g in set(groups[mine & (groups != noise)].tolist()):
            inside = (groups == g) & has
            hits = int((inside & mine).sum())
            p, r = hits / int(inside.sum()), hits / n
            f = _f1(p, r)
            if f > best[3]:
                best = (int(g), p, r, f)
        out.append(LabelRow(str(label), int(n), *best[:1], float(best[1]), float(best[2]),
                            float(best[3])))
    out.sort(key=lambda row: (-row.f1, -row.genes))
    return out


def learn_from_groups(groups, labels, noise: int = -1) -> Learned | None:
    """For each gene of known label in a group, do its group-mates (itself excluded) name it?

    The honest version of "can I learn this label from my clustering": a gene is never its
    own evidence, genes left unclustered get no call, and the result is set against the
    commonest label's share.
    """
    groups = np.asarray(groups)
    has, values = _labelled(labels)
    if int(has.sum()) == 0:
        return None
    right = called = 0
    for g in set(groups[groups != noise].tolist()):
        members = np.flatnonzero((groups == g) & has)
        if len(members) < 2:
            continue
        counts = Counter(values[members])
        for i in members:
            counts[values[i]] -= 1
            label, n = counts.most_common(1)[0]
            counts[values[i]] += 1
            if n == 0:
                continue
            called += 1
            right += int(label == values[i])
    commonest = Counter(values[has]).most_common(1)[0][1] / int(has.sum())
    candidates = sum(row.unknown for row in group_rows(groups, labels, noise)
                     if row.labelled >= CANDIDATE_MIN_LABELLED
                     and row.precision >= CANDIDATE_PRECISION)
    return Learned(right / called if called else float("nan"), float(commonest), int(called),
                   int(has.sum()), int(candidates))
