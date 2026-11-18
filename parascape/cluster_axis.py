"""Which cluster does an external per-gene measurement reach?

This is the engine behind the bradyzoite drug analysis, made reusable. A single-cell experiment gives
clusters and their marker genes. A bulk experiment -- a drug treatment, a knockout, a time course --
gives one value per gene. The two were never meant to meet, and the only thing they share is the gene
space. The question they can answer together is: **does any cluster have its markers moved more than
genes like them would be moved anyway?**

Three rules are built in because they are the ones that are easy to forget and expensive to get wrong.

**The null is a matched random set, never the size of the category.** A small marker set wins a
comparison by shape alone, so the baseline is what the same procedure gives on random gene sets drawn
to match the marker set in SIZE and in EXPRESSION LEVEL. Dividing an effect by how many genes are in
the category is the mistake this replaces.

**A negative is a limit until the detection floor says otherwise.** Every cluster carries the smallest
effect its own comparison could have told from its null, so "no effect" and "could not have seen one"
are different answers and are reported differently.

**A cluster the axis barely measured says UNKNOWN, not zero.** Where fewer than `min_measured` of a
cluster's markers carry a value in the external axis, no test is run. An empty cell that means
"nobody looked" must not read as an empty cell that means "nothing there".
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

#: Fewer measured markers than this and a cluster is reported as unknown rather than tested.
MIN_MEASURED = 15

#: How many matched random sets stand behind each cluster.
N_NULL = 2000

#: The window a matched gene is drawn from, as a fraction of all genes, around the member it stands
#: in for. Each member is matched against its own neighbourhood in the expression ordering rather
#: than against a wide bin, because a set at the EXTREME of a bin is not controlled by that bin.
MATCH_WINDOW = 0.10

#: How far, in standard deviations of the marker set itself, the matched sets may sit from the
#: markers in expression before an effect stops being distinguishable from expression.
MAX_EXPRESSION_GAP = 0.25

#: Kept for callers that passed it; the matching no longer bins.
N_BINS = 10


@dataclass(frozen=True)
class AxisResult:
    """One external axis tested against every cluster."""

    table: pd.DataFrame
    axis_name: str
    n_genes_in_axis: int
    n_null: int = N_NULL
    notes: tuple = field(default_factory=tuple)

    @property
    def reached(self) -> list:
        """Clusters the axis moves more than matched random genes, after correction."""
        t = self.table
        return list(t.index[t["reaches"].to_numpy(dtype=bool)])

    @property
    def unknown(self) -> list:
        """Clusters the axis did not measure well enough to test."""
        t = self.table
        return list(t.index[~t["testable"].to_numpy(dtype=bool)])

    def summary(self) -> str:
        reached, unknown = self.reached, self.unknown
        untouched = [c for c in self.table.index if c not in reached and c not in unknown]
        parts = [f"{self.axis_name}: reaches {', '.join(reached) if reached else 'no cluster'}"]
        if untouched:
            parts.append(f"no effect above the limit of detection in {', '.join(untouched)}")
        if unknown:
            parts.append(f"not measured well enough to test in {', '.join(unknown)}")
        return "; ".join(parts) + "."


def _rank_window(expression: pd.Series, width: float) -> tuple:
    """Order the genes by expression and give each one a window of comparable neighbours.

    Coarse bins do not control a marker set that sits at the EXTREME of its bin: deciles let the
    twenty brightest genes be compared against the fifty brightest, which are dimmer on average, and
    a monotone axis then reads as an effect. Matching each gene against its own neighbourhood in the
    expression ordering removes that. Found by the test that asks for exactly this (2026-11-18).

    :returns: the genes in expression order, and the half-width of the window in gene counts.
    """
    order = expression.sort_values().index.to_numpy()
    half = max(1, int(round(len(order) * width / 2)))
    return order, half


def _matched_draw(rng, members, order, position, half, taken):
    """A random gene set matched gene by gene on where each member sits in the expression ordering."""
    n = len(order)
    out = []
    for gene in members.index:
        centre = position[gene]
        lo, hi = max(0, centre - half), min(n, centre + half + 1)
        for _ in range(8):                      # a few tries to avoid the set drawing itself
            pick = order[rng.integers(lo, hi)]
            if pick not in taken:
                break
        out.append(pick)
    return out


def test_clusters_against_axis(markers: dict, axis: pd.Series, expression: pd.Series, *,
                               axis_name: str = "axis", n_null: int = N_NULL,
                               min_measured: int = MIN_MEASURED, n_bins: int = N_BINS,
                               seed: int = 0, absolute: bool = True) -> AxisResult:
    """Ask, for every cluster, whether the axis moves its markers more than matched random genes.

    :param markers: cluster name -> the genes that mark it.
    :param axis: one value per gene from the external experiment, for example a log fold change.
    :param expression: one value per gene used only for matching, for example mean expression in the
        single-cell data. Matching on it is what stops a highly expressed marker set being scored
        against genes nobody could detect.
    :param absolute: score the size of the move and not its direction. A drug that pushes half a
        programme up and half down has still reached it; set False to ask about direction instead.
    :returns: an :class:`AxisResult` whose table has one row per cluster.

    No cluster is ever reported as unaffected when it was not measured: see `testable`.
    """
    axis = pd.Series(axis).dropna()
    axis.index = axis.index.astype(str)
    expression = pd.Series(expression).dropna()
    expression.index = expression.index.astype(str)
    shared = axis.index.intersection(expression.index)
    if len(shared) == 0:
        raise ValueError("the axis and the expression table share no genes; check the identifiers")
    axis, expression = axis.loc[shared], expression.loc[shared]
    values = axis.abs() if absolute else axis

    order, half = _rank_window(expression, MATCH_WINDOW)
    position = {g: i for i, g in enumerate(order)}
    rng = np.random.default_rng(seed)

    rows = {}
    for cluster, genes in markers.items():
        wanted = pd.Index([str(g) for g in dict.fromkeys(genes)])
        measured = wanted.intersection(values.index)
        row = {"n_markers": len(wanted), "n_measured": len(measured)}
        if len(measured) < min_measured:
            # NOT a zero: the axis did not measure enough of this cluster to say anything
            row.update(observed=np.nan, null_mean=np.nan, null_sd=np.nan, z=np.nan, p=np.nan,
                       detection_limit=np.nan, testable=False, reaches=False,
                       expression_gap=np.nan, expression_matched=False,
                       reason=f"only {len(measured)} of {len(wanted)} markers are in the axis")
            rows[cluster] = row
            continue
        members = values.loc[measured]
        observed = float(members.mean())
        taken = set(measured)
        sets = [_matched_draw(rng, members, order, position, half, taken) for _ in range(n_null)]
        draws = np.array([float(values.loc[g].mean()) for g in sets])
        # How well the matching actually worked. A set at the EDGE of the expression range cannot be
        # matched, because there are no brighter genes to draw: the null is then systematically
        # dimmer and a monotone axis reads as an effect that is really expression. Measure it rather
        # than assume it (2026-11-18).
        member_expr = float(expression.loc[measured].mean())
        null_expr = float(np.mean([float(expression.loc[g].mean()) for g in sets]))
        spread = float(expression.loc[measured].std(ddof=1)) or 1.0
        expression_gap = abs(member_expr - null_expr) / spread
        null_mean, null_sd = float(draws.mean()), float(draws.std(ddof=1))
        z = (observed - null_mean) / null_sd if null_sd > 0 else np.nan
        # two-sided, from the null itself rather than from an assumed shape
        p = (1.0 + np.sum(np.abs(draws - null_mean) >= abs(observed - null_mean))) / (len(draws) + 1.0)
        matched = expression_gap <= MAX_EXPRESSION_GAP
        row.update(observed=observed, null_mean=null_mean, null_sd=null_sd, z=z, p=float(p),
                   # the smallest effect this comparison could have told from its own null
                   detection_limit=float(1.96 * null_sd), testable=True, reaches=False,
                   expression_gap=expression_gap, expression_matched=matched,
                   reason="" if matched else
                   (f"the matched sets are {expression_gap:.2f} standard deviations away in "
                    "expression, so an effect here cannot be told from expression itself"))
        rows[cluster] = row

    table = pd.DataFrame.from_dict(rows, orient="index")
    table["q"] = _benjamini_hochberg(table["p"])
    moved = (table["observed"] - table["null_mean"]).abs()
    table["over_detection_limit"] = moved / table["detection_limit"]
    table["reaches"] = (table["testable"].to_numpy(dtype=bool)
                        & table["expression_matched"].fillna(False).to_numpy(dtype=bool)
                        & (table["q"] < 0.05).fillna(False).to_numpy(dtype=bool)
                        & (table["over_detection_limit"] >= 1.0).fillna(False).to_numpy(dtype=bool))
    notes = []
    if not table["testable"].all():
        notes.append("a cluster with too few measured markers is reported unknown, never as no effect")
    unmatched = table.index[table["testable"].fillna(False).to_numpy(dtype=bool)
                            & ~table["expression_matched"].fillna(True).to_numpy(dtype=bool)]
    if len(unmatched):
        notes.append("the matching could not be centred for " + ", ".join(map(str, unmatched))
                     + "; an effect there cannot be told from expression itself")
    return AxisResult(table=table, axis_name=axis_name, n_genes_in_axis=len(axis),
                      n_null=n_null, notes=tuple(notes))


def _benjamini_hochberg(p: pd.Series) -> pd.Series:
    """Corrected probabilities across the clusters tested; a row with no test keeps no value."""
    out = pd.Series(np.nan, index=p.index, dtype=float)
    have = p.dropna().sort_values()
    if have.empty:
        return out
    n = len(have)
    adjusted = (have.to_numpy() * n / np.arange(1, n + 1))
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0, 1)
    out.loc[have.index] = adjusted
    return out
