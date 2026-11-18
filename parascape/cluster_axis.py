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
        """One sentence naming all three outcomes: reached, not moved, not measured.

        All three are named even when a group is empty of them, because a summary that lists
        only what was found reads as though nothing else was looked at.
        """
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


def _matched_draw(rng, centres, n, half, taken):
    """A random gene set matched gene by gene on where each member sits in the expression ordering.

    Positions, not names: a draw of two thousand sets against five clusters is ten thousand of
    these, and looking each gene up by label made a run of the whole contrast catalogue take
    minutes. The rng is called in the same order as before, so the sets drawn are the same ones.

    :param centres: the position in the expression ordering of each member of the real set.
    :param taken: the positions of the real set, which a matched set avoids where it can.
    :returns: one position per member.
    """
    out = np.empty(len(centres), dtype=int)
    for i, centre in enumerate(centres):
        lo, hi = max(0, centre - half), min(n, centre + half + 1)
        for _ in range(8):                      # a few tries to avoid the set drawing itself
            pick = int(rng.integers(lo, hi))
            if pick not in taken:
                break
        out[i] = pick
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
        # Everything below works on positions in the expression ordering, so a draw is an array
        # index rather than a label lookup.
        by_value = values.reindex(order).to_numpy(dtype=float)
        by_expression = expression.reindex(order).to_numpy(dtype=float)
        centres = np.array([position[g] for g in measured], dtype=int)
        taken = set(centres.tolist())
        sets = [_matched_draw(rng, centres, len(order), half, taken) for _ in range(n_null)]
        picks = np.vstack(sets)
        draws = np.nanmean(by_value[picks], axis=1)
        # How well the matching actually worked. A set at the EDGE of the expression range cannot be
        # matched, because there are no brighter genes to draw: the null is then systematically
        # dimmer and a monotone axis reads as an effect that is really expression. Measure it rather
        # than assume it (2026-11-18).
        member_expr = float(expression.loc[measured].mean())
        null_expr = float(np.nanmean(by_expression[picks]))
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


#: A cluster keeps at most this many markers when they are derived from an expression block.
TOP_MARKERS = 150

#: A gene is a marker of its best cluster only if it is this many times higher there than in the
#: runner-up. 1.0 would make every gene a marker of whichever cluster happens to be highest.
MIN_MARGIN = 1.25


def markers_from_blocks(frame: pd.DataFrame, columns: dict, *, top: int = TOP_MARKERS,
                        min_margin: float = MIN_MARGIN, gene_column: str = "gene_id") -> dict:
    """Marker sets for each cluster, from a block of one average-expression column per cluster.

    This is how the bradyzoite clusters were turned into gene sets, and it is here rather than in a
    script because the alternative is every user inventing their own rule and none of them saying
    which. A published single-cell paper usually ships exactly this: one average expression per
    gene per cluster, five or ten columns side by side.

    A gene marks the cluster where it is highest, and only when it is at least `min_margin` times
    the runner-up, so a gene expressed evenly everywhere marks nothing. The `top` kept are the ones
    with the largest margin, because a long list of weak markers is dominated by genes that are
    merely slightly higher.

    :param frame: one row per gene, with `gene_column` and the cluster columns.
    :param columns: cluster name -> the column holding its average expression.
    :param min_margin: the ratio over the runner-up a marker must reach.
    :returns: cluster name -> its marker genes, longest margin first. A cluster with no gene of its
        own is present with an empty list rather than absent, because "nothing marks it" is a result.
    :raises ValueError: when a named column is not in the frame, or fewer than two clusters are
        given -- a single column has no runner-up to be higher than.
    """
    if len(columns) < 2:
        raise ValueError("two or more cluster columns are needed; a margin needs a runner-up")
    missing = [c for c in columns.values() if c not in frame.columns]
    if missing:
        raise ValueError("not in the table: " + ", ".join(map(str, missing)))
    if gene_column not in frame.columns:
        raise ValueError(f"no {gene_column!r} column to name the genes by")

    names = list(columns)
    block = frame[[columns[n] for n in names]].apply(pd.to_numeric, errors="coerce")
    block.index = frame[gene_column].astype(str).to_numpy()
    block = block.dropna(how="any")
    values = block.to_numpy(dtype=float)
    out = {n: [] for n in names}
    if len(block) == 0:
        return out

    order = np.argsort(-values, axis=1)
    best, second = order[:, 0], order[:, 1]
    rows = np.arange(len(values))
    top_value = values[rows, best]
    runner_up = values[rows, second]
    # A margin over a runner-up of zero is unbounded, and a gene measured in one cluster only is
    # not evidence of specificity; both are held to the same ratio by a small floor.
    margin = top_value / np.maximum(runner_up, 1e-9)
    keep = (margin >= min_margin) & (top_value > 0)
    for i, name in enumerate(names):
        mine = rows[keep & (best == i)]
        mine = mine[np.argsort(-margin[mine])][:top]
        out[name] = [str(block.index[j]) for j in mine]
    return out


def axis_catalogue(release: str | None = None) -> pd.DataFrame:
    """Every within-study contrast ParaScape has stored, as candidate axes.

    The supply is the expression traits release: one shrunken log2 fold change per gene per
    contrast, each computed inside one study against that study's own control. A drug treatment, a
    knockout, a pH switch and a time point are all the same shape of thing here, which is the point
    -- the engine does not care which, and the person choosing does.

    :returns: one row per contrast with its study, its label and how many genes it measures,
        indexed by the column name. Empty when no release is present.
    """
    from . import perturbation_atlas as PA
    try:
        store = PA.load(release)
    except (FileNotFoundError, IndexError, OSError, KeyError):
        return pd.DataFrame()
    registry, lfc = store["registry"].copy(), store["lfc"]
    registry["n_genes"] = [int(lfc[c].notna().sum()) for c in registry.index]
    registry["release"] = store["release"]
    keep = [c for c in ("study", "contrast_label", "contrast", "n_treatment", "n_control",
                        "n_genes", "release") if c in registry.columns]
    return registry[keep]


def axis_values(column: str, release: str | None = None) -> pd.Series:
    """The per-gene fold changes of one stored contrast, by gene identifier.

    :raises KeyError: when no stored contrast is named `column`.
    """
    from . import perturbation_atlas as PA
    store = PA.load(release)
    if column not in store["lfc"].columns:
        raise KeyError(f"no stored contrast named {column!r}")
    return store["lfc"][column].dropna()


#: Tags that name a QUANTITY rather than a group, so a set of columns carrying one is not a
#: block of clusters. Each of these was a false positive in the chooser before it was listed.
QUANTITY_WORDS = frozenset({"margin", "enriched", "padj", "pvalue", "qvalue", "se", "score",
                            "rank", "sd", "mean", "n", "count", "fdr", "p"})


def _block_tags(columns) -> list:
    """The group tags of a block of columns, or [] when the columns are not one block.

    A block is a set of columns that differ only in the middle: ``bzsub_A_expr`` to
    ``bzsub_E_expr`` is one average expression per subtype, while
    ``iron_depletion_protein_log2fc`` beside ``iron_depletion_protein_padj`` is two different
    quantities of one contrast and is NOT a block. The difference is whether a common quantity
    word is left at the end once the differing part is removed, and it is tested here rather
    than guessed by name, so a dataset added later is classified by the same rule.
    """
    parts = [str(c).split("_") for c in columns]
    if len(parts) < 2:
        return []
    pre = 0
    while all(len(p) > pre + 1 for p in parts) and len({p[pre] for p in parts}) == 1:
        pre += 1
    suf = 0
    while all(len(p) > pre + suf + 1 for p in parts) and len({p[-1 - suf] for p in parts}) == 1:
        suf += 1
    if pre == 0 or suf == 0:
        return []
    tags = ["_".join(p[pre:len(p) - suf]) for p in parts]
    if any(t == "" for t in tags) or len(set(tags)) != len(tags):
        return []
    # Two rules that each removed a real false positive from the chooser. A tag that NAMES a
    # quantity means the columns are two quantities of one measurement, not one measurement of
    # two groups: `stage_enriched_derived` beside `stage_margin_derived` is which stage and how
    # far ahead, not two stages. A tag holding `vs` is already a contrast, which belongs on the
    # axis side of this module and never on the cluster side.
    prefix = {p[i] for p in parts[:1] for i in range(pre)}
    if any(t in QUANTITY_WORDS or "vs" in t.split("_") for t in tags) or "vs" in prefix:
        return []
    return tags


def cluster_blocks(frame: pd.DataFrame) -> dict:
    """The groups-of-one-value-per-cluster already in the gene table, ready to mark.

    Found from the declared datasets rather than from a list kept here, so a dataset added to
    the catalogue appears in the chooser without anyone remembering to add it twice.

    :param frame: the gene table, used to check every column of a block is present and numeric.
    :returns: dataset key -> ``{"name", "provides", "columns": {cluster tag: column}}``.
    """
    from . import datasets as DS
    out = {}
    for dataset in DS.REGISTRY:
        tags = _block_tags(dataset.columns)
        if not tags:
            continue
        if any(c not in frame.columns for c in dataset.columns):
            continue
        if not all(pd.api.types.is_numeric_dtype(frame[c]) for c in dataset.columns):
            continue
        out[dataset.key] = {"name": dataset.name, "provides": dataset.provides,
                            "columns": dict(zip(tags, dataset.columns))}
    return out
