"""Two rates on every location call (Einar 2026-11-26, from the Infer benchmark E7 outside check).

A call's stated rate is learned on LABELLED genes -- for hyperLOPIT, the marker proteins, which are the clearest
examples of each niche. On the other genes it overstates: in E7, calls stated at 98 % agreed with hyperLOPIT's own
classifier on non-marker genes only 47 % of the time. So a location call carries both:

- `rate on markers`: the stated rate, as the method reports it;
- `rate on other genes`: among the calls stated in the same band, the share whose niche matches the hyperLOPIT
  classifier (`lopit_map`) on genes that are not markers. This is the one to read.

Applied only when the answer's classes are hyperLOPIT niches (most predicted classes are in the classifier's
vocabulary); any other answer is returned unchanged.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

CLASSIFIER = "lopit_map"
COARSE = "compartment_best"
UNASSIGNED = {"unassigned", "nan", "None", "", "NaN"}
BANDS = (0.0, 0.3, 0.5, 0.7, 0.9, 1.0001)
MIN_BAND = 10


def predicted_class(calls: pd.DataFrame, column: str) -> pd.Series:
    """The class a call names: the column itself, or the class in a harness wording '<column> = <class>'."""
    s = calls[column].astype(str)
    return s.str.split(" = ", n=1).str[-1] if s.str.contains(" = ", regex=False).any() else s


def add_two_rates(calls: pd.DataFrame, nodes: pd.DataFrame, predicted: str, stated: str):
    """(calls with `rate on markers` and `rate on other genes`, band table) -- or (calls, None) when the answer
    is not about hyperLOPIT niches or the classifier is not in the table."""
    if not len(calls) or CLASSIFIER not in nodes or predicted not in calls or stated not in calls:
        return calls, None
    raw = nodes.set_index(nodes["gene_id"].astype(str))[CLASSIFIER]
    clf = raw.where(raw.notna(), None).map(lambda v: None if v is None else str(v))
    pred = predicted_class(calls, predicted)
    vocab = set(clf.dropna().unique()) - UNASSIGNED
    if pred.isin(vocab).mean() < 0.5 and COARSE in nodes:
        # An answer in coarse compartments (`compartment_best` is an exact coarsening of the classifier's niches:
        # every niche falls in one compartment) is checked against the classifier's niche coarsened the same way.
        both = nodes[[CLASSIFIER, COARSE]].astype(str)
        both = both[~both[CLASSIFIER].isin(UNASSIGNED) & ~both[COARSE].isin(UNASSIGNED)]
        coarse = both.groupby(CLASSIFIER)[COARSE].agg(lambda v: v.value_counts().index[0])
        clf = clf.map(lambda v: None if v is None else coarse.get(v, v))
        vocab = set(coarse.unique())
    if pred.isin(vocab).mean() < 0.5:
        return calls, None
    got = clf.reindex(calls["gene_id"].astype(str))
    # A gene hyperLOPIT never measured has no classifier call; it cannot be checked, which is not a disagreement.
    missing = (got.isna() | got.astype(str).isin(UNASSIGNED)).to_numpy()
    assigned = got.astype(str).to_numpy()
    defined = pred.isin(vocab).to_numpy() & ~missing
    agrees = pred.to_numpy() == assigned
    rate = calls[stated].astype(float).to_numpy()
    band = np.clip(np.digitize(rate, BANDS) - 1, 0, len(BANDS) - 2)
    rows, other = [], {}
    for b in range(len(BANDS) - 1):
        m = (band == b) & defined
        rows.append({"stated band": f"{BANDS[b]:.1f}-{min(BANDS[b + 1], 1):.1f}", "calls": int((band == b).sum()),
                     "checkable": int(m.sum()), "mean stated": float(rate[band == b].mean()) if (band == b).any()
                     else np.nan, "agree with hyperLOPIT classifier": float(agrees[m].mean()) if m.sum() else np.nan})
        other[b] = float(agrees[m].mean()) if m.sum() >= MIN_BAND else np.nan
    out = calls.copy()
    out["rate on markers"] = rate
    out["rate on other genes"] = [other[b] for b in band]
    out["hyperLOPIT classifier says"] = np.where(missing, "", assigned)
    return out, pd.DataFrame(rows)


def report(table: pd.DataFrame) -> str:
    """The band table as HTML for the report."""
    if table is None or not len(table):
        return ""
    if not table["checkable"].sum():
        return ("<br><b>Two rates.</b> The stated rate was learned on marker proteins and overstates on other genes "
                "(Infer benchmark: calls stated at 98 % agreed with hyperLOPIT's own classifier only 47 % of the time "
                "on genes that are not markers). None of these genes was measured by hyperLOPIT, so their rate on "
                "other genes cannot be checked here: read the stated rate as an upper bound.")
    body = "".join(f"<tr><td>{r['stated band']}</td><td align=right>{r['calls']:,}</td>"
                   f"<td align=right>{r['checkable']:,}</td><td align=right>"
                   f"{'' if pd.isna(r['agree with hyperLOPIT classifier']) else format(r['agree with hyperLOPIT classifier'], '.0%')}"
                   f"</td></tr>" for _, r in table.iterrows())
    return ("<br><b>Two rates.</b> The stated rate was learned on marker proteins, the clearest example of each "
            "niche; on other genes it overstates. Read <i>rate on other genes</i>: how often calls stated in that band "
            "agree with hyperLOPIT's own classifier on genes that are not markers."
            "<table cellspacing=4><tr><th align=left>stated</th><th>calls</th><th>checkable</th>"
            "<th>agree on other genes</th></tr>" + body + "</table>")
