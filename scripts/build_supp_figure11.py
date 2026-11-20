"""Supplementary figure 11a: the design choices behind figure 11, each tested with data.

Figure 11 scores every strategy as accuracy when answered minus the rate of guessing the
COMMONEST answer, on held-out genes. Here:

  A  three baselines per target: one over the number of classes, the chance agreement of two
     draws from the label frequencies, and the commonest answer (the strictest; chosen)
  B  accuracy over guessing under each baseline, by kind of target
  C  answering only some genes: correct calls per held-out gene (abstaining counts as wrong)
     against always giving the commonest answer
  D  why mutant phenotype is hard to beat: share of the commonest class per target
  E  held-out calls per target
  F  best lift against held-out calls: small targets do not carry the result
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402
import build_figure11_v2 as F11                         # noqa: E402
import paper_panels as PP                               # noqa: E402

F = "supplementary_figure_11a"


def rows_fig(names, series, xlabel, legend, ref=None, ref_label=None, log=False):
    FS.use()
    nr = len(names); rh = 0.2
    w, h = 1.6 + 1.6 + 1.8, 0.15 + nr * rh + 0.6
    fig = plt.figure(figsize=(w, h)); ax = fig.add_axes([1.6 / w, 0.6 / h, 1.6 / w, nr * rh / h])
    y = np.arange(nr)
    cols = [FS.CONTROL, FS.DARK_GREY, FS.ACCENT]
    for i, (vals, lab) in enumerate(zip(series, legend)):
        ax.plot(vals, y, "o", color=cols[-len(series):][i], ms=3.5, alpha=FS.ALPHA, label=lab)
    if ref is not None:
        ax.axvline(ref, color=FS.INK, lw=0.6, ls="--", label=ref_label)
    if log:
        ax.set_xscale("log")
    ax.set_yticks(y); ax.set_yticklabels(names); ax.set_ylim(nr - 0.5, -0.5)
    ax.tick_params(axis="y", length=0); ax.set_xlabel(xlabel)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return fig


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    t = pd.read_csv(os.path.join(F11.FM, "by_target_and_strategy.tsv"), sep="\t")
    best = pd.read_csv(os.path.join(F11.FM, "best_per_target.tsv"), sep="\t")
    tg = pd.read_parquet(os.path.join(F11.CACHE, "nodes.parquet"))
    pf = pd.read_parquet(os.path.join(F11.CACHE, "pf_nodes.parquet"))
    srcs = [os.path.join(F11.FM, f) for f in ("best_per_target.tsv", "by_target_and_strategy.tsv")]
    rows = []
    for (org, tgt), (name, kind) in F11.TARGETS.items():
        nodes = tg if org == "Tg" else pf
        if tgt not in nodes.columns:
            continue
        v = nodes[tgt].dropna().astype(str)
        v = v[~v.str.lower().isin(["nan", "none", "", "unassigned", "unknown", "<na>"])]
        p = v.value_counts(normalize=True)
        b = best[(best.organism == org) & (best.target == tgt)].iloc[0]
        rows.append({"target": name, "kind": kind, "classes": len(p), "uniform": 1 / len(p),
                     "frequency": float((p ** 2).sum()), "commonest": float(p.max()),
                     "best_accuracy": b.right_when_answered, "answer_rate": b.answer_rate,
                     "held_out_calls": b.held_out_calls, "best_lift": b.lift_over_baseline})
    d = pd.DataFrame(rows)
    order = list(d.sort_values(["kind", "best_lift"], ascending=[True, False]).target)
    d = d.set_index("target").loc[order]

    PP.panel(F, "a", fig=rows_fig(order, [d.uniform * 100, d.frequency * 100, d.commonest * 100],
                                  "Baseline accuracy (%)", ["1 / classes", "two random draws",
                                                            "commonest answer (chosen)"]),
             data=d[["classes", "uniform", "frequency", "commonest"]],
             stats={"targets": len(d)}, sources=srcs + [(os.path.join(F11.CACHE, "nodes.parquet"),
                                                         [k[1] for k in F11.TARGETS])],
             what="three baselines per target")
    lifts = pd.DataFrame({b_: (d.best_accuracy - d[b_]) * 100 for b_ in ("uniform", "frequency",
                                                                        "commonest")})
    lifts["kind"] = d.kind
    byk = lifts.groupby("kind").median()
    kinds = [k for k in ("where", "when", "sequence", "phenotype") if k in byk.index]
    PP.panel(F, "b", fig=rows_fig([F11.KIND_NAMES[k] for k in kinds],
                                  [byk.loc[kinds, c] for c in ("uniform", "frequency", "commonest")],
                                  "Best accuracy over baseline (points)",
                                  ["1 / classes", "two random draws", "commonest answer (chosen)"],
                                  ref=0, ref_label="no gain"),
             data=lifts, stats=byk.reset_index(), sources=srcs,
             what="lift under each baseline, by kind of target")
    d["correct_per_gene"] = d.best_accuracy * d.answer_rate
    PP.panel(F, "c", fig=rows_fig(order, [d.commonest * 100, d.correct_per_gene * 100],
                                  "Correct per held-out gene (%)",
                                  ["always the commonest answer (ref.)", "best strategy, "
                                   "abstaining counted wrong"]),
             data=d[["answer_rate", "best_accuracy", "correct_per_gene", "commonest"]],
             stats={"targets where best strategy beats always-commonest per gene":
                    int((d.correct_per_gene > d.commonest).sum())}, sources=srcs,
             what="answering only some genes, abstentions counted wrong")
    PP.panel(F, "d", fig=rows_fig(order, [d.commonest * 100], "Commonest class (%)",
                                  ["share of labelled genes"]),
             data=d[["classes", "commonest"]], stats={"median commonest share, phenotype":
                                                      float(d[d.kind == "phenotype"].commonest
                                                            .median())}, sources=srcs,
             what="share of the commonest class per target")
    PP.panel(F, "e", fig=rows_fig(order, [d.held_out_calls], "Held-out calls",
                                  ["per target"], log=True),
             data=d[["held_out_calls"]], stats={"min": int(d.held_out_calls.min()),
                                                "max": int(d.held_out_calls.max())}, sources=srcs,
             what="held-out calls per target")
    r = spearmanr(d.held_out_calls, d.best_lift)
    FS.use()
    ff, ax = FS.square_figure()
    sel = d.kind.isin(["where", "when"])
    ax.plot(d.held_out_calls[~sel], d.best_lift[~sel] * 100, "o", color=FS.CONTROL, ms=4,
            alpha=FS.ALPHA, lw=0, label="phenotype or sequence (ref.)")
    ax.plot(d.held_out_calls[sel], d.best_lift[sel] * 100, "o", color=FS.ACCENT, ms=4,
            alpha=FS.ALPHA, lw=0, label="location or phase")
    ax.axhline(0, color=FS.INK, lw=0.6); ax.set_xscale("log")
    ax.set_xlabel("Held-out calls"); ax.set_ylabel("Best accuracy over guessing (points)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "f", fig=ff, data=d[["held_out_calls", "best_lift", "kind"]],
             stats={"Spearman rho": r.statistic, "p": r.pvalue}, sources=srcs,
             what="best lift against held-out calls")
    print(d.round(3).to_string())
    PP.compose(F, "supplementary", "Every design choice behind figure 11, varied or checked.")


if __name__ == "__main__":
    main()
