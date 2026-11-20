"""Figure 11, one claim: the map fills in where a protein is and when in the cell cycle it is made at
a validated rate well above guessing; mutant phenotypes it predicts no better than the commonest
answer.

Einar, 2026-11-19: placeholders unhelpful. All panels real, from the held-out benchmark
results/filling_the_map_2026_11_18 (every target x strategy: answer rate, accuracy when
answered, the commonest-answer baseline), and the node tables for how many genes still lack
each label.

  A  lift over guessing, every target x strategy (map)                          (the overview)
  B  per target, best strategy: accuracy when answered vs commonest answer      (validated rate)
  C  accuracy over baseline vs answer rate, every target x strategy             (the trade-off)
  D  lift by kind of target, all strategies                                    (place yes, phenotype no)
  E  per target: genes still unlabelled, and how many the best strategy fills   (what it adds)
  F  per strategy: median lift across targets                                  (what generalises)
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402

OUT = os.path.join(ROOT, "results", "figure11_v2")
PNL = os.path.join(OUT, "panels")
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
FM = os.path.join(ROOT, "results", "filling_the_map_2026_11_18")
TARGETS = {  # (organism, target) -> (display name, kind)
    ("Tg", "compartment"): ("Tg compartment", "where"),
    ("Tg", "lopit_unified"): ("Tg hyperLOPIT", "where"),
    ("Pf", "lopit_pf_location"): ("Pf hyperLOPIT", "where"),
    ("Tg", "cellcycle_phase"): ("Tg cell-cycle phase", "when"),
    ("Tg", "dtm_class"): ("Tg membrane topology", "sequence"),
    ("Pf", "is_exported"): ("Pf exported", "sequence"),
    ("Tg", "screen_any_phenotype"): ("Tg any screen phenotype", "phenotype"),
    ("Tg", "screen_actin_phenotype"): ("Tg actin phenotype", "phenotype"),
    ("Tg", "screen_apicoplast_phenotype"): ("Tg apicoplast phenotype", "phenotype"),
    ("Tg", "screen_egress_phenotype"): ("Tg egress phenotype", "phenotype"),
    ("Tg", "screen_replication_phenotype"): ("Tg replication phenotype", "phenotype"),
    ("Pf", "pb_transferred_phenotype"): ("Pf essentiality (from Pb)", "phenotype"),
}
KIND_NAMES = {"where": "location", "when": "cell-cycle phase", "sequence": "sequence class",
              "phenotype": "mutant phenotype"}

TSHORT = {"Tg compartment": "Tg compartment", "Tg hyperLOPIT": "Tg LOPIT",
          "Pf hyperLOPIT": "Pf LOPIT", "Tg cell-cycle phase": "Tg phase",
          "Tg membrane topology": "Tg topology", "Pf exported": "Pf exported",
          "Tg any screen phenotype": "Tg any pheno.", "Tg actin phenotype": "Tg actin",
          "Tg apicoplast phenotype": "Tg apicoplast", "Tg egress phenotype": "Tg egress",
          "Tg replication phenotype": "Tg replication",
          "Pf essentiality (from Pb)": "Pf essential"}


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def rows_fig(nr, left=1.5, value=1.6, right=1.6, rh=0.2):
    w, h = left + value + right, 0.15 + nr * rh + 0.6
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([left / w, 0.6 / h, value / w, nr * rh / h])
    ax.set_ylim(nr - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    return fig, ax


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    t = pd.read_csv(os.path.join(FM, "by_target_and_strategy.tsv"), sep="\t")
    t["key"] = list(zip(t.organism, t.target))
    t = t[t.key.isin(TARGETS)]
    t["name"] = t.key.map(lambda k: TARGETS[k][0])
    t["kind"] = t.key.map(lambda k: TARGETS[k][1])
    best = pd.read_csv(os.path.join(FM, "best_per_target.tsv"), sep="\t")
    best["key"] = list(zip(best.organism, best.target))
    best = best[best.key.isin(TARGETS)].copy()
    best["name"] = best.key.map(lambda k: TARGETS[k][0])
    best["kind"] = best.key.map(lambda k: TARGETS[k][1])
    kind_rank = {"where": 0, "when": 1, "sequence": 2, "phenotype": 3}
    best = best.sort_values(["kind", "lift_over_baseline"], key=lambda s: s.map(kind_rank)
                            if s.name == "kind" else -s)
    order = list(best["name"])

    # A: lift map, targets x strategies
    strat = list(t.groupby("strategy")["lift_over_baseline"].median().sort_values(
        ascending=False).index)
    m = t.pivot_table(index="name", columns="strategy", values="lift_over_baseline").reindex(
        index=order, columns=strat) * 100
    FS.use()
    nr, nc = m.shape; cell = 0.16
    w, h = 1.0 + nc * cell + 0.9, 0.15 + nr * cell + 1.25
    figa = plt.figure(figsize=(w, h))
    ax = figa.add_axes([1.0 / w, 1.25 / h, nc * cell / w, nr * cell / h])
    cmap = matplotlib.colormaps["RdBu_r"].copy(); cmap.set_bad(FS.BACKGROUND)
    im = ax.imshow(np.ma.masked_invalid(m.to_numpy()), cmap=cmap, vmin=-50, vmax=50,
                   interpolation="nearest")
    ax.set_xticks(range(nc)); ax.set_xticklabels([s.replace("_", " ") for s in strat],
                                                 rotation=90)
    ax.set_yticks(range(nr)); ax.set_yticklabels([TSHORT.get(o, o) for o in order])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = figa.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    cb.set_label("Accuracy over guessing (points)"); cb.outline.set_linewidth(0.6)

    # B: best strategy per target (bars), commonest answer (tick)
    b = best.set_index("name").loc[order]
    FS.use()
    nb = [TSHORT.get(o, o) for o in order]
    figb, ax = FS.ref_bars(b["right_when_answered"].to_numpy() * 100, nb,
                           refs=b["commonest_answer_baseline"].to_numpy() * 100,
                           highlight=[TSHORT.get(o, o) for o, k in zip(order, b["kind"])
                                      if k in ("where", "when")],
                           ref_label="commonest answer", ylabel="Held-out calls correct (%)")
    ax.set_ylim(0, 105)

    # C: lift vs answer rate
    FS.use()
    figc, ax = FS.square_figure()
    ax.axhline(0, color=FS.INK, lw=0.6)
    sel = t["kind"].isin(["where", "when"])
    ax.plot(t.loc[~sel, "answer_rate"] * 100, t.loc[~sel, "lift_over_baseline"] * 100, "o",
            color=FS.CONTROL, ms=3, alpha=FS.ALPHA, lw=0, label="phenotype or sequence (ref.)")
    ax.plot(t.loc[sel, "answer_rate"] * 100, t.loc[sel, "lift_over_baseline"] * 100, "o",
            color=FS.ACCENT, ms=3, alpha=FS.ALPHA, lw=0, label="location or phase")
    ax.set_xlabel("Held-out genes answered (%)")
    ax.set_ylabel("Accuracy over guessing (points)")
    ax.set_xlim(0, 103)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    # D: lift by kind
    kinds = ["where", "when", "sequence", "phenotype"]
    # A strategy that answered nothing has no lift; it is left out, not scored as zero.
    groups = [t.loc[t.kind == k, "lift_over_baseline"].dropna().to_numpy() * 100 for k in kinds]
    names = [KIND_NAMES[k] for k in kinds]
    FS.use()
    figd, ax = FS.figure(4, vertical=True)
    FS.violin(ax, groups, names, control=names[2:], vertical=True, draw_points=True)
    ax.axhline(0, color=FS.INK, lw=0.6)
    ax.set_xticks(range(4)); ax.set_xticklabels(names, rotation=45, ha="right",
                                                rotation_mode="anchor")
    ax.set_xlim(-0.6, 3.4); ax.tick_params(axis="x", length=0)
    ax.set_ylabel("Accuracy over guessing (points)")
    FS.pad_limits(ax, np.concatenate(groups))
    for k, g in zip(kinds, groups):
        print(f"  {k:10s} median lift {np.nanmedian(g):+.1f} points over {len(g)} strategies")

    # E: genes still unlabelled and what the best strategy would fill
    tg = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    pf = pd.read_parquet(os.path.join(CACHE, "pf_nodes.parquet"))
    rows = []
    for _, r in best.iterrows():
        nodes = tg if r.organism == "Tg" else pf
        if r.target not in nodes.columns:
            continue
        # isna() first: a pandas NA prints as "<NA>", so a text test alone counted none missing.
        raw = nodes[r.target]
        col = raw.astype(str).str.strip().str.lower()
        missing = int((raw.isna() | col.isin(["nan", "none", "", "unassigned", "unknown",
                                              "<na>"])).sum())
        filled = missing * r.answer_rate
        rows.append((r["name"], missing, filled, filled * r.right_when_answered, r.kind,
                     r.lift_over_baseline))
    e = pd.DataFrame(rows, columns=["name", "missing", "filled", "right", "kind", "lift"])
    e = e[(e.missing > 0)].set_index("name")
    e = e.reindex([o for o in order if o in e.index])   # fully labelled targets drop out
    FS.use()
    ne = [TSHORT.get(o, o) for o in e.index]
    good = (e["lift"] > 0.05).to_numpy()
    fige, ax = FS.grouped_bars([e["missing"].to_numpy(), np.where(good, e["right"], np.nan)],
                               ["genes unlabelled", "expected correct new calls"], ne,
                               ylabel="Genes", log=True)
    ax.set_yticks([10, 100, 1000, 10000]); ax.set_yticklabels(["10", "100", "1,000", "10,000"])
    ax.minorticks_off()
    print(e.round(0).to_string())

    # F: per strategy, accuracy over guessing across targets (median bar, one point per target)
    s = t.groupby("strategy")["lift_over_baseline"].agg(["median", lambda v: (v > 0).mean()])
    s.columns = ["median", "share_pos"]
    s = s.sort_values("median", ascending=False)
    FS.use()
    figf, ax = FS.dist_bars([t.loc[t.strategy == k, "lift_over_baseline"].dropna().to_numpy() * 100
                             for k in s.index], [k.replace("_", " ") for k in s.index],
                            highlight=[k.replace("_", " ") for k in s.index if s.loc[k, "median"] > 0],
                            ylabel="Accuracy over guessing (points)")
    ax.axhline(0, color=FS.INK, lw=0.6)

    import paper_panels as PP
    F = "figure_11"
    from scipy.stats import kruskal
    FMF = [os.path.join(FM, f) for f in ("best_per_target.tsv", "by_target_and_strategy.tsv",
                                         "summary.json")]
    PP.panel(F, "a", fig=figa, data=m.rename_axis("target"), stats={
             "targets": m.shape[0], "strategies": m.shape[1]}, sources=FMF,
             what="lift over guessing, every target x strategy")
    PP.panel(F, "b", fig=figb, data=b, stats=b[["strategy", "answer_rate",
             "right_when_answered", "commonest_answer_baseline", "lift_over_baseline"]]
             .reset_index(), sources=FMF, what="best strategy per target vs commonest answer")
    PP.panel(F, "c", fig=figc, data=t.drop(columns=["key"]), stats={"target x strategy": len(t)},
             sources=FMF, what="accuracy over guessing vs answer rate")
    kw = kruskal(*[g_ for g_ in groups if len(g_)])
    PP.panel(F, "d", fig=figd, data=pd.concat([pd.DataFrame({"kind": n_, "lift_points": g_})
             for n_, g_ in zip(names, groups)]), stats=pd.DataFrame(
             [{"kind": n_, "strategies": len(g_), "median_lift_points": float(np.median(g_))}
              for n_, g_ in zip(names, groups)] + [{"kind": "Kruskal-Wallis", "H": kw.statistic,
                                                     "p": kw.pvalue}]),
             sources=FMF, what="lift by kind of target")
    PP.panel(F, "e", fig=fige, data=e, stats={"expected correct = unlabelled x answer rate x "
             "accuracy when answered": "shown only where lift > 5 points"},
             sources=FMF + [os.path.join(CACHE, "nodes.parquet"),
                            os.path.join(CACHE, "pf_nodes.parquet")],
             what="unlabelled genes vs expected correct new calls")
    PP.panel(F, "f", fig=figf, data=s.rename_axis("strategy"), stats={"targets": t.name.nunique()},
             sources=FMF, what="median lift per strategy")
    out = PP.compose(F, "main", "The map fills in where a protein is and when it is made at a validated rate; mutant phenotypes it predicts no better than the commonest answer.", rows=[3, 3])
    print("figure 11 ->", out)


if __name__ == "__main__":
    main()
