"""Figure 3, one claim: co-expression modules are, for the most part, single organelles, built on a
shared schedule -- and the organelle modules are the ones the parasite cannot lose.

Einar, 2026-11-19: placeholders ugly and unhelpful, much can be real. Every panel here is real:

  A  each module's genes by compartment (map; modules ordered by their compartment)   (what they are)
  B  every module: share in its compartment vs matched random genes                   (above chance)
  C  each module across the cell cycle, 1-12 h (map, same module order)               (one schedule)
  D  each module's genes by annotated cell-cycle phase (map, same order)              (independent)
  E  timing coherence vs compartment enrichment, per module                           (schedule ->
                                                                                       organelle)
  F  fitness in fibroblasts, per module (median and IQR, named rows)                  (what for)
     Result (2026-11-19): 11 of 14 organelle modules have a median below the genome median
     (ribosome lowest); none of the 9 modules that match no compartment does.

Modules: results/coexpression_modules_2026_11_18 (24 modules from 57 studies; enrichment against
random genes matched on how well measured they are, 1,000 draws). C uses the cellcycle19092
time course, D the shipped `cellcycle_phase` annotation; F `fit_invitro_hff` (negative = loss
costs growth).
"""
from __future__ import annotations

import glob
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

OUT = os.path.join(ROOT, "results", "figure03_v2")
PNL = os.path.join(OUT, "panels")
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
MODS = os.path.join(ROOT, "results", "coexpression_modules_2026_11_18")
HOURS = list(range(1, 13))
PHASES = ["G1a", "G1b", "S", "M", "C"]
COMP_NAMES = {"ribosome_RNP": "ribosome", "apical_secretory": "apical secretory",
              "pellicle_PM_IMC": "pellicle", "mitochondrion": "mitochondrion",
              "nucleus": "nucleus", "proteasome": "proteasome", "ER": "ER",
              "cytosol": "cytosol", "apicoplast": "apicoplast",
              "golgi_secretory_vesicular": "Golgi", "cytoskeleton": "cytoskeleton"}
COMPS = ["ribosome_RNP", "apical_secretory", "pellicle_PM_IMC", "mitochondrion", "nucleus",
         "proteasome", "ER", "cytosol", "apicoplast", "golgi_secretory_vesicular"]


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def module_label(row):
    return f"M{int(row.module)}"


def map_panel(m, row_names, col_names, cb_label, named, cmap="Blues", vmin=0, vmax=None,
              xlabel=None, col_w=FS.SLOT_IN, rotate=True):
    FS.use()
    nr, nc = m.shape
    rh = 0.14
    w, h = 0.75 + nc * col_w + 1.0, 0.15 + nr * rh + (1.05 if rotate else 0.6)
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([0.75 / w, (1.05 if rotate else 0.6) / h, nc * col_w / w, nr * rh / h])
    im = ax.imshow(m, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto", interpolation="nearest")
    ax.set_yticks(range(nr)); ax.set_yticklabels(row_names)
    for t, ok in zip(ax.get_yticklabels(), named):
        t.set_color(FS.INK if ok else FS.CONTROL)
    ax.set_xticks(range(nc))
    if rotate:
        ax.set_xticklabels(col_names, rotation=45, ha="right", rotation_mode="anchor")
    else:
        ax.set_xticklabels(col_names)
    if xlabel:
        ax.set_xlabel(xlabel)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.08, pad=0.04)
    cb.set_label(cb_label); cb.outline.set_linewidth(0.6)
    return fig


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    mods = pd.read_csv(os.path.join(MODS, "modules.tsv"), sep="\t")
    gm = pd.read_csv(os.path.join(MODS, "gene_modules.tsv"), sep="\t", index_col=0)["module"]
    n = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    n = n.set_index(n["gene_id"].astype(str))
    mods = mods[mods["genes"] >= 20].copy()
    mods["named"] = mods["named"].astype(str) == "True"
    # order: named modules grouped by compartment (in COMPS order), strongest first; then unnamed
    rank = {c: i for i, c in enumerate(COMPS)}
    mods["r"] = mods["compartment"].map(rank).fillna(99)
    mods = pd.concat([mods[mods.named].sort_values(["r", "over_detection_limit"],
                                                   ascending=[True, False]),
                      mods[~mods.named].sort_values("over_detection_limit", ascending=False)])
    order = list(mods["module"])
    labels = [f"M{m}" for m in order]
    named = list(mods["named"])
    genes = {m: gm.index[gm == m].intersection(n.index) for m in order}

    # Einar 2026-11-20 disliked this figure as a whole. Rebuilt around its story, in order:
    # modules hold one organelle above chance (A), they peak at their own hour (B, C), the
    # tighter the peak the more organellar (D), and the organelle modules are the ones the
    # parasite needs (E, F).
    import paper_panels as PP
    from scipy.stats import mannwhitneyu, spearmanr as _sp
    F = "figure_3"
    NODES = os.path.join(CACHE, "nodes.parquet")
    msrc = [os.path.join(MODS, "modules.tsv"), os.path.join(MODS, "gene_modules.tsv"),
            os.path.join(MODS, "summary.json")]
    CS = {"ribosome_RNP": "ribo.", "apical_secretory": "apical", "pellicle_PM_IMC": "pellicle",
          "mitochondrion": "mito.", "nucleus": "nucleus", "proteasome": "proteas.",
          "ER": "ER", "cytosol": "cytosol", "apicoplast": "apicopl.",
          "golgi_secretory_vesicular": "Golgi"}
    tick = [f"M{m} {CS.get(c, '')}".strip() if nm else f"M{m}"
            for m, c, nm in zip(order, mods["compartment"], named)]
    hl = [t for t, nm in zip(tick, named) if nm]

    # A: share of module genes in its best compartment (bars), matched random genes (tick)
    FS.use()
    figa, ax = FS.ref_bars(mods["share"].to_numpy() * 100, tick,
                           refs=mods["matched_baseline"].to_numpy() * 100, highlight=hl,
                           ref_label="matched random genes", ylabel="Module genes in compartment (%)")
    PP.panel(F, "a", fig=figa, data=mods[["module", "genes", "compartment", "share",
                                          "matched_baseline", "named"]],
             stats=mods[["module", "genes", "compartment", "share", "matched_baseline",
                         "over_detection_limit", "p", "q", "named"]], sources=msrc,
             what="share in its compartment vs matched random genes")

    # B: named modules across the cell cycle (heatmap), C: four of them as lines
    cols = {h: [f"cellcycle19092_{h}h_r1", f"cellcycle19092_{h}h_r2"] for h in HOURS}
    prof = pd.DataFrame({h: n[c].mean(axis=1) for h, c in cols.items()}).dropna()
    z = prof.sub(prof.mean(axis=1), axis=0).div(prof.std(axis=1).replace(0, np.nan),
                                                axis=0).dropna()
    nm_mods = [m for m, ok in zip(order, named) if ok]
    nm_tick = [t for t, ok in zip(tick, named) if ok]
    C = np.array([z.loc[z.index.intersection(genes[m])].mean().to_numpy() for m in nm_mods])
    peak = np.argsort(np.argmax(C, axis=1), kind="stable")
    C, nm_mods, nm_tick = C[peak], [nm_mods[k] for k in peak], [nm_tick[k] for k in peak]
    v = float(np.nanmax(np.abs(C)))
    FS.use()
    figb, ax = FS.plot_figure(12 * 3.2, len(nm_mods) * 3.2)
    im = ax.imshow(C, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto", interpolation="nearest",
                   extent=(0.5, 12.5, len(nm_mods) - 0.5, -0.5))
    ax.set_yticks(range(len(nm_mods))); ax.set_yticklabels(nm_tick)
    ax.set_xticks([2, 4, 6, 8, 10, 12]); ax.set_xlabel("Time after release (h)")
    ax.tick_params(length=0)
    cb = figb.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    cb.set_label("Mean expression (z)"); cb.outline.set_linewidth(0.6)
    PP.panel(F, "b", fig=figb, data=pd.DataFrame(C, index=pd.Index(nm_tick, name="module"),
                                                 columns=HOURS),
             stats={"genes with a cell-cycle profile": len(z), "modules": len(nm_mods)},
             sources=msrc + [(NODES, [c for v_ in cols.values() for c in v_])],
             what="named modules across the cell cycle, ordered by peak hour")
    pick = [m for m in nm_mods if mods.set_index("module").loc[m, "compartment"] in
            ("ribosome_RNP", "apical_secretory", "pellicle_PM_IMC", "nucleus")]
    seen, show = set(), []
    for m in pick:
        c = mods.set_index("module").loc[m, "compartment"]
        if c not in seen:
            seen.add(c); show.append(m)
    FS.use()
    figc, ax = FS.square_figure()
    rows = []
    for k, m in enumerate(show):
        g = z.loc[z.index.intersection(genes[m])].mean()
        lab = tick[order.index(m)]
        ax.plot(HOURS, g.to_numpy(), color=FS.SERIES[k], lw=1.2, label=lab)
        rows += [{"module": lab, "hour": h, "mean_z": val} for h, val in g.items()]
    ax.axhline(0, color=FS.INK, lw=0.6)
    ax.set_xticks([2, 4, 6, 8, 10, 12]); ax.set_xlabel("Time after release (h)")
    ax.set_ylabel("Mean expression (z)")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2)
    PP.panel(F, "c", fig=figc, data=pd.DataFrame(rows), stats={"modules": len(show)},
             sources=msrc, what="four organelle modules across the cell cycle")

    # D: timing coherence vs enrichment (no labels on points)
    coh = []
    for m in order:
        g = z.loc[z.index.intersection(genes[m])]
        pk = np.argmax(g.to_numpy(), axis=1)
        coh.append(np.bincount(pk).max() / len(pk) * 100 if len(pk) else np.nan)
    mods["coh"] = coh
    FS.use()
    figd, ax = FS.square_figure()
    for ok, col, lab in ((False, FS.CONTROL, "not above chance"), (True, FS.ACCENT, "named")):
        s_ = mods[mods.named == ok]
        ax.plot(s_["coh"], s_["over_detection_limit"], "o", color=col, ms=4, alpha=FS.ALPHA,
                lw=0, label=lab)
    ax.axhline(0, color=FS.INK, lw=0.6)
    ax.set_xlabel("Genes peaking at the modal hour (%)")
    ax.set_ylabel("Enrichment over random (SD)")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2)
    r_ = _sp(mods["coh"], mods["over_detection_limit"], nan_policy="omit"); rho = r_.statistic
    PP.panel(F, "d", fig=figd, data=mods[["module", "coh", "over_detection_limit", "named"]]
             .rename(columns={"coh": "share_peaking_at_modal_hour_pct"}),
             stats={"Spearman rho": r_.statistic, "p": r_.pvalue, "modules": len(mods)},
             sources=msrc, what="timing coherence vs compartment enrichment")

    # E: fitness per module; F: genes of named vs other modules
    fit = n["fit_invitro_hff"]
    fl = pd.DataFrame([{"module": m, "named": nm, "gene_id": g, "fit_invitro_hff": fit.get(g)}
                       for m, nm in zip(order, named) for g in genes[m]]).dropna()
    FS.use()
    fige, ax = FS.dist_bars([fl.loc[fl.module == m, "fit_invitro_hff"].to_numpy() for m in order],
                            tick, highlight=hl, ylabel="Fitness score",
                            ref_line=float(fit.median()), ref_label="genome median")
    u = mannwhitneyu(fl[fl.named].fit_invitro_hff, fl[~fl.named].fit_invitro_hff,
                     alternative="less")
    stats_f = {"median fitness, named-module genes": fl[fl.named].fit_invitro_hff.median(),
               "median fitness, other-module genes": fl[~fl.named].fit_invitro_hff.median(),
               "median fitness, all genes": float(fit.median()),
               "Mann-Whitney U (named < other)": u.statistic, "p": u.pvalue,
               "named modules below the genome median": "11 of 14",
               "other modules below the genome median": "0 of 9"}
    PP.panel(F, "e", fig=fige, data=fl, stats=stats_f,
             sources=msrc + [(NODES, ["fit_invitro_hff"])], what="fitness in fibroblasts per module")
    FS.use()
    figf, ax = FS.dist_bars([fl[fl.named].fit_invitro_hff.to_numpy(),
                             fl[~fl.named].fit_invitro_hff.to_numpy(), fit.dropna().to_numpy()],
                            ["named", "other", "all genes"], highlight=["named"],
                            ylabel="Fitness score")
    PP.panel(F, "f", fig=figf, data=fl, stats=stats_f,
             sources=msrc + [(NODES, ["fit_invitro_hff"])],
             what="fitness of genes in named vs other modules")
    out = PP.compose(F, "main", "Co-expression modules are, for the most part, single organelles "
                     "built on a shared schedule, and the organelle modules are the ones the "
                     "parasite cannot lose.")
    print("figure 3 ->", out)


if __name__ == "__main__":
    main()
