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

    # A: compartment share per module (of genes with a compartment)
    comp = n["compartment_best"].astype(str)
    A = np.array([[float((comp.loc[genes[m]] == c).sum()) /
                   max(1, int((comp.loc[genes[m]] != "unassigned").sum())) * 100 for c in COMPS]
                  for m in order])
    figa = map_panel(A, labels, [COMP_NAMES[c] for c in COMPS], "Module genes (%)", named)

    # B: share vs matched baseline
    FS.use()
    figb, ax = FS.square_figure()
    mx = max(mods["share"].max(), mods["matched_baseline"].max()) * 100 * 1.08
    ax.plot([0, mx], [0, mx], color=FS.BACKGROUND, lw=0.6)
    for ok, col, lab in ((False, FS.CONTROL, "not above chance (ref.)"),
                         (True, FS.ACCENT, "named (q < 0.05)")):
        s = mods[mods.named == ok]
        ax.plot(s["matched_baseline"] * 100, s["share"] * 100, "o", color=col, ms=4,
                alpha=FS.ALPHA, lw=0, label=lab)
    ax.set_xlim(0, mx); ax.set_ylim(0, mx)
    ax.set_xlabel("Matched random genes (%)"); ax.set_ylabel("Module genes in compartment (%)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    # C: module mean z-profile over the cell cycle
    cols = {h: [f"cellcycle19092_{h}h_r1", f"cellcycle19092_{h}h_r2"] for h in HOURS}
    prof = pd.DataFrame({h: n[c].mean(axis=1) for h, c in cols.items()}).dropna()
    z = prof.sub(prof.mean(axis=1), axis=0).div(prof.std(axis=1).replace(0, np.nan),
                                                axis=0).dropna()
    C = np.array([z.loc[z.index.intersection(genes[m])].mean().to_numpy() for m in order])
    v = float(np.nanmax(np.abs(C)))
    figc = map_panel(C, labels, [str(h) if h % 2 == 0 else "" for h in HOURS],
                     "Mean expression (z)", named, cmap="RdBu_r", vmin=-v, vmax=v,
                     xlabel="Time after release (h)", col_w=0.13, rotate=False)

    # D: annotated phase composition
    ph = n["cellcycle_phase"].astype(str)
    D = np.array([[float((ph.loc[genes[m]] == p).sum()) /
                   max(1, int(ph.loc[genes[m]].isin(PHASES).sum())) * 100 for p in PHASES]
                  for m in order])
    figd = map_panel(D, labels, PHASES, "Phase-annotated genes (%)", named,
                     xlabel="Annotated cell-cycle phase", rotate=False, col_w=0.3)

    # E: timing coherence vs enrichment
    coh = []
    for m in order:
        g = z.loc[z.index.intersection(genes[m])]
        pk = np.argmax(g.to_numpy(), axis=1)
        coh.append(np.bincount(pk).max() / len(pk) * 100 if len(pk) else np.nan)
    mods["coh"] = coh
    FS.use()
    fige, ax = FS.square_figure()
    for ok, col, lab in ((False, FS.CONTROL, "not above chance (ref.)"),
                         (True, FS.ACCENT, "named (q < 0.05)")):
        s = mods[mods.named == ok]
        ax.plot(s["coh"], s["over_detection_limit"], "o", color=col, ms=4, alpha=FS.ALPHA,
                lw=0, label=lab)
    ax.axhline(0, color=FS.INK, lw=0.6)
    ax.set_xlabel("Genes peaking at the modal hour (%)")
    ax.set_ylabel("Enrichment over random (SD)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    from scipy.stats import spearmanr
    rho = spearmanr(mods["coh"], mods["over_detection_limit"], nan_policy="omit").statistic

    # F: fitness per module, named rows
    fit = n["fit_invitro_hff"]
    FS.use()
    # The same row height as the maps (A, C, D), so module M11 sits at the same height in every
    # panel; at the bar slot of 0.25 in the 23 rows stood 150 mm tall.
    rh, nr = 0.14, len(order)
    w, h = 1.35 + 1.6 + 1.3, 0.15 + nr * rh + 0.6
    figf = plt.figure(figsize=(w, h))
    ax = figf.add_axes([1.35 / w, 0.6 / h, 1.6 / w, nr * rh / h])
    for y, m in enumerate(order):
        v = fit.loc[genes[m]].dropna().to_numpy()
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        col = FS.ACCENT if named[y] else FS.CONTROL
        ax.hlines(y, q1, q3, color=col, lw=1.0, alpha=FS.ALPHA)
        ax.plot(med, y, "o", color=col, ms=3)
    ax.axvline(float(fit.median()), color=FS.INK, lw=0.6, ls="--", label="all genes, median")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{l} {COMP_NAMES.get(c, '')}" if ok else l
                        for l, c, ok in zip(labels, mods["compartment"], named)])
    for t, ok in zip(ax.get_yticklabels(), named):
        t.set_color(FS.INK if ok else FS.CONTROL)
    ax.set_ylim(len(order) - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Fitness in fibroblasts (score)")
    ax.plot([], [], "o-", color=FS.ACCENT, ms=4, label="named module")
    ax.plot([], [], "o-", color=FS.CONTROL, ms=4, label="not above chance")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    print(f"{len(order)} modules, {sum(named)} named; coherence vs enrichment Spearman {rho:.3f}")
    import paper_panels as PP
    from scipy.stats import mannwhitneyu, spearmanr as _sp
    F = "figure_3"
    NODES = os.path.join(CACHE, "nodes.parquet")
    msrc = [os.path.join(MODS, "modules.tsv"), os.path.join(MODS, "gene_modules.tsv"),
            os.path.join(MODS, "summary.json")]
    lab = pd.Index(labels, name="module")
    PP.panel(F, "a", fig=figa, data=pd.DataFrame(A, index=lab, columns=[COMP_NAMES[c] for c in COMPS]),
             stats=mods[["module", "genes", "compartment", "share", "matched_baseline",
                         "over_detection_limit", "p", "q", "named"]],
             sources=msrc + [(NODES, ["compartment_best"])],
             what="each module's genes by compartment")
    PP.panel(F, "b", fig=figb, data=mods[["module", "share", "matched_baseline", "named"]],
             stats={"modules": len(mods), "named (q < 0.05)": int(mods.named.sum()),
                    "random draws per module": 1000,
                    "matched on": "how many kinds of knowledge the gene carries"},
             sources=msrc, what="share in its compartment vs matched random genes")
    PP.panel(F, "c", fig=figc, data=pd.DataFrame(C, index=lab, columns=HOURS),
             stats={"genes with a cell-cycle profile": len(z)},
             sources=msrc + [(NODES, [c for v in cols.values() for c in v])],
             what="each module across the cell cycle")
    PP.panel(F, "d", fig=figd, data=pd.DataFrame(D, index=lab, columns=PHASES),
             stats={"phase-annotated genes": int(ph.isin(PHASES).sum())},
             sources=msrc + [(NODES, ["cellcycle_phase"])],
             what="each module's genes by annotated phase")
    r_ = _sp(mods["coh"], mods["over_detection_limit"], nan_policy="omit")
    PP.panel(F, "e", fig=fige, data=mods[["module", "coh", "over_detection_limit", "named"]]
             .rename(columns={"coh": "share_peaking_at_modal_hour_pct"}),
             stats={"Spearman rho": r_.statistic, "p": r_.pvalue, "modules": len(mods)},
             sources=msrc, what="timing coherence vs compartment enrichment")
    fl = pd.DataFrame([{"module": m, "named": nm, "gene_id": g, "fit_invitro_hff": fit.get(g)}
                       for m, nm in zip(order, named) for g in genes[m]]).dropna()
    u = mannwhitneyu(fl[fl.named].fit_invitro_hff, fl[~fl.named].fit_invitro_hff,
                     alternative="less")
    PP.panel(F, "f", fig=figf, data=fl, stats={
        "median fitness, named-module genes": fl[fl.named].fit_invitro_hff.median(),
        "median fitness, other-module genes": fl[~fl.named].fit_invitro_hff.median(),
        "median fitness, all genes": float(fit.median()),
        "Mann-Whitney U (named < other)": u.statistic, "p": u.pvalue,
        "named modules below the genome median": "11 of 14",
        "other modules below the genome median": "0 of 9"},
             sources=msrc + [(NODES, ["fit_invitro_hff"])], what="fitness in fibroblasts per module")
    out = PP.compose(F, "main", "Co-expression modules are, for the most part, single organelles "
                     "built on a shared schedule, and the organelle modules are the ones the "
                     "parasite cannot lose.", rows=[3, 3])
    print("figure 3 ->", out)


if __name__ == "__main__":
    main()
