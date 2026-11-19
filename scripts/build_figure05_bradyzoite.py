"""Figure 5: a tissue cyst holds several bradyzoite states, and published treatments reach them
unequally.

The bradyzoite single-cell analysis (v150f: 15,177 cells, ten regions M1-M9 and G10, published
subtypes A-E) mapped against the bulk treatment corpus. Twelve panels, all from data, per the spec
in `paper_figure_spec.py`.

Run with the toxo-sc environment, which reads the h5ad:
    /home/carruthers/anaconda3/envs/toxo-sc/bin/python scripts/build_figure05_bradyzoite.py

Identity colours are imported from the bradyzoite paper (`pubfig.REGION_COLOUR`,
`SUBTYPE_COLOUR`), because the skill requires one colour per identity across every figure, and the
two papers describe the same regions.
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
BRADY = "/media/carruthers/mnt3/claude/brady_metabolic_classes"
sys.path.insert(0, os.path.join(BRADY, "notebooks"))

from parascape import figstyle as FS            # noqa: E402
import pubfig as PF                             # noqa: E402

H5AD = f"{BRADY}/data/processed/v150f_clustered.h5ad"
RES = f"{BRADY}/results_v150f"
OUT = os.path.join(ROOT, "results", "figure05_bradyzoite")

#: Region order on every axis: grouped by the published subtype each region overlaps most.
ORDER = ["M5", "M8", "M3", "M4", "M7", "M1", "M2", "M6", "M9", "G10"]
REGION_COLOUR = dict(PF.REGION_COLOUR)
REGION_COLOUR.setdefault("G10", "#C8C8C8")
SUBTYPE_COLOUR = dict(PF.SUBTYPE_COLOUR)
TARGETS = {"CDPK1": "TGME49_301440", "PKG": "TGME49_311360", "DHFR-TS": "TGME49_249180"}


def newest(pattern: str) -> str:
    found = sorted(glob.glob(os.path.join(ROOT, "results", pattern)))
    if not found:
        raise SystemExit(f"missing: {pattern}")
    return found[-1]


def save(fig, name: str) -> None:
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=600 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("  wrote", name)


def load():
    import anndata as ad
    a = ad.read_h5ad(H5AD)
    um = a.obsm["X_umap"]
    obs = a.obs.copy()
    obs["ux"], obs["uy"] = um[:, 0], um[:, 1]
    return a, obs


def umap_panel(obs, colour_of, legend_items, letter, cmap_values=None, cmap_label=""):
    FS.use()
    fig, ax = FS.square_figure(side_mm=52.0)
    order = np.random.default_rng(0).permutation(len(obs))
    o = obs.iloc[order]
    if cmap_values is not None:
        v = cmap_values.iloc[order] if hasattr(cmap_values, "iloc") else cmap_values[order]
        sc = ax.scatter(o["ux"], o["uy"], c=v, cmap="Blues", s=0.6, lw=0, rasterized=True,
                        alpha=FS.ALPHA)
        cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
        cb.set_label(cmap_label)
        cb.outline.set_linewidth(0.6)
    else:
        ax.scatter(o["ux"], o["uy"], c=[colour_of(r) for r in o.itertuples()], s=0.6, lw=0,
                   rasterized=True, alpha=FS.ALPHA)
        for lab, col in legend_items:
            ax.plot([], [], "o", color=col, markersize=4, label=lab, linestyle="none")
        # OUTSIDE the axes, to the right: inside, it sat on top of the cells (skill: a legend
        # never goes over data).
        FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0), ncol=1, markerscale=1.0,
                  handletextpad=0.2)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    FS.letter(fig, letter)
    return fig


def panel_a(obs):
    items = [(r, REGION_COLOUR[r]) for r in ORDER if r in set(obs["cluster"].astype(str))]
    return umap_panel(obs, lambda r: REGION_COLOUR.get(str(r.cluster), FS.BACKGROUND), items, "a")


def panel_b(obs):
    items = [(k, c) for k, c in SUBTYPE_COLOUR.items()] + [("unassigned", FS.BACKGROUND)]
    return umap_panel(obs, lambda r: SUBTYPE_COLOUR.get(str(r.published_subtype), FS.BACKGROUND),
                      items, "b")


def panel_c(obs):
    return umap_panel(obs, None, [], "c", cmap_values=obs["cc_SM"].astype(float),
                      cmap_label="S/M score")


def panel_d():
    """Marker genes per region as a dot plot: size = share of cells, colour = log2 fold change."""
    m = pd.read_csv(f"{RES}/cluster_markers.tsv", sep="\t")
    m = m[m["cluster"].astype(str).isin(ORDER)]
    top = (m.sort_values("log2_fold_change", ascending=False)
             .groupby("cluster").head(2))
    genes = list(dict.fromkeys(top.sort_values(
        "cluster", key=lambda s: s.map({r: i for i, r in enumerate(ORDER)}))["gene_id"]))
    regions = [r for r in ORDER if r in set(m["cluster"])]
    FS.use()
    fig = plt.figure(figsize=(FS.MARGIN_MM[0] / 25.4 + len(genes) * 0.20 + 0.7,
                              FS.MARGIN_MM[1] / 25.4 + len(regions) * 0.20 + 0.9))
    ax = fig.add_axes([0.24, 0.30, 0.62, 0.62])
    for gi, g in enumerate(genes):
        for ri, r in enumerate(regions):
            row = m[(m["gene_id"] == g) & (m["cluster"] == r)]
            if row.empty:
                continue
            size = float(row["frac_cells_in_cluster"].iloc[0]) * 60
            lfc = float(row["log2_fold_change"].iloc[0])
            ax.scatter(gi, ri, s=size, c=[[0.04, 0.47, 0.72, min(1.0, 0.25 + lfc / 4)]], lw=0)
    ax.set_xticks(range(len(genes)))
    # A product name, shortened; a bare accession number tells the reader nothing.
    names = (m.drop_duplicates("gene_id").set_index("gene_id")["product"].astype(str)
               .str.replace(", putative", "", regex=False).str.slice(0, 22))
    ax.set_xticklabels([names.get(g, g.replace("TGME49_", "")) for g in genes], rotation=90)
    ax.set_yticks(range(len(regions))); ax.set_yticklabels(regions)
    ax.set_xlim(-0.6, len(genes) - 0.4); ax.set_ylim(-0.6, len(regions) - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("Marker gene"); ax.set_ylabel("Region")
    FS.letter(fig, "d")
    return fig


def panel_e():
    c = pd.read_csv(newest("bradyzoite_reach_*/contrasts.tsv"), sep="\t")
    subs = [s for s in "ABCDE" if f"q_{s}" in c.columns]
    counts = [int((c[f"q_{s}"] < 0.05).sum()) for s in subs]
    FS.use()
    fig, ax = FS.figure(len(subs), vertical=True)
    x = range(len(subs))
    ax.bar(x, counts, width=FS.bar_width(), lw=0, alpha=FS.ALPHA,
           color=[SUBTYPE_COLOUR[s] for s in subs])
    ax.set_xticks(list(x)); ax.set_xticklabels(subs)
    ax.set_xlabel("Published subtype"); ax.set_ylabel("Contrasts reaching it")
    ax.set_xlim(-0.6, len(subs) - 0.4); ax.set_ylim(0, max(counts) * 1.12)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "e")
    return fig


def panel_f():
    c = pd.read_csv(newest("bradyzoite_reach_*/contrasts.tsv"), sep="\t")
    subs = [s for s in "ABCDE" if f"q_{s}" in c.columns]
    FS.use()
    fig, ax = FS.figure(len(subs), vertical=True)
    groups = [c[f"q_{s}"].dropna().to_numpy() for s in subs]
    parts = ax.violinplot(groups, positions=list(range(len(subs))), widths=FS.bar_width(),
                          showextrema=False)
    for body, s in zip(parts["bodies"], subs):
        body.set_facecolor(SUBTYPE_COLOUR[s]); body.set_alpha(FS.ALPHA); body.set_linewidth(0)
    for i, g in enumerate(groups):
        ax.plot([i - 0.2, i + 0.2], [np.median(g)] * 2, color=FS.INK, lw=0.9)
    ax.axhline(0.05, color=FS.INK, lw=0.6, ls=(0, (2, 2)), label="q = 0.05")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.set_xticks(range(len(subs))); ax.set_xticklabels(subs)
    ax.set_xlabel("Published subtype"); ax.set_ylabel("q")
    ax.set_xlim(-0.6, len(subs) - 0.4); ax.set_ylim(-0.04, 1.06)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "f")
    return fig


def panel_g():
    """Treatment by region: signed z of the treatment effect on each region's markers."""
    d = pd.read_csv(f"{RES}/drug_by_cluster.tsv", sep="\t")
    d = d[d["cluster_kind"] == "our region"]
    pivot = d.pivot_table(index="treatment_label", columns="cluster", values="signed_z")
    pivot = pivot[[r for r in ORDER if r in pivot.columns]].dropna(how="all")
    FS.use()
    nrow, ncol = pivot.shape
    cell = 0.17
    fig = plt.figure(figsize=(2.6 + ncol * cell + 0.6, 0.9 + nrow * cell))
    ax = fig.add_axes([2.4 / (2.6 + ncol * cell + 0.6), 0.6 / (0.9 + nrow * cell),
                       ncol * cell / (2.6 + ncol * cell + 0.6), nrow * cell / (0.9 + nrow * cell)])
    lim = float(np.nanmax(np.abs(pivot.to_numpy()))) or 1.0
    im = ax.imshow(pivot.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="equal")
    for (i, j), v in np.ndenumerate(pivot.to_numpy()):
        if np.isnan(v):
            ax.text(j, i, "×", ha="center", va="center", fontsize=FS.FONT_PT, color=FS.INK)
    ax.set_xticks(range(ncol)); ax.set_xticklabels(pivot.columns, rotation=90)
    ax.set_yticks(range(nrow)); ax.set_yticklabels(pivot.index)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cb.set_label("Signed z"); cb.outline.set_linewidth(0.6)
    FS.letter(fig, "g")
    return fig


def panel_h(a, obs):
    """Drug-target expression by region: CDPK1, the target with the widest regional range."""
    gid = TARGETS["CDPK1"]
    if gid not in a.var_names:
        raise SystemExit(f"{gid} not in the h5ad")
    x = a[:, gid].X
    x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    regions = [r for r in ORDER if r in set(obs["cluster"].astype(str))]
    groups = [x[(obs["cluster"].astype(str) == r).to_numpy()] for r in regions]
    FS.use()
    fig, ax = FS.figure(len(regions), vertical=True)
    parts = ax.violinplot(groups, positions=list(range(len(regions))), widths=FS.bar_width(),
                          showextrema=False)
    for body, r in zip(parts["bodies"], regions):
        body.set_facecolor(REGION_COLOUR[r]); body.set_alpha(FS.ALPHA); body.set_linewidth(0)
    for i, g in enumerate(groups):
        ax.plot([i - 0.2, i + 0.2], [np.median(g)] * 2, color=FS.INK, lw=0.9)
    ax.set_xticks(range(len(regions))); ax.set_xticklabels(regions)
    ax.set_xlabel("Region"); ax.set_ylabel("CDPK1 expression (log1p)")
    ax.set_xlim(-0.6, len(regions) - 0.4)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "h")
    return fig


def bar_by_region(values: dict, ylabel: str, letter: str, reference=None, ref_label=""):
    regions = list(ORDER)
    FS.use()
    fig, ax = FS.figure(len(regions), vertical=True)
    ax.bar(range(len(regions)), [values.get(r, 0) or 0 for r in regions], width=FS.bar_width(),
           lw=0,
           alpha=FS.ALPHA, color=[REGION_COLOUR[r] for r in regions])
    if reference is not None:
        ax.axhline(reference, color=FS.INK, lw=0.6, ls=(0, (2, 2)), label=ref_label)
        FS.legend(ax, loc="best")
    ax.axhline(0, color=FS.INK, lw=0.6)
    # A region this measurement does not cover is marked x on the axis, never left blank:
    # a blank reads as zero, and zero is a value.
    missing = [i for i, r in enumerate(regions) if values.get(r) is None or
               (isinstance(values.get(r), float) and np.isnan(values.get(r)))]
    for i in missing:
        ax.plot(i, 0, marker="x", color=FS.INK, markersize=4, markeredgewidth=0.8,
                linestyle="none", zorder=6, label="not measured" if i == missing[0] else None)
    if missing:
        FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.set_xticks(range(len(regions))); ax.set_xticklabels(regions)
    ax.set_xlabel("Region"); ax.set_ylabel(ylabel)
    ax.set_xlim(-0.6, len(regions) - 0.4)
    vals = [values[r] for r in regions if values.get(r) is not None] + \
        ([reference] if reference is not None else []) + [0]
    FS.pad_limits(ax, vals)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, letter)
    return fig


def panel_i():
    d = pd.read_csv(f"{RES}/druggability.tsv", sep="\t")
    sig = d[d["q_bh"] < 0.05]
    counts = sig["region_highest"].astype(str).value_counts().to_dict()
    return bar_by_region({r: counts.get(r, 0) for r in ORDER}, "Druggable candidates", "i")


def panel_j():
    d = pd.read_csv(f"{RES}/gse229548_invivo_screen_by_region.tsv", sep="\t")
    d = d[d["gene_set"].astype(str).str.startswith("candidates highest in ")]
    d = d.assign(region=d["gene_set"].str.replace("candidates highest in ", "", regex=False))
    d = d.groupby("region")["z_vs_null"].mean().to_dict()
    return bar_by_region(d, "In vivo fitness cost (z vs null)", "j")


def panel_k():
    d = pd.read_csv(f"{RES}/cyst_wall_and_bradyzoite_sources.tsv", sep="\t")
    d = d[d["region"].astype(str).isin(ORDER) & d["z"].notna()]
    vals = d.groupby("region")["z"].mean().to_dict()
    return bar_by_region(vals, "Cyst-wall enrichment (z vs null)", "k")


def panel_l():
    d = pd.read_csv(f"{RES}/differentiation_timecourse_mapping.tsv", sep="\t")
    d = d[d["region"].astype(str).isin(ORDER)]
    vals = d.groupby("region")["trend_r"].mean().to_dict()
    return bar_by_region(vals, "Trend with differentiation (r)", "l")


def main() -> None:
    a, obs = load()
    print(f"loaded {a.shape[0]:,} cells")
    for name, fn in (("a_region_map", lambda: panel_a(obs)),
                     ("b_subtype_map", lambda: panel_b(obs)),
                     ("c_cellcycle_map", lambda: panel_c(obs)),
                     ("d_markers", panel_d), ("e_reach_counts", panel_e),
                     ("f_reach_q", panel_f), ("g_treatment_by_region", panel_g),
                     ("h_cdpk1_by_region", lambda: panel_h(a, obs)),
                     ("i_druggability", panel_i), ("j_in_vivo_cost", panel_j),
                     ("k_cyst_wall", panel_k), ("l_differentiation", panel_l)):
        try:
            save(fn(), name)
        except Exception as e:
            print(f"  FAILED {name}: {type(e).__name__}: {e}")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
