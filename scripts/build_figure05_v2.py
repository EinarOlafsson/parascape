"""Figure 5, in the order Einar set out on 2026-11-19.

    1 old clustering   2 new clustering (+ how new maps to old)   3 confounders
    4 defining traits  5 metabolic traits   6 drug treatments mapped to clusters
    7 enrichment       8 drug UMAPs: reached / n.d., candidates for n.d., autophagy-PLVAC

Most panels are the bradyzoite paper's own (report/figures_v2), reused as drawn because Einar
asked for them by name ("like the annotated figure you made", "like you did before"). Three are
new and drawn here: how the new regions map onto the published subtypes, a UMAP labelled with
the treatment that reaches each region (n.d. where none does), and the same UMAP with each n.d.
region given a candidate from its own enriched druggable genes.

Run with the toxo-sc env (reads the h5ad).
"""
from __future__ import annotations

import glob
import os
import shutil
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
V2 = f"{BRADY}/report/figures_v2"
OUT = os.path.join(ROOT, "results", "figure05_v2")
PNL = os.path.join(OUT, "panels")

ORDER = ["M5", "M8", "M3", "M4", "M7", "M1", "M2", "M6", "M9", "G10"]
SUBTYPES = list("ABCDE")
RC = dict(PF.REGION_COLOUR); RC.setdefault("G10", "#C8C8C8")

#: A candidate for a region no treatment reaches, ONLY where the inhibitor-target pair is
#: established. LHVS for cathepsin L (CPL) is Einar's and is standard in Toxoplasma. Anything
#: else is shown as the enriched target GENE, not as an invented drug.
KNOWN_INHIBITOR = {"CPL": "LHVS"}


def brady(rel):
    return os.path.join(V2, rel)


def new_to_old(obs):
    """How the new regions map onto the published subtypes: share of each region's cells."""
    t = pd.crosstab(obs["cluster"].astype(str), obs["published_subtype"].astype(str))
    t = t.reindex([r for r in ORDER if r in t.index])
    t = t[[s for s in SUBTYPES if s in t.columns]]
    share = t.div(t.sum(axis=1), axis=0)
    FS.use()
    nr, nc = share.shape
    cell = FS.SLOT_IN
    w, h = 0.8 + nc * cell + 0.9, 0.5 + nr * cell + 0.7
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([0.8 / w, 0.7 / h, nc * cell / w, nr * cell / h])
    im = ax.imshow(share.to_numpy(), cmap="Blues", vmin=0, vmax=1, aspect="equal")
    ax.set_xticks(range(nc)); ax.set_xticklabels(share.columns)
    ax.set_yticks(range(nr)); ax.set_yticklabels(share.index)
    ax.set_xlabel("Published subtype"); ax.set_ylabel("Region")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
    cb.set_label("Share of the region's cells"); cb.outline.set_linewidth(0.6)
    return fig


def reached():
    """Per region: the treatment reaching it hardest, and its signed effect on the markers."""
    d = pd.read_csv(f"{RES}/drug_by_cluster.tsv", sep="\t")
    d = d[(d["cluster_kind"] == "our region") & (d["reaches_cluster"].astype(str) == "True")]
    out = {}
    for reg, g in d.groupby("cluster"):
        b = g.reindex(g["signed_z"].abs().sort_values(ascending=False).index).iloc[0]
        out[reg] = (str(b["treatment_label"]).split(",")[0], float(b["signed_effect"]))
    return out


def candidates():
    """Per region: its top druggable genes, for regions no treatment reaches."""
    d = pd.read_csv(f"{RES}/druggability.tsv", sep="\t")
    d = d[d["q_bh"] < 0.05].sort_values("rank_in_region")
    # A symbol where there is one, else a short product name. Skipping symbol-less genes dropped
    # M3's most significant candidate (nucleoside diphosphate kinase) and printed 'nan'.
    def name(sym, prod):
        if isinstance(sym, str) and sym.strip() and sym.lower() != "nan":
            return sym
        prod = str(prod).replace(", putative", "").strip()
        return prod if prod and prod.lower() != "nan" else None
    out = {}
    for r, g in d.groupby("region_highest"):
        out[r] = [n for n in (name(sy, pr) for sy, pr in zip(g["symbol"], g["product"])) if n]
    return out


def drug_umap(obs, label_of, legend_items):
    FS.use()
    fig, ax = FS.square_figure(side_mm=58.0)
    order = np.random.default_rng(0).permutation(len(obs))
    o = obs.iloc[order]
    ax.scatter(o["ux"], o["uy"], c=[RC.get(str(c), FS.BACKGROUND) for c in o["cluster"]],
               s=0.5, lw=0, rasterized=True, alpha=0.35)
    for reg in ORDER:
        m = obs["cluster"].astype(str) == reg
        if not m.any():
            continue
        cx, cy = float(obs.loc[m, "ux"].median()), float(obs.loc[m, "uy"].median())
        text, colour, style = label_of(reg)
        ax.text(cx, cy, f"{reg}\n{text}", ha="center", va="center", color=colour,
                fontstyle=style, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
    for lab, col in legend_items:
        ax.plot([], [], "s", color=col, ms=5, label=lab, linestyle="none")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    return fig


def main():
    import anndata as ad
    a = ad.read_h5ad(H5AD, backed="r")
    obs = a.obs.copy()
    obs["ux"], obs["uy"] = a.obsm["X_umap"][:, 0], a.obsm["X_umap"][:, 1]
    os.makedirs(PNL, exist_ok=True)
    R, C = reached(), candidates()

    def save(fig, name):
        p = os.path.join(PNL, f"{name}.png")
        fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)
        return p

    # 8a: which treatment reaches each region, one number (signed marker effect); n.d. otherwise
    def lab_reached(reg):
        if reg in R:
            drug, eff = R[reg]
            return f"{drug}\n{eff:+.2f}", FS.INK, "normal"
        return "n.d.", FS.CONTROL, "normal"

    # 8b: the same, n.d. replaced by a candidate from the region's own enriched druggable genes
    def lab_candidate(reg):
        if reg in R:
            return lab_reached(reg)
        genes = C.get(reg, [])
        for g in genes:
            if g in KNOWN_INHIBITOR:
                return f"{KNOWN_INHIBITOR[g]}\n({g})", FS.ACCENT, "normal"
        return (f"{genes[0]}" if genes else "n.d."), FS.ACCENT, "italic"

    # The two drug maps use the bradyzoite paper's own leader-line map (pubfig.annotated_map):
    # labels outside the map, one line to each cluster. Labels placed at cluster centres piled
    # on top of each other wherever two clusters touch, which is most of them.
    PF.use("v150f")
    d = PF.load(with_expression=False)

    def text_reached(reg):
        if reg in R:
            drug, eff = R[reg]
            return f"{PF.with_group(reg)}: {drug} {eff:+.2f}"
        return f"{PF.with_group(reg)}: n.d."

    def text_candidate(reg):
        if reg in R:
            return text_reached(reg)
        genes = C.get(reg, [])
        for g in genes:
            if g in KNOWN_INHIBITOR:
                return f"{PF.with_group(reg)}: {KNOWN_INHIBITOR[g]} (targets {g})"
        return f"{PF.with_group(reg)}: target {genes[0]}" if genes else \
            f"{PF.with_group(reg)}: n.d."

    def leader_map(label_fn, w_mm=180.0, h_mm=50.0):
        """The paper leader-line map. 180 x 50 mm (was 200 x 62) so three fit an A4 page."""
        labels = {r: label_fn(r) for r in PF.REGIONS}
        fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4))
        g = fig.add_gridspec(1, 3, width_ratios=[1.35, 0.62, 1.35], wspace=0.03)
        PF.annotated_map(fig, g[0], g[1], g[2], d, labels, fontsize=FS.FONT_PT)
        return fig

    def old_map():
        """The published subtypes A-E on the UMAP, house style, no cell counts on the map."""
        FS.use()
        fig, ax = FS.square_figure(side_mm=40.0)
        o = obs.iloc[np.random.default_rng(0).permutation(len(obs))]
        sub = o["published_subtype"].astype(str)
        for i, k in enumerate(SUBTYPES):
            m = sub == k
            ax.scatter(o.loc[m, "ux"], o.loc[m, "uy"], s=0.5, lw=0, rasterized=True,
                       color=FS.SERIES[i], alpha=0.5, label=f"subtype {k}")
        m = ~sub.isin(SUBTYPES)
        ax.scatter(o.loc[m, "ux"], o.loc[m, "uy"], s=0.5, lw=0, rasterized=True,
                   color=FS.BACKGROUND, alpha=0.5, label="unassigned")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0), markerscale=8)
        return fig

    def score_map(col, label):
        """A per-cell score on the UMAP, house style, square."""
        FS.use()
        fig, ax = FS.square_figure(side_mm=40.0)
        o = obs.iloc[np.argsort(obs[col].to_numpy())]
        sc = ax.scatter(o["ux"], o["uy"], c=o[col], s=0.4, lw=0, rasterized=True, cmap="magma_r")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        cb = fig.colorbar(sc, ax=ax, fraction=0.06, pad=0.04)
        cb.set_label(label); cb.outline.set_linewidth(0.6)
        return fig

    def enrichment_map():
        """The four enrichments in one region x enrichment map, regions in the paper order.

        Einar: the enrichment plots in a fixed order. Four bar panels in four units became one map;
        each row is divided by its own largest absolute value, so colour compares regions within a
        row, never rows with each other. Grey = not measured.
        """
        dz = pd.read_csv(f"{RES}/differentiation_timecourse_mapping.tsv", sep="\t")
        cw = pd.read_csv(f"{RES}/cyst_wall_and_bradyzoite_sources.tsv", sep="\t")
        dr = pd.read_csv(f"{RES}/druggability.tsv", sep="\t")
        iv = pd.read_csv(f"{RES}/gse229548_invivo_screen_by_region.tsv", sep="\t")
        iv = iv[iv["gene_set"].astype(str).str.startswith("candidates highest in ")]
        iv = iv.assign(region=iv["gene_set"].str.replace("candidates highest in ", "", regex=False))
        rows = {
            "trend with differentiation (r)": dz.groupby("region")["trend_r"].mean(),
            "cyst-wall genes (z)": cw[cw["z"].notna()].groupby("region")["z"].mean(),
            "druggable candidates (n)": dr[dr["q_bh"] < 0.05]["region_highest"].astype(str)
                                        .value_counts(),
            "in vivo fitness cost (z)": iv.groupby("region")["z_vs_null"].mean(),
        }
        m = pd.DataFrame({k: v.reindex(ORDER) for k, v in rows.items()}).T
        if "druggable candidates (n)" in m.index:
            m.loc["druggable candidates (n)"] = m.loc["druggable candidates (n)"].fillna(0)
        scaled = m.div(m.abs().max(axis=1), axis=0)
        FS.use()
        nr, nc = scaled.shape; cell = 0.27   # 0.2 in let "G10" run into "M9"
        w, h = 1.75 + nc * cell + 0.9, 0.15 + nr * cell + 0.6
        fig = plt.figure(figsize=(w, h))
        ax = fig.add_axes([1.75 / w, 0.6 / h, nc * cell / w, nr * cell / h])
        cmap = matplotlib.colormaps["RdBu_r"].copy(); cmap.set_bad(FS.BACKGROUND)
        im = ax.imshow(np.ma.masked_invalid(scaled.to_numpy(float)), cmap=cmap, vmin=-1, vmax=1,
                       interpolation="nearest")
        ax.set_xticks(range(nc)); ax.set_xticklabels(ORDER)
        ax.set_yticks(range(nr)); ax.set_yticklabels(list(scaled.index))
        for t, r in zip(ax.get_xticklabels(), ORDER):
            t.set_color(RC.get(r, FS.INK))
        ax.set_xlabel("Region")
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(length=0)
        cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
        cb.set_label("Value / row maximum"); cb.outline.set_linewidth(0.6)
        print(m.round(2).to_string())
        return fig

    built = {
        "1_old_map": save(old_map(), "1_old_map"),
        "2_new_map": save(leader_map(lambda r: PF.with_group(r), w_mm=110.0, h_mm=50.0),
                          "2_new_map"),
        "3_g1a": save(score_map("cc_G1a", "G1a score"), "3_g1a"),
        "3_sm": save(score_map("cc_SM", "S/M score"), "3_sm"),
        "7_enrichment": save(enrichment_map(), "7_enrichment"),
        "q_autophagy": save(leader_map(lambda r: PF.autophagy_map_labels(d)[r]), "q_autophagy"),
        "2b_new_to_old": save(new_to_old(obs), "2b_new_to_old"),
        "8a_drug_reached": save(leader_map(text_reached), "8a_drug_reached"),
        "8b_drug_candidates": save(leader_map(text_candidate), "8b_drug_candidates"),
    }
    for k, v in built.items():
        print("  built", k)

    # Einar, 2026-11-19: split into figure 5 (the map, its confounders, what reaches each region)
    # and supplementary figure S5 (what defines each region and its enrichments), each in his
    # order, three panels to a row.
    # A4 at true size (Einar, 2026-11-19) does not hold 17 panels: three pages, in his order.
    # 5: the map and what reaches each region; S5: the confounders; S6: what defines regions.
    main_seq = [
        ("old clustering", built["1_old_map"]),
        ("new clustering", built["2_new_map"]),
        ("new to old", built["2b_new_to_old"]),
        ("drugs to clusters", brady("Figure8_drug_by_cluster/panels/c_treatment_by_our_region_signed/c_treatment_by_our_region_signed.png")),
        ("drug reached", built["8a_drug_reached"]),
        ("drug candidates", built["8b_drug_candidates"]),
        ("autophagy PLVAC", built["q_autophagy"]),
    ]
    s5_seq = [
        ("cell cycle G1a", built["3_g1a"]),
        ("cell cycle S/M", built["3_sm"]),
        ("recrudescence", brady("Figure4_reactivation/panels/a_reactivation_bar_jitter/a_reactivation_bar_jitter.png")),
        ("tachy vs brady", brady("Figure6_drugs_and_perturbations/panels/sa_perturbation_map_waldman_wt_alkaline/sa_perturbation_map_waldman_wt_alkaline.png")),
    ]
    s6_seq = [
        ("defining traits", brady("Figure5_what_distinguishes_regions/panels/b_region_property_enrichment/b_region_property_enrichment.png")),
        ("metabolic traits", brady("Figure2_metabolism/panels/b_programme_by_region/b_programme_by_region.png")),
        ("enrichment", built["7_enrichment"]),
    ]
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import assemble_paper_figures as A
    A.PER_ROW = 4                  # compose() then splits any row wider than A4
    letters = "abcdefghijklmnopqrstuvwxyz"
    for fid, seq, rows in (("05", main_seq, None), ("S5", s5_seq, None), ("S6", s6_seq, None)):
        missing = [n for n, p in seq if not os.path.exists(p)]
        if missing:
            print("  MISSING:", missing)
        present = [(n, p) for n, p in seq if os.path.exists(p)]
        # Panels reused from the bradyzoite paper are drawn at its 6 pt; enlarged to this
        # paper's 8.
        panels = [{"png": p, "letter": letters[i],
                   "scale": FS.FONT_PT / 6.0 if p.startswith(V2) else 1.0}
                  for i, (n, p) in enumerate(present)]
        out = A.compose(fid, panels, OUT, rows=rows)
        print(f"  figure {fid}: {len(panels)} panels -> {out}")


if __name__ == "__main__":
    main()
