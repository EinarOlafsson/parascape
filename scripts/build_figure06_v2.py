"""Figure 6, one claim: organelles built in waves during division can be placed from transcription
alone, and those made continuously cannot -- because the wave is what transcription records.

Einar, 2026-11-19: "i dont see the point of this figure ... allot of the place holder data here
you can already make graphs from data we have". The old figure 6 mixed two claims -- localisation
from RNA, and what sets protein abundance -- which instruction 101 forbids (every panel the same
shape of claim). This keeps the localisation claim and SHOWS ITS MECHANISM from the synchronised
cell-cycle time course (cellcycle19092, 1-12 h, 6,504 genes) and hyperLOPIT (lopit_unified).

  A  placed from RNA alone vs from all evidence, per compartment        (the result)
  B  mean expression across the cell cycle, wave-built vs continuous    (the wave)
  C  share of genes peaking at each hour, per compartment (map)        (how tight the wave is)
  D  timing coherence against placeability from RNA                     (mechanism -> result)
  E  single genes: a wave-built and a continuous compartment            (what a reader sees)
  F  replicate agreement, median and IQR per compartment (named rows) (the wave is signal)
All six from data; no placeholders.
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
from parascape import figstyle as FS        # noqa: E402

OUT = os.path.join(ROOT, "results", "figure06_v2")
PNL = os.path.join(OUT, "panels")
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
HOURS = list(range(1, 13))

#: The fine hyperLOPIT labels (lopit_map) grouped into the niches panel A scores, by EXACT
#: name. Substring matching on the coarse lopit_unified put apicoplast in the wave group and
#: found no rhoptries, micronemes or plasma membrane at all, which emptied panel E.
NICHE = {
    "rhoptries": ["rhoptries 1", "rhoptries 2"], "micronemes": ["micronemes"],
    "IMC": ["IMC"], "apical": ["apical 1", "apical 2"],
    "ribosome": ["40S ribosome", "60S ribosome"],
    "plasma membrane": ["PM - integral", "PM - peripheral 1", "PM - peripheral 2"],
    "ER": ["ER", "ER 2"], "mitochondrion": ["mitochondrion - soluble", "mitochondrion - membranes"],
    "nucleus": ["nucleus - chromatin", "nucleus - non-chromatin"], "cytosol": ["cytosol"],
    "apicoplast": ["apicoplast"], "dense granules": ["dense granules"],
}
#: From panel A's own result: placed as well or better by RNA alone (wave-built) vs not.
WAVE = ["rhoptries", "micronemes", "IMC", "apical", "ribosome"]
CONTINUOUS = ["plasma membrane", "ER", "mitochondrion", "nucleus", "cytosol"]


def newest(pat):
    f = sorted(glob.glob(os.path.join(ROOT, "results", pat)))
    return f[-1] if f else None


def profiles():
    n = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    cols = {h: [f"cellcycle19092_{h}h_r1", f"cellcycle19092_{h}h_r2"] for h in HOURS}
    m = pd.DataFrame({h: n[c].mean(axis=1) for h, c in cols.items()})
    m.index = n["gene_id"].astype(str)
    m = m.dropna()
    z = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0).dropna()
    comp = n.set_index(n["gene_id"].astype(str))["lopit_map"]
    z = z.join(comp.rename("comp")).dropna(subset=["comp"])
    z["comp"] = z["comp"].astype(str)
    # Replicate agreement per gene: Pearson r between the r1 and r2 12-hour profiles. The z range
    # used before is near-constant by construction (every gene is z-scored), so it said nothing.
    r1 = pd.DataFrame({h: n[f"cellcycle19092_{h}h_r1"] for h in HOURS}).set_axis(n["gene_id"].astype(str), axis=0)
    r2 = pd.DataFrame({h: n[f"cellcycle19092_{h}h_r2"] for h in HOURS}).set_axis(r1.index, axis=0)
    a1, a2 = r1.sub(r1.mean(axis=1), axis=0), r2.sub(r2.mean(axis=1), axis=0)
    rep = (a1 * a2).sum(axis=1) / np.sqrt((a1 ** 2).sum(axis=1) * (a2 ** 2).sum(axis=1))
    z["rep_r"] = rep.reindex(z.index)
    return z


def group(z, key):
    return z[z["comp"].isin(NICHE.get(key, [key]))][HOURS]


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def panel_b(z):
    FS.use()
    fig, ax = FS.square_figure()
    for key in CONTINUOUS:
        g = group(z, key)
        if len(g) >= 5:
            ax.plot(HOURS, g.mean(), color=FS.CONTROL, lw=1.0, alpha=FS.ALPHA)
    for i, key in enumerate(WAVE):
        g = group(z, key)
        if len(g) >= 5:
            ax.plot(HOURS, g.mean(), color=FS.SERIES[i % len(FS.SERIES)], lw=1.0,
                    alpha=FS.ALPHA, label=key)
    ax.plot([], [], color=FS.CONTROL, lw=1.0, label="other niches")
    ax.axhline(0, color=FS.INK, lw=0.6)
    ax.set_xlabel("Time after release (h)"); ax.set_ylabel("Mean expression (z)")
    ax.set_xticks([2, 4, 6, 8, 10, 12])
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    rows = []
    for key in WAVE + CONTINUOUS:
        g = group(z, key)
        if len(g) >= 5:
            for h, v in g.mean().items():
                rows.append({"compartment": key, "group": "S/M subtranscriptome" if key in WAVE
                             else "continuous", "hour": h, "mean_z": v, "genes": len(g)})
    fig._data = pd.DataFrame(rows)
    fig._stats = {f"genes, {k}": len(group(z, k)) for k in WAVE + CONTINUOUS}
    return fig


def peak_hours(z, key):
    g = group(z, key)
    return np.array(HOURS)[np.argmax(g.to_numpy(), axis=1)] if len(g) else np.array([])


def panel_c(z):
    """Share of each compartment's genes peaking at each hour: a map, compartments as rows.

    Ten violins in a row came out wide with tilted names (Einar, 2026-11-19: wide graphs need
    another form). The map reads names flat and shows the wave as a bright column.
    """
    keys = [k for k in WAVE + CONTINUOUS if len(group(z, k)) >= 5]
    m = np.array([np.bincount(peak_hours(z, k), minlength=13)[1:13] / len(group(z, k)) * 100
                  for k in keys])
    FS.use()
    fig, ax = FS.square_figure()
    im = ax.imshow(m, cmap="Blues", vmin=0, aspect="auto", interpolation="nearest",
                   extent=(0.5, 12.5, len(keys) - 0.5, -0.5))
    ax.set_yticks(range(len(keys))); ax.set_yticklabels(keys)
    for t, k in zip(ax.get_yticklabels(), keys):
        t.set_color(FS.ACCENT if k in WAVE else FS.INK)
    ax.axhline(len([k for k in keys if k in WAVE]) - 0.5, color=FS.INK, lw=0.6)
    ax.set_xticks([2, 4, 6, 8, 10, 12]); ax.set_xlabel("Time of peak expression (h)")
    ax.tick_params(axis="y", length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    cb.set_label("Genes peaking (%)"); cb.outline.set_linewidth(0.6)
    fig._data = pd.DataFrame(m, index=keys, columns=HOURS).rename_axis("compartment")
    fig._stats = pd.DataFrame([{"compartment": k, "genes": len(group(z, k)),
                                "modal_hour": int(np.argmax(r)) + 1,
                                "share_at_modal_hour_pct": float(r.max())}
                               for k, r in zip(keys, m)])
    return fig


def coherence(z, key):
    """How tightly a compartment's genes peak together: the share at the modal hour."""
    h = peak_hours(z, key)
    if len(h) < 5:
        return np.nan
    return float(np.bincount(h).max() / len(h))


def panel_d(z):
    s = pd.read_csv(newest("coexpression_components_localization_*/scores.tsv"), sep="\t")
    rows = []
    for _, r in s.iterrows():
        key = str(r["niche"])
        c = coherence(z, key)
        if not np.isnan(c):
            rows.append((key, c, float(r["coexp_cv_auc"])))
    d = pd.DataFrame(rows, columns=["comp", "coh", "auc"])
    FS.use()
    # Twelve labelled compartments do not fit a 46 mm square without labels landing between two
    # points; larger, so every label sits beside its own point.
    fig, ax = FS.square_figure()
    wave = d["comp"].isin(WAVE)
    ax.plot(d.loc[~wave, "coh"] * 100, d.loc[~wave, "auc"], "o", color=FS.CONTROL, ms=4,
            alpha=FS.ALPHA, lw=0, label="other niches")
    ax.plot(d.loc[wave, "coh"] * 100, d.loc[wave, "auc"], "o", color=FS.ACCENT, ms=4,
            alpha=FS.ALPHA, lw=0, label="S/M subtranscriptome")
    FS.pad_limits(ax, d["coh"] * 100, axis="x", frac=0.3)
    FS.pad_limits(ax, d["auc"], axis="y", frac=0.12)
    ax.set_xlabel("Genes peaking at the modal hour (%)")
    ax.set_ylabel("AUC, RNA only")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02))
    from scipy.stats import spearmanr
    r = spearmanr(d["coh"], d["auc"])
    fig._data = d.rename(columns={"coh": "share_peaking_at_modal_hour", "auc": "auc_rna_only"})
    fig._stats = {"Spearman rho, coherence vs AUC": r.statistic, "p": r.pvalue,
                  "compartments": len(d)}
    return fig


def panel_e(z):
    FS.use()
    fig, ax = FS.square_figure()
    for key, col in (("plasma membrane", FS.CONTROL), ("rhoptries", FS.ACCENT)):
        g = group(z, key)
        for _, row in g.head(60).iterrows():
            ax.plot(HOURS, row.to_numpy(), color=col, lw=0.4, alpha=0.25)
        ax.plot(HOURS, g.mean(), color=col, lw=1.4, alpha=1.0,
                label=key)
    ax.set_xlabel("Time (h)"); ax.set_ylabel("mRNA level (z-score)")
    ax.set_xticks([2, 4, 6, 8, 10, 12])
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02))
    rows = []
    for key in ("plasma membrane", "rhoptries"):
        g = group(z, key)
        for gid, row in g.head(60).iterrows():
            rows += [{"compartment": key, "gene_id": gid, "hour": h, "z": v}
                     for h, v in row.items()]
    fig._data = pd.DataFrame(rows)
    fig._stats = {"genes drawn per compartment": 60,
                  "genes, rhoptries": len(group(z, "rhoptries")),
                  "genes, plasma membrane": len(group(z, "plasma membrane"))}
    return fig


def panel_f(z):
    """Replicate agreement per gene, per compartment: violins, wave-built in the accent."""
    keys = [k for k in WAVE + CONTINUOUS if len(group(z, k)) >= 5]
    FS.use()
    short = {"plasma membrane": "plasma memb.", "mitochondrion": "mito."}
    fig, ax = FS.dist_bars([z.loc[group(z, k).index, "rep_r"].dropna().to_numpy() for k in keys],
                           [short.get(k, k) for k in keys], highlight=WAVE,
                           ylabel="Replicate agreement (Pearson r)")
    ax.set_ylim(-1.05, 1.05)
    from scipy.stats import mannwhitneyu
    rows = [{"compartment": k, "group": "S/M subtranscriptome" if k in WAVE else "continuous",
             "gene_id": gid, "replicate_r": v}
            for k in keys for gid, v in z.loc[group(z, k).index, "rep_r"].dropna().items()]
    fig._data = pd.DataFrame(rows)
    w = fig._data[fig._data.group == "S/M subtranscriptome"]["replicate_r"]
    c = fig._data[fig._data.group == "continuous"]["replicate_r"]
    u = mannwhitneyu(w, c, alternative="greater")
    fig._stats = {"median r, wave-built": w.median(), "median r, continuous": c.median(),
                  "Mann-Whitney U (wave > continuous)": u.statistic, "p": u.pvalue,
                  "genes, wave-built": len(w), "genes, continuous": len(c)}
    return fig


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    z = profiles()
    print(f"{len(z):,} genes with a cell-cycle profile and a compartment")
    import paper_panels as PP
    F = "figure_5__final"
    a = newest("coexpression_components_localization_*/panel_which_evidence_places_a_protein.png")
    sc = newest("coexpression_components_localization_*/scores.tsv")
    sd = pd.read_csv(sc, sep="\t")
    from scipy.stats import wilcoxon
    pair = sd[["coexp_cv_auc"] + [c for c in sd.columns if c != "coexp_cv_auc"
                                   and "auc" in c][:1]].dropna()
    st = {"compartments": len(sd)}
    if pair.shape[1] == 2:
        w = wilcoxon(pair.iloc[:, 0], pair.iloc[:, 1])
        st.update({f"Wilcoxon, {pair.columns[0]} vs {pair.columns[1]}": w.statistic,
                   "p": w.pvalue})
    nodes = os.path.join(CACHE, "nodes.parquet")
    cc_cols = [f"cellcycle19092_{h}h_r{r}" for h in HOURS for r in (1, 2)] + ["lopit_map"]
    # Storyboard 2026-11-20: A when each niche is transcribed (heatmap, grouped by
    # subtranscriptome); B what that looks like gene by gene; C RNA alone places the S/M niches as
    # well as all data; D why: the tighter the timing, the better RNA places the niche.
    for letter, fn, what in (("a", panel_c, "hour of peak expression per hyperLOPIT niche"),
                             ("b", panel_e, "single-gene profiles, rhoptries vs plasma membrane")):
        f = fn(z)
        PP.panel(F, letter, fig=f, data=f._data, stats=f._stats, sources=[(nodes, cc_cols)],
                 what=what)
    FS.use()
    figc, ax = FS.square_figure()
    sm = sd["niche"].isin(WAVE)
    ax.plot([0.5, 1], [0.5, 1], color=FS.BACKGROUND, lw=0.6, label="RNA = all data")
    ax.scatter(sd.loc[~sm, "coexp_cv_auc"], sd.loc[~sm, "all_evidence_cv_auc"], s=20,
               color=FS.CONTROL, alpha=0.8, lw=0, label="other niches")
    ax.scatter(sd.loc[sm, "coexp_cv_auc"], sd.loc[sm, "all_evidence_cv_auc"], s=20,
               color=FS.ACCENT, alpha=0.8, lw=0, label="S/M subtranscriptome")
    ax.set_xlim(0.6, 1.02); ax.set_ylim(0.6, 1.02)
    ax.set_xlabel("AUC, RNA only"); ax.set_ylabel("AUC, all data")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02))
    PP.panel(F, "c", fig=figc, data=sd, stats=st, sources=[sc],
             what="placement from RNA alone vs all data, per niche")
    f = panel_d(z)
    PP.panel(F, "d", fig=f, data=f._data, stats=f._stats, sources=[(nodes, cc_cols), sc],
             what="timing coherence vs placement from RNA alone")
    out = PP.compose(F, "main", "Proteins of organelles transcribed in the S/M subtranscriptome "
                     "are placed from RNA alone.")
    print("figure 6 ->", out)


if __name__ == "__main__":
    main()
