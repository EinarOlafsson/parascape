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
    ax.plot([], [], color=FS.CONTROL, lw=1.0, label="continuous (ref.)")
    ax.axhline(0, color=FS.INK, lw=0.6)
    ax.set_xlabel("Time after release (h)"); ax.set_ylabel("Mean expression (z)")
    ax.set_xticks([2, 4, 6, 8, 10, 12])
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
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
    ax.set_xticks([2, 4, 6, 8, 10, 12]); ax.set_xlabel("Hour of peak expression")
    ax.tick_params(axis="y", length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    cb.set_label("Genes peaking (%)"); cb.outline.set_linewidth(0.6)
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
    fig, ax = FS.square_figure(side_mm=66.0)
    wave = d["comp"].isin(WAVE)
    ax.plot(d.loc[~wave, "coh"] * 100, d.loc[~wave, "auc"], "o", color=FS.CONTROL, ms=4,
            alpha=FS.ALPHA, lw=0, label="continuous (ref.)")
    ax.plot(d.loc[wave, "coh"] * 100, d.loc[wave, "auc"], "o", color=FS.ACCENT, ms=4,
            alpha=FS.ALPHA, lw=0, label="wave-built")
    FS.pad_limits(ax, d["coh"] * 100, axis="x", frac=0.3)
    FS.pad_limits(ax, d["auc"], axis="y", frac=0.12)
    FS.place_labels(ax, fig, list(zip(d["coh"] * 100, d["auc"])), list(d["comp"]),
                    [FS.ACCENT if w else FS.CONTROL for w in wave])
    ax.set_xlabel("Genes peaking at the modal hour (%)")
    ax.set_ylabel("AUC, RNA only")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return fig


def panel_e(z):
    FS.use()
    fig, ax = FS.square_figure()
    for key, col in (("plasma membrane", FS.CONTROL), ("rhoptries", FS.ACCENT)):
        g = group(z, key)
        for _, row in g.head(60).iterrows():
            ax.plot(HOURS, row.to_numpy(), color=col, lw=0.4, alpha=0.25)
        ax.plot(HOURS, g.mean(), color=col, lw=1.4, alpha=1.0,
                label=f"{key}{' (ref.)' if col == FS.CONTROL else ''}")
    ax.set_xlabel("Time after release (h)"); ax.set_ylabel("Expression (z)")
    ax.set_xticks([2, 4, 6, 8, 10, 12])
    FS.legend(ax, loc="upper left")
    return fig


def panel_f(z):
    """Replicate agreement per compartment: median and interquartile range, as named rows."""
    keys = [k for k in WAVE + CONTINUOUS if len(group(z, k)) >= 5]
    FS.use()
    fig, ax = FS.figure(len(keys), vertical=False, value_mm=40.0)
    for y, k in enumerate(keys):
        v = z.loc[group(z, k).index, "rep_r"].dropna().to_numpy()
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        col = FS.ACCENT if k in WAVE else FS.CONTROL
        ax.hlines(y, q1, q3, color=col, lw=1.2, alpha=FS.ALPHA)
        ax.plot(med, y, "o", color=col, ms=4, alpha=1.0)
    ax.set_yticks(range(len(keys))); ax.set_yticklabels(keys)
    ax.set_ylim(len(keys) - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1); ax.set_xticks([0, 0.5, 1.0])
    ax.set_xlabel("Replicate agreement (Pearson r)")
    ax.plot([], [], "o-", color=FS.ACCENT, ms=4, label="wave-built")
    ax.plot([], [], "o-", color=FS.CONTROL, ms=4, label="continuous (ref.)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return fig


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    z = profiles()
    print(f"{len(z):,} genes with a cell-cycle profile and a compartment")
    a = newest("coexpression_components_localization_*/panel_which_evidence_places_a_protein.png")
    panels = [("a", a), ("b", save(panel_b(z), "b_wave")), ("c", save(panel_c(z), "c_peak")),
              ("d", save(panel_d(z), "d_coherence")), ("e", save(panel_e(z), "e_genes")),
              ("f", save(panel_f(z), "f_amplitude"))]
    import assemble_paper_figures as A
    A.PER_ROW = 3
    out = A.compose("06", [{"png": p, "letter": l} for l, p in panels], OUT)
    print("figure 6 ->", out)


if __name__ == "__main__":
    main()
