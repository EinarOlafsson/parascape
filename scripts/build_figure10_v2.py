"""Figure 10, one claim: most of the host response to losing a secreted effector is shared across
effectors -- largely the response to losing export through the translocon -- and the parasite
itself changes in the same cell for almost none of them.

Einar, 2026-11-19: placeholders unhelpful. All panels real: host transcription after effector loss
(results/hosttx_per_host_gene_2026_11_18, 21 effectors, 7,625 host genes) and the dual
host-parasite perturb-seq (results/dual_perturb_single_cell_2026_11_18).

  A  per effector: host genes moved, and those it alone moves                    (mostly shared)
  B  how many effectors move each host gene                                      (the distribution)
  C  overlap of the host genes moved, effector x effector (map)                  (translocon block)
  D  the shared host genes: direction and size of change per effector (map)      (one response)
  E  per knockout in single cells: parasite vs host read change (log2)           (parasite unmoved)
  F  the host genes moved by the most effectors                                  (what is shared)
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

OUT = os.path.join(ROOT, "results", "figure10_v2")
PNL = os.path.join(OUT, "panels")
HR = os.path.join(ROOT, "results", "host_response_shared_2026_11_18")
HT = os.path.join(ROOT, "results", "hosttx_per_host_gene_2026_11_18",
                  "hosttx_host_gene_effects.parquet")
DP = os.path.join(ROOT, "results", "dual_perturb_single_cell_2026_11_18")


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def rows_fig(nr, left=0.9, value=1.6, right=1.6, rh=0.14):
    w, h = left + value + right, 0.15 + nr * rh + 0.6
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([left / w, 0.6 / h, value / w, nr * rh / h])
    ax.set_ylim(nr - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    return fig, ax


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    pe = pd.read_csv(os.path.join(HR, "per_effector.tsv"), sep="\t")
    pe["translocon"] = pe["translocon"].astype(str) == "True"
    ht = pd.read_parquet(HT)
    ht["effector"] = ht["effector"].astype(str)
    moved = {e: set(g["host_id"]) for e, g in ht.groupby("effector")}
    eff_order = list(pe.sort_values(["translocon", "host_genes_moved"],
                                    ascending=[False, False])["effector"])
    eff_order = [e for e in eff_order if e in moved]

    # A
    a = pe.set_index("effector").loc[eff_order]
    FS.use()
    figa, ax = rows_fig(len(a))
    y = np.arange(len(a))
    ax.hlines(y, a["its_own_only"].clip(lower=1), a["host_genes_moved"], color=FS.CONTROL, lw=0.6)
    ax.plot(a["host_genes_moved"], y, "o", color=FS.CONTROL, ms=3, alpha=FS.ALPHA,
            label="all moved (ref.)")
    ax.plot(a["its_own_only"].clip(lower=1), y, "o", color=FS.ACCENT, ms=3, alpha=FS.ALPHA,
            label="moved by it alone")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{e}*" if t else e for e, t in zip(a.index, a["translocon"])])
    ax.set_xscale("log"); ax.set_xticks([10, 100, 1000]); ax.minorticks_off()
    ax.set_xticklabels(["10", "100", "1,000"])
    ax.set_xlabel("Host genes moved")
    ax.plot([], [], " ", label="* translocon component")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    own = a["its_own_only"].sum() / a["host_genes_moved"].sum() * 100
    print(f"host genes moved by one effector alone: {own:.1f} % of effector-gene pairs")

    # B
    epg = pd.read_csv(os.path.join(HR, "effectors_per_host_gene.tsv"), sep="\t")
    counts = epg["effectors_moving_it"].value_counts().sort_index()
    FS.use()
    figb, ax = FS.square_figure()
    ax.bar(counts.index, counts.values, width=0.8, color=FS.CONTROL, alpha=FS.ALPHA, lw=0)
    ax.set_yscale("log"); ax.set_yticks([1, 10, 100, 1000]); ax.minorticks_off()
    ax.set_yticklabels(["1", "10", "100", "1,000"])
    ax.set_xlabel("Effectors moving the host gene"); ax.set_ylabel("Host genes")
    print(f"host genes moved by >= 2 effectors: "
          f"{(epg.effectors_moving_it >= 2).mean() * 100:.1f} % of {len(epg)}")

    # C
    J = np.array([[len(moved[x] & moved[z]) / max(1, len(moved[x] | moved[z])) * 100
                   for z in eff_order] for x in eff_order])
    FS.use()
    nl = len(eff_order); cell = 0.14
    w, h = 0.75 + nl * cell + 0.9, 0.15 + nl * cell + 0.65
    figc = plt.figure(figsize=(w, h))
    ax = figc.add_axes([0.75 / w, 0.65 / h, nl * cell / w, nl * cell / h])
    im = ax.imshow(J, cmap="Blues", vmin=0, vmax=60, interpolation="nearest")
    labs = [f"{e}*" if e in set(pe[pe.translocon].effector) else e for e in eff_order]
    ax.set_xticks(range(nl)); ax.set_xticklabels(labs, rotation=90)
    ax.set_yticks(range(nl)); ax.set_yticklabels(labs)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = figc.colorbar(im, ax=ax, fraction=0.06, pad=0.04, extend="max")
    cb.set_label("Host genes in common (%)"); cb.outline.set_linewidth(0.6)

    # D: the shared host genes (moved by >= 10 effectors), log2 FC per effector
    shared = set(epg.loc[epg.effectors_moving_it >= 10, "host_id"])
    sub = ht[ht.host_id.isin(shared)]
    m = sub.pivot_table(index="host_symbol", columns="effector", values="hosttx_log2fc",
                        aggfunc="mean").reindex(columns=eff_order)
    m = m.loc[m.mean(axis=1).sort_values().index]
    FS.use()
    nr, nc = m.shape; rw, cw = 0.045, 0.14
    w, h = 0.75 + nc * cw + 0.9, 0.15 + nr * rw + 0.65
    figd = plt.figure(figsize=(w, h))
    ax = figd.add_axes([0.75 / w, 0.65 / h, nc * cw / w, nr * rw / h])
    v = float(np.nanpercentile(np.abs(m.to_numpy()), 95))
    cmap = matplotlib.colormaps["RdBu_r"].copy(); cmap.set_bad(FS.BACKGROUND)
    im = ax.imshow(np.ma.masked_invalid(m.to_numpy()), cmap=cmap, vmin=-v, vmax=v,
                   aspect="auto", interpolation="nearest")
    ax.set_xticks(range(nc)); ax.set_xticklabels(labs, rotation=90)
    ax.set_yticks([]); ax.set_ylabel(f"Shared host genes")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = figd.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
    cb.set_label("Host change (log$_2$ FC)"); cb.outline.set_linewidth(0.6)
    ax.fill_between([], [], color=FS.BACKGROUND, label="not moved")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.55, 1.0))
    up = (m > 0).sum().sum() / m.notna().sum().sum() * 100
    print(f"shared host genes: {nr}; entries up {up:.1f} %")

    # E: dual perturb-seq, parasite vs host evidence per knockout
    frames = [pd.read_csv(f, sep="\t") for f in glob.glob(os.path.join(DP, "*_targets.tsv"))]
    dp = pd.concat(frames, ignore_index=True)
    # Effect sizes, not q: 80 of 86 knockouts have q = 1 on both axes and sat on one point.
    lp = np.log2(dp["parasite_reads_median"] / dp["parasite_reads_median_others"])
    lh = np.log2(dp["host_reads_median"] / dp["host_reads_median_others"])
    sig = (dp["parasite_q"] < 0.05) | (dp["host_q"] < 0.05)
    FS.use()
    fige, ax = FS.square_figure()
    ax.axhline(0, color=FS.INK, lw=0.6); ax.axvline(0, color=FS.INK, lw=0.6)
    isM = dp["target"].astype(str) == "MYR1"
    ax.plot(lh[~sig], lp[~sig], "o", color=FS.CONTROL, ms=3.5, alpha=FS.ALPHA, lw=0,
            label="knockout, q $\\geq$ 0.05 (ref.)")
    ax.plot(lh[sig & ~isM], lp[sig & ~isM], "o", color=FS.DARK_GREY, ms=3.5, lw=0,
            label="knockout, q < 0.05")
    ax.plot(lh[isM], lp[isM], "o", color=FS.ACCENT, ms=4.5, lw=0, label="MYR1")
    ax.set_xlabel("Host reads, knockout / others (log$_2$)")
    ax.set_ylabel("Parasite reads, knockout / others (log$_2$)")
    lim = float(max(np.abs(lp).max(), np.abs(lh).max())) * 1.1
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    print(f"knockouts {len(dp)}; host q<0.05 {(dp.host_q < .05).sum()}, parasite q<0.05 "
          f"{(dp.parasite_q < .05).sum()}")

    # F: host genes moved by the most effectors
    top = epg.sort_values("effectors_moving_it", ascending=False).head(15)
    sym = ht.drop_duplicates("host_id").set_index("host_id")["host_symbol"]
    FS.use()
    figf, ax = rows_fig(len(top), left=0.8, right=0.3, rh=0.2)
    yy = np.arange(len(top))
    ax.hlines(yy, 0, top["effectors_moving_it"], color=FS.ACCENT, lw=0.6, alpha=FS.ALPHA)
    ax.plot(top["effectors_moving_it"], yy, "o", color=FS.ACCENT, ms=3.5, alpha=FS.ALPHA)
    ax.set_yticks(yy); ax.set_yticklabels([sym.get(h, h) for h in top["host_id"]],
                                          fontstyle="italic")
    ax.set_xlim(0, len(eff_order) + 0.5); ax.set_xlabel("Effectors moving it")

    panels = [("a", save(figa, "a_per_effector")), ("b", save(figb, "b_effectors_per_gene")),
              ("c", save(figc, "c_overlap")), ("d", save(figd, "d_shared_genes")),
              ("e", save(fige, "e_dual_perturb")), ("f", save(figf, "f_top_host_genes"))]
    import assemble_paper_figures as A_
    out = A_.compose("10", [{"png": p, "letter": l} for l, p in panels], OUT, rows=[3, 3])
    print("figure 10 ->", out)


if __name__ == "__main__":
    main()
