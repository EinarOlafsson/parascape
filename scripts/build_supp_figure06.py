"""Supplementary figure 6a: the design choices behind figure 6, each tested with data.

  A  compartment labels: exact fine hyperLOPIT names (chosen) vs substring matching on the
     coarse labels -- genes assigned to each niche under each
  B  the two replicate time courses: share of genes whose peak hour agrees within one hour
  C  timing coherence: share at the modal hour (chosen) vs concentration of peak phase on the
     12-hour circle (mean resultant length); the two agree
  D  compartments and their gene counts against the five-gene minimum
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
import build_figure06_v2 as F6                          # noqa: E402
import paper_panels as PP                               # noqa: E402

F = "supplementary_figure_6a"


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    n = pd.read_parquet(os.path.join(F6.CACHE, "nodes.parquet"))
    n = n.set_index(n["gene_id"].astype(str))
    keys = F6.WAVE + F6.CONTINUOUS
    cols = [f"cellcycle19092_{h}h_r{r}" for h in F6.HOURS for r in (1, 2)]
    srcs = [(os.path.join(F6.CACHE, "nodes.parquet"), cols + ["lopit_map", "lopit_unified"])]
    z = F6.profiles()

    # A: exact vs substring
    fine = n["lopit_map"].astype(str); coarse = n["lopit_unified"].astype(str)
    rows = []
    for k in keys:
        rows.append({"niche": k, "exact_fine_names": int(fine.isin(F6.NICHE.get(k, [k])).sum()),
                     "substring_on_coarse": int(coarse.str.lower().str.contains(k.lower(),
                                                                                regex=False).sum())})
    a = pd.DataFrame(rows)
    FS.use()
    fa, ax = FS.grouped_bars([a.substring_on_coarse.to_numpy(), a.exact_fine_names.to_numpy()],
                             ["substring, coarse labels", "exact, fine labels (used)"],
                             [k.replace("plasma membrane", "plasma memb.").replace(
                                 "mitochondrion", "mito.") for k in a.niche],
                             ylabel="Genes assigned")
    PP.panel(F, "a", fig=fa, data=a, stats={"niches with no genes under substring matching":
             int((a.substring_on_coarse == 0).sum())}, sources=srcs,
             what="exact fine labels vs substring matching")

    # B: replicate peak agreement
    r1 = pd.DataFrame({h: n[f"cellcycle19092_{h}h_r1"] for h in F6.HOURS}).loc[z.index]
    r2 = pd.DataFrame({h: n[f"cellcycle19092_{h}h_r2"] for h in F6.HOURS}).loc[z.index]
    ok = r1.notna().all(axis=1) & r2.notna().all(axis=1)
    p1 = np.argmax(r1[ok].to_numpy(), axis=1); p2 = np.argmax(r2[ok].to_numpy(), axis=1)
    dist = np.minimum(np.abs(p1 - p2), 12 - np.abs(p1 - p2))
    pk = pd.DataFrame({"comp": z.loc[ok, "comp"], "within1": dist <= 1})
    bb = []
    for k in keys:
        m = pk.comp.isin(F6.NICHE.get(k, [k]))
        if m.sum() >= 5:
            bb.append({"niche": k, "genes": int(m.sum()), "agree_within_1h_pct":
                       float(pk.within1[m].mean() * 100)})
    b = pd.DataFrame(bb)
    FS.use()
    fb, ax = FS.ref_bars(b.agree_within_1h_pct.to_numpy(),
                         [k.replace("plasma membrane", "plasma memb.").replace("mitochondrion",
                                                                               "mito.")
                          for k in b.niche], highlight=[k for k in b.niche if k in F6.WAVE],
                         ylabel="Replicates agree on peak (%)")
    ax.axhline(3 / 12 * 100, color=FS.INK, lw=0.6, ls=(0, (1.5, 1.5)),
               label="chance (within 1 h)")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02)); ax.set_ylim(0, 100)
    PP.panel(F, "b", fig=fb, data=b, stats={"genes": int(ok.sum())}, sources=srcs,
             what="replicate agreement on the peak hour")

    # C: coherence measures
    cc = []
    for k in keys:
        g = F6.group(z, k)
        if len(g) < 5:
            continue
        ph = np.argmax(g.to_numpy(), axis=1)
        ang = 2 * np.pi * ph / 12
        cc.append({"niche": k, "modal_share_pct": np.bincount(ph).max() / len(ph) * 100,
                   "resultant_length": float(np.hypot(np.cos(ang).mean(), np.sin(ang).mean()))})
    c = pd.DataFrame(cc)
    r = spearmanr(c.modal_share_pct, c.resultant_length)
    FS.use()
    fc, ax = FS.square_figure()
    wv = c.niche.isin(F6.WAVE)
    ax.plot(c.modal_share_pct[~wv], c.resultant_length[~wv], "o", color=FS.CONTROL, ms=4,
            alpha=FS.ALPHA, lw=0, label="continuous (ref.)")
    ax.plot(c.modal_share_pct[wv], c.resultant_length[wv], "o", color=FS.ACCENT, ms=4,
            alpha=FS.ALPHA, lw=0, label="wave-built")
    ax.plot([], [], " ", label=f"Spearman $\\rho$ = {r.statistic:.2f}")
    ax.set_xlabel("Genes peaking at the modal hour (%)")
    ax.set_ylabel("Peak-phase concentration (resultant length)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "c", fig=fc, data=c, stats={"Spearman rho": r.statistic, "p": r.pvalue},
             sources=srcs, what="two measures of timing coherence")

    # D: gene counts vs minimum
    dd = pd.DataFrame([{"niche": k, "genes": len(F6.group(z, k))}
                       for k in list(F6.NICHE)]).sort_values("genes", ascending=False)
    FS.use()
    fd, ax = FS.figure(len(dd), vertical=True)
    FS.bars(ax, dd.genes.to_numpy(), list(dd.niche), highlight=[k for k in dd.niche
                                                                if k in F6.WAVE])
    ax.axhline(5, color=FS.INK, lw=0.6, ls="--", label="minimum (5 genes)")
    ax.set_yscale("log"); ax.set_ylabel("Genes with a profile"); ax.set_xlim(-0.6, len(dd) - 0.4)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "d", fig=fd, data=dd, stats={"compartments": len(dd)}, sources=srcs,
             what="genes per compartment against the minimum")
    print(a.to_string(index=False)); print(b.round(1).to_string(index=False))
    PP.compose(F, "supplementary", "Every design choice behind figure 6, varied or checked.")


if __name__ == "__main__":
    main()
