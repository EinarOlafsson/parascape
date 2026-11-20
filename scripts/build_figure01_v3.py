"""Figure 1 (storyboard 2026-11-20): published data cover most Toxoplasma genes for expression and
fitness, almost none for host response or mutant phenotype.

  A  flowchart: published data -> ParaScape -> uses               (what ParaScape is)
  B  cumulative datasets per year, by data type                  (what was published, when)
  C  presence matrix, data type x gene                           (which genes are measured)
  D  data types measured per gene                                (how much is known per gene)
  E  % of genes measured per data type, secreted vs all          (are host-facing proteins better)
"""
from __future__ import annotations

import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import datasets as DS                    # noqa: E402
from parascape import figstyle as FS                    # noqa: E402
import build_figure01_v2 as V2                          # noqa: E402
import draw_figure01a_flowchart as FC                   # noqa: E402
import paper_panels as PP                               # noqa: E402

CACHE = V2.CACHE
F = "figure_1"
LEVEL = {"post_translation": "protein", "transcription": "RNA", "translation": "translation",
         "DNA": "DNA", "reference": "reference"}


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    NODES = os.path.join(CACHE, "nodes.parquet")
    REG = os.path.join(ROOT, "parascape", "datasets.py")
    nodes = pd.read_parquet(NODES)
    cols = V2.knowledge_columns(nodes)
    G = len(nodes)

    # A
    PP.panel(F, "a", fig=FC.draw(), data=pd.DataFrame({"input": FC.INPUTS[:len(FC.USES)],
             "use": FC.USES}), stats={"datasets": len(DS.REGISTRY), "genes": G},
             sources=[REG], what="ParaScape: published data in, what it does, what it is used for")

    # B: cumulative datasets per year, stacked by level of biology
    rows = []
    for d in DS.REGISTRY:
        y = re.findall(r"\b(19[5-9]\d|20[0-2]\d)\b", str(d.citation or ""))
        if y:
            rows.append({"key": d.key, "year": int(y[0]), "level": LEVEL.get(d.level, d.level)})
    yr = pd.DataFrame(rows)
    years = np.arange(2005, 2027)
    order = list(yr.level.value_counts().index)
    cum = np.array([[int(((yr.level == l) & (yr.year <= t)).sum()) for t in years] for l in order])
    FS.use()
    figb, ax = FS.square_figure()
    colors = [FS.ACCENT, FS.TEAL, FS.NAVY, FS.ACCENT_LIGHT, FS.MAGENTA][:len(order)]
    ax.stackplot(years, cum, colors=colors, alpha=FS.ALPHA, lw=0, labels=order)
    ax.set_xlim(years[0], years[-1]); ax.set_xticks([2005, 2015, 2025])
    ax.set_xlabel("Year of publication"); ax.set_ylabel("Datasets (cumulative)")
    FS.legend(ax, loc="upper left")
    PP.panel(F, "b", fig=figb, data=pd.DataFrame(cum.T, index=pd.Index(years, name="year"),
             columns=order), stats={"datasets with a year": len(yr),
             "datasets without a year": len(DS.REGISTRY) - len(yr)}, sources=[REG],
             what="cumulative datasets per year by data type")

    # C: presence matrix, data types as rows, genes as columns; % covered as a right margin
    share = {c: nodes[v].notna().any(axis=1).mean() * 100 for c, v in cols.items()}
    ordc = sorted(share, key=share.get, reverse=True)
    names = [V2.kind_name(c) for c in ordc]
    m = np.column_stack([nodes[cols[c]].notna().any(axis=1).to_numpy() for c in ordc]).T
    m = m[:, np.lexsort(m[::-1])[::-1]]
    m = m[:, np.argsort(-m.sum(axis=0), kind="stable")]
    FS.use()
    figc = plt.figure(figsize=(120 / 25.4, 56 / 25.4))
    ax = figc.add_axes([0.17, 0.14, 0.60, 0.83])
    rgb = np.ones(m.shape + (3,))
    thin = np.array([c in V2.THIN for c in ordc])
    rgb[m & ~thin[:, None]] = matplotlib.colors.to_rgb(FS.DARK_GREY)
    rgb[m & thin[:, None]] = matplotlib.colors.to_rgb(FS.ACCENT)
    ax.imshow(rgb, aspect="auto", interpolation="nearest", rasterized=True)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names); ax.tick_params(axis="y", length=0)
    ax.set_xticks([0, G - 1]); ax.set_xticklabels(["1", f"{G:,}"])
    ax.set_xlabel("Genes, most to least measured")
    mx = figc.add_axes([0.82, 0.14, 0.15, 0.83], sharey=ax)
    mx.barh(range(len(names)), [share[c] for c in ordc], height=0.6,
            color=[FS.ACCENT if t else FS.CONTROL for t in thin], alpha=FS.ALPHA, lw=0)
    mx.set_xlim(0, 100); mx.set_xticks([0, 100]); mx.set_xlabel("Genes (%)")
    plt.setp(mx.get_yticklabels(), visible=False); mx.tick_params(axis="y", length=0)
    pres = pd.DataFrame({V2.kind_name(c): nodes[cols[c]].notna().any(axis=1) for c in ordc})
    pres.insert(0, "gene_id", nodes["gene_id"].astype(str).to_numpy())
    PP.panel(F, "c", fig=figc, data=pres, stats=pd.DataFrame(
        [{"data type": V2.kind_name(c), "genes_measured_pct": share[c]} for c in ordc]),
             sources=[(NODES, [c for v in cols.values() for c in v])],
             what="which genes are measured by which data type")

    # D: data types per gene
    per_gene = m.sum(axis=0)
    FS.use()
    figd, ax = FS.square_figure()
    ax.hist(per_gene, bins=np.arange(-0.5, len(ordc) + 1.5), color=FS.CONTROL, alpha=FS.ALPHA,
            lw=0)
    ax.axvline(np.median(per_gene), color=FS.INK, lw=0.6, ls=(0, (1.5, 1.5)), label="median")
    ax.set_xlabel("Data types per gene"); ax.set_ylabel("Genes")
    FS.legend(ax, loc="upper left")
    PP.panel(F, "d", fig=figd, data=pd.DataFrame({"gene_id": nodes["gene_id"].astype(str),
             "data_types": per_gene}), stats={"median": float(np.median(per_gene)),
             "genes with <= 5 data types": int((per_gene <= 5).sum())}, sources=[NODES],
             what="data types measured per gene")

    # E: heat strip, % measured per data type in secreted proteins and all genes
    sec = nodes["lopit_map"].astype(str).isin(V2.SECRETED).to_numpy()
    ce = [c for c in ordc if c != "place"]
    strip = np.array([[nodes[cols[c]].notna().any(axis=1).to_numpy().mean() * 100 for c in ce],
                      [nodes[cols[c]].notna().any(axis=1).to_numpy()[sec].mean() * 100
                       for c in ce]])
    FS.use()
    fige, ax = FS.plot_figure(len(ce) * 5.4, 2 * 6.0, extra_bottom_mm=10)
    im = ax.imshow(strip, cmap="Blues", vmin=0, vmax=100, aspect="auto")
    ax.set_yticks([0, 1]); ax.set_yticklabels(["all genes", "secreted"])
    ax.set_xticks(range(len(ce)))
    ax.set_xticklabels([V2.kind_name(c) for c in ce], rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fige.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Genes (%)"); cb.outline.set_linewidth(0.6)
    PP.panel(F, "e", fig=fige, data=pd.DataFrame(strip, index=["all genes", "secreted"],
             columns=[V2.kind_name(c) for c in ce]), stats={"secreted proteins": int(sec.sum())},
             sources=[(NODES, ["lopit_map"])], what="% measured per data type, secreted vs all")
    print(PP.compose(F, "main", "Published data cover most Toxoplasma genes for expression and "
                     "fitness, almost none for host response or mutant phenotype."))


if __name__ == "__main__":
    main()
