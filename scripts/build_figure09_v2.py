"""Figure 9, one claim: where a protein is, and whether losing it costs fitness, transfer between
Toxoplasma and Plasmodium; how its message is written does not.

Einar, 2026-11-19: placeholders unhelpful. All panels real, on the 2,221 one-to-one orthologue
pairs of build_what_transfers_between_species.py (many-to-many orthogroups set aside):

  A  agreement beyond chance per layer: Cohen kappa (location) or Spearman (rest)   (the claim)
  B  location in Toxoplasma vs Plasmodium, row-normalised (map)                      (what moves)
  C  fitness: Toxoplasma score per fifth of the Plasmodium mutagenesis index       (transfers)
  D  codon bias: CAI in each species                                                 (does not)
  E  protein length in each species                                                  (pairing control)
  F  location agreement per Toxoplasma compartment vs guessing from Plasmodium shares (where)
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
import build_what_transfers_between_species as WT       # noqa: E402

OUT = os.path.join(ROOT, "results", "figure09_v2")
PNL = os.path.join(OUT, "panels")
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
COMP_SHORT = {"endoplasmic reticulum": "ER", "plasma membrane": "plasma memb.",
              "mitochondrion": "mito.",
              "golgi": "Golgi", "imc": "IMC", "eif3": "eIF3", "maurer's cleft": "Maurer's cleft",
              "pvm": "PVM"}


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def kappa(a, b):
    labs = sorted(set(a) | set(b))
    po = float((a == b).mean())
    pe = sum(float((a == l).mean()) * float((b == l).mean()) for l in labs)
    return (po - pe) / (1 - pe)


def scatter(x, y, xlabel, ylabel, log=False):
    both = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"),
                         "y": pd.to_numeric(y, errors="coerce")}).dropna()
    if log:
        both = both[(both.x > 0) & (both.y > 0)]
    r = spearmanr(both.x, both.y).statistic
    FS.use()
    fig, ax = FS.square_figure()
    ax.plot(both.x, both.y, "o", color=FS.ACCENT, ms=1.3, alpha=0.35, lw=0, rasterized=True,
            label="orthologue pair")
    ax.plot([], [], " ", label=f"Spearman $\\rho$ = {r:.2f}")
    if log:
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xticks([100, 1000, 10000]); ax.set_xticklabels(["100", "1,000", "10,000"])
        ax.set_yticks([100, 1000, 10000]); ax.set_yticklabels(["100", "1,000", "10,000"])
        ax.minorticks_off()
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return fig, r, len(both)


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    tg = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    pf = pd.read_parquet(os.path.join(CACHE, "pf_nodes.parquet"))
    keep_tg = ["orthogroup", "gene_id", "lopit_map", "fit_invitro_hff", "codon_cai_ribosomal",
               "length"]
    keep_pf = ["orthogroup", "gene_id", "lopit_pf_location", "piggybac_mis",
               "codon_cai_ribosomal", "length"]
    tg = tg[[c for c in keep_tg if c in tg.columns]].dropna(subset=["orthogroup"])
    pf = pf[[c for c in keep_pf if c in pf.columns]].dropna(subset=["orthogroup"])
    pairs, counted = WT.one_to_one(tg, pf)
    print(f"{counted['one_to_one']:,} one-to-one pairs")

    loc = pd.DataFrame({"tg": WT._normalise(pairs["lopit_map"]),
                        "pf": WT._normalise(pairs["lopit_pf_location"])}).dropna()
    k = kappa(loc.tg, loc.pf)

    # Fitness as violins of the Toxoplasma score per fifth of the Plasmodium index: the index
    # piles up at its two ends, which made the scatter a pair of stripes.
    fb = pd.DataFrame({"tg": pd.to_numeric(pairs["fit_invitro_hff"], errors="coerce"),
                       "pf": pd.to_numeric(pairs["piggybac_mis"], errors="coerce")}).dropna()
    r_fit, n_fit = spearmanr(fb.tg, fb.pf).statistic, len(fb)
    fb["q"] = pd.qcut(fb["pf"].rank(method="first"), 5, labels=False)
    edges = [f"{fb.loc[fb.q == i, 'pf'].min():.2f}-{fb.loc[fb.q == i, 'pf'].max():.2f}"
             for i in range(5)]
    FS.use()
    figc, ax = FS.figure(5, vertical=True)
    FS.violin(ax, [fb.loc[fb.q == i, "tg"].to_numpy() for i in range(5)], edges,
              control=edges[1:], vertical=True, draw_points=False)
    ax.set_xticks(range(5)); ax.set_xticklabels(edges, rotation=45, ha="right",
                                                 rotation_mode="anchor")
    ax.set_xlim(-0.6, 4.4); ax.tick_params(axis="x", length=0)
    ax.set_xlabel("P. falciparum MIS (quintile)")
    ax.set_ylabel("T. gondii phenotype score")
    FS.pad_limits(ax, fb["tg"])
    figd, r_cai, n_cai = scatter(pairs["codon_cai_ribosomal_tg"], pairs["codon_cai_ribosomal_pf"],
                                 "CAI, T. gondii",
                                 "CAI, P. falciparum")
    fige, r_len, n_len = scatter(pairs["length_tg"], pairs["length_pf"],
                                 "Toxoplasma protein length (aa)",
                                 "Plasmodium protein length (aa)", log=True)

    # A
    layers = [("location", k, len(loc)), ("fitness", r_fit, n_fit),
              ("codon bias", r_cai, n_cai), ("length (control)", r_len, n_len)]
    print({a: round(b, 3) for a, b, _ in layers})
    FS.use()
    figa, ax = FS.ref_bars([l[1] for l in layers], [l[0].replace(" (control)", "") for l in layers],
                           highlight=["location", "fitness", "codon bias"],
                           ylabel="Agreement ($\\kappa$ or $\\rho$)")
    ax.set_ylim(0, 1); ax.axhline(0, color=FS.INK, lw=0.6)

    # B
    rows_ = loc.tg.value_counts()
    rows_ = list(rows_[rows_ >= 10].index)
    cols_ = list(loc.pf.value_counts().index)
    cols_ = rows_ + [c for c in cols_ if c not in rows_]
    m = pd.crosstab(loc.tg, loc.pf).reindex(index=rows_, columns=cols_, fill_value=0)
    m = m.div(m.sum(axis=1), axis=0) * 100
    m = m.loc[:, m.max(axis=0) >= 2]      # Plasmodium niches no Toxoplasma orthologue reaches
    FS.use()
    nr, nc = m.shape; cell = 0.16
    w, h = 1.3 + nc * cell + 0.9, 0.15 + nr * cell + 1.25
    figb = plt.figure(figsize=(w, h))
    ax = figb.add_axes([1.3 / w, 1.25 / h, nc * cell / w, nr * cell / h])
    im = ax.imshow(m.to_numpy(), cmap="Blues", vmin=0, vmax=100, interpolation="nearest")
    ax.set_xticks(range(nc)); ax.set_xticklabels([COMP_SHORT.get(c, c) for c in m.columns],
                                                 rotation=90)
    ax.set_yticks(range(nr)); ax.set_yticklabels([COMP_SHORT.get(c, c) for c in rows_])
    ax.set_xlabel("P. falciparum niche"); ax.set_ylabel("T. gondii niche")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = figb.colorbar(im, ax=ax, fraction=0.06, pad=0.04)
    cb.set_label("Orthologues (% of row)"); cb.outline.set_linewidth(0.6)

    # F: per Toxoplasma compartment, agreement vs the Plasmodium share of that label
    share_pf = loc.pf.value_counts(normalize=True) * 100
    f = []
    for c in rows_:
        s = loc[loc.tg == c]
        f.append((c, float((s.pf == c).mean() * 100), float(share_pf.get(c, 0.0)), len(s)))
    f = pd.DataFrame(f, columns=["comp", "agree", "chance", "n"]).sort_values("agree",
                                                                           ascending=False)
    FS.use()
    figf, ax = FS.ref_bars(f["agree"].to_numpy(), [COMP_SHORT.get(c, c) for c in f["comp"]],
                           refs=f["chance"].to_numpy(),
                           highlight=[COMP_SHORT.get(c, c) for c in f["comp"]],
                           ref_label="share of the compartment in Plasmodium",
                           ylabel="Orthologues in the same compartment (%)")
    ax.set_ylim(0, 105)
    print(f.round(1).to_string(index=False))

    import paper_panels as PP
    F = "figure_7__final"
    from scipy.stats import spearmanr as _sp
    NOD = (os.path.join(CACHE, "nodes.parquet"), keep_tg)
    PFN = os.path.join(CACHE, "pf_nodes.parquet")
    def sp(a, b):
        x = pd.DataFrame({"a": pd.to_numeric(a, errors="coerce"),
                          "b": pd.to_numeric(b, errors="coerce")}).dropna()
        r = _sp(x.a, x.b); return r.statistic, r.pvalue, len(x)
    rf, pf_, nf = sp(pairs["fit_invitro_hff"], pairs["piggybac_mis"])
    rc, pc, nc_ = sp(pairs["codon_cai_ribosomal_tg"], pairs["codon_cai_ribosomal_pf"])
    rl, pl_, nl_ = sp(pairs["length_tg"], pairs["length_pf"])
    agree = float((loc.tg == loc.pf).mean())
    PP.panel(F, "a", fig=figb, data=m.rename_axis("toxoplasma"), stats={
             "Cohen kappa": k, "raw agreement": agree, "pairs": len(loc),
             "one-to-one orthologue pairs": counted["one_to_one"]}, sources=[NOD, PFN],
             what="location in Toxoplasma vs Plasmodium")
    PP.panel(F, "b", fig=figa, data=pd.DataFrame(layers, columns=["layer", "agreement", "pairs"]),
             stats=pd.DataFrame([{"layer": "location", "measure": "Cohen kappa", "value": k,
                                  "pairs": len(loc), "raw agreement": agree},
                                 {"layer": "fitness", "measure": "Spearman", "value": rf,
                                  "p": pf_, "pairs": nf},
                                 {"layer": "codon bias", "measure": "Spearman", "value": rc,
                                  "p": pc, "pairs": nc_},
                                 {"layer": "length (control)", "measure": "Spearman",
                                  "value": rl, "p": pl_, "pairs": nl_}]),
             sources=[NOD, PFN], what="agreement beyond chance per layer")
    PP.panel(F, "c", fig=figc, data=fb, stats={"Spearman rho": rf, "p": pf_, "pairs": nf},
             sources=[NOD, PFN], what="fitness across the two species")
    PP.panel(F, "d", fig=figd, data=pairs[["codon_cai_ribosomal_tg", "codon_cai_ribosomal_pf"]],
             stats={"Spearman rho": rc, "p": pc, "pairs": nc_}, sources=[NOD, PFN],
             what="codon adaptation in the two species")
    out = PP.compose(F, "main", "Localisation and essentiality transfer between Toxoplasma and "
                     "Plasmodium orthologues; codon usage does not.")
    print("figure 9 ->", out)


if __name__ == "__main__":
    main()
