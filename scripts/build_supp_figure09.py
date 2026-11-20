"""Supplementary figure 9a: the design choices behind figure 9, each tested with data.

  A  one shared vocabulary for the two atlases: location agreement with the raw names vs mapped
  B  one-to-one orthologues only vs every pairing of many-to-many orthogroups
  C  location agreement against 2,000 shuffles of the pairing (the chance level kappa removes)
  D  the pairing control: protein length ratio of orthologues vs random pairs
  E  fitness: how well the Toxoplasma score separates Plasmodium essential from dispensable (AUC)
  F  codon bias: observed correlation against 2,000 shuffled pairings
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
from sklearn.metrics import roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402
import build_what_transfers_between_species as WT       # noqa: E402
import paper_panels as PP                               # noqa: E402

F = "supplementary_figure_9a"
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")


def two_bars(vals, names, highlight, ylabel, ref=None, ref_label=None):
    FS.use()
    fig, ax = FS.figure(len(vals), vertical=True)
    FS.bars(ax, np.asarray(vals, float), names, highlight=[highlight])
    if ref is not None:
        ax.axhline(ref, color=FS.INK, lw=0.6, ls="--", label=ref_label)
        FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.set_ylabel(ylabel); ax.set_xlim(-0.6, len(vals) - 0.4)
    return fig


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    tg0 = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    pf0 = pd.read_parquet(os.path.join(CACHE, "pf_nodes.parquet"))
    keep_tg = ["orthogroup", "gene_id", "lopit_map", "fit_invitro_hff", "codon_cai_ribosomal",
               "length"]
    keep_pf = ["orthogroup", "gene_id", "lopit_pf_location", "piggybac_mis",
               "codon_cai_ribosomal", "length"]
    tg = tg0[keep_tg].dropna(subset=["orthogroup"]); pf = pf0[keep_pf].dropna(subset=["orthogroup"])
    pairs, counted = WT.one_to_one(tg, pf)
    srcs = [(os.path.join(CACHE, "nodes.parquet"), keep_tg),
            (os.path.join(CACHE, "pf_nodes.parquet"), keep_pf)]
    rng = np.random.default_rng(0)

    # A: vocabulary
    raw = pairs[["lopit_map", "lopit_pf_location"]].dropna()
    raw_agree = float((raw.lopit_map.astype(str).str.lower() ==
                       raw.lopit_pf_location.astype(str).str.lower()).mean() * 100)
    loc = pd.DataFrame({"tg": WT._normalise(pairs["lopit_map"]),
                        "pf": WT._normalise(pairs["lopit_pf_location"])}).dropna()
    mapped = float((loc.tg == loc.pf).mean() * 100)
    PP.panel(F, "a", fig=two_bars([raw_agree, mapped], ["raw names", "one vocabulary"],
                                  "one vocabulary", "Same compartment (%)"),
             data=pd.DataFrame({"vocabulary": ["raw names", "one vocabulary"],
                                "same_compartment_pct": [raw_agree, mapped]}),
             stats={"pairs, raw": len(raw), "pairs, mapped": len(loc)}, sources=srcs,
             what="raw atlas names vs one shared vocabulary")

    # B: one-to-one vs all pairings
    t2 = tg.dropna(subset=["lopit_map"]); p2 = pf.dropna(subset=["lopit_pf_location"])
    allp = t2.merge(p2, on="orthogroup", suffixes=("_tg", "_pf"))
    allp_agree = float((WT._normalise(allp.lopit_map) == WT._normalise(allp.lopit_pf_location))
                       .mean() * 100)
    PP.panel(F, "b", fig=two_bars([mapped, allp_agree], ["one-to-one", "all pairings"],
                                  "one-to-one", "Same compartment (%)"),
             data=pd.DataFrame({"pairs": ["one-to-one", "all pairings"],
                                "same_compartment_pct": [mapped, allp_agree],
                                "n": [len(loc), len(allp)]}),
             stats={"orthogroups set aside as many-to-many": counted["many_to_many_set_aside"],
                    "one-to-one pairs": counted["one_to_one"]}, sources=srcs,
             what="one-to-one orthologues vs every pairing")

    # C: shuffled pairing
    null = np.array([float((rng.permutation(loc.tg.to_numpy()) == loc.pf.to_numpy()).mean() * 100)
                     for _ in range(2000)])
    FS.use()
    fc, ax = FS.square_figure()
    ax.hist(null, bins=30, color=FS.CONTROL, alpha=FS.ALPHA, lw=0, label="shuffled pairing")
    ax.axvline(mapped, color=FS.ACCENT, lw=1.4, label="orthologues")
    ax.set_xlabel("Same compartment (%)"); ax.set_ylabel("Shuffles")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "c", fig=fc, data=pd.DataFrame({"shuffled_same_pct": null}),
             stats={"observed": mapped, "null mean": null.mean(), "null sd": null.std(ddof=1),
                    "permutation p": (1 + int((null >= mapped).sum())) / 2001}, sources=srcs,
             what="location agreement against shuffled pairings")

    # D: length ratio
    lr = np.log2(pairs.length_tg / pairs.length_pf).dropna()
    rp = np.log2(pairs.length_tg.sample(frac=1, random_state=0).to_numpy() /
                 pairs.length_pf.to_numpy())
    rp = rp[np.isfinite(rp)]
    FS.use()
    fd, ax = FS.figure(2, vertical=True)
    FS.violin(ax, [rp, lr.to_numpy()], ["random pairs", "orthologues"], control=["random pairs"],
              vertical=True, draw_points=False)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["random pairs", "orthologues"], rotation=45,
                                              ha="right", rotation_mode="anchor")
    ax.set_ylabel("Length, Toxoplasma / Plasmodium (log$_2$)"); ax.set_xlim(-0.6, 1.4)
    ax.axhline(0, color=FS.INK, lw=0.6)
    PP.panel(F, "d", fig=fd, data=pd.DataFrame({"orthologue_log2_ratio": lr.reset_index(drop=True)}),
             stats={"IQR orthologues": float(np.subtract(*np.percentile(lr, [75, 25]))),
                    "IQR random pairs": float(np.subtract(*np.percentile(rp, [75, 25])))},
             sources=srcs, what="length ratio: orthologues vs random pairs")

    # E: fitness AUC
    fb = pd.DataFrame({"tg": pd.to_numeric(pairs.fit_invitro_hff, errors="coerce"),
                       "pf": pd.to_numeric(pairs.piggybac_mis, errors="coerce")}).dropna()
    ess = (fb.pf < 0.5).astype(int)
    auc = roc_auc_score(ess, -fb.tg)
    nul = np.array([roc_auc_score(rng.permutation(ess.to_numpy()), -fb.tg) for _ in range(500)])
    PP.panel(F, "e", fig=two_bars([auc, nul.mean()], ["orthologues", "shuffled"], "orthologues",
                                  "AUC, Plasmodium essential", ref=0.5, ref_label="chance"),
             data=fb.assign(pf_essential=ess), stats={"AUC": auc, "shuffled mean": nul.mean(),
             "essential cut": "Plasmodium mutagenesis index < 0.5", "pairs": len(fb)},
             sources=srcs, what="Toxoplasma fitness separating Plasmodium essential genes")

    # F: codon null
    cc = pairs[["codon_cai_ribosomal_tg", "codon_cai_ribosomal_pf"]].dropna()
    obs = spearmanr(cc.iloc[:, 0], cc.iloc[:, 1]).statistic
    nc = np.array([spearmanr(rng.permutation(cc.iloc[:, 0].to_numpy()), cc.iloc[:, 1]).statistic
                   for _ in range(2000)])
    FS.use()
    ff, ax = FS.square_figure()
    ax.hist(nc, bins=30, color=FS.CONTROL, alpha=FS.ALPHA, lw=0, label="shuffled pairing")
    ax.axvline(obs, color=FS.ACCENT, lw=1.4, label="orthologues")
    ax.set_xlabel("Spearman $\\rho$, codon adaptation"); ax.set_ylabel("Shuffles")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "f", fig=ff, data=pd.DataFrame({"shuffled_rho": nc}),
             stats={"observed rho": obs, "null sd": nc.std(ddof=1),
                    "permutation p": (1 + int((nc >= obs).sum())) / 2001, "pairs": len(cc)},
             sources=srcs, what="codon adaptation against shuffled pairings")
    print(raw_agree, mapped, allp_agree, auc, obs)
    PP.compose(F, "supplementary", "Every design choice behind figure 9, varied or checked.")


if __name__ == "__main__":
    main()
