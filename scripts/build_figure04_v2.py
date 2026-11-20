"""Figure 4, one claim: published treatment responses group by laboratory before biology -- two
contrasts from one study resemble each other more than the same perturbation done in two labs.

Einar, 2026-11-19: point unclear, placeholders unhelpful. All panels real, from the perturbation
atlas (201 log2 fold-change contrasts, 28 studies) and results/what_structures_contrasts_2026_11_18:

  A  similarity of every pair of contrasts, by same/other study x same/other kind   (the claim)
  B  all 201 x 201 similarities, contrasts grouped by study                         (the blocks)
  C  how well clusters of contrasts recover the study vs the kind of perturbation   (lab first)
  D  per study: similarity within vs to other studies (named rows)                  (every lab)
  E  one perturbation, two laboratories: BFD1 loss vs wild type, per gene           (the cost)
  F  the same pairs after removing each study's mean response                       (no fix)

Results (2026-11-19): same-study pairs median |rho| 0.303 / 0.164 (same / other kind) against
0.092 / 0.084 across studies; clusters recover the study (ARI 0.418 at 19 clusters) better than
the kind (0.283). BFD1 loss in two labs agrees at rho 0.226: above a typical cross-study pair,
below a typical same-study pair. Removing each study's mean response does not raise same-kind
agreement across studies (0.097 against 0.111 for other kind): the laboratory effect is not an
offset, which is why figure 2 matches the controls instead.
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402
from parascape import perturbation_atlas as PA          # noqa: E402

OUT = os.path.join(ROOT, "results", "figure04_v2")
PNL = os.path.join(OUT, "panels")
#: BFD1 loss vs wild type under bradyzoite-inducing stress, in two laboratories.
BFD1_A = "expr_lfc_GSE132238_genotype_delta_bfd1_vs_wt_condition_stress_time_72_hr"
BFD1_B = "expr_lfc_GSE218921_genotype_bfd1_ko_vs_wild_type_treatment_high_ph_low_c02_low"


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def pair_groups(rho, study, axis):
    cols = list(rho.columns)
    iu = np.triu_indices(len(cols), 1)
    s = np.array([study[c] for c in cols]); a = np.array([axis[c] for c in cols])
    same_s = (s[:, None] == s[None, :])[iu]
    same_a = (a[:, None] == a[None, :])[iu]
    v = np.abs(rho.to_numpy()[iu])
    ok = ~np.isnan(v)
    return {"same study\nsame kind": v[ok & same_s & same_a],
            "same study\nother kind": v[ok & same_s & ~same_a],
            "other study\nsame kind": v[ok & ~same_s & same_a],
            "other study\nother kind": v[ok & ~same_s & ~same_a]}


def violin_panel(values, ticks, grey, ylabel, legend):
    """Four violins with one-line tick names; which pair a violin belongs to is its colour.

    Two-line names ("same study / same kind") collided under 45-degree tilt; colour plus a legend
    carries the second factor instead.
    """
    FS.use()
    fig, ax = FS.figure(len(values), vertical=True)
    names = [f"{t}{i}" for i, t in enumerate(ticks)]
    FS.violin(ax, values, names, control=[names[i] for i in grey], vertical=True,
              draw_points=False)
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(ticks, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.6, len(values) - 0.4); ax.tick_params(axis="x", length=0)
    ax.set_ylabel(ylabel); ax.set_ylim(0, 1)
    for lab, col in legend:
        ax.fill_between([], [], color=col, alpha=FS.ALPHA, lw=0, label=lab)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return fig


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    store = PA.load()
    lfc, reg = store["lfc"], store["registry"]
    study = reg["study"].astype(str).to_dict()
    axis = reg["ladder_axis"].astype(str).to_dict()
    cols = sorted(lfc.columns, key=lambda c: (study[c], axis[c], c))
    lfc = lfc[cols]
    rho = lfc.corr(method="spearman", min_periods=1000)
    print(f"{len(cols)} contrasts, {len(set(study.values()))} studies")

    # A
    g = pair_groups(rho, study, axis)
    for k, v in g.items():
        print(f"  {k.replace(chr(10), ' '):24s} n={len(v):6d} median |rho| {np.median(v):.3f}")
    figa = violin_panel(list(g.values()), ["same kind", "other kind"] * 2, [2, 3],
                        "Similarity (|Spearman $\\rho$|)",
                        [("same study", FS.ACCENT), ("other study (ref.)", FS.CONTROL)])

    # B: heatmap grouped by study, block outlines
    FS.use()
    figb, ax = FS.square_figure()
    im = ax.imshow(rho.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
    edges = np.flatnonzero(np.array([study[c] for c in cols][1:]) !=
                           np.array([study[c] for c in cols][:-1])) + 0.5
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("Contrasts, grouped by study"); ax.set_ylabel("Contrasts, grouped by study")
    cb = figb.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Spearman $\\rho$"); cb.outline.set_linewidth(0.6)

    # C: cluster recovery of study vs kind
    d = 1 - rho.abs().fillna(0).to_numpy()
    np.fill_diagonal(d, 0); d = (d + d.T) / 2
    link = hierarchy.linkage(squareform(d, checks=False), method="average")
    ks = np.arange(2, 41)
    s_lab = [study[c] for c in cols]; a_lab = [axis[c] for c in cols]
    ari_s, ari_a = [], []
    for k in ks:
        lab = hierarchy.fcluster(link, t=k, criterion="maxclust")
        ari_s.append(adjusted_rand_score(s_lab, lab)); ari_a.append(adjusted_rand_score(a_lab, lab))
    FS.use()
    figc, ax = FS.square_figure()
    ax.plot(ks, ari_a, color=FS.CONTROL, lw=1.2, label="kind of perturbation (ref.)")
    ax.plot(ks, ari_s, color=FS.ACCENT, lw=1.2, label="study")
    ax.axhline(0, color=FS.INK, lw=0.6)
    ax.set_xlabel("Clusters of contrasts"); ax.set_ylabel("Adjusted Rand index")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    print(f"  ARI at 19 clusters: study {ari_s[17]:.3f}, kind {ari_a[17]:.3f}")

    # D: per study within vs across
    studies = sorted(set(s_lab))
    rows = []
    R = rho.abs()
    for st in studies:
        ins = [c for c in cols if study[c] == st]
        out = [c for c in cols if study[c] != st]
        if len(ins) < 2:
            continue
        w = R.loc[ins, ins].to_numpy()[np.triu_indices(len(ins), 1)]
        x = R.loc[ins, out].to_numpy().ravel()
        rows.append((st, np.nanmedian(w), np.nanmedian(x), len(ins)))
    dd = pd.DataFrame(rows, columns=["study", "within", "across", "n"]).sort_values("within",
                                                                                   ascending=False)
    FS.use()
    figd, ax = FS.ref_bars(dd["within"].to_numpy(), list(dd["study"]),
                           refs=dd["across"].to_numpy(), highlight=list(dd["study"]),
                           ref_label="to other studies",
                           ylabel="Median similarity within study (|$\\rho$|)")
    ax.set_ylim(0, 1)
    print(f"  studies with within > across: {(dd.within > dd.across).sum()} of {len(dd)}")

    # E: BFD1 loss in two labs
    both = lfc[[BFD1_A, BFD1_B]].dropna()
    r_e = spearmanr(both[BFD1_A], both[BFD1_B]).statistic
    FS.use()
    fige, ax = FS.square_figure()
    ax.plot(both[BFD1_A], both[BFD1_B], "o", color=FS.ACCENT, ms=1.2, alpha=0.3, lw=0,
            rasterized=True, label="gene")
    # Each axis on its own range: laboratory 1 moved far fewer genes far, and a shared range
    # squashed its response into a vertical stripe.
    lx = float(np.nanpercentile(np.abs(both[BFD1_A]), 99.8))
    ly = float(np.nanpercentile(np.abs(both[BFD1_B]), 99.8))
    ax.set_xlim(-lx, lx); ax.set_ylim(-ly, ly)
    ax.axhline(0, color=FS.INK, lw=0.4); ax.axvline(0, color=FS.INK, lw=0.4)
    ax.plot([], [], " ", label=f"Spearman $\\rho$ = {r_e:.2f}")
    ax.set_xlabel("BFD1 loss, laboratory 1 (log$_2$ FC)")
    ax.set_ylabel("BFD1 loss, laboratory 2 (log$_2$ FC)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    within_a = R.loc[BFD1_A, [c for c in cols if study[c] == study[BFD1_A] and c != BFD1_A]]
    print(f"  BFD1 across labs rho {r_e:.3f} over {len(both)} genes; median |rho| of lab-1 "
          f"BFD1 contrast to its own study's other contrasts {within_a.median():.3f}")

    # F: remove each study's mean response, recompute
    resid = lfc.copy()
    for st in studies:
        ins = [c for c in cols if study[c] == st]
        if len(ins) >= 2:
            resid[ins] = lfc[ins].sub(lfc[ins].mean(axis=1), axis=0)
    rho2 = resid.corr(method="spearman", min_periods=1000)
    g2 = pair_groups(rho2, study, axis)
    gf = {"before\nsame kind": g["other study\nsame kind"],
          "before\nother kind": g["other study\nother kind"],
          "study mean\nremoved,\nsame kind": g2["other study\nsame kind"],
          "study mean\nremoved,\nother kind": g2["other study\nother kind"]}
    for k, v in gf.items():
        print(f"  F {k.replace(chr(10), ' '):34s} median {np.median(v):.3f}")
    figf = violin_panel(list(gf.values()), ["same kind", "other kind"] * 2, [0, 1],
                        "Similarity across studies (|$\\rho$|)",
                        [("as published (ref.)", FS.CONTROL),
                         ("study mean removed", FS.ACCENT)])

    import paper_panels as PP
    from scipy.stats import mannwhitneyu
    F = "figure_4"
    ATLAS = [os.path.join(ROOT, "results", "what_structures_contrasts_2026_11_18", f)
             for f in ("pairs.tsv", "summary.json")]
    atlas_dir = f"perturbation atlas release {store.get('release', '')}".strip()
    def long(groups):
        return pd.concat([pd.DataFrame({"group": k.replace(chr(10), " "), "abs_rho": v})
                          for k, v in groups.items()], ignore_index=True)
    ga = long(g)
    u1 = mannwhitneyu(g["same study\nsame kind"], g["other study\nsame kind"], alternative="greater")
    u2 = mannwhitneyu(g["other study\nsame kind"], g["other study\nother kind"], alternative="greater")
    PP.panel(F, "a", fig=figa, data=ga, stats=pd.DataFrame(
        [{"group": k.replace(chr(10), " "), "pairs": len(v), "median_abs_rho": float(np.median(v))}
         for k, v in g.items()] +
        [{"group": "Mann-Whitney same study/same kind > other study/same kind",
          "pairs": None, "median_abs_rho": None, "U": u1.statistic, "p": u1.pvalue},
         {"group": "Mann-Whitney other study: same kind > other kind",
          "pairs": None, "median_abs_rho": None, "U": u2.statistic, "p": u2.pvalue}]),
             sources=ATLAS, what="similarity of every pair of contrasts by study and kind")
    PP.panel(F, "b", fig=figb, data=rho, stats={"contrasts": len(cols),
             "studies": len(set(s_lab))}, sources=ATLAS,
             what="all contrast similarities, grouped by study")
    PP.panel(F, "c", fig=figc, data=pd.DataFrame({"clusters": ks, "ARI_study": ari_s,
             "ARI_kind": ari_a}), stats={"ARI study at 19 clusters": ari_s[17],
             "ARI kind at 19 clusters": ari_a[17], "linkage": "average, 1 - |rho|"},
             sources=ATLAS, what="cluster recovery of study vs kind of perturbation")
    PP.panel(F, "d", fig=figd, data=dd, stats={"studies with within > across": int(
        (dd.within > dd.across).sum()), "studies": len(dd)}, sources=ATLAS,
             what="per study: similarity within vs to other studies")
    PP.panel(F, "e", fig=fige, data=both.rename(columns={BFD1_A: "laboratory_1_log2fc",
             BFD1_B: "laboratory_2_log2fc"}).rename_axis("gene_id"),
             stats={"Spearman rho": r_e, "genes": len(both), "laboratory 1": BFD1_A,
                    "laboratory 2": BFD1_B,
                    "median |rho| of lab-1 contrast to its own study": float(within_a.median())},
             sources=ATLAS, what="BFD1 loss in two laboratories")
    PP.panel(F, "f", fig=figf, data=long(gf), stats=pd.DataFrame(
        [{"group": k.replace(chr(10), " "), "pairs": len(v), "median_abs_rho": float(np.median(v))}
         for k, v in gf.items()]), sources=ATLAS,
             what="cross-study similarity after removing each study's mean response")
    out = PP.compose(F, "main", "Published treatment responses group by laboratory before "
                     "biology.", rows=[3, 3])
    print("figure 4 ->", out)


if __name__ == "__main__":
    main()
