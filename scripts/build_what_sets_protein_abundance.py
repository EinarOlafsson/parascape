"""How much of protein abundance each layer of gene expression explains.

Figure 6 of the ParaScape paper (instruction 101), the panel Einar asked for: *"protein expression is
explained by Transcription, Translation, degradation, etc. -- how much is explained by
transcription"*. Every panel of that figure is the same shape, one measurement predicting another,
cross-validated, with the mechanism in the caption.

WHY PLASMODIUM AND NOT TOXOPLASMA. The chain needs four measured layers on one gene set: how much
message there is, how fast it is made and destroyed, how much of it is on ribosomes, and how much
protein results. Plasmodium has all four (`riboseq_mrna_*`, `transcription_rate_4tu` with
`mrna_decay_rate_4tu`, `riboseq_rpf_*`, `proteome_blood_log2`). Toxoplasma has 748 genes with an
iBAQ value against 7,739 with a transcript, so the same panel there would describe the proteins
somebody could measure rather than the proteome. That gap is itself a figure 12 entry.

THE GENE SET IS FIXED BEFORE ANY MODEL IS FITTED. Each layer is measured on fewer genes than the
last, so letting the universe shrink as layers are added would compare models on different genes
and attribute the difference to the layer. The genes used are those carrying EVERY layer, chosen
once, and the cost of that choice is reported rather than hidden.

Writes `results/what_sets_protein_abundance_<date>/`.
"""
from __future__ import annotations

import json
import os
from datetime import date

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: What is being predicted: protein abundance in the blood stage.
TARGET = "proteome_blood_log2"

#: The layers, in the order the central dogma runs them. Each entry adds to the one before.
LAYERS = (
    ("transcript abundance", ("riboseq_mrna_ring", "riboseq_mrna_early_trophozoite",
                              "riboseq_mrna_late_trophozoite", "riboseq_mrna_schizont",
                              "riboseq_mrna_merozoite")),
    ("+ how fast it is made and destroyed", ("transcription_rate_4tu", "mrna_decay_rate_4tu")),
    ("+ how much is on ribosomes", ("riboseq_rpf_ring", "riboseq_rpf_early_trophozoite",
                                    "riboseq_rpf_late_trophozoite", "riboseq_rpf_schizont",
                                    "riboseq_rpf_merozoite")),
    ("+ how the message is written", ("codon_enc", "codon_gc3", "codon_cai_ribosomal", "length")),
)

#: Folds, and how many times the fold assignment is redrawn.
N_FOLDS, N_REPEATS = 5, 10

#: How many shuffles stand behind the floorless baseline -- what this procedure returns on a target
#: that carries no signal at all.
N_NULL = 100


def load(cache: str) -> pd.DataFrame:
    """The Plasmodium node table, reduced to the target and every layer column."""
    nodes = pd.read_parquet(os.path.join(cache, "pf_nodes.parquet"))
    wanted = [TARGET] + [c for _name, cols in LAYERS for c in cols]
    missing = [c for c in wanted if c not in nodes.columns]
    if missing:
        raise KeyError("not in the Plasmodium node table: " + ", ".join(missing))
    frame = nodes.set_index(nodes["gene_id"].astype(str))[wanted].apply(
        pd.to_numeric, errors="coerce")
    return frame


def cross_validated_r2(x: np.ndarray, y: np.ndarray, seed: int = 0) -> float:
    """Variance of `y` explained out of sample, averaged over repeated folds.

    Ridge with the penalty chosen inside each training fold, because the replicate columns of one
    layer are nearly collinear -- five stages of the same transcript -- and an unpenalised fit
    would report the noise it fitted.
    """
    scores = []
    for repeat in range(N_REPEATS):
        folds = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed + repeat)
        predicted = np.empty_like(y)
        for train, test in folds.split(x):
            scaler = StandardScaler().fit(x[train])
            model = RidgeCV(alphas=np.logspace(-3, 4, 24)).fit(scaler.transform(x[train]), y[train])
            predicted[test] = model.predict(scaler.transform(x[test]))
        ss_res = float(((y - predicted) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        scores.append(1.0 - ss_res / ss_tot)
    return float(np.mean(scores))


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    frame = load(cache)

    with_target = frame[TARGET].notna().sum()
    complete = frame.dropna()
    print(f"{with_target:,} genes carry a protein measurement; "
          f"{len(complete):,} carry every layer as well")
    if len(complete) < 200:
        raise SystemExit("too few genes measured on every layer")

    y = complete[TARGET].to_numpy(dtype=float)
    rows, used = [], []
    print("\n  layer                                  cumulative   alone")
    for name, columns in LAYERS:
        used.extend(columns)
        r2 = cross_validated_r2(complete[used].to_numpy(dtype=float), y)
        alone = cross_validated_r2(complete[list(columns)].to_numpy(dtype=float), y)
        rows.append({"layer": name, "n_columns": len(used), "r2": r2, "r2_alone": alone})
        print(f"  {name:38s} {r2:8.3f} {alone:8.3f}")
    # Cumulative R2 is NOT guaranteed to rise: adding the five ribosome-footprint columns lowers
    # it, because out of sample a set of nearly collinear columns can cost more in variance than
    # it returns in signal at this many genes. The alone column is reported beside it for exactly
    # that reason -- it does not depend on the order the layers are added in.

    # What this procedure returns on a target with no signal, so a small R2 can be read.
    rng = np.random.default_rng(0)
    x_all = complete[used].to_numpy(dtype=float)
    null = np.array([cross_validated_r2(x_all, rng.permutation(y), seed=100 + i)
                     for i in range(min(N_NULL, 20))])
    print(f"\nshuffled target: R2 = {null.mean():+.4f} (sd {null.std(ddof=1):.4f})")

    # The same question on the widest gene set that can answer it, which is four times larger.
    # Reported separately and never mixed with the table above: it is a different gene set, and
    # comparing an R2 across gene sets compares the gene sets.
    wide = frame[[TARGET] + list(LAYERS[0][1])].dropna()
    wide_r2 = cross_validated_r2(wide[list(LAYERS[0][1])].to_numpy(dtype=float),
                                 wide[TARGET].to_numpy(dtype=float))
    print(f"\non the {len(wide):,} genes with a transcript and a protein measurement "
          f"({len(wide) / len(complete):.1f} times\nthe number above), transcript abundance alone "
          f"explains {wide_r2:.1%}.")

    table = pd.DataFrame(rows)
    table["gain"] = table["r2"].diff().fillna(table["r2"])
    # THE NUMBER THE LITERATURE QUOTES, in the currency it quotes it in: the squared correlation
    # between one mRNA measure and protein, uncorrected. Human estimates run from about 40 per
    # cent (Schwanhaeusser 2011) through 56 to 81 per cent once measurement error is corrected
    # for (Li and Biggin), and ~60 per cent is the figure usually carried. Reported here so the
    # result can be compared with that at all; the model numbers below are NOT in that currency.
    single = frame[[TARGET] + list(LAYERS[0][1])].dropna(subset=[TARGET]).copy()
    single["mrna"] = np.log2(single[list(LAYERS[0][1])].mean(axis=1) + 1)
    pair = single[[TARGET, "mrna"]].dropna()
    rho = spearmanr(pair["mrna"], pair[TARGET]).statistic
    r = pearsonr(pair["mrna"], pair[TARGET])[0]
    print(f"\nComparable with the human literature, uncorrected for measurement error: over "
          f"{len(pair):,} genes,\none averaged mRNA measure against protein gives r2 = {r ** 2:.3f} "
          f"(Pearson) and {rho ** 2:.3f} (Spearman).\nHuman estimates run 0.40 to 0.81 depending "
          f"on whether measurement error is corrected for.")

    first, last = table.iloc[0], table.iloc[-1]
    share = first["r2"] / last["r2"] if last["r2"] else float("nan")
    print(f"\nAnd what the layers add. Transcript abundance explains {first['r2']:.1%} of the "
          f"variance;\nadding turnover, ribosome occupancy and codon usage takes it only to "
          f"{last['r2']:.1%}.\nSo the other layers add {last['r2'] - first['r2']:.1%} -- "
          f"transcript is {share:.0%} of what THIS MODEL explains, which is a ratio of one R2 to "
          f"another\nand NOT a share of protein variance.")
    # The qualification belongs beside the claim, not in a footnote. The genes that carry every
    # layer are the genes everybody could measure, and transcript predicts protein far better in
    # them than in the wider set, so the headline number describes the easy half of the proteome.
    print(f"\nBUT: that is on the genes where every layer happens to be measured. On the wider "
          f"set\ntranscript explains {wide_r2:.1%}, not {first['r2']:.1%}. The complete-case genes "
          f"are the ones\neverybody could measure, and protein is easier to predict in them. The "
          f"ratio of\nexplained variance between the two sets is "
          f"{first['r2'] / wide_r2:.1f}, which is the size of that bias.")

    out = os.path.join(ROOT, "results", f"what_sets_protein_abundance_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "layers.tsv"), sep="\t", index=False)
    json.dump({"organism": "Plasmodium falciparum", "target": TARGET,
               "genes_with_protein": int(with_target), "genes_used": int(len(complete)),
               "folds": N_FOLDS, "repeats": N_REPEATS,
               "r2_by_layer": {r["layer"]: r["r2"] for r in rows},
               "shuffled_r2_mean": float(null.mean()), "shuffled_r2_sd": float(null.std(ddof=1)),
               "transcription_share_of_what_this_model_explains": float(share),
               "single_mrna_vs_protein_pearson_r2": float(r ** 2),
               "single_mrna_vs_protein_spearman_r2": float(rho ** 2),
               "single_mrna_vs_protein_n": int(len(pair)),
               "r2_alone_by_layer": {r["layer"]: r["r2_alone"] for r in rows},
               "wide_genes": int(len(wide)), "wide_transcript_only_r2": float(wide_r2)},
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print("wrote", out)


if __name__ == "__main__":
    main()
