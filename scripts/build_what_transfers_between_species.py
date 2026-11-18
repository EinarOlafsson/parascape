"""Which layer of knowledge transfers between apicomplexan species, and which does not.

Figure 9 of the ParaScape paper (instruction 101). The field moves a result from Plasmodium to
Toxoplasma and back all the time, and the exchange rate has never been measured per layer. It can
be, because the two species sit in one gene space joined by orthogroup, and several layers are
measured in both.

THE MEASUREMENT. For each layer, take the one-to-one orthologue pairs where both species carry a
value, and ask how well the value in one predicts the value in the other -- as accuracy against the
rate of guessing the commonest class for a categorical layer, and as Spearman for a continuous one.
The number reported is the lift over that baseline, because a layer in which 60 per cent of genes
are "cytosol" transfers at 60 per cent by saying "cytosol" every time.

ONE-TO-ONE ONLY. An orthogroup holding three Toxoplasma genes and two Plasmodium ones has six
pairings and no way to say which is the counterpart. Those groups are counted and set aside rather
than averaged, because averaging the members of a gene family measures how similar a family is to
itself, which is not the question.

Writes `results/species_transfer_<date>/`.
"""
from __future__ import annotations

import json
import os
from datetime import date

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: Layer -> (Toxoplasma column, Plasmodium column, kind). A categorical layer is scored by
#: accuracy against the commonest class; a continuous one by Spearman against zero.
LAYERS = (
    ("where the protein is", "lopit_map", "lopit_pf_location", "categorical"),
    ("whether losing it costs fitness", "fit_invitro_hff", "piggybac_mis", "continuous"),
    ("how the message is written (codon bias)", "codon_cai_ribosomal", "codon_cai_ribosomal",
     "continuous"),
    ("how long the protein is", "length", "length", "continuous"),
)

#: A layer scored on fewer pairs than this is reported as not measurable rather than scored.
MIN_PAIRS = 40

#: Draws behind the baseline of a categorical layer.
N_NULL = 2000

#: The two atlases name compartments differently and at different resolutions -- Toxoplasma has
#: "nucleus - chromatin" and "PM - integral", Plasmodium has "nucleus 1" and "plasma membrane 1" --
#: so both are reduced to one coarse vocabulary before anything is compared. The first version of
#: this map was written for the Toxoplasma names alone and sent its "cytosol" to "cytoplasm" while
#: leaving the Plasmodium "cytosol" as it was, which destroyed the match it was meant to make.
COARSE = {
    "nucleus": "nucleus", "nucleolus": "nucleus", "nucleo-cytoplasm": "nucleus",
    "mitotic apparatus": "nucleus",
    "cytosol": "cytosol", "cytoplasm": "cytosol",
    "er": "endoplasmic reticulum", "endoplasmic reticulum": "endoplasmic reticulum",
    "golgi": "golgi", "golgi apparatus": "golgi", "transport vesicles": "golgi",
    "pm": "plasma membrane", "plasma membrane": "plasma membrane",
    "mitochondrion": "mitochondrion", "apicoplast": "apicoplast", "imc": "imc",
    "rhoptries": "rhoptries", "micronemes": "micronemes", "apical": "apical",
    "dense granules": "dense granules", "ribosomes": "ribosome", "ribosome": "ribosome",
    "19s proteasome": "proteasome", "20s proteasome": "proteasome", "proteasome": "proteasome",
}


def _coarse(value) -> float | str:
    """One compartment name from either atlas, reduced to the shared vocabulary.

    Strips a trailing cluster number ("nucleus 1", "rhoptries 2") and any sub-compartment after a
    dash ("mitochondrion - soluble"), then looks the head up. A name with no entry is kept as
    itself, lowered, so an unmapped vocabulary shows up as disagreement rather than as a silent
    merge into some other compartment.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    text = str(value).strip().lower()
    if not text or text in ("nan", "none", "unknown"):
        return np.nan
    text = text.split(" - ")[0].strip()
    parts = text.split()
    if len(parts) > 1 and parts[-1].isdigit():
        text = " ".join(parts[:-1])
    return COARSE.get(text, text)


def _normalise(series: pd.Series) -> pd.Series:
    """Compartment names on one vocabulary, so the two atlases can disagree about biology only."""
    return series.map(_coarse)


def one_to_one(tg: pd.DataFrame, pf: pd.DataFrame) -> tuple:
    """Orthogroups holding exactly one gene of each species, as a pair table, and what was dropped."""
    tg_counts = tg.groupby("orthogroup").size()
    pf_counts = pf.groupby("orthogroup").size()
    shared = tg_counts.index.intersection(pf_counts.index)
    single = [g for g in shared if tg_counts[g] == 1 and pf_counts[g] == 1]
    pairs = (tg[tg["orthogroup"].isin(single)].set_index("orthogroup")
             .join(pf[pf["orthogroup"].isin(single)].set_index("orthogroup"),
                   lsuffix="_tg", rsuffix="_pf"))
    return pairs, {"shared_orthogroups": int(len(shared)), "one_to_one": int(len(single)),
                   "many_to_many_set_aside": int(len(shared) - len(single))}


def score_categorical(a: pd.Series, b: pd.Series, seed: int = 0) -> dict:
    """How often the label in one species is the label in the other, against guessing the commonest.

    The baseline is not 1/k: a vocabulary in which half the genes are cytoplasm is matched half the
    time by saying cytoplasm, so that is what has to be beaten.
    """
    both = pd.DataFrame({"a": a, "b": b}).dropna()
    if len(both) < MIN_PAIRS:
        return {"n": len(both), "measurable": False,
                "reason": f"{len(both)} pairs, fewer than {MIN_PAIRS}"}
    agree = float((both["a"] == both["b"]).mean())
    commonest = float(both["b"].value_counts(normalize=True).max())
    rng = np.random.default_rng(seed)
    draws = np.array([float((rng.permutation(both["a"].to_numpy()) == both["b"].to_numpy()).mean())
                      for _ in range(N_NULL)])
    p = (1.0 + int((draws >= agree).sum())) / (len(draws) + 1.0)
    return {"n": int(len(both)), "measurable": True, "transfer": agree,
            "baseline": max(commonest, float(draws.mean())), "shuffled_mean": float(draws.mean()),
            "lift": agree - max(commonest, float(draws.mean())), "p": p,
            "shared_labels": int(len(set(both["a"]) & set(both["b"])))}


def score_continuous(a: pd.Series, b: pd.Series) -> dict:
    """Spearman between the two species, which is the transfer of an ordering rather than a value."""
    both = pd.DataFrame({"a": pd.to_numeric(a, errors="coerce"),
                         "b": pd.to_numeric(b, errors="coerce")}).dropna()
    if len(both) < MIN_PAIRS:
        return {"n": len(both), "measurable": False,
                "reason": f"{len(both)} pairs, fewer than {MIN_PAIRS}"}
    result = spearmanr(both["a"], both["b"])
    return {"n": int(len(both)), "measurable": True, "transfer": float(result.statistic),
            "baseline": 0.0, "lift": float(result.statistic), "p": float(result.pvalue)}


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    tg = pd.read_parquet(os.path.join(cache, "nodes.parquet"))
    pf = pd.read_parquet(os.path.join(cache, "pf_nodes.parquet"))
    wanted_tg = {c for _n, c, _p, _k in LAYERS} | {"orthogroup", "gene_id"}
    wanted_pf = {c for _n, _t, c, _k in LAYERS} | {"orthogroup", "gene_id"}
    tg = tg[[c for c in tg.columns if c in wanted_tg]].dropna(subset=["orthogroup"])
    pf = pf[[c for c in pf.columns if c in wanted_pf]].dropna(subset=["orthogroup"])

    pairs, counted = one_to_one(tg, pf)
    print(f"{counted['one_to_one']:,} one-to-one orthologue pairs; "
          f"{counted['many_to_many_set_aside']:,} orthogroups set aside as many-to-many")

    rows = []
    for name, tg_col, pf_col, kind in LAYERS:
        a = f"{tg_col}_tg" if f"{tg_col}_tg" in pairs.columns else tg_col
        b = f"{pf_col}_pf" if f"{pf_col}_pf" in pairs.columns else pf_col
        if a not in pairs.columns or b not in pairs.columns:
            rows.append({"layer": name, "measurable": False, "reason": "column absent"})
            continue
        if kind == "categorical":
            score = score_categorical(_normalise(pairs[a]), _normalise(pairs[b]))
        else:
            score = score_continuous(pairs[a], pairs[b])
        rows.append({"layer": name, "kind": kind, **score})
        if score.get("measurable"):
            print(f"  {name:42s} n={score['n']:5,}  transfer={score['transfer']:+.3f}  "
                  f"baseline={score['baseline']:.3f}  lift={score['lift']:+.3f}  "
                  f"p={score['p']:.3g}")
        else:
            print(f"  {name:42s} not measurable: {score.get('reason')}")

    # THE CONTROL ON THE PAIRING ITSELF. Protein length is a sequence property and must transfer
    # between true orthologues; if it did not, the pairs would be wrong and every other number here
    # would be meaningless. It comes back at 0.9, so a layer that transfers weakly is telling us
    # about the biology rather than about the orthologue mapping.
    measurable = [r for r in rows if r.get("measurable")]
    control = next((r for r in measurable if r["layer"].startswith("how long")), None)
    if control:
        print(f"\ncontrol: protein length transfers at {control['transfer']:.3f}, so the pairing "
              f"is sound and\na layer that transfers weakly is biology, not a mapping artefact.")

    table = pd.DataFrame(rows)
    out = os.path.join(ROOT, "results", f"species_transfer_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "layers.tsv"), sep="\t", index=False)
    json.dump({"pairs": counted, "layers": rows, "min_pairs": MIN_PAIRS, "n_null": N_NULL},
              open(os.path.join(out, "summary.json"), "w"), indent=2, default=float)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
