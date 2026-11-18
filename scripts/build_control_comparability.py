"""Which studies can be compared at all, and what normalisation makes more of them comparable.

Figure 2 of the ParaScape paper (instruction 101). Einar's design, 2026-11-18:

    "cluster the studies by how well their controls (non treated, wildtype, etc.) match onto each
    other. some controls in their current state might not agree, but if we say add some
    normalization steps that make biological sense, like normalize for cell cycle stage, they might
    agree more, the conditions in these clusters will be a fairer comparison than studies where the
    controls dont agree."

THE ARGUMENT. A perturbation measured in study A is only comparable with one measured in study B if
the two untreated arms describe the same starting parasite. Nobody checks this, because no single
laboratory can: it needs every study on one scale. So the first result is a map of which published
studies may be set against each other, and the second is which biological variable, once removed,
moves studies into agreement.

WHAT IS CONTROLLED, AND WHY IT HAS TO BE. Removing any covariate estimated FROM the data raises
agreement whether or not the covariate means anything, so the cell-cycle result is reported against
a null in which the same procedure is run on random gene sets of the same size and shape. The
number that matters is not how much agreement rises; it is how much it rises ABOVE the null.

Writes `results/control_comparability_<date>/`.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import expression_contrast as EC   # noqa: E402

#: The ladder, which names the baseline (control) group of every contrast it could parse.
LADDER = os.path.join(ROOT, "results", "contrast_ladder_2026_11_06", "contrasts.tsv")

#: Cell-cycle phase per gene, from Xue et al. 2020: 964 genes over G1a, G1b, S, M and C.
PHASE_FILE = "cellcycle_phase"

#: A study needs this many control runs to give a profile worth correlating.
MIN_CONTROL_RUNS = 2

#: Genes measured in fewer than this many studies are dropped before anything is correlated.
MIN_STUDIES_PER_GENE = 0.8

#: How many random gene sets stand behind the cell-cycle correction.
N_NULL = 200


def _tail(label: str) -> str:
    """The group name at the end of a ladder label: 'treatment: none | Pru' -> 'pru'."""
    return str(label).rsplit("|", 1)[-1].strip().lower()


def control_groups_per_study(path: str = LADDER) -> dict:
    """Study -> the set of condition names its contrasts treat as the baseline.

    Read from the ladder rather than guessed from words like 'wild type', because the ladder
    baselines were reviewed against the GEO text, and 'untreated' means different things in
    different studies. A study whose baseline the ladder could not parse is simply absent.
    """
    ladder = pd.read_csv(path, sep="\t")
    baseline, treated = {}, {}
    for row in ladder.itertuples():
        for m in str(row.baseline_members).split("||"):
            baseline.setdefault(row.study, set()).add(_tail(m))
        for m in str(row.treated_members).split("||"):
            treated.setdefault(row.study, set()).add(_tail(m))
    # A GROUP THAT IS EVER THE TREATED ARM IS NOT A CONTROL. The ladder pairs groups on one axis
    # at a time, so a drug-treated group is the baseline of the contrast against a higher dose and
    # the treated arm of the contrast against vehicle. Taking the union of baselines made every
    # group of GSE103113 a "control", including gra18 and pru_lps, which would have compared
    # treated parasites between studies and called the result a control profile.
    return {s: b - treated.get(s, set()) for s, b in baseline.items() if b - treated.get(s, set())}


def control_profile(study: str, controls: set, root: str, quantifier: str = "salmon",
                    reference: str = "ME49") -> tuple:
    """The mean expression of the control runs of one study, as log2 TPM, and how many runs.

    The parasite-only reference is used, not ME49_host: a sample that is 95 per cent host has far
    fewer parasite reads than its depth suggests, and a profile built on the combined index would
    be comparing host content between studies rather than parasite expression.
    """
    try:
        _counts, tpm, _lengths, samples = EC.load(study, quantifier=quantifier,
                                                  reference=reference, root=root)
    except (FileNotFoundError, OSError, ValueError, KeyError):
        return None, 0
    if "condition" not in samples.columns:
        return None, 0
    keep = samples[samples["condition"].astype(str).str.strip().str.lower().isin(controls)]
    runs = [r for r in keep["run_accession"].astype(str) if r in tpm.columns]
    if len(runs) < MIN_CONTROL_RUNS:
        return None, len(runs)
    mean = tpm[runs].mean(axis=1)
    return np.log2(mean + 1.0), len(runs)


def build_profiles(root: str) -> tuple:
    """One control profile per study, as a gene x study frame, and a record of what was used."""
    controls = control_groups_per_study()
    available = EC.available(root=root)
    have = set(available["study"].astype(str))
    profiles, record = {}, []
    for study in sorted(controls):
        if study not in have:
            record.append({"study": study, "used": False, "reason": "no built matrix"})
            continue
        profile, n = control_profile(study, controls[study], root)
        if profile is None:
            record.append({"study": study, "used": False, "n_control_runs": n,
                           "reason": f"{n} control runs, fewer than {MIN_CONTROL_RUNS}"})
            continue
        profiles[study] = profile
        record.append({"study": study, "used": True, "n_control_runs": n,
                       "controls": "; ".join(sorted(controls[study])), "reason": ""})
    frame = pd.DataFrame(profiles)
    return frame, pd.DataFrame(record)


def shared_genes(frame: pd.DataFrame) -> pd.DataFrame:
    """Genes measured in most studies, each profile centred on its own median.

    Centring is what makes two studies comparable at all: an absolute log2 TPM carries the depth
    and the library of its own experiment. What is compared is the SHAPE of the profile.
    """
    enough = frame.notna().mean(axis=1) >= MIN_STUDIES_PER_GENE
    out = frame.loc[enough].copy()
    return out.sub(out.median(axis=0), axis=1)


def agreement(frame: pd.DataFrame) -> pd.DataFrame:
    """Spearman between every pair of control profiles, over the genes both measure."""
    studies = list(frame.columns)
    out = pd.DataFrame(np.eye(len(studies)), index=studies, columns=studies)
    for i, a in enumerate(studies):
        for b in studies[i + 1:]:
            both = frame[[a, b]].dropna()
            rho = spearmanr(both[a], both[b]).statistic if len(both) > 200 else np.nan
            out.loc[a, b] = out.loc[b, a] = rho
    return out


def cluster(scores: pd.DataFrame, threshold: float = 0.7) -> pd.Series:
    """Group the studies whose controls agree, by average linkage on 1 - Spearman.

    `threshold` is a correlation, so 0.7 means "studies whose controls agree at least this well end
    up together". It is stated here rather than tuned to make a pretty picture.
    """
    distance = (1.0 - scores).to_numpy(dtype=float).copy()
    np.fill_diagonal(distance, 0.0)
    distance = np.nan_to_num(distance, nan=1.0)
    distance = (distance + distance.T) / 2.0
    link = hierarchy.linkage(squareform(distance, checks=False), method="average")
    labels = hierarchy.fcluster(link, t=1.0 - threshold, criterion="distance")
    return pd.Series(labels, index=scores.index, name="cluster")


def phase_table(cache: str) -> pd.Series:
    """Gene -> cell-cycle phase, from the shipped node table."""
    path = os.path.join(cache, "nodes.parquet")
    if not os.path.exists(path):
        for name in ("toxoplasma_nodes.parquet", "tg_nodes.parquet"):
            if os.path.exists(os.path.join(cache, name)):
                path = os.path.join(cache, name)
                break
    nodes = pd.read_parquet(path)
    column = next((c for c in nodes.columns if "cellcycle_phase" in c), None)
    if column is None:
        raise KeyError("no cell-cycle phase column in the node table")
    out = nodes.set_index(nodes["gene_id"].astype(str))[column]
    return out.dropna().astype(str)


def _covariate(frame: pd.DataFrame, groups: dict) -> np.ndarray:
    """One number per study: where its controls sit on the axis these gene groups define.

    The groups are the five cell-cycle phases. Each study gets the mean profile value of each
    group, and the studies are placed on the first principal component of those five numbers --
    the single axis along which studies differ most in the balance of their phase genes.
    """
    rows = []
    for _name, genes in groups.items():
        members = frame.index.intersection(pd.Index(genes))
        rows.append(frame.loc[members].mean(axis=0).to_numpy() if len(members) else
                    np.zeros(frame.shape[1]))
    block = np.vstack(rows)                              # groups x studies
    block = block - block.mean(axis=1, keepdims=True)
    _u, _s, vt = np.linalg.svd(block, full_matrices=False)
    axis = vt[0]
    return (axis - axis.mean()) / (axis.std() or 1.0)


def _regress_out(frame: pd.DataFrame, covariate: np.ndarray) -> pd.DataFrame:
    """Remove, from every gene, the part of its across-study variation the covariate explains."""
    x = np.vstack([np.ones_like(covariate), covariate]).T
    values = frame.to_numpy(dtype=float)
    out = values.copy()
    ok = ~np.isnan(values)
    for i in range(values.shape[0]):
        m = ok[i]
        if m.sum() < 4:
            continue
        beta, *_ = np.linalg.lstsq(x[m], values[i, m], rcond=None)
        out[i, m] = values[i, m] - x[m] @ beta + beta[0]
    return pd.DataFrame(out, index=frame.index, columns=frame.columns)


def correct_for(frame: pd.DataFrame, groups: dict) -> pd.DataFrame:
    """The profiles with the axis those gene groups define removed."""
    return _regress_out(frame, _covariate(frame, groups))


def median_offdiagonal(scores: pd.DataFrame) -> float:
    """The median agreement between different studies, which is the number the figure reports."""
    values = scores.to_numpy(dtype=float).copy()
    np.fill_diagonal(values, np.nan)
    return float(np.nanmedian(values))


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE to the release whose rnaseq matrices are to be read")
    root = os.path.join(cache, "rnaseq")

    print("reading control runs...")
    raw, record = build_profiles(root)
    used = record[record["used"]]
    print(f"{len(used)} studies with a control profile, of {len(record)} the ladder names")
    if raw.shape[1] < 4:
        raise SystemExit("too few studies with control profiles to cluster")
    frame = shared_genes(raw)
    print(f"{len(frame)} genes measured in at least {MIN_STUDIES_PER_GENE:.0%} of them")

    before = agreement(frame)
    phases = phase_table(cache)
    groups = {p: list(g.index) for p, g in phases.groupby(phases)}
    n_phase_genes = sum(len(frame.index.intersection(pd.Index(g))) for g in groups.values())
    print(f"cell-cycle phases: {len(groups)} groups, {n_phase_genes} of their genes in the matrix")
    after = agreement(correct_for(frame, groups))

    # THE CONTROL ON THE TEST ITSELF. Removing any data-derived axis raises agreement; the question
    # is whether the cell-cycle one raises it more than a random axis of the same shape would.
    rng = np.random.default_rng(0)
    pool = list(frame.index)
    sizes = [len(frame.index.intersection(pd.Index(g))) for g in groups.values()]
    null = []
    for _ in range(N_NULL):
        drawn, taken = {}, set()
        for i, size in enumerate(sizes):
            pick = [g for g in rng.choice(pool, size=min(size, len(pool)), replace=False)
                    if g not in taken]
            taken.update(pick)
            drawn[f"random_{i}"] = pick
        null.append(median_offdiagonal(agreement(correct_for(frame, drawn))))
    null = np.array(null)

    base, corrected = median_offdiagonal(before), median_offdiagonal(after)
    gain = corrected - base
    null_gain = null - base
    p = (1.0 + int((null_gain >= gain).sum())) / (len(null_gain) + 1.0)
    print(f"\nmedian agreement between studies: {base:.3f} -> {corrected:.3f} "
          f"(gain {gain:+.3f})")
    print(f"a random axis of the same shape gains {null_gain.mean():+.3f} "
          f"(sd {null_gain.std(ddof=1):.3f}); p = {p:.4f}")

    groups_before, groups_after = cluster(before), cluster(after)
    moved = [s for s in before.index
             if (before.loc[s].drop(s).median() < 0.7) <= (after.loc[s].drop(s).median() >= 0.7)
             and after.loc[s].drop(s).median() >= 0.7 > before.loc[s].drop(s).median()]

    out = os.path.join(ROOT, "results", f"control_comparability_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    before.to_csv(os.path.join(out, "agreement_before.tsv"), sep="\t")
    after.to_csv(os.path.join(out, "agreement_after_cellcycle.tsv"), sep="\t")
    record.to_csv(os.path.join(out, "studies.tsv"), sep="\t", index=False)
    pd.DataFrame({"cluster_before": groups_before,
                  "cluster_after": groups_after,
                  "median_agreement_before": before.apply(lambda c: c.drop(c.name).median()),
                  "median_agreement_after": after.apply(lambda c: c.drop(c.name).median())}
                 ).to_csv(os.path.join(out, "clusters.tsv"), sep="\t")
    json.dump({"studies": int(before.shape[0]), "genes": int(len(frame)),
               "median_agreement_before": base, "median_agreement_after": corrected,
               "gain": gain, "null_gain_mean": float(null_gain.mean()),
               "null_gain_sd": float(null_gain.std(ddof=1)), "p": p, "n_null": N_NULL,
               "clusters_before": int(groups_before.nunique()),
               "clusters_after": int(groups_after.nunique()),
               "studies_crossing_0_7": moved},
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print(f"\nclusters: {groups_before.nunique()} -> {groups_after.nunique()}")
    print("wrote", out)


if __name__ == "__main__":
    main()
