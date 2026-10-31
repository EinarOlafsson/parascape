"""Compare expression between two groups of samples: strains, knockouts, drugs, stages (instruction 69).

Everything a user may want to choose is a parameter, because the conventions differ and each is defended
in the literature:

* **quantifier** -- ``salmon`` (selective alignment with genome decoys, GC and sequence bias corrected)
  or ``kallisto`` (pseudoalignment);
* **reference** -- ``ME49`` for every study (shared gene identifiers), or a strain's own genome (``GT1``
  Type I, ``VEG`` Type III) where reads came from that strain: mapping Type I reads to ME49 penalises the
  genes that differ most between strains, which is exactly what a strain comparison is about;
* **normalisation** -- ``tmm`` (Robinson & Oshlack 2010), ``median_ratio`` (Anders & Huber 2010),
  ``upper_quartile`` (Bullard et al. 2010), ``cpm`` (library size only) or ``tpm``;
* **within-study scale** -- every sample against its own study's control group, so studies run on
  different machines, depths and protocols can be read on one scale (:func:`control_relative`).

The test is **limma-trend** (Law et al. 2014): log2-CPM with a prior count, a two-group linear model, a
mean-variance trend, and empirical-Bayes moderation of each gene's variance towards the trend (Smyth
2004), giving a moderated t with residual + prior degrees of freedom and Benjamini-Hochberg FDR.
Implemented here in numpy/scipy because neither R nor DESeq2 is available; checked by
`tests/test_expression_contrast.py` against simulations (null p-values uniform, FDR controlled) before
any result is shipped.
"""
from __future__ import annotations

from . import paths as _data_paths

import gzip
import json
import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import special, stats

NAS = os.environ.get("PARASCAPE_OMICS", os.path.join(_data_paths.data_dir(), "rnaseq"))
QUANTIFIERS = ("salmon", "kallisto")
#: References a study may be quantified against; strain references only for strain studies.
REFERENCES = ("ME49", "ME49_host", "GT1", "VEG")
#: Gene identifiers the competitive reference gives host transcripts (`scripts/omics_quant.py`).
HOST_PREFIX = "HOST:"
NORMALISATIONS = ("tmm", "median_ratio", "upper_quartile", "cpm", "tpm")


# --------------------------------------------------------------------------- loading
def tx2gene(reference: str = "ME49", root: str = NAS) -> pd.Series:
    """Transcript -> gene for a reference, from its cDNA headers (host genes prefixed HOST:)."""
    t = pd.read_csv(os.path.join(root, "references", reference, "tx2gene.tsv"), sep="\t",
                    header=None, names=["tx", "gene"])
    return t.set_index("tx")["gene"]


def read_run(path: str, quantifier: str) -> pd.DataFrame:
    """One run's transcript table as (tx, count, tpm, length)."""
    if quantifier == "salmon":
        f = pd.read_csv(os.path.join(path, "quant.sf.gz"), sep="\t")
        return pd.DataFrame({"tx": f["Name"], "count": f["NumReads"], "tpm": f["TPM"],
                             "length": f["EffectiveLength"]})
    f = pd.read_csv(os.path.join(path, "abundance.tsv.gz"), sep="\t")
    return pd.DataFrame({"tx": f["target_id"], "count": f["est_counts"], "tpm": f["tpm"],
                         "length": f["eff_length"]})


def study_runs(study: str, root: str = NAS) -> pd.DataFrame:
    """The quantified runs of one study, with their sample metadata (from ENA, recorded at quantification)."""
    d = os.path.join(root, "rnaseq", study)
    rows = []
    for run in sorted(os.listdir(d)) if os.path.isdir(d) else []:
        done = os.path.join(d, run, "done.json")
        if os.path.exists(done):
            with open(done) as fh:
                rows.append(json.load(fh))
    out = pd.DataFrame(rows)
    if len(out):
        out["condition"] = out.apply(_condition, axis=1)
    return out


def _condition(r) -> str:
    """A first guess at the condition, from the sample title with its replicate suffix removed.

    Only a default: `contrast` takes the groups explicitly, and the UI lets a person regroup samples,
    because sample titles are free text and no rule reads every lab's naming correctly.
    """
    import re
    t = str(r.get("sample_title") or r.get("sample_alias") or r.get("run_accession"))
    # A replicate number is set off by a separator or named as one ("_1", " rep2", "replicate 3");
    # a number glued to a word is part of the name -- "gra18" is the GRA18 knockout, not "gra", rep 18.
    t = re.sub(r"(?i)(?:[\s_\-]+(?:biological\s+)?(?:rep(?:licate)?|br|r)?[\s_\-]*|"
               r"(?:biological\s+)?(?:rep(?:licate)?)[\s_\-]*)\d+$", "", t).strip()
    return t or str(r.get("run_accession"))


def matrices(study: str, quantifier: str = "salmon", reference: str = "ME49",
             root: str = NAS) -> tuple:
    """(counts, tpm, lengths, samples): gene x run, summed over each gene's transcripts."""
    if quantifier not in QUANTIFIERS:
        raise ValueError(f"quantifier must be one of {QUANTIFIERS}")
    samples = study_runs(study, root)
    if not len(samples):
        raise ValueError(f"no quantified runs for {study}")
    t2g = tx2gene(reference, root)
    counts, tpm, length = {}, {}, {}
    for run in samples["run_accession"]:
        path = os.path.join(root, "rnaseq", study, run, f"{quantifier}_{reference}")
        if not os.path.isdir(path):
            continue
        f = read_run(path, quantifier)
        f["gene"] = f["tx"].map(t2g).fillna(f["tx"])
        g = f.groupby("gene")
        counts[run] = g["count"].sum()
        tpm[run] = g["tpm"].sum()
        # A gene's effective length: its transcripts' lengths weighted by their abundance.
        w = f.assign(w=f["tpm"] * f["length"]).groupby("gene")["w"].sum() / g["tpm"].sum().replace(0, np.nan)
        length[run] = w.fillna(g["length"].mean())
    samples = samples[samples["run_accession"].isin(counts)].reset_index(drop=True)
    counts, tpm, length = (pd.DataFrame(counts).fillna(0.0), pd.DataFrame(tpm).fillna(0.0),
                           pd.DataFrame(length))
    host = counts.index.astype(str).str.startswith(HOST_PREFIX)
    if host.any():
        # The competitive reference: host genes are what host reads were assigned to. They are dropped
        # from the parasite matrices, and their share of each sample is kept as a quality figure --
        # a sample that is 95% host has far fewer parasite reads than its depth suggests.
        share = counts[host].sum() / counts.sum().replace(0, np.nan)
        samples["host_share"] = samples["run_accession"].map(share).to_numpy()
        counts, length = counts[~host], length[~host]
        tpm = tpm[~host]
        tpm = tpm / tpm.sum() * 1e6                 # TPM re-closed over parasite genes only
    return counts, tpm, length, samples



def matrix_dir(root: str = NAS) -> str:
    """Where built gene-level matrices live under an omics root."""
    return os.path.join(root, "matrices")


def build_matrices(study: str, root: str = NAS) -> list:
    """Write one parquet per (quantifier, reference) for a study: counts, TPM and lengths, gene x run,
    with the sample sheet beside them. These -- not the per-run files -- are what the pack ships."""
    out = []
    os.makedirs(os.path.join(matrix_dir(root), study), exist_ok=True)
    for q in QUANTIFIERS:
        for ref in REFERENCES:
            try:
                counts, tpm, length, samples = matrices(study, q, ref, root)
            except (ValueError, FileNotFoundError):
                continue
            if counts.shape[1] == 0:
                continue
            base = os.path.join(matrix_dir(root), study, f"{q}_{ref}")
            pd.concat({"count": counts, "tpm": tpm, "length": length}, axis=1).to_parquet(
                base + ".parquet", compression="zstd")
            samples.to_csv(base + "_samples.tsv", sep="\t", index=False)
            out.append(base)
    return out


def available(root: str = NAS) -> pd.DataFrame:
    """Every study with built matrices, and the conventions each was quantified under."""
    d = matrix_dir(root)
    rows = []
    for study in sorted(os.listdir(d)) if os.path.isdir(d) else []:
        for f in sorted(os.listdir(os.path.join(d, study))):
            if f.endswith(".parquet"):
                q, ref = f[:-len(".parquet")].split("_", 1)
                rows.append({"study": study, "quantifier": q, "reference": ref})
    return pd.DataFrame(rows, columns=["study", "quantifier", "reference"])


def load(study: str, quantifier: str = "salmon", reference: str = "ME49_host",
         root: str = NAS) -> tuple:
    """Read counts, TPM, lengths and samples; absent compact-pack lengths remain NaN."""
    base = os.path.join(matrix_dir(root), study, f"{quantifier}_{reference}")
    if not os.path.exists(base + ".parquet"):
        raise FileNotFoundError(f"{study} has no {quantifier}/{reference} matrix in {matrix_dir(root)}")
    m = pd.read_parquet(base + ".parquet")
    samples = pd.read_csv(base + "_samples.tsv", sep="\t")
    length = m["length"] if "length" in m.columns.get_level_values(0) else m["tpm"] * np.nan
    return m["count"], m["tpm"], length, samples

# --------------------------------------------------------------------------- normalisation
def norm_factors(counts: pd.DataFrame, method: str = "tmm") -> pd.Series:
    """Per-sample scaling relative to library size (1.0 = library size alone), as edgeR defines them."""
    X = counts.to_numpy(float)
    lib = X.sum(axis=0)
    if method in ("cpm", "tpm"):
        return pd.Series(1.0, index=counts.columns)
    if method == "upper_quartile":
        q = np.array([np.percentile(c[c > 0], 75) if (c > 0).any() else np.nan for c in X.T]) / lib
        f = q / np.exp(np.nanmean(np.log(q)))
    elif method == "median_ratio":
        # DESeq2 size factors, expressed as normalisation factors on top of library size.
        keep = (X > 0).all(axis=1)
        logs = np.log(X[keep])
        geo = logs.mean(axis=1, keepdims=True)
        sf = np.exp(np.median(logs - geo, axis=0))
        f = sf / lib
        f = f / np.exp(np.mean(np.log(f)))
    elif method == "tmm":
        f = _tmm(X, lib)
    else:
        raise ValueError(f"normalisation must be one of {NORMALISATIONS}")
    return pd.Series(f, index=counts.columns)


def _tmm(X, lib, logratio_trim: float = 0.3, sum_trim: float = 0.05) -> np.ndarray:
    """edgeR's calcNormFactors(method="TMM"), reference = sample whose upper quartile is closest to mean."""
    uq = np.array([np.percentile(c, 75) for c in X.T]) / lib
    ref = int(np.argmin(np.abs(uq - uq.mean())))
    f = np.ones(X.shape[1])
    for j in range(X.shape[1]):
        obs, rf, nO, nR = X[:, j], X[:, ref], lib[j], lib[ref]
        ok = (obs > 0) & (rf > 0)
        logR = np.log2((obs[ok] / nO) / (rf[ok] / nR))
        absE = (np.log2(obs[ok] / nO) + np.log2(rf[ok] / nR)) / 2
        v = (nO - obs[ok]) / nO / obs[ok] + (nR - rf[ok]) / nR / rf[ok]
        n = len(logR)
        if n == 0:
            continue
        lo_l, hi_l = np.floor(n * logratio_trim) + 1, n + 1 - np.floor(n * logratio_trim)
        lo_s, hi_s = np.floor(n * sum_trim) + 1, n + 1 - np.floor(n * sum_trim)
        rl, rs = stats.rankdata(logR), stats.rankdata(absE)
        keep = (rl >= lo_l) & (rl <= hi_l) & (rs >= lo_s) & (rs <= hi_s)
        f[j] = 2 ** (np.sum(logR[keep] / v[keep]) / np.sum(1 / v[keep])) if keep.any() else 1.0
    return f / np.exp(np.mean(np.log(f)))


def log_cpm(counts: pd.DataFrame, factors: pd.Series, prior_count: float = 2.0) -> pd.DataFrame:
    """log2 counts per million on effective library sizes, with edgeR's library-scaled prior count."""
    lib = counts.sum(axis=0) * factors
    prior = prior_count * lib / lib.mean()
    return np.log2((counts + prior) / (lib + 2 * prior) * 1e6)


def normalised(counts: pd.DataFrame, tpm: pd.DataFrame | None = None,
               method: str = "tmm") -> pd.DataFrame:
    """log2 expression on the chosen scale -- the matrix a profile comparison reads."""
    if method == "tpm":
        if tpm is None:
            raise ValueError("TPM needs the TPM matrix")
        return np.log2(tpm + 1.0)
    return log_cpm(counts, norm_factors(counts, method))


def control_relative(expression: pd.DataFrame, control: list) -> pd.DataFrame:
    """Every sample against the mean of its own study's control samples, in log2.

    The scale on which studies from different labs, machines and depths can be laid side by side: what
    is compared is each perturbation's effect, never absolute levels measured under different protocols.
    """
    ctrl = [c for c in control if c in expression.columns]
    if not ctrl:
        raise ValueError("no control samples in this matrix")
    return expression.sub(expression[ctrl].mean(axis=1), axis=0)


# --------------------------------------------------------------------------- the test
@dataclass
class Contrast:
    """The result of comparing group B with group A (positive log2FC = higher in B)."""
    table: pd.DataFrame          # gene, log2fc, mean_expression, t, p, fdr, ...
    group_a: list
    group_b: list
    quantifier: str
    reference: str
    normalisation: str
    prior_df: float
    min_expression: float
    processing_record: dict = field(default_factory=dict)

    def significant(self, fdr: float = 0.05, min_abs_log2fc: float = 1.0) -> pd.DataFrame:
        """Genes at or below `fdr` with an absolute log2 fold change of at least `min_abs_log2fc`."""
        t = self.table
        return t[(t["fdr"] <= fdr) & (t["log2fc"].abs() >= min_abs_log2fc)]


def _trigamma_inverse(x: np.ndarray) -> np.ndarray:
    """limma's trigammaInverse: Newton iteration on 1/trigamma."""
    x = np.asarray(x, float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.all(np.abs(dif / y) < 1e-8):
            break
    return y


def _trend(a: np.ndarray, logs2: np.ndarray, span: float = 0.3) -> np.ndarray:
    """A smooth trend of log variance on mean expression: local linear regression (tricube weights)."""
    order = np.argsort(a)
    xs, ys = a[order], logs2[order]
    n = len(xs)
    k = max(int(span * n), 10)
    grid = np.quantile(xs, np.linspace(0, 1, 60))
    fit = np.empty(len(grid))
    for i, g in enumerate(grid):
        d = np.abs(xs - g)
        h = np.partition(d, min(k, n - 1))[min(k, n - 1)] or 1e-9
        w = np.clip(1 - (d / h) ** 3, 0, None) ** 3
        W = w.sum()
        if W <= 0:
            fit[i] = ys.mean()
            continue
        xm, ym = (w * xs).sum() / W, (w * ys).sum() / W
        sxx = (w * (xs - xm) ** 2).sum()
        b = (w * (xs - xm) * (ys - ym)).sum() / sxx if sxx > 0 else 0.0
        fit[i] = ym + b * (g - xm)
    out = np.empty(n)
    out[order] = np.interp(xs, grid, fit)
    return out


def squeeze_variances(s2: np.ndarray, df: float, trend: np.ndarray | None = None) -> tuple:
    """(posterior variances, prior df, prior variances): Smyth (2004) moderation, with a trend if given."""
    s2 = np.maximum(np.asarray(s2, float), 1e-12)
    z = np.log(s2)
    e = z - special.digamma(df / 2) + np.log(df / 2)
    if trend is None:
        emean = np.mean(e)
        resid = e - emean
    else:
        emean = trend
        resid = e - trend
    n = len(e)
    evar = np.sum(resid ** 2) / (n - 1) - special.polygamma(1, df / 2)
    if evar > 0:
        d0 = 2 * float(_trigamma_inverse(np.array([evar]))[0])
        s0 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0 = np.inf
        s0 = np.exp(emean)
    if np.isinf(d0):
        return np.broadcast_to(s0, s2.shape).astype(float), d0, s0
    post = (d0 * s0 + df * s2) / (d0 + df)
    return post, d0, s0


def contrast(counts: pd.DataFrame, group_a: list, group_b: list, normalisation: str = "tmm",
             tpm: pd.DataFrame | None = None, min_expression: float = 1.0,
             quantifier: str = "salmon", reference: str = "ME49", *, sample_metadata: pd.DataFrame | None = None) -> Contrast:
    """Moderated two-group comparison of B against A, limma-trend style.

    `min_expression` is the filter edgeR's filterByExpr approximates: genes kept only where at least as
    many samples as the smaller group reach this CPM -- genes too lowly expressed to test are left out
    rather than reported as unchanged.
    """
    from . import comparison_provenance as P
    design = P.selection(counts, group_a, group_b, sample_metadata)
    a, b = design['group_a'], design['group_b']
    if normalisation not in NORMALISATIONS:
        raise ValueError('Unknown comparison normalization')
    if not np.isfinite(min_expression) or min_expression < 0:
        raise ValueError('Minimum expression must be finite and nonnegative')
    sub = counts[a + b]
    if any(not pd.api.types.is_numeric_dtype(dtype) or pd.api.types.is_complex_dtype(dtype) for dtype in sub.dtypes):
        raise ValueError('Comparison counts must have real numeric storage; text values are not coerced')
    try:
        values = sub.to_numpy(dtype=float)
    except (ValueError, TypeError) as error:
        raise ValueError('Comparison counts must be numeric') from error
    if not len(sub) or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Comparison counts must be nonempty, finite and nonnegative; missing values are not zero')
    inputs = dict(counts=P.matrix_record(sub))
    lib = sub.sum(axis=0)
    if not np.isfinite(lib.to_numpy(dtype=float)).all() or (lib <= 0).any():
        raise ValueError('Every selected run needs a positive finite count library')
    cpm = sub / lib * 1e6
    keep = (cpm >= min_expression).sum(axis=1) >= min(len(a), len(b))
    sub = sub[keep]
    if sub.empty:
        raise ValueError('No genes pass the selected expression filter')
    if normalisation == "tpm":
        factors = None
        if (tpm is None or not tpm.index.is_unique or not tpm.columns.is_unique
                or not counts.index.equals(tpm.index) or not counts.columns.equals(tpm.columns)):
            raise ValueError('TPM must match the exact count matrix gene and run axes')
        if any(not pd.api.types.is_numeric_dtype(dtype) or pd.api.types.is_complex_dtype(dtype) for dtype in tpm[a+b].dtypes):
            raise ValueError('TPM must have real numeric storage; text values are not coerced')
        numeric = tpm[a + b].to_numpy(dtype=float)
        if not np.isfinite(numeric).all() or (numeric < 0).any():
            raise ValueError('TPM must be finite and nonnegative; missing values are not zero')
        inputs['tpm'] = P.matrix_record(tpm[a + b])
        E = np.log2(tpm.loc[sub.index, a + b] + 1.0)
    else:
        if normalisation == 'median_ratio' and not (sub > 0).all(axis=1).any():
            raise ValueError('Median-ratio normalization needs genes positive in every selected run')
        factors = norm_factors(sub, normalisation)
        E = log_cpm(sub, factors)
    if not np.isfinite(E.to_numpy(dtype=float)).all():
        raise ValueError('The selected normalization does not yield finite expression for these runs')
    ea, eb = E[a].to_numpy(), E[b].to_numpy()
    ma, mb = ea.mean(axis=1), eb.mean(axis=1)
    df = len(a) + len(b) - 2
    s2 = (((ea - ma[:, None]) ** 2).sum(axis=1) + ((eb - mb[:, None]) ** 2).sum(axis=1)) / df
    amean = E.to_numpy().mean(axis=1)
    trend = _trend(amean, np.log(np.maximum(s2, 1e-12)))
    post, d0, _s0 = squeeze_variances(s2, df, trend - special.digamma(df / 2) + np.log(df / 2))
    se = np.sqrt(post * (1 / len(a) + 1 / len(b)))
    lfc = mb - ma
    t = lfc / se
    dft = df + (d0 if np.isfinite(d0) else 1e6)
    p = 2 * stats.t.sf(np.abs(t), dft)
    table = pd.DataFrame({"gene": sub.index, "log2fc": lfc, "mean_expression": amean,
                          "mean_a": ma, "mean_b": mb, "t": t, "p": p, "fdr": _bh(p),
                          "variance": s2, "moderated_variance": post}).sort_values("p")
    import hashlib
    import scipy
    from datetime import datetime, timezone
    from pathlib import Path
    record = dict(schema_version=1, created_utc=datetime.now(timezone.utc).isoformat(), design=design, model_inputs=inputs,
                  quantifier=quantifier, reference=reference, normalization=normalisation,
                  min_expression_cpm=float(min_expression), tested_genes=len(sub),
                  excluded_genes=len(counts)-len(sub),
                  filtering_library_sizes={str(k): float(v) for k, v in lib.items()},
                  normalization_library_sizes={str(k): float(v) for k, v in sub.sum(axis=0).items()},
                  normalization_factors=None if factors is None else {str(k): float(v) for k, v in factors.items()},
                  prior_degrees_of_freedom=float(d0) if np.isfinite(d0) else str(float(d0)),
                  normalization_scope='Selected runs only, within the caller-selected study; no cross-study batch correction',
                  expression_scale='log2(TPM+1)' if normalisation == 'tpm' else 'log2 CPM with prior count2',
                  algorithm='Existing limma-trend-style moderated two-group model; B minus A',
                  software=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
                                model_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                provenance_source_sha256=hashlib.sha256(Path(P.__file__).read_bytes()).hexdigest()),
                  operations=['Validate exact disjoint selected run identities and submitted sample links',
                              'Refuse nonfinite, negative, empty or zero-library selected counts',
                              'Filter by minimum CPM in at least the smaller group sample count',
                              'Apply selected within-comparison normalization, moderated variance model and BH FDR'],
                  biological_quality_score=None, control_pairing_verified=False)
    return Contrast(table.reset_index(drop=True), a, b, quantifier, reference, normalisation,
                    float(d0), float(min_expression), record)


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out
