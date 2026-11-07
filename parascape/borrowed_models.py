"""The borrowed profile-inference strategies (instruction 77), as estimators the harness can score.

Every one of them is an sklearn-style classifier (`fit(X, y)`, `predict_proba(X)`), registered in
`harness.MODELS` under its own name. That is the whole point of the shape: a strategy written this
way inherits the harness unchanged -- four partitions grouped by orthogroup, nested cross-validation
inside train + tune, `guard_fit` at every stage, isotonic calibration, the Clopper-Pearson bound,
the baselines, permutation importance by family and the once-only frozen test. Nothing here can
skip a check, because nothing here runs one.

Two of the nine are not estimators, and are said so rather than faked:

* **Selective prediction with a reject option** is how the harness already grades every strategy:
  the fill threshold is the lowest calibrated probability whose lower bound clears the bar, and
  `coverage_at_bar` is coverage at fixed precision. It is a property of the harness, not a model.
* **Multi-task learning** and **empirical-Bayes shrinkage** are listed in `UNTESTED` with the reason.

**Changing an estimator here changes a strategy without changing its hash** (`harness.Strategy` hashes the
model name and settings, not this code). Bump the strategy `version` with any change, or the changed model
would get a second look at a frozen partition already spent.

The evidence family of each column is needed by stacking and propensity weighting; the caller sets
it once with `set_families` (it is already part of the harness evidence identity, so the hash of a
run covers it).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

#: column -> evidence family, set by the caller before a run (`set_families`).
_FAMILIES: dict = {}
#: column -> information categories, for the attention columns propensity weighting models.
_CATEGORIES: dict = {}


def set_families(families, categories=None) -> None:
    """Tell the family-aware strategies which family (and category) each column belongs to."""
    _FAMILIES.clear()
    _FAMILIES.update({str(k): str(v) for k, v in pd.Series(families).items()})
    _CATEGORIES.clear()
    if categories is not None:
        _CATEGORIES.update({str(k): tuple(v) for k, v in dict(categories).items()})


def _gbt(**settings):
    from sklearn.ensemble import HistGradientBoostingClassifier

    params = dict(max_iter=150, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=10,
                  l2_regularization=1.0, early_stopping=False, random_state=0)
    params.update(settings)
    return HistGradientBoostingClassifier(**params)


class _GBT:
    """Boosted trees fitted on a subset of rows, dropping columns that are empty in that subset.

    The harness removes columns empty across all training genes; a bag or a fold sees fewer genes,
    and a column can be empty there alone, which the trees cannot bin.
    """

    def __init__(self, **gbt):
        self.gbt = gbt

    def fit(self, X, y, sample_weight=None):
        """Fit on the columns that carry at least one value in these rows."""
        self.columns = [c for c in X.columns if X[c].notna().any()]
        self.constant = not self.columns or len(np.unique(y)) < 2
        if self.constant:
            self.prior = float(np.mean(y))
            return self
        self.model = _gbt(**self.gbt).fit(X[self.columns], y, sample_weight=sample_weight)
        return self

    def predict_proba(self, X):
        """Probabilities on the kept columns (the prior when nothing could be fitted)."""
        if self.constant:
            return np.column_stack([np.full(len(X), 1 - self.prior), np.full(len(X), self.prior)])
        return self.model.predict_proba(X[self.columns])


class _Base:
    """Shared sklearn plumbing: `classes_`, a two-column `predict_proba`, `get_params`."""

    classes_ = np.array([0, 1])

    def get_params(self, deep: bool = True) -> dict:
        """The constructor settings, as sklearn expects."""
        return dict(self._settings)

    def _proba(self, score) -> np.ndarray:
        score = np.clip(np.asarray(score, dtype=float), 0.0, 1.0)
        return np.column_stack([1.0 - score, score])


class _Standardise:
    """Median-impute and z-score with statistics from the fitting rows only."""

    def fit(self, X: pd.DataFrame):
        values = X.to_numpy(dtype=float)
        self.median = np.nanmedian(values, axis=0)
        self.median = np.where(np.isfinite(self.median), self.median, 0.0)
        filled = np.where(np.isnan(values), self.median, values)
        self.mean, self.sd = filled.mean(axis=0), filled.std(axis=0)
        self.sd = np.where(self.sd > 0, self.sd, 1.0)
        return self

    def transform(self, X: pd.DataFrame) -> np.ndarray:
        values = X.to_numpy(dtype=float)
        filled = np.where(np.isnan(values), self.median, values)
        return (filled - self.mean) / self.sd


# --------------------------------------------------------------------------- 1 positive-unlabelled
class PUBagging(_Base):
    """Positive-unlabelled bagging (Mordelet and Vert 2014): "not recorded" is not "not true".

    Each bag trains boosted trees on every positive against a random subsample of the unlabelled
    genes the size of the positives; a gene is scored by the bags that did not draw it. The
    unlabelled genes are therefore never treated as a solid negative class, only as a background
    that is mostly negative.
    """

    def __init__(self, n_bags: int = 10, seed: int = 0, **gbt):
        self._settings = dict(n_bags=n_bags, seed=seed, **gbt)
        self.n_bags, self.seed, self.gbt = n_bags, seed, gbt

    def fit(self, X, y):
        """Fit one model per bag."""
        y = np.asarray(y).astype(bool)
        pos, unl = np.flatnonzero(y), np.flatnonzero(~y)
        rng = np.random.default_rng(self.seed)
        size = min(len(unl), max(len(pos), 1))
        self.models = []
        for _ in range(self.n_bags):
            draw = rng.choice(unl, size=size, replace=False)
            rows = np.concatenate([pos, draw])
            target = np.concatenate([np.ones(len(pos), int), np.zeros(len(draw), int)])
            self.models.append(_GBT(**self.gbt).fit(X.iloc[rows], target))
        return self

    def predict_proba(self, X):
        """Average of the bag probabilities (all bags: scoring rows are never fitting rows)."""
        score = np.mean([m.predict_proba(X)[:, 1] for m in self.models], axis=0)
        return self._proba(score)


# --------------------------------------------------------------------------- 2 propensity weighting
class PropensityWeighted(_Base):
    """Inverse-propensity weighting for labels missing not at random.

    A well-studied gene is more likely to carry a positive label at all. A small model of the label
    from the **attention** columns alone (publication counts, study effort) estimates that tilt;
    each training gene is weighted by the inverse of the probability attention alone gives its
    label, so the weighted sample carries no association between label and attention. The main
    model then sees every allowed column but cannot profit from attention. The principled version
    of deleting attention columns: it removes the bias from the labels, not the columns from the
    model. Falls back to unweighted when no attention column is present.

    It fails where attention nearly decides the label (no overlap: weights hit `clip`, and the
    weighted sample keeps the tilt). `weight_summary` reports the range; a maximum at the clip
    is that failure showing.
    """

    def __init__(self, clip: float = 20.0, **gbt):
        self._settings = dict(clip=clip, **gbt)
        self.clip, self.gbt = clip, gbt

    def _attention(self, X) -> list:
        return [c for c in X.columns if "attention" in _CATEGORIES.get(str(c), ())]

    def fit(self, X, y):
        """Estimate the attention propensity, then fit the weighted model."""
        from sklearn.linear_model import LogisticRegression

        y = np.asarray(y).astype(int)
        self.attention_columns = self._attention(X)
        weight = np.ones(len(y))
        if self.attention_columns and len(np.unique(y)) == 2:
            scaler = _Standardise().fit(X[self.attention_columns])
            propensity = LogisticRegression(max_iter=2000).fit(
                scaler.transform(X[self.attention_columns]), y).predict_proba(
                scaler.transform(X[self.attention_columns]))[:, 1]
            p_label = np.where(y == 1, propensity, 1.0 - propensity)
            weight = np.clip(1.0 / np.clip(p_label, 1e-6, None), 0, self.clip)
            weight = weight / weight.mean()
        self.sample_weight_ = weight
        self.weight_summary = (float(weight.min()), float(np.median(weight)), float(weight.max()))
        self.model = _GBT(**self.gbt).fit(X, y, sample_weight=weight)
        return self

    def predict_proba(self, X):
        """Probabilities from the weighted model."""
        return self.model.predict_proba(X)


# --------------------------------------------------------------------------- 3 low-rank completion
class LowRankCompletion(_Base):
    """Low-rank matrix completion: predict from the whole co-occurrence pattern of the evidence.

    The gene x evidence matrix is sparse. Soft-impute (iterated truncated SVD on the standardised
    matrix, observed cells held fixed) fills it at rank `rank`; a logistic model is fitted on the
    gene factors. A gene is predicted from where it sits in the shared low-dimensional structure
    rather than from any chosen column.
    """

    def __init__(self, rank: int = 30, iterations: int = 15, C: float = 1.0):
        self._settings = dict(rank=rank, iterations=iterations, C=C)
        self.rank, self.iterations, self.C = rank, iterations, C

    def _complete(self, Z, mask):
        filled = np.where(mask, Z, 0.0)
        for _ in range(self.iterations):
            approx = (filled @ self.basis.T) @ self.basis
            filled = np.where(mask, Z, approx)
        return filled @ self.basis.T

    def fit(self, X, y):
        """Learn the low-rank basis on the fitting rows, then the logistic model on the factors."""
        from sklearn.linear_model import LogisticRegression

        values = X.to_numpy(dtype=float)
        mask = ~np.isnan(values)
        self.mu = np.nanmean(values, axis=0)
        self.mu = np.where(np.isfinite(self.mu), self.mu, 0.0)
        self.sd = np.nanstd(values, axis=0)
        self.sd = np.where(np.isfinite(self.sd) & (self.sd > 0), self.sd, 1.0)
        Z = np.where(mask, (values - self.mu) / self.sd, 0.0)
        filled = Z.copy()
        rank = int(min(self.rank, min(Z.shape) - 1))
        for _ in range(self.iterations):
            _u, _s, vt = np.linalg.svd(filled, full_matrices=False)
            self.basis = vt[:rank]
            filled = np.where(mask, Z, (filled @ self.basis.T) @ self.basis)
        factors = filled @ self.basis.T
        self.model = LogisticRegression(max_iter=3000, C=self.C).fit(factors, np.asarray(y).astype(int))
        return self

    def predict_proba(self, X):
        """Complete each scoring row against the learnt basis, then score its factors."""
        values = X.to_numpy(dtype=float)
        mask = ~np.isnan(values)
        Z = np.where(mask, (values - self.mu) / self.sd, 0.0)
        return self.model.predict_proba(self._complete(Z, mask))


# --------------------------------------------------------------------------- 4 per-family stacking
class FamilyStacking(_Base):
    """A model per evidence family, combined by a meta-learner on out-of-fold predictions.

    Each family's contribution becomes one measurable weight instead of being hidden among dozens.
    The meta-learner is a non-negative logistic combination fitted on predictions each family model
    made for genes it was not fitted on (internal 3-fold), so a family that only memorises gets no
    weight. `weights_` names what each family is worth.
    """

    def __init__(self, folds: int = 3, min_columns: int = 1, seed: int = 0, **gbt):
        self._settings = dict(folds=folds, min_columns=min_columns, seed=seed, **gbt)
        self.folds, self.min_columns, self.seed, self.gbt = folds, min_columns, seed, gbt

    def fit(self, X, y):
        """Out-of-fold family scores, then the meta-learner, then the family models on everything."""
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import StratifiedKFold

        y = np.asarray(y).astype(int)
        groups: dict = {}
        for c in X.columns:
            groups.setdefault(_FAMILIES.get(str(c), "unregistered"), []).append(c)
        self.families = {f: cols for f, cols in sorted(groups.items()) if len(cols) >= self.min_columns}
        oof = np.full((len(y), len(self.families)), 0.5)
        folds = StratifiedKFold(self.folds, shuffle=True, random_state=self.seed)
        for j, cols in enumerate(self.families.values()):
            for tr, te in folds.split(np.zeros(len(y)), y):
                oof[te, j] = _GBT(**self.gbt).fit(X.iloc[tr][cols], y[tr]).predict_proba(X.iloc[te][cols])[:, 1]
        logit = np.log(np.clip(oof, 1e-4, 1 - 1e-4) / np.clip(1 - oof, 1e-4, 1))
        self.meta = LogisticRegression(max_iter=3000, C=1.0).fit(logit, y)
        coef = np.clip(self.meta.coef_[0], 0, None)
        self.meta.coef_ = coef.reshape(1, -1)
        self.weights_ = pd.Series(coef, index=list(self.families)).sort_values(ascending=False)
        self.models = {f: _GBT(**self.gbt).fit(X[cols], y) for f, cols in self.families.items()}
        return self

    def predict_proba(self, X):
        """Family scores for the scoring rows, combined by the meta-learner."""
        scores = np.column_stack([self.models[f].predict_proba(X[cols])[:, 1]
                                  for f, cols in self.families.items()])
        logit = np.log(np.clip(scores, 1e-4, 1 - 1e-4) / np.clip(1 - scores, 1e-4, 1))
        return self.meta.predict_proba(logit)


# --------------------------------------------------------------------------- 5 lookalike ranking
class Lookalike(_Base):
    """Rank every gene by how much it looks like the known positives (a seed set).

    The score is the mean cosine similarity to the `k` most similar positives minus the same to the
    `k` most similar negatives, on standardised evidence. No classifier boundary is learnt, which is
    the better-posed question when positives are few and negatives are invented.
    """

    def __init__(self, k: int = 10):
        self._settings = dict(k=k)
        self.k = k

    def fit(self, X, y):
        """Remember the standardised positives and negatives."""
        y = np.asarray(y).astype(bool)
        self.scaler = _Standardise().fit(X)
        Z = self.scaler.transform(X)
        Z = Z / np.clip(np.linalg.norm(Z, axis=1, keepdims=True), 1e-9, None)
        self.positives, self.negatives = Z[y], Z[~y]
        train_score = self._raw(Z, exclude_self=True)
        self.lo, self.hi = np.percentile(train_score, [1, 99])
        return self

    def _raw(self, Z, exclude_self: bool = False):
        def top(ref):
            sims = Z @ ref.T
            if exclude_self:
                sims = np.where(sims > 0.999999, -np.inf, sims)
            k = min(self.k, ref.shape[0])
            return np.sort(sims, axis=1)[:, -k:].mean(axis=1)
        return top(self.positives) - top(self.negatives)

    def predict_proba(self, X):
        """The similarity contrast, rescaled to [0, 1] by its range on the fitting genes."""
        Z = self.scaler.transform(X)
        Z = Z / np.clip(np.linalg.norm(Z, axis=1, keepdims=True), 1e-9, None)
        score = (self._raw(Z) - self.lo) / max(self.hi - self.lo, 1e-9)
        return self._proba(score)


# --------------------------------------------------------------------------- 7 network propagation
class GraphPropagation(_Base):
    """Personalised PageRank over a similarity graph of the evidence.

    Fitting and scoring genes are joined in one k-nearest-neighbour graph on standardised evidence
    (only the evidence of a scoring gene is used, never its label). A random walk with restart from
    the known positives, minus one from the known negatives, is the score: a gene inherits from
    neighbours of neighbours, not only from its nearest ones -- which is what the kNN strategy
    already built cannot do.
    """

    def __init__(self, k: int = 15, restart: float = 0.3, steps: int = 30):
        self._settings = dict(k=k, restart=restart, steps=steps)
        self.k, self.restart, self.steps = k, restart, steps

    def fit(self, X, y):
        """Remember the standardised fitting genes and their labels."""
        self.scaler = _Standardise().fit(X)
        self.Z = self.scaler.transform(X)
        self.y = np.asarray(y).astype(bool)
        return self

    def predict_proba(self, X):
        """Walk the joint graph and read the scoring genes."""
        from sklearn.neighbors import NearestNeighbors
        import scipy.sparse as sp

        Q = self.scaler.transform(X)
        A = np.vstack([self.Z, Q])
        n = A.shape[0]
        k = min(self.k, n - 1)
        index = NearestNeighbors(n_neighbors=k + 1).fit(A)
        _d, nbr = index.kneighbors(A)
        rows = np.repeat(np.arange(n), k)
        W = sp.csr_matrix((np.ones(n * k), (rows, nbr[:, 1:].ravel())), shape=(n, n))
        W = W + W.T
        degree = np.asarray(W.sum(axis=1)).ravel()
        P = sp.diags(1.0 / np.clip(degree, 1, None)) @ W

        def walk(seed_mask):
            seeds = np.zeros(n)
            seeds[:len(self.y)][seed_mask] = 1.0
            if seeds.sum() == 0:
                return seeds
            seeds /= seeds.sum()
            r = seeds.copy()
            for _ in range(self.steps):
                r = self.restart * seeds + (1 - self.restart) * (P.T @ r)
            return r * n

        score = walk(self.y) - walk(~self.y) * (self.y.mean() / max(1 - self.y.mean(), 1e-9))
        q = score[len(self.y):]
        lo, hi = np.percentile(score[:len(self.y)], [1, 99])
        return self._proba((q - lo) / max(hi - lo, 1e-9))


BORROWED = {
    "pu_bagging": PUBagging,
    "propensity_weighted": PropensityWeighted,
    "low_rank_completion": LowRankCompletion,
    "family_stacking": FamilyStacking,
    "lookalike": Lookalike,
    "graph_propagation": GraphPropagation,
}

#: The strategies of instruction 77 that are not scored as models here, and why.
UNTESTED = {
    "selective_prediction": (
        "Built into the harness rather than a model: every strategy is graded by coverage at fixed "
        "precision (`coverage_at_bar`), with the fill threshold as the reject option."),
    "multi_task": (
        "Needs auxiliary labels that may be used beside the target. The harness excludes every "
        "column of the target's own information category, and an auxiliary label of that category "
        "would bring it back in by another door. No pair of registry questions has yet been "
        "checked for which labels may stand beside which, so it is untested until the registry "
        "records that."),
    "empirical_bayes_shrinkage": (
        "Needs a per-cell standard error. The only columns that carry one are the expression fold "
        "changes, which are already shrunk when they are built (`expression_traits`); the node "
        "table columns carry no uncertainty to shrink by, and rescaling a column changes nothing "
        "for a tree or a standardised linear model."),
}
