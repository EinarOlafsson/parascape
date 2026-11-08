"""One scoring harness, written once, for every question in the registry.

This has now been written three times by hand -- host-facing v1 and v2, function, partners -- and
each copy repeated the same omission. All three used 5-fold cross-validation grouped by orthogroup,
which is a held-out measurement and not circular, but **no partition was held back from
everything**: the same folds that produced the reported number also chose the feature set, 14
variants for host-facing and 10 for function. Quoting the cross-validated score of the variant that
cross-validation picked is optimistic, so 0.727, 0.617 and 0.891 are upper bounds on what a frozen
test would give. Worse, `parascape.splits` already provided the right primitive and the runs did not
call it.

What this module enforces automatically, for every question:

1. **Four partitions** from `splits.make_split`, whole orthogroups to one side, seed recorded; and
   `splits.cold_node_pairs` for a ranked-list question, so both endpoints of a held-out pair are
   unseen.
2. **Nested cross-validation inside train + tune.** All searching -- model, settings, feature
   families, the fill threshold -- happens there. It is where strategies are compared.
3. **`split.guard_fit` at every stage**, through `StageGuard`, which also records the calls. A leak
   raises here rather than being noticed in review, and a question that cannot show its guard calls
   is not allowed to fill.
4. **A holdout budget.** The `test` partition is touched once per question per strategy. The use is
   appended to `results/harness_holdout_log.jsonl` with a hash of the strategy definition; a second
   use of the same partition by a *changed* strategy is refused, because it would make the frozen
   test a third tuning set. Re-running the identical strategy is allowed and recorded as a repeat.
5. **Both numbers, never one.** The nested cross-validated estimate and the frozen-test result, each
   with a Clopper-Pearson lower 90 % bound, coverage at the bar, calibration error, the baselines and
   permutation importance by family.
6. **A question with too few known answers says so.** The plain version of the rule: to grade a
   question honestly we hide some genes whose answer is already known, predict them, and see how
   often we were right. That needs enough known answers to hide a meaningful number and still have
   enough left to learn from -- `question.min_labels_for_frozen_test`, 400 by default. Below it the
   held-back set is skipped, the question is flagged `no_frozen_test` ("not enough known answers to
   grade"), and it may **rank genes but not write anything into the map**. The flag is not a failure
   of the method; it is a statement about how much is published.
7. **The exclusion is computed from information categories, not from a hand-written list.** Inferring
   a label of category X holds out every column of category X, whatever its name, family or study
   (`information_categories`, instruction 73), and `identity` columns always. The hand-written
   `evidence_forbidden` list is applied in addition and recorded as the override it is: it was the
   mechanism, it missed three pooled fitness screens whose family names say nothing about fitness,
   and it is now the thing that expresses what a category cannot. `evidence_selection` returns the
   requested and the effective selection with the reason for every removal, and
   `HarnessResult.evidence_selection` carries it into the run record.

**Filling a cell requires the frozen-test lower bound to clear the bar**, which is what `may_fill`
reports -- never the cross-validated one.

    result = run_question(question, labels=y, features=X, families=fam, groups=og,
                          strategies=[Strategy.of('all families', 'hist_gbt')])
    result.summary()      # the two numbers side by side
    may_fill(result)
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import information_categories as IC
from . import organisms as O
from . import splits as S
from .claims import CertaintyModel, _ece, _isotonic

BAR = 0.90                       # the fill bar: held-out precision, as a lower bound
CONFIDENCE = 0.90                # of the one-sided Clopper-Pearson bound
HOLDOUT_LOG = "results/harness_holdout_log.jsonl"
STAGES_REQUIRED = ("setting_selection", "refit", "calibration")


def clopper_pearson_lower(k: int, n: int, confidence: float = CONFIDENCE) -> float:
    """Exact one-sided lower confidence bound on a proportion -- the number the fill rule is about.

    The point estimate of precision is not the fill rule; this is. With no calls there is no bound,
    reported as NaN rather than as zero, because "made no calls" and "was wrong" differ.
    """
    from scipy.stats import beta

    k, n = int(k), int(n)
    if n <= 0 or k < 0 or k > n:
        return float("nan")
    if k == 0:
        return 0.0
    if k == n:
        return float(beta.ppf(1 - confidence, k, 1))
    return float(beta.ppf(1 - confidence, k, n - k + 1))


def _hist_gbt(**settings):
    """Histogram gradient boosting with missing values passed through, no imputation anywhere."""
    from sklearn.ensemble import HistGradientBoostingClassifier

    params = dict(max_iter=200, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=10,
                  l2_regularization=1.0, early_stopping=False, random_state=0)
    params.update(settings)
    return HistGradientBoostingClassifier(**params)


def _logistic(**settings):
    """Median-imputed, scaled logistic regression -- the linear comparison for the boosted trees."""
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    params = dict(max_iter=2000, C=1.0)
    params.update(settings)
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                         LogisticRegression(**params))


MODELS = {"hist_gbt": _hist_gbt, "logistic": _logistic}
# The borrowed profile-inference strategies (instruction 77), scored by this same harness.
from .borrowed_models import BORROWED as _BORROWED  # noqa: E402
MODELS.update(_BORROWED)


@dataclass(frozen=True)
class Strategy:
    """Everything that would change a result, named and hashable, so the holdout budget can be kept.

    The hash covers the model, its settings, the evidence families used and an explicit version
    string. Two runs with the same hash are the same strategy and may look at the frozen partition
    again; anything else is a changed strategy and is refused.
    """

    name: str
    model: str
    families: tuple = ()             # empty = every family the question allows
    settings: str = "{}"             # canonical JSON, so the dataclass stays hashable
    version: str = "1"

    def __post_init__(self):
        if not self.name or self.model not in MODELS:
            raise ValueError(f"Unknown model {self.model!r}; declare one of {sorted(MODELS)}")
        json.loads(self.settings)

    @classmethod
    def of(cls, name: str, model: str, families=(), version: str = "1", **settings) -> "Strategy":
        """Build a strategy from keyword model settings, canonicalising them for the hash."""
        return cls(name, model, tuple(sorted(families)),
                   json.dumps(settings, sort_keys=True, separators=(",", ":")), version)

    @property
    def definition(self) -> dict:
        """The strategy as a plain dictionary -- what the hash is taken over."""
        return {"model": self.model, "families": list(self.families),
                "settings": json.loads(self.settings), "version": self.version}

    @property
    def identity(self) -> str:
        """Content hash of the definition; the key of the holdout budget together with the question."""
        return hashlib.sha256(json.dumps(self.definition, sort_keys=True,
                                         separators=(",", ":")).encode()).hexdigest()

    def build(self):
        """A fresh unfitted estimator for this strategy."""
        return MODELS[self.model](**json.loads(self.settings))

    def columns(self, features: pd.DataFrame, families: pd.Series) -> list:
        """The feature columns this strategy uses, by family."""
        if not self.families:
            return list(features.columns)
        keep = set(self.families)
        return [c for c in features.columns if families.get(c, "unregistered") in keep]


class StageGuard:
    """Calls `SplitManifest.guard_fit` and records every call, so the audit is in the result.

    `FIT_ACCESS` is the permission table: `model_fit` and `feature_selection` see only `train`,
    `setting_selection` and `refit` see `train` + `tune`, `calibration` sees only `calibration`. Any
    other reach raises. Nested cross-validation inside train + tune therefore fits under `refit`,
    not `model_fit` -- the distinction `splits.FIT_ACCESS` already draws.
    """

    def __init__(self, split):
        self.split = split
        self._calls: list[tuple] = []

    def guard(self, stage: str, entities) -> None:
        """Check a stage's access to these entities and record that the check happened."""
        entities = list(entities)
        self.split.guard_fit(stage, entities)
        self._calls.append((stage, len(entities)))

    @property
    def calls(self) -> tuple:
        """Every guarded stage, in order, with the number of entities it touched."""
        return tuple(self._calls)

    @property
    def stages(self) -> frozenset:
        """The distinct stages that have been guarded."""
        return frozenset(stage for stage, _n in self._calls)


def partition_question(question, entities, groups, *, organism: str = O.TOXOPLASMA, seed: int = 0,
                       benchmark_id: str | None = None, pairs=None,
                       fractions=(0.55, 0.15, 0.15, 0.15)):
    """Four orthogroup-grouped partitions for a question, and a cold-node pair partition if needed.

    Whole orthogroups go to one side, because paralogs share a fate and the GRA and ROP families are
    full of them. A ranked-list question additionally gets `cold_node_pairs`, so a held-out pair has
    *both* endpoints outside training rather than one.
    """
    split = S.make_split(list(entities), list(groups), organism=organism,
                         benchmark_id=benchmark_id or question.id, group_kind="homology",
                         seed=seed, fractions=fractions)
    if question.answer_type == "ranked_list":
        if pairs is None:
            raise ValueError(f"{question.id}: a ranked-list question needs its pairs for a cold-node split")
        return split, S.cold_node_pairs(pairs, split)
    return split, None


def evidence_selection(question, columns, families, *, force_include=None) -> dict:
    """What this question may see, and why anything was removed: the instruction 73 run record.

    The exclusion is computed from **information categories**, not from the hand-written family
    list. A column is removed when any of its categories matches the category of the question's own
    labels (`information_categories.label_categories`), whatever its name, family or study -- which
    is what catches a pooled fitness screen called `gra17_synthlethal`, where two passes of
    hand-written guessing did not. `identity` columns are removed for every question, because an
    accession is not evidence.

    The hand-written `evidence_forbidden` list is kept and applied **in addition**, and recorded
    separately: it expresses what a category cannot, such as holding `pubmed` out of the host-facing
    question. It is an override, never the mechanism.

    `force_include` is the escape hatch instruction 73 asks for: a mapping of column name **or**
    category name to the reason for including it anyway. Every use is returned in the record, so a
    deliberately circular run is visible on the scorecard rather than impossible.

    Keys: `label_categories`, `requested`, `effective`, `excluded_by_category`,
    `excluded_by_family`, `force_included`, `unannotated`.
    """
    columns = [str(c) for c in columns]
    families = pd.Series(families) if not isinstance(families, pd.Series) else families
    label = tuple(IC.label_categories(question))
    forced = {str(k): str(v) for k, v in dict(force_include or {}).items()}
    excluded_category, excluded_family, force_used, unannotated = {}, {}, {}, []
    effective = []
    for column in columns:
        carried = IC.categories_of(column)
        if not carried:
            unannotated.append(column)
        hit = tuple(sorted((set(carried) & set(label)) | (set(carried) & set(IC.NEVER_A_FEATURE))))
        reason = forced.get(column) or next((forced[c] for c in hit if c in forced), "")
        if hit and reason:
            force_used[column] = f"{','.join(hit)}: {reason}"
        elif hit:
            excluded_category[column] = ",".join(hit)
            continue
        family = families.get(column, "unregistered")
        if not question.allows(family):
            excluded_family[column] = family
            continue
        effective.append(column)
    return {"label_categories": label, "requested": tuple(columns), "effective": tuple(effective),
            "excluded_by_category": excluded_category, "excluded_by_family": excluded_family,
            "force_included": force_used, "unannotated": tuple(unannotated)}


def allowed_features(question, features: pd.DataFrame, families: pd.Series,
                     *, force_include=None) -> pd.DataFrame:
    """Drop every column the question may not see -- the caller's chance to filter first.

    Both exclusions apply: the information category of the question's labels (and `identity`), and
    the hand-written forbidden families. `run_question` refuses such a column rather than silently
    removing it, so that a leak is a failure and not a quiet correction; this helper is how a caller
    asks for the filtering on purpose.
    """
    record = evidence_selection(question, features.columns, families, force_include=force_include)
    return features[list(record["effective"])]


def _support_and_correctness(probability: pd.Series, truth: pd.Series):
    """Turn positive-class probabilities into the (prediction, support, correct) triple claims uses.

    `claims.CertaintyModel` maps a strategy's *support* to its held-out rate of being right, so the
    harness speaks that language: support is the probability of whichever answer was given, and for
    a positive call P(right) is precision, which is exactly what the fill bar is about.
    """
    probability = probability.astype(float)
    prediction = probability >= 0.5
    support = np.where(prediction, probability, 1.0 - probability)
    correct = (prediction == truth.reindex(probability.index).astype(bool))
    return prediction, pd.Series(support, index=probability.index), correct


def _curve(x, y):
    """An isotonic support -> P(right) curve, reusing `claims._isotonic` rather than refitting one."""
    return _isotonic(np.asarray(x, dtype=float), np.asarray(y, dtype=float))


def _apply(curve, support):
    """Evaluate an isotonic curve on new support values."""
    xs, ys = curve
    s = np.asarray(support, dtype=float)
    if not len(xs):
        return np.full(s.shape, np.nan)
    return np.interp(s, xs, ys)


def _threshold_for_bar(calibrated, correct, prediction, bar: float = BAR) -> float:
    """The lowest calibrated probability whose positive-call lower bound reaches the bar.

    Chosen on train + tune only. If nothing reaches the bar the threshold that gets closest is
    returned, so the result still reports how far away it is instead of refusing to say.
    """
    calibrated = np.asarray(calibrated, dtype=float)
    correct = np.asarray(correct, dtype=bool)
    prediction = np.asarray(prediction, dtype=bool)
    best, best_bound = 0.5, -1.0
    candidates = sorted({round(float(v), 4) for v in calibrated[prediction & np.isfinite(calibrated)]})
    for threshold in candidates or [0.5]:
        called = prediction & (calibrated >= threshold)
        n = int(called.sum())
        if n < 5:
            continue
        bound = clopper_pearson_lower(int(correct[called].sum()), n)
        if bound >= bar:
            return float(threshold)
        if bound > best_bound:
            best, best_bound = float(threshold), bound
    return best


def gene_rows(name: str, probability: pd.Series, truth: pd.Series, curve, threshold: float) -> pd.DataFrame:
    """One row per scored gene of a partition: its probability, the call, and whether it was right.

    ``prediction`` is the answer given (``positive``/``negative``); ``called`` is whether it is a positive call
    at the fill threshold (the calls `_block` counts). The same rows give every summary number of `_block`.
    """
    prediction, support, correct = _support_and_correctness(probability, truth)
    calibrated = pd.Series(_apply(curve, support), index=probability.index)
    return pd.DataFrame({"entity": probability.index.astype(str), "partition": name,
                         "probability": probability.astype(float).to_numpy(),
                         "prediction": np.where(prediction, "positive", "negative"),
                         "calibrated": calibrated.to_numpy(),
                         "called": (prediction & (calibrated >= threshold)).to_numpy(),
                         "truth": np.where(truth.reindex(probability.index).astype(bool), "positive", "negative"),
                         "correct": correct.to_numpy()})


def _block(name: str, probability: pd.Series, truth: pd.Series, curve, threshold: float,
           bar: float = BAR) -> dict:
    """The same numbers for every partition: precision, its lower bound, coverage, calibration, AUC."""
    from sklearn.metrics import average_precision_score, roc_auc_score

    prediction, support, correct = _support_and_correctness(probability, truth)
    calibrated = pd.Series(_apply(curve, support), index=probability.index)
    called = prediction & (calibrated >= threshold)
    k, n = int(correct[called].sum()), int(called.sum())
    y = truth.reindex(probability.index).astype(bool)
    two_classes = 0 < int(y.sum()) < len(y)
    ok = np.isfinite(calibrated)
    return {"partition": name,
            "n": int(len(probability)),
            "positives": int(y.sum()),
            "threshold": float(threshold),
            "calls": n,
            "correct": k,
            "precision": float(k / n) if n else float("nan"),
            "precision_lower90": clopper_pearson_lower(k, n),
            "reaches_bar": bool(n and clopper_pearson_lower(k, n) >= bar),
            "coverage_at_bar": float(n / len(probability)) if len(probability) else float("nan"),
            "calibration_error": float(_ece(calibrated[ok], correct[ok].astype(float)))
            if ok.sum() >= 10 else float("nan"),
            "accuracy": float(correct.mean()),
            "roc_auc": float(roc_auc_score(y, probability)) if two_classes else float("nan"),
            "average_precision": float(average_precision_score(y, probability)) if two_classes else float("nan")}


def by_stratum(probability: pd.Series, truth: pd.Series, curve, threshold: float, strata: pd.Series,
               bar: float = BAR) -> dict:
    """The frozen-test numbers again within each stratum (instruction 64.16): understudied against
    well-studied genes, common against rare classes. The threshold and calibration are the ones chosen on
    train + tune, applied unchanged; a stratum with fewer than 10 genes or one class reports its counts only."""
    out = {}
    strata = strata.reindex(probability.index)
    for name in sorted(strata.dropna().unique(), key=str):
        idx = strata.index[strata == name]
        y = truth.reindex(idx).astype(bool)
        if len(idx) < 10 or y.nunique() < 2:
            out[str(name)] = {"n": int(len(idx)), "positives": int(y.sum()), "too_few": True}
            continue
        block = _block(f"stratum {name}", probability.loc[idx], truth, curve, threshold, bar)
        out[str(name)] = {k: block[k] for k in ("n", "positives", "calls", "precision",
                                                "precision_lower90", "roc_auc", "average_precision",
                                                "calibration_error")}
    return out


def _grouped_folds(entities, labels, groups, n_folds: int):
    """Grouped folds over an index, falling back to fewer folds than asked when groups are scarce."""
    from sklearn.model_selection import GroupKFold

    entities = np.asarray(list(entities))
    g = np.asarray([groups[e] for e in entities])
    folds = int(min(n_folds, len(set(g.tolist()))))
    if folds < 2:
        raise ValueError("Fewer than two protected groups in this partition; no grouped folds")
    for train_index, test_index in GroupKFold(n_splits=folds).split(entities, None, g):
        yield entities[train_index], entities[test_index]


def _fit_predict(strategy, features, families, labels, train, score, guard, stage: str):
    """Fit a strategy on `train` under a guarded stage and return its probabilities on `score`."""
    guard.guard(stage, train)
    columns = strategy.columns(features, families)
    usable = [c for c in columns
              if str(features[c].dtype) == "category" or features.loc[train, c].notna().any()]
    model = strategy.build().fit(features.loc[train, usable], labels.loc[train].astype(int))
    probability = model.predict_proba(features.loc[score, usable])[:, 1]
    return model, pd.Series(probability, index=pd.Index(score)), usable


def _select_strategy(strategies, features, families, labels, train, groups, guard, inner_folds):
    """Inner cross-validation: pick the strategy with the best average precision on train folds only."""
    from sklearn.metrics import average_precision_score

    scores = {}
    for strategy in strategies:
        values = []
        for inner_train, inner_test in _grouped_folds(train, labels, groups, inner_folds):
            if labels.loc[inner_train].nunique() < 2 or labels.loc[inner_test].nunique() < 2:
                continue
            _m, probability, _u = _fit_predict(strategy, features, families, labels,
                                               inner_train, inner_test, guard, "refit")
            values.append(average_precision_score(labels.loc[inner_test].astype(bool), probability))
        scores[strategy.name] = float(np.mean(values)) if values else float("nan")
    ranked = sorted(strategies, key=lambda s: (-(scores[s.name] if np.isfinite(scores[s.name]) else -1), s.name))
    return ranked[0], scores


def _baselines(labels, features, groups, train, score, bar: float = BAR) -> pd.DataFrame:
    """Majority, prevalence-weighted random, best single feature and ortholog copy, on one partition.

    Nothing is reported except as improvement over the best of these. Ortholog copy has coverage 0
    under orthogroup-grouped partitions by construction, which is correct: a baseline that only works
    because paralogs leak is not a baseline.
    """
    train, score = list(train), list(score)
    y = labels.loc[score].astype(bool)
    prevalence = float(labels.loc[train].astype(bool).mean())
    rows = []

    def row(name, called, correct, coverage=1.0, auc=float("nan")):
        k, n = int(np.sum(correct & called)), int(np.sum(called))
        rows.append({"baseline": name, "calls": n, "precision": float(k / n) if n else float("nan"),
                     "precision_lower90": clopper_pearson_lower(k, n),
                     "reaches_bar": bool(n and clopper_pearson_lower(k, n) >= bar),
                     "coverage": float(coverage), "roc_auc": auc})

    majority = prevalence > 0.5
    row("training majority", np.full(len(y), majority), (y.values == majority))
    rng = np.random.default_rng(0)
    coin = rng.random(len(y)) < prevalence
    row("prevalence-weighted random", coin, coin == y.values)

    numeric = [c for c in features.columns if str(features[c].dtype) != "category"
               and features.loc[train, c].notna().sum() > 10 and labels.loc[train].nunique() > 1]
    best, best_auc, best_cut, best_sign = None, 0.5, 0.0, 1.0
    if numeric:
        from sklearn.metrics import roc_auc_score
        yt = labels.loc[train].astype(bool)
        for column in numeric:
            values = features.loc[train, column].astype(float)
            ok = values.notna()
            if ok.sum() < 10 or yt[ok].nunique() < 2:
                continue
            auc = roc_auc_score(yt[ok], values[ok].fillna(0.0))
            auc = max(auc, 1 - auc)
            if auc > best_auc:
                sign = 1.0 if roc_auc_score(yt[ok], values[ok]) >= 0.5 else -1.0
                signed = sign * values[ok]
                cuts = np.quantile(signed, np.linspace(0.05, 0.95, 19))
                bounds = [(clopper_pearson_lower(int((yt[ok] & (signed >= c)).sum()),
                                                 int((signed >= c).sum())), c) for c in cuts]
                best, best_auc = column, auc
                best_cut = float(max(bounds)[1]) * 1.0
                best_sign = sign
    if best is not None:
        signed = best_sign * features.loc[score, best].astype(float)
        called = (signed >= best_cut).fillna(False).values
        row(f"best single feature ({best})", called, y.values, auc=best_auc)

    # Only peers outside the scored set may be copied from, so a gene never copies itself or a
    # paralog that is being scored in the same breath. Under orthogroup-grouped partitions that
    # leaves nothing to copy and coverage is 0 -- which is correct, not a defect.
    outside = [e for e in train if e not in set(score)]
    by_group: dict = {}
    for entity, value in labels.loc[outside].astype(bool).items():
        by_group.setdefault(groups[entity], []).append(bool(value))
    called, correct = [], []
    for entity in score:
        peers = by_group.get(groups[entity], [])
        called.append(bool(peers) and float(np.mean(peers)) >= 0.5)
        correct.append(bool(peers) and (float(np.mean(peers)) >= 0.5) == bool(y[entity]))
    row("ortholog copy (abstains)", np.array(called), np.array(correct),
        coverage=float(np.mean([bool(by_group.get(groups[e], [])) for e in score])))
    return pd.DataFrame(rows)


def _family_importance(model, features, families, labels, columns, entities, repeats: int = 3,
                       seed: int = 0) -> pd.DataFrame:
    """Permutation importance by family, on a partition the model was not fitted on.

    A whole family is shuffled together, because 36 expression columns from one experiment are one
    piece of evidence and not 36, and shuffling them one at a time hides a family whose members are
    redundant with each other.
    """
    from sklearn.metrics import average_precision_score

    entities = list(entities)
    y = labels.loc[entities].astype(bool)
    if y.nunique() < 2:
        return pd.DataFrame(columns=["family", "n_features", "importance_mean", "importance_sd"])
    frame = features.loc[entities, columns]
    base = average_precision_score(y, model.predict_proba(frame)[:, 1])
    rng = np.random.default_rng(seed)
    rows: list = []
    grouped: dict = {}
    for column in columns:
        grouped.setdefault(families.get(column, "unregistered"), []).append(column)
    for family, members in sorted(grouped.items()):
        drops = []
        for _ in range(repeats):
            shuffled = frame.copy()
            order = rng.permutation(len(shuffled))
            for column in members:
                shuffled[column] = shuffled[column].to_numpy()[order]
            drops.append(base - average_precision_score(y, model.predict_proba(shuffled)[:, 1]))
        rows.append({"family": family, "n_features": len(members),
                     "importance_mean": float(np.mean(drops)), "importance_sd": float(np.std(drops))})
    return (pd.DataFrame(rows).sort_values("importance_mean", ascending=False)
            .reset_index(drop=True).assign(baseline_average_precision=base))


def evidence_identity(question, features: pd.DataFrame, selection=None) -> str:
    """Hash of the question's own definition and the exact feature set used to answer it.

    The strategy hash alone is not enough to key the holdout budget, and the first demo run showed
    why: widening or narrowing `evidence_forbidden` changes the evidence a model sees while leaving
    the model, its settings and its version untouched. Without this, re-running "the same strategy"
    on a changed feature table would pass the budget as a repeat and the frozen partition would
    quietly become a tuning set again.

    `selection` is the `evidence_selection` record, and the force-included columns are hashed with
    it: including one circular column on purpose changes what the model sees exactly as widening the
    forbidden list does, so it must spend the partition too.
    """
    payload = {"answer_type": question.answer_type, "label_source": question.label_source,
               "label_family": question.label_family,
               "evidence_allowed": sorted(question.evidence_allowed),
               "evidence_forbidden": sorted(question.evidence_forbidden),
               "label_categories": sorted((selection or {}).get("label_categories", ())),
               "force_included": sorted((selection or {}).get("force_included", {})),
               "min_labels_for_frozen_test": int(question.min_labels_for_frozen_test),
               "features": sorted(map(str, features.columns))}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _log_path(path=None, root=None) -> Path:
    """Where the holdout budget is kept; one persistent file per repository."""
    if path is not None:
        return Path(path)
    return Path(root or ".") / HOLDOUT_LOG


def holdout_log(path=None, root=None) -> pd.DataFrame:
    """Every recorded use of a frozen test partition, as a table."""
    file = _log_path(path, root)
    if not file.exists():
        return pd.DataFrame(columns=["question_id", "strategy", "strategy_hash", "split_identity",
                                     "seed", "n_test", "precision", "precision_lower90", "when"])
    rows = [json.loads(line) for line in file.read_text().splitlines() if line.strip()]
    return pd.DataFrame(rows)


def holdout_uses(question_id: str, split_identity: str, path=None, root=None) -> list:
    """Recorded uses of one question's one frozen partition, newest last."""
    frame = holdout_log(path, root)
    if not len(frame):
        return []
    hit = frame[(frame["question_id"] == question_id) & (frame["split_identity"] == split_identity)]
    return hit.to_dict("records")


def record_holdout_use(entry: dict, path=None, root=None) -> Path:
    """Append one frozen-test use to the persistent log, creating it if needed."""
    file = _log_path(path, root)
    file.parent.mkdir(parents=True, exist_ok=True)
    with file.open("a") as stream:
        stream.write(json.dumps(entry, sort_keys=True) + "\n")
    return file


@dataclass
class HarnessResult:
    """Both numbers for one question: the nested cross-validated estimate and the frozen-test result.

    `frozen_test` is None when the question cannot afford a held-back partition (`no_frozen_test`) or
    when the holdout budget refused a second look (`frozen_test_refused`); `refusal` says which.
    Filling cells depends on `frozen_test` alone -- see `may_fill`.
    """

    question_id: str
    answer_type: str
    question_class: str
    n_labels: int
    seed: int
    split_identity: str
    bar: float
    strategy: str
    strategy_hash: str
    threshold: float
    cross_validated: dict
    evidence_hash: str = ""
    evidence_selection: dict = field(default_factory=dict)
    frozen_test: dict | None = None
    inner_scores: dict = field(default_factory=dict)
    chosen_per_outer_fold: tuple = ()
    calibration: dict = field(default_factory=dict)
    baselines: pd.DataFrame = field(default_factory=pd.DataFrame)
    importance: pd.DataFrame = field(default_factory=pd.DataFrame)
    flags: tuple = ()
    guard_calls: tuple = ()
    refusal: str = ""
    rows: pd.DataFrame = field(default_factory=pd.DataFrame)
    #: The deployed fit (train + tune), the columns it reads and the calibration curve fitted on the calibration
    #: partition: what `deploy` applies to genes with no label. Not part of the serialised record.
    deployed: dict = field(default_factory=dict, repr=False)

    @property
    def dominance(self) -> dict:
        """Whether one evidence family carries most of the importance (`diagnostics.dominance`)."""
        from .diagnostics import dominance
        return dominance(self.importance)

    @property
    def may_fill(self) -> bool:
        """Whether this question may write cells as facts: the frozen-test lower bound clears the bar."""
        if self.frozen_test is None or "no_frozen_test" in self.flags:
            return False
        if not set(STAGES_REQUIRED) <= {stage for stage, _n in self.guard_calls}:
            return False
        bound = self.frozen_test.get("precision_lower90", float("nan"))
        return bool(np.isfinite(bound) and bound >= self.bar)

    def summary(self) -> pd.DataFrame:
        """The two numbers side by side, which is how a result from this harness is read."""
        rows = [dict(self.cross_validated)]
        if self.frozen_test is not None:
            rows.append(dict(self.frozen_test))
        frame = pd.DataFrame(rows)
        frame.insert(0, "question", self.question_id)
        frame.insert(1, "strategy", self.strategy)
        return frame

    def as_dict(self) -> dict:
        """A JSON-serialisable record of the run, for results directories and notebooks."""
        return {"question_id": self.question_id, "answer_type": self.answer_type,
                "question_class": self.question_class, "n_labels": self.n_labels,
                "seed": self.seed, "split_identity": self.split_identity, "bar": self.bar,
                "strategy": self.strategy, "strategy_hash": self.strategy_hash,
                "evidence_hash": self.evidence_hash, "threshold": self.threshold,
                "evidence_selection": {
                    key: (list(value) if isinstance(value, tuple) else value)
                    for key, value in self.evidence_selection.items()},
                "cross_validated": self.cross_validated,
                "frozen_test": self.frozen_test, "inner_scores": self.inner_scores,
                "chosen_per_outer_fold": list(self.chosen_per_outer_fold),
                "calibration": self.calibration, "flags": list(self.flags),
                "guard_calls": [list(call) for call in self.guard_calls],
                "refusal": self.refusal, "may_fill": self.may_fill}


def deploy(result: HarnessResult, features: pd.DataFrame) -> pd.DataFrame:
    """Apply the deployed fit to genes with no label: probability, calibrated probability, and the call.

    `called` is True only where the answer is yes AND the calibrated probability reaches the threshold chosen on
    train + tune -- the same rule that made the frozen-test precision -- and `fills` is True only where, in
    addition, the question may fill (`may_fill`). Every other row ranks a gene; none of them is a fact.
    """
    d = result.deployed or {}
    if not d:
        raise ValueError("this result has no deployed fit")
    probability = pd.Series(d["model"].predict_proba(features[d["columns"]])[:, 1], index=features.index)
    yes = probability >= 0.5
    support = np.where(yes, probability, 1.0 - probability)
    calibrated = pd.Series(_apply(d["curve"], support), index=features.index)
    called = yes & (calibrated >= result.threshold)
    return pd.DataFrame({"probability": probability, "calibrated": calibrated, "called": called,
                         "fills": called & bool(result.may_fill)})


def may_fill(result: HarnessResult) -> bool:
    """May this question fill cells? True only when the FROZEN-TEST lower bound clears the bar.

    The one function the rest of the matrix is allowed to ask. A cross-validated bound above the bar
    is not sufficient and never has been: the folds that produced it also chose the strategy.
    """
    return result.may_fill


def run_question(question, *, labels, features, families, groups, strategies,
                 organism: str = O.TOXOPLASMA, seed: int = 0, outer_folds: int = 3,
                 inner_folds: int = 3,
                 bar: float = BAR, pairs=None, log_path=None, root=None, record: bool = True,
                 importance_repeats: int = 3, force_include=None, strata=None) -> HarnessResult:
    """Score one question end to end: four partitions, nested cross-validation, one frozen test.

    `labels` is a boolean series indexed by entity, `features` a frame over at least those entities,
    `families` a mapping from feature to evidence family, `groups` from entity to protected group
    (orthogroup). A column of the labels' own **information category**, or one whose family the
    question forbids, raises here rather than being dropped quietly -- call `allowed_features` first
    if filtering is what was meant. `force_include={'column or category': 'reason'}` includes a
    circular column on purpose; it is recorded in `HarnessResult.evidence_selection`, shown on the
    scorecard and hashed into the evidence identity, so it also spends the holdout budget.

    For a ranked-list question the entities are the *nodes* and `pairs` the labelled pairs: the
    cold-node partition keeps only same-partition pairs, so a held-out pair has both endpoints
    unseen, and the guards are called with endpoint nodes, which is the unit `splits` partitions.
    """
    labels = pd.Series(labels).dropna()
    labels = labels.astype(bool)
    entities = list(labels.index)
    groups = pd.Series(groups).reindex(entities)
    groups = groups.where(groups.notna(), pd.Series({e: "self:" + str(e) for e in entities}))
    families = pd.Series(families)
    strategies = list(strategies)
    if not strategies:
        raise ValueError(f"{question.id}: no strategy to score")

    selection = evidence_selection(question, features.columns, families, force_include=force_include)
    if selection["excluded_by_category"]:
        circular = sorted(selection["excluded_by_category"])
        raise ValueError(
            f"{question.id}: columns of its own information category "
            f"{list(selection['label_categories'])} are present in the feature table: "
            f"{[(c, selection['excluded_by_category'][c]) for c in circular[:10]]} "
            f"({len(circular)} columns). Inferring a label of category X excludes every column of "
            f"category X (instruction 73). Call allowed_features first, or pass force_include="
            "{'<column or category>': '<reason>'} to include one on purpose and have it recorded.")
    forbidden = sorted(selection["excluded_by_family"])
    if forbidden:
        raise ValueError(f"{question.id}: forbidden evidence families present in the feature table: "
                         f"{forbidden[:10]} ({len(forbidden)} columns). The label family "
                         f"{question.label_family!r} and everything in evidence_forbidden must be held out.")
    for strategy in strategies:
        bad = sorted(f for f in strategy.families if not question.allows(f))
        if bad:
            raise ValueError(f"{question.id}: strategy {strategy.name!r} asks for forbidden families {bad}")

    # Measured redundancy (instruction 73, Einar 2026-11-07: report only): inputs that move with the label
    # column although their category differs are named on the result, never removed.
    from .redundancy import overlap as _overlap
    label_column = str(getattr(question, "label_source", "")).rsplit(":", 1)[-1].split(" ")[0]
    selection["measured_overlap"] = _overlap(label_column, features.columns, organism)
    evidence = evidence_identity(question, features, selection)
    split, pair_partition = partition_question(question, entities, groups.loc[entities],
                                               organism=organism, seed=seed, pairs=pairs)
    guard = StageGuard(split)
    train = list(split.entities("train"))
    tune = list(split.entities("tune"))
    calibration_entities = list(split.entities("calibration"))
    test = list(split.entities("test"))
    pool = train + tune
    flags: list[str] = []

    # --- nested cross-validation, inside train + tune only -------------------------------------
    out_of_fold = pd.Series(np.nan, index=pd.Index(pool), dtype=float)
    chosen, inner_scores = [], {}
    for outer_train, outer_test in _grouped_folds(pool, labels, groups, outer_folds):
        guard.guard("setting_selection", outer_train)
        best, scores = _select_strategy(strategies, features, families, labels, outer_train,
                                        groups, guard, inner_folds)
        chosen.append(best.name)
        inner_scores[f"fold {len(chosen)}"] = scores
        _model, probability, _used = _fit_predict(best, features, families, labels,
                                                  outer_train, outer_test, guard, "refit")
        out_of_fold.loc[outer_test] = probability.values
    counts = pd.Series(chosen).value_counts()
    strategy = next(s for s in strategies if s.name == counts.index[0])

    scored = out_of_fold.dropna()
    prediction, support, correct = _support_and_correctness(scored, labels)
    guard.guard("setting_selection", list(scored.index))
    pool_curve = _curve(support, correct.astype(float))
    threshold = _threshold_for_bar(_apply(pool_curve, support), correct, prediction, bar)
    cross_validated = _block("cross_validated (train+tune, nested)", scored, labels, pool_curve,
                             threshold, bar)
    cross_validated["folds"] = len(chosen)
    rows = [gene_rows("cross_validated", scored, labels, pool_curve, threshold)]

    # --- the deployed fit, then calibration on the calibration partition only -------------------
    model, calibration_probability, used = _fit_predict(strategy, features, families, labels,
                                                        pool, calibration_entities, guard, "refit")
    guard.guard("calibration", calibration_entities)
    _pred_c, support_c, correct_c = _support_and_correctness(calibration_probability, labels)
    out = np.full(len(support_c), np.nan)
    try:
        for inner_train, inner_test in _grouped_folds(calibration_entities, labels, groups, 3):
            fold_curve = _curve(support_c.loc[inner_train], correct_c.loc[inner_train].astype(float))
            position = [calibration_entities.index(e) for e in inner_test]
            out[position] = _apply(fold_curve, support_c.loc[inner_test])
    except ValueError:
        pass
    finite = np.isfinite(out)
    calibration_error = float(_ece(out[finite], correct_c.to_numpy(dtype=float)[finite])) \
        if finite.sum() >= 10 else float("nan")
    curve = _curve(support_c, correct_c.astype(float))
    certainty = CertaintyModel(strategy.name, question.id, int(len(support_c)), calibration_error,
                               curve[0], curve[1])
    calibration = {"n": int(len(support_c)), "calibration_error": calibration_error,
                   "cross_fitted": bool(finite.sum() >= 10), "calibrated": certainty.calibrated}

    importance = _family_importance(model, features, families, labels, used, calibration_entities,
                                    repeats=importance_repeats, seed=seed)

    # --- the frozen partition: once, under a recorded budget ------------------------------------
    frozen, refusal = None, ""
    if not question.frozen_test_possible(len(labels)):
        flags.append("no_frozen_test")
        refusal = (f"only {len(labels)} genes have a known answer for this question, and grading "
                   f"needs at least {question.min_labels_for_frozen_test}. With fewer, any set we "
                   f"hide is too small to tell a good strategy from a lucky one. The strategy can "
                   f"still rank genes by how likely the answer is, and can put nothing into the "
                   f"map. More published answers, not a better model, is what lifts this.")
    else:
        prior = holdout_uses(question.id, split.identity, log_path, root)
        changed = [use for use in prior
                   if use.get("strategy_hash") != strategy.identity
                   or use.get("evidence_hash") != evidence]
        if changed:
            flags += ["frozen_test_refused", "needs_fresh_partition"]
            why = ("a changed strategy" if any(use.get("strategy_hash") != strategy.identity
                                               for use in changed) else "changed evidence")
            refusal = (f"the frozen partition of {question.id} (split {split.identity[:12]}) has "
                       f"already been used by {len(changed)} other version(s) "
                       f"({', '.join(sorted({use.get('strategy', '?') for use in changed}))}). A "
                       f"second look with {why} makes it another tuning set, so it is "
                       f"refused: this question needs a fresh partition (change the seed).")
        else:
            if prior:
                flags.append("frozen_test_repeat")
            test_probability = pd.Series(model.predict_proba(features.loc[test, used])[:, 1],
                                         index=pd.Index(test))
            frozen = _block("frozen_test (touched once)", test_probability, labels, curve,
                            threshold, bar)
            rows.append(gene_rows("frozen_test", test_probability, labels, curve, threshold))
            if strata is not None:
                frozen["by_stratum"] = by_stratum(test_probability, labels, curve, threshold,
                                                  pd.Series(strata), bar)
            if record:
                record_holdout_use({"question_id": question.id, "strategy": strategy.name,
                                    "strategy_hash": strategy.identity,
                                    "evidence_hash": evidence,
                                    "split_identity": split.identity, "seed": int(seed),
                                    "n_test": frozen["n"], "calls": frozen["calls"],
                                    "precision": frozen["precision"],
                                    "precision_lower90": frozen["precision_lower90"],
                                    "when": time.strftime("%Y-%m-%dT%H:%M:%S")},
                                   log_path, root)
    if pair_partition is not None:
        flags.append(f"cold_node_pairs:{len(pair_partition.pairs)}"
                     f" kept, {len(pair_partition.excluded_cross_partition)} cross-partition dropped")

    baselines = pd.concat(
        [_baselines(labels, features, groups, pool, list(scored.index), bar).assign(partition="cross_validated")]
        + ([_baselines(labels, features, groups, pool, test, bar).assign(partition="frozen_test")]
           if frozen is not None else []), ignore_index=True)

    missing = set(STAGES_REQUIRED) - guard.stages
    if missing:
        flags.append("guard_calls_missing:" + ",".join(sorted(missing)))
    return HarnessResult(question_id=question.id, answer_type=question.answer_type,
                         question_class=question.question_class, n_labels=int(len(labels)),
                         seed=int(seed), split_identity=split.identity, bar=float(bar),
                         strategy=strategy.name, strategy_hash=strategy.identity,
                         evidence_hash=evidence, evidence_selection=selection,
                         threshold=float(threshold),
                         cross_validated=cross_validated,
                         frozen_test=frozen, inner_scores=inner_scores,
                         chosen_per_outer_fold=tuple(chosen), calibration=calibration,
                         baselines=baselines, importance=importance, flags=tuple(flags),
                         guard_calls=guard.calls, refusal=refusal,
                         rows=pd.concat(rows, ignore_index=True),
                         deployed={"model": model, "columns": list(used), "curve": curve})


# --- plain language -------------------------------------------------------------------------------
# The stored records keep their identifiers so an old run stays readable, but nothing a person reads
# should make them learn one. Every status a run can end in has a sentence here.
PLAIN = {
    'no_frozen_test': 'not enough known answers to grade',
    'frozen_test_refused': 'already graded once; a changed strategy needs a fresh set of hidden genes',
    'frozen_test_repeat': 'graded again with the identical strategy',
    'needs_fresh_partition': 'needs a new set of hidden genes',
}
PLAIN_TERMS = {
    'frozen_test': 'the hidden genes, checked once',
    'cross_validated': 'rotating the hidden genes through all of them',
    'precision_lower90': 'how often it is right, taking the unlucky end of the error bar',
    'may_fill': 'trusted enough to write into the map',
    'min_labels_for_frozen_test': 'the fewest known answers a question needs before it can be graded',
}


def plain(flag: str) -> str:
    """One flag or status in words, for a panel, a README or a sentence to Einar."""
    return PLAIN.get(flag, flag.replace('_', ' '))


def plain_status(flags=(), graded: bool = False, fills: bool = False) -> str:
    """The whole outcome of a run in one sentence, with no vocabulary to look up."""
    if not graded:
        reasons = [plain(f) for f in flags if f in PLAIN] or ['not graded yet']
        return f"ranks genes only -- {'; '.join(reasons)}"
    if fills:
        return 'graded on hidden genes and accurate enough to write into the map'
    return 'graded on hidden genes, but not accurate enough to write into the map'
