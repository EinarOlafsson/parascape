"""Answer an Infer question with a harness method instead of map groups (instruction 74, Infer step 3).

The question keeps its recipe: the same label (the recipe holdout), the same evidence after closure (circular
columns and the kinds the user left out removed). What changes is the method. For each yes/no form of the label
(`binary_labels`), the chosen method is graded through `harness.run_question` -- orthogroup-grouped nested
cross-validation and a frozen test touched once -- and then fitted on every labelled gene to rank the genes that
have no label. Every call carries the grade of the question it answers; a question whose frozen test does not reach
the fill bar ranks genes and never writes them as facts.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

EXTREME = 0.2
MIN_CLASS = 50

#: model key -> what the user reads. Map groups (the recipe) are the default and are not a harness model.
METHODS = {
    "hist_gbt": "Boosted trees",
    "logistic": "Logistic regression",
    "pu_bagging": "Positive-unlabelled bagging",
    "propensity_weighted": "Study-effort weighting",
    "low_rank_completion": "Low-rank completion",
    "family_stacking": "One model per evidence family, combined",
    "lookalike": "Lookalike ranking",
    "graph_propagation": "Network propagation",
}


def _safe(text: str) -> str:
    import re
    return re.sub(r"[^A-Za-z0-9]+", "_", str(text)).strip("_")[:80]


def binary_labels(column, kind, category, series):
    """[(suffix, wording, boolean labels over the labelled genes)] for one label column.

    A value column: the extreme 20 % in the direction it means (lowest for a fitness screen, largest absolute
    change for a two-signed contrast, otherwise highest). A class column: its two values, or each class with at
    least 50 genes against the rest.
    """
    s = series.dropna()
    out = []
    if kind == "value":
        v = s.astype(float)
        if category == "fitness":
            cut = v.quantile(EXTREME)
            out.append(("low", f"among the most depleted 20 % ({column} <= {cut:.3g})", v <= cut))
        elif (v < 0).any() and (v > 0).any():
            cut = v.abs().quantile(1 - EXTREME)
            out.append(("change", f"among the largest 20 % of changes (|{column}| >= {cut:.3g})", v.abs() >= cut))
        else:
            cut = v.quantile(1 - EXTREME)
            out.append(("high", f"among the highest 20 % ({column} >= {cut:.3g})", v >= cut))
        return out
    values = s.astype(str)
    counts = values.value_counts()
    if len(counts) == 2:
        pos = "True" if "True" in counts.index else sorted(counts.index)[-1]
        out.append((_safe(pos), f"{column} = {pos}", values == pos))
        return out
    for cls, n in counts.items():
        if MIN_CLASS <= n <= len(values) - MIN_CLASS:
            out.append((_safe(cls), f"{column} = {cls}", values == cls))
    return out


def label_kind(series) -> str:
    """'value' for a numeric column with many values, else 'class'."""
    numeric = pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series)
    return "value" if numeric and series.nunique() > 12 else "class"


def answer(nodes: pd.DataFrame, recipe, model: str = "hist_gbt", organism: str | None = None, seed: int = 0,
           log=None, record: bool = True) -> dict:
    """Grade the method on the recipe question, then rank the unlabelled genes.

    Returns ``{"refusal", "grades": DataFrame, "calls": DataFrame}``: one grade row per yes/no form of the label
    (frozen-test AUC, average precision, the lower bound, may it fill), and the unlabelled genes with the
    probability the method gives each form, best first.
    """
    from . import borrowed_models as B, datasets as ds, harness as H, information_categories as IC
    from . import organisms as O
    from . import recipes as R
    from .question_registry import Question
    say = log or (lambda _m: None)
    organism = organism or getattr(recipe, "organism", O.TOXOPLASMA)
    if model not in METHODS:
        raise ValueError(f"unknown method {model!r}")
    closure = R.close(nodes, recipe)
    if not closure.ok:
        return {"refusal": closure.refusal, "grades": pd.DataFrame(), "calls": pd.DataFrame()}
    table = nodes.set_index("gene_id")
    column = closure.holdout_column
    series = table[column]
    # "unassigned", "unknown" and the like say a gene has no label -- they are the genes to answer, not a class.
    from .search import ABSENCE_LABELS
    if not pd.api.types.is_numeric_dtype(series):
        series = series.where(~series.astype(str).str.lower().isin({str(a).lower() for a in ABSENCE_LABELS}))
    try:
        category = next(iter(IC.categories_of(column, organism)), "")
    except Exception:
        category = ""
    numeric = [c for c in closure.columns if c in table.columns
               and (pd.api.types.is_numeric_dtype(table[c]) or pd.api.types.is_bool_dtype(table[c]))]
    if not numeric:
        return {"refusal": "the evidence the question may use holds no numeric column a model can read",
                "grades": pd.DataFrame(), "calls": pd.DataFrame()}

    def family_of(c):
        found = ds.provenance(c, organism)
        return found.key if found is not None else "unregistered"
    families = pd.Series({c: family_of(c) for c in numeric})
    B.set_families(families)
    groups = table["orthogroup"].astype(object) if "orthogroup" in table else pd.Series(table.index, index=table.index)
    groups = groups.where(groups.notna(), pd.Series({g: "self:" + g for g in table.index}))
    label_family = family_of(column) if column in families.index else column
    grades, calls = [], []
    unlabelled = table.index[series.isna()]
    for suffix, wording, labels in binary_labels(column, label_kind(series), category, series):
        qid = f"infer_{_safe(column)}_{suffix}_{model}"
        q = Question(id=qid, text=f"{recipe.question} -- yes when {wording}", answer_type="boolean",
                     question_class="B", label_source=f"nodes.parquet:{column}", label_family=label_family,
                     evidence_forbidden=(label_family,), label_tier="observed-screen", label_category=category,
                     status="scored")
        t0 = time.monotonic()
        full = table.loc[labels.index, numeric].astype("float64")
        features = H.allowed_features(q, full, families)
        strategy = H.Strategy.of(METHODS[model], model)
        say(f"{METHODS[model]} on '{wording}': grading on held-out genes")
        r = H.run_question(q, labels=labels.astype(bool), features=features, families=families,
                           groups=groups.reindex(labels.index), strategies=[strategy], seed=seed,
                           importance_repeats=1, record=record, log_path=H.HOLDOUT_LOG)
        ft = r.frozen_test or {}
        grades.append(dict(label=wording, labelled=int(len(labels)), positives=int(labels.sum()),
                           frozen_auc=ft.get("roc_auc", np.nan), average_precision=ft.get("average_precision", np.nan),
                           share_of_yes=float(labels.mean()), precision_lower90=ft.get("precision_lower90", np.nan),
                           may_fill=bool(r.may_fill), flags=", ".join(r.flags), seconds=round(time.monotonic() - t0, 1)))
        if len(unlabelled):
            # Fitted exactly as the harness fits: the allowed features, columns with any value, integer labels.
            used = [c for c in strategy.columns(features, families) if features[c].notna().any()]
            est = strategy.build().fit(features[used], labels.astype(int))
            prob = est.predict_proba(table.loc[unlabelled, used].astype("float64"))[:, 1]
            calls.append(pd.DataFrame({"gene_id": unlabelled.astype(str), "answer": wording,
                                       "probability": prob, "may_fill": bool(r.may_fill)}))
    calls = pd.concat(calls, ignore_index=True).sort_values("probability", ascending=False) if calls else pd.DataFrame()
    return {"refusal": "", "grades": pd.DataFrame(grades), "calls": calls.reset_index(drop=True)}
