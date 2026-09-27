"""M7 - Classifier: can genus composition predict healthy vs periodontitis?

Leakage control: the model receives raw genus counts. Prevalence filtering,
the CLR transform and scaling all live inside a scikit-learn Pipeline, so
they are re-fitted on the training part of every fold and never see the
test samples. L1 logistic regression tunes its penalty in an inner CV loop
(nested CV); the random forest uses fixed, standard settings.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from . import config


class PrevalenceCLR(BaseEstimator, TransformerMixin):
    """Prevalence filter learned on training data, then CLR per sample.

    ``fit`` decides which taxa to keep using only the training samples;
    ``transform`` applies that same column set to any samples.
    """

    def __init__(self, min_prevalence: float = config.MIN_PREVALENCE,
                 min_mean_rel: float = config.MIN_MEAN_REL_ABUND,
                 pseudocount: float = config.PSEUDOCOUNT):
        self.min_prevalence = min_prevalence
        self.min_mean_rel = min_mean_rel
        self.pseudocount = pseudocount

    def fit(self, X, y=None):
        X = np.asarray(X, dtype=float)
        prevalence = (X > 0).mean(axis=0)
        rel = X / X.sum(axis=1, keepdims=True)
        self.keep_ = (prevalence >= self.min_prevalence) & (rel.mean(axis=0) >= self.min_mean_rel)
        return self

    def transform(self, X):
        logged = np.log(np.asarray(X, dtype=float)[:, self.keep_] + self.pseudocount)
        return logged - logged.mean(axis=1, keepdims=True)


def make_models(seed: int = config.SEED) -> dict[str, Pipeline | GridSearchCV]:
    """The two model pipelines (unfitted)."""
    lasso = Pipeline([
        ("clr", PrevalenceCLR()),
        ("scale", StandardScaler()),
        ("model", LogisticRegression(l1_ratio=1.0, solver="liblinear", class_weight="balanced",
                                     max_iter=5000, random_state=seed)),
    ])
    lasso_cv = GridSearchCV(lasso, {"model__C": np.logspace(-2, 1.5, 8)}, scoring="roc_auc",
                            cv=StratifiedKFold(3, shuffle=True, random_state=seed), n_jobs=1)
    forest = Pipeline([
        ("clr", PrevalenceCLR()),
        ("model", RandomForestClassifier(n_estimators=config.RF_TREES, max_features="sqrt",
                                         class_weight="balanced", random_state=seed, n_jobs=1)),
    ])
    return {"L1 logistic regression": lasso_cv, "Random forest": forest}


def _importance(fitted, feature_names: np.ndarray) -> pd.Series:
    """Per-feature importance of one fitted model (coef for LR, impurity for RF)."""
    pipe = fitted.best_estimator_ if hasattr(fitted, "best_estimator_") else fitted
    names = feature_names[pipe.named_steps["clr"].keep_]
    model = pipe.named_steps["model"]
    values = model.coef_[0] if hasattr(model, "coef_") else model.feature_importances_
    return pd.Series(values, index=names)


def _fit_fold(model, X, y, train, test, names, want_importance):
    fitted = clone(model).fit(X[train], y[train])
    prob = fitted.predict_proba(X[test])[:, 1]
    return test, prob, (_importance(fitted, names) if want_importance else None)


def repeated_cv(model, X: np.ndarray, y: np.ndarray, names: np.ndarray,
                n_repeats: int = config.CV_REPEATS, n_folds: int = config.CV_FOLDS,
                seed: int = config.SEED, n_jobs: int = -1, want_importance: bool = True) -> dict:
    """Repeated stratified k-fold CV; returns out-of-fold probabilities and AUCs."""
    jobs = []
    for r in range(n_repeats):
        splitter = StratifiedKFold(n_folds, shuffle=True, random_state=seed + r)
        for train, test in splitter.split(X, y):
            jobs.append((r, train, test))
    out = Parallel(n_jobs=n_jobs)(
        delayed(_fit_fold)(model, X, y, tr, te, names, want_importance) for _, tr, te in jobs)
    probs = np.zeros((n_repeats, len(y)))
    importances = []
    for (r, _, _), (test, prob, imp) in zip(jobs, out):
        probs[r, test] = prob
        if imp is not None:
            importances.append(imp)
    aucs = np.array([roc_auc_score(y, probs[r]) for r in range(n_repeats)])
    return {"probs": probs, "aucs": aucs, "importances": importances}


def bootstrap_auc(y: np.ndarray, score: np.ndarray, n_boot: int = config.N_BOOTSTRAP,
                  seed: int = config.SEED, grid: np.ndarray | None = None) -> dict:
    """Bootstrap 95% CI for AUC and a pointwise band for the ROC curve.

    Samples are resampled with replacement within each class so every
    resample contains both classes.
    """
    rng = np.random.default_rng(seed)
    grid = np.linspace(0, 1, 101) if grid is None else grid
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    aucs, tprs = [], []
    for _ in range(n_boot):
        idx = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        aucs.append(roc_auc_score(y[idx], score[idx]))
        fpr, tpr, _ = roc_curve(y[idx], score[idx])
        tprs.append(np.interp(grid, fpr, tpr))
    fpr, tpr, _ = roc_curve(y, score)
    return {
        "auc": float(roc_auc_score(y, score)),
        "ci_low": float(np.percentile(aucs, 2.5)), "ci_high": float(np.percentile(aucs, 97.5)),
        "grid": grid, "tpr_mean": np.interp(grid, fpr, tpr),
        "tpr_low": np.percentile(tprs, 2.5, axis=0), "tpr_high": np.percentile(tprs, 97.5, axis=0),
    }


def _null_one(model, X, y, names, n_repeats, seed):
    rng = np.random.default_rng(seed)
    y_perm = rng.permutation(y)
    res = repeated_cv(model, X, y_perm, names, n_repeats=n_repeats, seed=seed, n_jobs=1,
                      want_importance=False)
    return float(res["aucs"].mean())


def permutation_null(model, X: np.ndarray, y: np.ndarray, names: np.ndarray,
                     n_perm: int = config.N_LABEL_PERMUTATIONS, n_repeats: int = 3,
                     seed: int = config.SEED) -> np.ndarray:
    """Mean CV AUC for ``n_perm`` label-shuffled datasets (same pipeline)."""
    return np.array(Parallel(n_jobs=-1)(
        delayed(_null_one)(model, X, y, names, n_repeats, seed + 1000 + i) for i in range(n_perm)))


def summarise_importance(importances: list[pd.Series], kind: str) -> pd.DataFrame:
    """Aggregate fold-level importances (missing = 0 because the taxon was filtered)."""
    table = pd.concat(importances, axis=1).fillna(0.0)
    out = pd.DataFrame({"mean": table.mean(axis=1), "sd": table.std(axis=1)})
    if kind == "coef":
        out["selection_frequency"] = (table != 0).mean(axis=1)
        out = out.reindex(out["mean"].abs().sort_values(ascending=False).index)
    else:
        out = out.sort_values("mean", ascending=False)
    return out
