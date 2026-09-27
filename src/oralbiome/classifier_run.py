"""M7 runner: cross-validated AUCs, permutation null, importances and figures."""
from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import classifier as clf
from . import config, style
from .preprocessing import Preprocessed

MODEL_COLORS = {"L1 logistic regression": "#4a3aa7", "Random forest": "#1baf7a"}


def plot_roc(boots: dict, n: dict[str, int]) -> None:
    """ROC curves with bootstrap 95% bands for both models."""
    fig, ax = plt.subplots(figsize=(6.4, 6))
    for name, b in boots.items():
        color = MODEL_COLORS[name]
        ax.fill_between(b["grid"], b["tpr_low"], b["tpr_high"], color=color, alpha=0.15, lw=0)
        ax.plot(b["grid"], b["tpr_mean"], color=color, lw=2,
                label=f"{name}: AUC {b['auc']:.2f} [{b['ci_low']:.2f}, {b['ci_high']:.2f}]")
    ax.plot([0, 1], [0, 1], ls="--", color=style.NEUTRAL_COLOR, lw=1, label="chance (AUC 0.50)")
    ax.set_xlabel("False positive rate (healthy called periodontitis)")
    ax.set_ylabel("True positive rate (periodontitis detected)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.01)
    ax.set_title(f"{style.title_prefix()}Out-of-fold ROC, {config.CV_REPEATS}x repeated "
                 f"{config.CV_FOLDS}-fold CV\nHealthy n={n['H']}, periodontitis n={n['P']}; "
                 f"shaded = bootstrap 95% band", fontsize=10.5)
    ax.legend(loc="lower right", fontsize=8.5)
    fig.tight_layout()
    style.save(fig, "12_roc_curves")


def plot_importance(lasso: pd.DataFrame, forest: pd.DataFrame, n_show: int = 15) -> None:
    """L1 coefficients (with selection frequency) and random-forest importances."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6.2))
    top = lasso.head(n_show).iloc[::-1]
    colors = [style.GROUP_COLORS["Periodontitis"] if v > 0 else style.GROUP_COLORS["Healthy"]
              for v in top["mean"]]
    ax1.barh(range(len(top)), top["mean"], xerr=top["sd"], color=colors,
             error_kw={"ecolor": style.INK_SECONDARY, "lw": 0.8})
    ax1.set_yticks(range(len(top)))
    ax1.set_yticklabels([f"{t}  ({f:.0%})" for t, f in zip(top.index, top["selection_frequency"])],
                        fontsize=8.5)
    ax1.axvline(0, color=style.INK_SECONDARY, lw=0.8)
    ax1.set_xlabel("Mean standardised coefficient across folds (+/- SD)\n"
                   "orange = associated with periodontitis, blue = with healthy")
    ax1.set_title("L1 logistic regression\n(% = folds in which the genus was selected)",
                  fontsize=10.5)

    top = forest.head(n_show).iloc[::-1]
    ax2.barh(range(len(top)), top["mean"], xerr=top["sd"], color=style.NEUTRAL_COLOR,
             error_kw={"ecolor": style.INK_SECONDARY, "lw": 0.8})
    ax2.set_yticks(range(len(top)))
    ax2.set_yticklabels(top.index, fontsize=8.5)
    ax2.set_xlabel("Mean decrease in impurity across folds (+/- SD)")
    ax2.set_title("Random forest\n(importance has no direction)", fontsize=10.5)
    fig.suptitle(f"{style.title_prefix()}Genera the classifiers rely on (fold-level models, "
                 f"{config.CV_REPEATS * config.CV_FOLDS} fits each); associational, not causal",
                 fontweight="bold")
    fig.tight_layout()
    style.save(fig, "13_feature_importance")


def plot_null(results: dict) -> None:
    """Histogram of permuted-label AUCs with the observed AUC marked."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, (name, r) in zip(axes, results.items()):
        ax.hist(r["null"], bins=20, color=style.OTHER_COLOR, edgecolor="white")
        ax.axvline(r["observed"], color=MODEL_COLORS[name], lw=2.5)
        ax.axvline(0.5, color=style.INK_SECONDARY, ls="--", lw=1)
        ax.text(r["observed"], ax.get_ylim()[1] * 0.95, f" observed\n {r['observed']:.2f}",
                color=MODEL_COLORS[name], va="top", fontsize=9, fontweight="bold")
        ax.set_title(f"{name}\npermutation p = {r['p']:.3f} ({len(r['null'])} shuffles)",
                     fontsize=10.5)
        ax.set_xlabel("Mean cross-validated AUC")
    axes[0].set_ylabel("Number of label shuffles")
    fig.suptitle(f"{style.title_prefix()}Negative control: same pipeline on shuffled labels",
                 fontweight="bold")
    fig.tight_layout()
    style.save(fig, "14_permutation_null")


def run(pre: Preprocessed, disease: str = "P", reference: str = "H",
        n_perm: int = config.N_LABEL_PERMUTATIONS, figures: bool = True, tag: str = "primary") -> dict:
    """Run M7 on one comparison. Raw genus counts go in; everything else is in-fold."""
    meta = pre.metadata
    X = pre.genus_counts[meta.index].T.to_numpy(dtype=float)
    y = (meta["group_code"] == disease).to_numpy(dtype=int)
    names = pre.genus_counts.index.to_numpy()
    flagged = meta.get("flagged_atypical", pd.Series(False, index=meta.index)).to_numpy(bool)

    summary, boots, nulls, importances = {}, {}, {}, {}
    for name, model in clf.make_models().items():
        cv = clf.repeated_cv(model, X, y, names)
        mean_prob = cv["probs"].mean(axis=0)
        boots[name] = clf.bootstrap_auc(y, mean_prob)
        null = clf.permutation_null(model, X, y, names, n_perm=n_perm)
        observed = float(cv["aucs"].mean())
        p = (1 + np.sum(null >= observed)) / (1 + len(null))
        nulls[name] = {"null": null, "observed": observed, "p": float(p)}
        sens = clf.repeated_cv(model, X[~flagged], y[~flagged], names, n_repeats=5,
                               want_importance=False)
        kind = "coef" if "logistic" in name else "impurity"
        importances[name] = clf.summarise_importance(cv["importances"], kind)
        summary[name] = {
            "auc_mean_over_repeats": observed,
            "auc_repeat_range": [float(cv["aucs"].min()), float(cv["aucs"].max())],
            "auc_pooled_oof": boots[name]["auc"],
            "auc_bootstrap_ci95": [boots[name]["ci_low"], boots[name]["ci_high"]],
            "null_auc_mean": float(null.mean()),
            "null_auc_95th_percentile": float(np.percentile(null, 95)),
            "permutation_p": float(p),
            "n_permutations": len(null),
            "auc_excluding_flagged_mean_5_repeats": float(sens["aucs"].mean()),
            "n_excluding_flagged": int((~flagged).sum()),
        }
        pd.DataFrame({"sample_id": meta.index, "label": y, "mean_oof_probability": mean_prob}) \
            .to_csv(config.RESULTS / f"m7_oof_predictions_{kind}_{tag}.csv", index=False)
        importances[name].to_csv(config.RESULTS / f"m7_importance_{kind}_{tag}.csv")
        np.savetxt(config.RESULTS / f"m7_null_aucs_{kind}_{tag}.csv", null, delimiter=",")

    # Ablation: does the signal survive without the single most-used genus?
    top = importances["L1 logistic regression"].index[0]
    keep = names != top
    summary["ablation_without_top_genus"] = {"removed": str(top)}
    for name, model in clf.make_models().items():
        abl = clf.repeated_cv(model, X[:, keep], y, names[keep], n_repeats=5, want_importance=False)
        summary["ablation_without_top_genus"][name] = float(abl["aucs"].mean())

    summary["n"] = {"H": int((y == 0).sum()), "P": int((y == 1).sum())}
    summary["cv"] = {"folds": config.CV_FOLDS, "repeats": config.CV_REPEATS,
                     "inner_folds_for_C": 3, "null_repeats": 3,
                     "note": "auc_pooled_oof = AUC of per-sample probabilities averaged over repeats "
                             "(shown in ROC figure, bootstrap CI); auc_mean_over_repeats = mean "
                             "of per-repeat AUCs (compared with the permutation null)"}
    (config.RESULTS / f"m7_classifier_summary_{tag}.json").write_text(json.dumps(summary, indent=2))
    if figures:
        plot_roc(boots, summary["n"])
        plot_importance(importances["L1 logistic regression"], importances["Random forest"])
        plot_null(nulls)
    return summary
