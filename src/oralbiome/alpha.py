"""M3 - Alpha diversity: Shannon, Simpson, observed richness, Pielou's evenness,
Chao1 and ACE.

Chao1 and ACE were added so the results can be compared directly with the
source study, which reports those two richness estimators.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config, style
from .preprocessing import Preprocessed
from .stats import benjamini_hochberg, compare_two_groups

METRICS: dict[str, str] = {
    "shannon": "Shannon index (H')",
    "simpson": "Gini-Simpson index (1 - D)",
    "observed": "Observed richness (taxa)",
    "pielou": "Pielou's evenness (J')",
    "chao1": "Chao1 richness estimate",
    "ace": "ACE richness estimate",
}
ACE_RARE_THRESHOLD: int = 10


def shannon(counts: np.ndarray) -> float:
    """H' = -sum p_i ln p_i over taxa with p_i > 0 (natural log)."""
    p = counts[counts > 0] / counts.sum()
    return float(-(p * np.log(p)).sum())


def simpson(counts: np.ndarray) -> float:
    """Gini-Simpson 1 - sum p_i^2: chance two random reads are different taxa."""
    p = counts / counts.sum()
    return float(1.0 - (p**2).sum())


def observed(counts: np.ndarray) -> int:
    """Number of taxa with at least one read."""
    return int((counts > 0).sum())


def pielou(counts: np.ndarray) -> float:
    """J' = H' / ln(S): Shannon relative to its maximum for S taxa."""
    s = observed(counts)
    return float(shannon(counts) / np.log(s)) if s > 1 else 0.0


def chao1(counts: np.ndarray) -> float:
    """Bias-corrected Chao1: S_obs + F1 (F1 - 1) / (2 (F2 + 1)).

    F1 and F2 are the numbers of taxa seen exactly once and twice. Many
    singletons relative to doubletons suggest many taxa were missed.
    """
    c = np.asarray(counts)
    f1, f2 = int((c == 1).sum()), int((c == 2).sum())
    return float(observed(c) + f1 * (f1 - 1) / (2 * (f2 + 1)))


def ace(counts: np.ndarray, rare_threshold: int = ACE_RARE_THRESHOLD) -> float:
    """Abundance-based Coverage Estimator (Chao & Lee 1992).

    Taxa with <= ``rare_threshold`` reads are "rare". With N_rare reads in rare
    taxa and F_i taxa seen i times: coverage C = 1 - F1 / N_rare,
    gamma^2 = max(S_rare / C * sum_i i(i-1) F_i / (N_rare (N_rare - 1)) - 1, 0),
    ACE = S_abund + S_rare / C + F1 / C * gamma^2.
    If every rare read is a singleton (C = 0) ACE is undefined; bias-corrected
    Chao1 is returned instead (logged in DECISIONS.md).
    """
    c = np.asarray(counts)
    c = c[c > 0]
    rare = c[c <= rare_threshold]
    s_abund = int((c > rare_threshold).sum())
    s_rare, n_rare = len(rare), float(rare.sum())
    if s_rare == 0:
        return float(s_abund)
    f1 = int((rare == 1).sum())
    coverage = 1 - f1 / n_rare
    if coverage <= 0:
        return chao1(counts)
    i = np.arange(1, rare_threshold + 1)
    f_i = np.array([(rare == k).sum() for k in i])
    denom = n_rare * (n_rare - 1)
    gamma2 = max(s_rare / coverage * (i * (i - 1) * f_i).sum() / denom - 1, 0.0) if denom else 0.0
    return float(s_abund + s_rare / coverage + f1 / coverage * gamma2)


def alpha_table(counts: pd.DataFrame) -> pd.DataFrame:
    """All six metrics for every sample (samples are columns of ``counts``)."""
    rows = {}
    for sample in counts.columns:
        c = counts[sample].to_numpy(dtype=float)
        rows[sample] = {"shannon": shannon(c), "simpson": simpson(c), "observed": observed(c),
                        "pielou": pielou(c), "chao1": chao1(c), "ace": ace(c)}
    return pd.DataFrame.from_dict(rows, orient="index")


def compare_alpha(alpha: pd.DataFrame, meta: pd.DataFrame, disease: str, reference: str,
                  label: str) -> pd.DataFrame:
    """Mann-Whitney U per metric with rank-biserial effect size and BH q-values."""
    groups = meta.loc[alpha.index, "group_code"]
    rows = []
    for metric in METRICS:
        res = compare_two_groups(alpha.loc[groups == disease, metric],
                                 alpha.loc[groups == reference, metric])
        rows.append({"analysis": label, "metric": metric, **res})
    table = pd.DataFrame(rows)
    table["q_bh_within_analysis"] = benjamini_hochberg(table["p"].to_numpy())
    return table


def high_richness_check(genus_alpha: pd.DataFrame, asv_alpha: pd.DataFrame,
                        meta: pd.DataFrame, disease: str, reference: str,
                        min_features: int = 1000) -> pd.DataFrame:
    """Post-hoc diagnostic: richness tests without the high-richness healthy samples.

    Healthy samples with > ``min_features`` observed features at rarefaction
    depth (the subset described in M5) are removed and the richness
    estimators re-tested. Not part of any BH family; descriptive only.
    """
    high = asv_alpha["observed"] > min_features
    keep = ~(high & (meta.loc[asv_alpha.index, "group_code"] == reference))
    rows = []
    for level, table in (("genus", genus_alpha), ("feature", asv_alpha)):
        for metric in ("observed", "chao1", "ace"):
            for subset, mask in (("all samples", slice(None)), ("without high-richness H", keep)):
                sub = table.loc[mask]
                g = meta.loc[sub.index, "group_code"]
                res = compare_two_groups(sub.loc[g == disease, metric], sub.loc[g == reference, metric],
                                         with_ci=False)
                rows.append({"level": level, "metric": metric, "subset": subset, **res})
    return pd.DataFrame(rows)


def plot_alpha(alpha: pd.DataFrame, tests: pd.DataFrame, meta: pd.DataFrame,
               disease: str, reference: str) -> None:
    """Six-panel boxplot with individual points, n, p, q and effect size."""
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.6))
    labels = {code: config.GROUP_LABELS[code] for code in (reference, disease)}
    flagged = meta.get("flagged_atypical", pd.Series(False, index=meta.index))
    for ax, metric in zip(axes.ravel(), METRICS):
        values, hollow = {}, {}
        for code, lab in labels.items():
            ids = meta.index[meta["group_code"] == code]
            values[lab] = alpha.loc[ids, metric].to_numpy()
            hollow[lab] = flagged.loc[ids].to_numpy()
        style.box_with_points(ax, values, hollow)
        row = tests.set_index("metric").loc[metric]
        ax.set_title(METRICS[metric], fontsize=10.5)
        ax.set_xlabel(f"p = {row['p']:.3f}, q = {row['q_bh_within_analysis']:.3f}\n"
                      f"rank-biserial r = {row['rank_biserial']:+.2f} "
                      f"[{row['rank_biserial_ci_low']:+.2f}, {row['rank_biserial_ci_high']:+.2f}]",
                      fontsize=8.5, color=style.INK_SECONDARY)
    fig.suptitle(f"{style.title_prefix()}Alpha diversity (genus level, rarefied to equal depth; "
                 f"q = BH across 6 metrics); open circles = flagged atypical samples",
                 fontweight="bold")
    fig.tight_layout()
    style.save(fig, "04_alpha_diversity")


def run(pre: Preprocessed, disease: str = "P", reference: str = "H") -> dict:
    """Compute alpha diversity, the primary test, sensitivity analyses and a diagnostic."""
    meta = pre.metadata
    genus_alpha = alpha_table(pre.genus_rarefied)
    asv_alpha = alpha_table(pre.asv_rarefied)
    flagged = meta.get("flagged_atypical", pd.Series(False, index=meta.index)).astype(bool)

    tests = pd.concat([
        compare_alpha(genus_alpha, meta, disease, reference, "primary_genus"),
        compare_alpha(asv_alpha, meta, disease, reference, "sensitivity_asv_level"),
        compare_alpha(genus_alpha.loc[~flagged], meta, disease, reference,
                      "sensitivity_genus_excluding_flagged"),
    ], ignore_index=True)

    per_sample = genus_alpha.add_prefix("genus_").join(asv_alpha.add_prefix("asv_"))
    per_sample = per_sample.join(meta[["group", "flagged_atypical"]] if "flagged_atypical"
                                 in meta else meta[["group"]])
    per_sample.to_csv(config.RESULTS / "m3_alpha_per_sample.csv")
    tests.to_csv(config.RESULTS / "m3_alpha_tests.csv", index=False)
    high_richness_check(genus_alpha, asv_alpha, meta, disease, reference).to_csv(
        config.RESULTS / "m3_richness_high_richness_check.csv", index=False)
    plot_alpha(genus_alpha, tests[tests["analysis"] == "primary_genus"], meta, disease, reference)

    primary = tests[tests["analysis"] == "primary_genus"].set_index("metric")
    return {m: {k: primary.loc[m, k] for k in ("median_reference", "median_disease", "p",
                                               "q_bh_within_analysis", "rank_biserial",
                                               "rank_biserial_ci_low", "rank_biserial_ci_high")}
            for m in METRICS}
