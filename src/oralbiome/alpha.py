"""M3 - Alpha diversity: Shannon, Simpson, observed richness, Pielou's evenness."""
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
}


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


def alpha_table(counts: pd.DataFrame) -> pd.DataFrame:
    """All four metrics for every sample (samples are columns of ``counts``)."""
    rows = {}
    for sample in counts.columns:
        c = counts[sample].to_numpy(dtype=float)
        rows[sample] = {"shannon": shannon(c), "simpson": simpson(c),
                        "observed": observed(c), "pielou": pielou(c)}
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


def plot_alpha(alpha: pd.DataFrame, tests: pd.DataFrame, meta: pd.DataFrame,
               disease: str, reference: str) -> None:
    """Four-panel boxplot with individual points, n, p and effect size."""
    fig, axes = plt.subplots(1, 4, figsize=(13, 4.4))
    labels = {code: config.GROUP_LABELS[code] for code in (reference, disease)}
    flagged = meta.get("flagged_atypical", pd.Series(False, index=meta.index))
    for ax, metric in zip(axes, METRICS):
        values, hollow = {}, {}
        for code, lab in labels.items():
            ids = meta.index[meta["group_code"] == code]
            values[lab] = alpha.loc[ids, metric].to_numpy()
            hollow[lab] = flagged.loc[ids].to_numpy()
        style.box_with_points(ax, values, hollow)
        row = tests.set_index("metric").loc[metric]
        ax.set_title(METRICS[metric], fontsize=10.5)
        ax.set_xlabel(f"Mann-Whitney p = {row['p']:.3f}\nrank-biserial r = {row['rank_biserial']:+.2f} "
                      f"[{row['rank_biserial_ci_low']:+.2f}, {row['rank_biserial_ci_high']:+.2f}]",
                      fontsize=8.5, color=style.INK_SECONDARY)
    fig.suptitle(f"{style.title_prefix()}Alpha diversity (genus level, rarefied to "
                 f"equal depth); open circles = flagged atypical samples",
                 fontweight="bold")
    fig.tight_layout()
    style.save(fig, "04_alpha_diversity")


def run(pre: Preprocessed, disease: str = "P", reference: str = "H") -> dict:
    """Compute alpha diversity, the primary test and two sensitivity analyses."""
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
    plot_alpha(genus_alpha, tests[tests["analysis"] == "primary_genus"], meta, disease, reference)

    primary = tests[tests["analysis"] == "primary_genus"].set_index("metric")
    return {m: {k: primary.loc[m, k] for k in ("median_reference", "median_disease", "p",
                                               "q_bh_within_analysis", "rank_biserial",
                                               "rank_biserial_ci_low", "rank_biserial_ci_high")}
            for m in METRICS}
