"""Data quality report: read depth, sparsity, group sizes and rarefaction curves."""
from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import gammaln

from . import config, style
from .data import Dataset, quality_report


def expected_richness(counts: np.ndarray, depth: int) -> float:
    """Expected number of taxa observed when subsampling ``depth`` reads.

    Uses the exact hypergeometric formula (Hurlbert 1971):
    E[S_n] = sum_i [1 - C(N - N_i, n) / C(N, n)].
    """
    counts = counts[counts > 0]
    total = counts.sum()
    if depth >= total:
        return float(len(counts))
    log_denominator = gammaln(total + 1) - gammaln(depth + 1) - gammaln(total - depth + 1)
    remaining = total - counts
    ok = remaining >= depth
    log_numerator = np.full(len(counts), -np.inf)
    r = remaining[ok]
    log_numerator[ok] = gammaln(r + 1) - gammaln(depth + 1) - gammaln(r - depth + 1)
    return float(np.sum(1.0 - np.exp(log_numerator - log_denominator)))


def rarefaction_curves(ds: Dataset, n_points: int = 15) -> pd.DataFrame:
    """Expected ASV richness at increasing depths for every sample."""
    max_depth = int(ds.counts.sum(axis=0).max())
    depths = np.unique(np.linspace(100, max_depth, n_points).astype(int))
    rows = []
    for sample in ds.counts.columns:
        col = ds.counts[sample].to_numpy()
        for d in depths[depths <= col.sum()]:
            rows.append((sample, int(d), expected_richness(col, int(d))))
    return pd.DataFrame(rows, columns=["sample_id", "depth", "expected_richness"])


def plot_qc(ds: Dataset, curves: pd.DataFrame) -> None:
    """Two-panel figure: read depth by group and rarefaction curves."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    meta = ds.metadata
    order = [g for g in config.GROUP_LABELS.values() if g in set(meta["group"])]
    rng = np.random.default_rng(config.SEED)
    for i, label in enumerate(order):
        depths = meta.loc[meta["group"] == label, "read_depth"] / 1000
        color = style.GROUP_COLORS.get(label, style.NEUTRAL_COLOR)
        axes[0].boxplot(depths, positions=[i], widths=0.5, showfliers=False,
                        medianprops={"color": style.INK})
        axes[0].scatter(i + rng.uniform(-0.12, 0.12, len(depths)), depths, s=22,
                        color=color, alpha=0.8, edgecolor="white", linewidth=0.5, zorder=3)
    axes[0].set_xticks(range(len(order)))
    short = {"Hypertension + periodontitis": "Hypertension +\nperiodontitis"}
    axes[0].set_xticklabels([f"{short.get(g, g)}\n(n={int((meta['group'] == g).sum())})"
                             for g in order], fontsize=8)
    axes[0].set_ylabel("Reads per sample (thousands)")
    axes[0].set_title("Sequencing depth by group")

    for sample, sub in curves.groupby("sample_id"):
        label = meta.loc[sample, "group"]
        axes[1].plot(sub["depth"] / 1000, sub["expected_richness"], lw=0.9, alpha=0.6,
                     color=style.GROUP_COLORS.get(label, style.NEUTRAL_COLOR))
    for label in order:
        axes[1].plot([], [], color=style.GROUP_COLORS.get(label, style.NEUTRAL_COLOR), lw=2,
                     label=f"{label} (n={int((meta['group'] == label).sum())})")
    axes[1].legend(loc="center right", fontsize=8)
    axes[1].axvline(meta["read_depth"].min() / 1000, color=style.INK_SECONDARY, ls="--", lw=1)
    axes[1].text(meta["read_depth"].min() / 1000, axes[1].get_ylim()[1] * 0.02,
                 " rarefaction depth", fontsize=8, color=style.INK_SECONDARY)
    axes[1].set_xlabel("Reads subsampled (thousands)")
    axes[1].set_ylabel("Expected ASVs observed")
    axes[1].set_title("Rarefaction curves (all 67 samples)" if len(meta) == 67 else "Rarefaction curves")
    fig.suptitle(f"{style.title_prefix()}Data quality overview", fontweight="bold")
    fig.tight_layout()
    style.save(fig, "00_data_quality")


def run(ds: Dataset) -> dict:
    """Write results/data_quality.json and the QC figure."""
    report = quality_report(ds)
    curves = rarefaction_curves(ds)
    # How close each curve is to flat at the rarefaction depth: slope of the
    # last segment in new ASVs per 1,000 reads.
    last = curves.sort_values("depth").groupby("sample_id").tail(2)
    slopes = last.groupby("sample_id").apply(
        lambda s: np.diff(s["expected_richness"]).item() / max(np.diff(s["depth"]).item(), 1) * 1000,
        include_groups=False,
    )
    report["rarefaction_end_slope_asvs_per_1000_reads_median"] = float(slopes.median())
    curves.to_csv(config.RESULTS / "qc_rarefaction_curves.csv", index=False)
    (config.RESULTS / "data_quality.json").write_text(json.dumps(report, indent=2))
    plot_qc(ds, curves)
    return report
