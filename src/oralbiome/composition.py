"""M2 - Composition: stacked relative-abundance bars and the core microbiome."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config, style
from .preprocessing import Preprocessed, relative_abundance

N_SHOWN = len(style.TAXON_COLORS)  # 8 named taxa + "Other"


def top_taxa_table(rel: pd.DataFrame, n: int = N_SHOWN) -> pd.DataFrame:
    """Keep the ``n`` taxa with highest mean abundance; sum the rest as 'Other'."""
    order = rel.mean(axis=1).sort_values(ascending=False).index
    top = rel.loc[order[:n]]
    other = rel.loc[order[n:]].sum(axis=0).rename("Other")
    return pd.concat([top, other.to_frame().T])


def group_means(rel: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Mean relative abundance of each taxon per group (taxa x groups)."""
    return rel.T.groupby(meta.loc[rel.columns, "group"]).mean().T


def core_microbiome(counts: pd.DataFrame, meta: pd.DataFrame,
                    threshold: float = config.CORE_PREVALENCE) -> pd.DataFrame:
    """Per-group prevalence of every taxon and core membership (prevalence >= threshold)."""
    present = counts > 0
    prevalence = present.T.groupby(meta.loc[counts.columns, "group"]).mean().T
    rel = relative_abundance(counts)
    table = prevalence.add_prefix("prevalence_")
    for group in prevalence.columns:
        table[f"core_{group}"] = prevalence[group] >= threshold
    table["mean_rel_abundance"] = rel.mean(axis=1)
    core_cols = [c for c in table.columns if c.startswith("core_")]
    return table.loc[table[core_cols].any(axis=1)].sort_values("mean_rel_abundance", ascending=False)


def _stack(ax: plt.Axes, table: pd.DataFrame, labels: list[str], width: float) -> None:
    bottom = np.zeros(table.shape[1])
    colors = style.TAXON_COLORS + [style.OTHER_COLOR]
    for i, taxon in enumerate(table.index):
        values = table.loc[taxon].to_numpy() * 100
        ax.bar(range(len(labels)), values, bottom=bottom, width=width, color=colors[i],
               edgecolor="white", linewidth=0.4, label=taxon)
        bottom += values
    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.set_ylim(0, 100)
    ax.grid(axis="x", visible=False)


def plot_composition(rel: pd.DataFrame, meta: pd.DataFrame, level: str, name: str) -> None:
    """Two-panel stacked bar chart: every sample, then the group means."""
    table = top_taxa_table(rel)
    groups = [g for g in config.GROUP_LABELS.values() if g in set(meta["group"])]
    first = table.index[0]
    ordered = []
    for g in groups:
        ids = meta.index[meta["group"] == g]
        ordered += list(table.loc[first, ids].sort_values(ascending=False).index)
    sample_table = table[ordered]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), width_ratios=[4.2, 1],
                                   sharey=True)
    _stack(ax1, sample_table, ordered, width=0.9)
    ax1.set_xticks(range(len(ordered)))
    flagged = meta.get("flagged_atypical", pd.Series(False, index=meta.index))
    ax1.set_xticklabels([f"{s}*" if flagged.get(s, False) else s for s in ordered],
                        rotation=90, fontsize=6.5)
    ax1.set_ylabel("Relative abundance (%)")
    ax1.set_title(f"Per sample ({level} level; * = flagged atypical sample)", pad=22)
    boundary = 0
    for g in groups[:-1]:
        boundary += int((meta["group"] == g).sum())
        ax1.axvline(boundary - 0.5, color=style.INK, lw=1.2)
    start = 0
    for g in groups:
        n = int((meta["group"] == g).sum())
        ax1.text(start + n / 2 - 0.5, 101.5, style.group_label_with_n(g, n), ha="center",
                 fontsize=9, color=style.GROUP_COLORS.get(g, style.INK), fontweight="bold")
        start += n

    means = group_means(table, meta)[groups]
    _stack(ax2, means, groups, width=0.7)
    ax2.set_xticks(range(len(groups)))
    ax2.set_xticklabels([f"{g}\n(n={int((meta['group'] == g).sum())})" for g in groups], fontsize=8)
    ax2.set_title("Group mean", pad=22)
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[::-1], labels[::-1], loc="center left", bbox_to_anchor=(1.02, 0.5),
                title=f"Top {N_SHOWN} {'genera' if level == 'genus' else 'phyla'}", fontsize=8)
    fig.suptitle(f"{style.title_prefix()}Salivary bacterial composition, {level} level",
                 fontweight="bold", y=1.0)
    fig.tight_layout()
    style.save(fig, name)


def plot_core(core: pd.DataFrame, meta: pd.DataFrame, n_show: int = 25) -> None:
    """Dot plot of per-group prevalence for the most abundant core genera."""
    groups = [g for g in config.GROUP_LABELS.values() if f"prevalence_{g}" in core.columns]
    show = core.head(n_show).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 0.28 * len(show) + 1.5))
    offsets = np.linspace(-0.15, 0.15, len(groups))
    for off, g in zip(offsets, groups):
        n = int((meta["group"] == g).sum())
        ax.scatter(show[f"prevalence_{g}"] * 100, np.arange(len(show)) + off, s=36,
                   color=style.GROUP_COLORS[g], label=style.group_label_with_n(g, n), zorder=3)
    ax.axvline(config.CORE_PREVALENCE * 100, color=style.INK_SECONDARY, ls="--", lw=1)
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels(show.index, fontsize=8)
    ax.set_xlabel("Prevalence: % of samples where the genus is detected")
    ax.set_title(f"{style.title_prefix()}Core genera (detected in >= "
                 f"{config.CORE_PREVALENCE:.0%} of a group)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    ax.set_title(ax.get_title(), pad=30)
    fig.tight_layout()
    style.save(fig, "03_core_microbiome")


def run(pre: Preprocessed) -> dict:
    """Write composition tables and figures; return core-microbiome counts."""
    meta = pre.metadata
    genus_rel_all = relative_abundance(pre.genus_counts)
    for level, rel in (("phylum", pre.phylum_rel), ("genus", genus_rel_all)):
        group_means(rel, meta).sort_values(meta["group"].iloc[0], ascending=False).to_csv(
            config.RESULTS / f"m2_{level}_group_mean_rel_abundance.csv")
    plot_composition(pre.phylum_rel, meta, "phylum", "01_composition_phylum")
    plot_composition(genus_rel_all, meta, "genus", "02_composition_genus")

    core = core_microbiome(pre.genus_counts, meta)
    core.to_csv(config.RESULTS / "m2_core_microbiome.csv")
    plot_core(core, meta)
    flags = [c for c in core.columns if c.startswith("core_")]
    summary = {c.replace("core_", "core_genera_"): int(core[c].sum()) for c in flags}
    summary["core_genera_shared_by_all_groups"] = int(core[flags].all(axis=1).sum())
    return summary
