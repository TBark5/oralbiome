"""M5 runner: tables, null control, pathogen cross-check and figures."""
from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config, style
from .differential import differential_table, permutation_null, targeted_species_test
from .preprocessing import Preprocessed, relative_abundance


def plot_volcano(table: pd.DataFrame, n: dict[str, int], null_mean: float) -> None:
    """Volcano plot: CLR difference vs -log10 p, FDR hits coloured by direction."""
    fig, ax = plt.subplots(figsize=(8.4, 6))
    sig = table["q_bh"] < config.ALPHA
    up = sig & (table["clr_mean_diff"] > 0)
    down = sig & (table["clr_mean_diff"] < 0)
    y = -np.log10(table["p"])
    ax.scatter(table.loc[~sig, "clr_mean_diff"], y[~sig], s=14, color=style.OTHER_COLOR,
               label=f"not significant (n={int((~sig).sum())})")
    ax.scatter(table.loc[up, "clr_mean_diff"], y[up], s=30, color=style.GROUP_COLORS["Periodontitis"],
               edgecolor="white", linewidth=0.5,
               label=f"higher in periodontitis, q < {config.ALPHA} (n={int(up.sum())})")
    ax.scatter(table.loc[down, "clr_mean_diff"], y[down], s=30, color=style.GROUP_COLORS["Healthy"],
               edgecolor="white", linewidth=0.5,
               label=f"higher in healthy, q < {config.ALPHA} (n={int(down.sum())})")
    if sig.any():
        threshold = table.loc[sig, "p"].max()
        ax.axhline(-np.log10(threshold), color=style.INK_SECONDARY, ls="--", lw=0.9)
        ax.text(0.99, -np.log10(threshold), f"FDR {config.ALPHA:.0%} cut-off ",
                transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=8,
                color=style.INK_SECONDARY)
    ax.axvline(0, color=style.GRID, lw=1)
    to_label = list(table.index[:8]) + [g for g in config.RED_COMPLEX_GENERA
                                        if g in table.index and g not in table.index[:8]]
    ax.set_xlabel("Difference in mean CLR abundance (periodontitis - healthy)")
    ax.set_ylabel("-log10(Mann-Whitney p)")
    ax.set_title(f"{style.title_prefix()}Differential abundance of {len(table)} genera\n"
                 f"Healthy n={n['H']}, periodontitis n={n['P']}; label-shuffled null: "
                 f"{null_mean:.2f} FDR hits on average", fontsize=11)
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    style.place_labels(ax, [(table.loc[t, "clr_mean_diff"], y[t],
                             t + (" (red complex)" if t in config.RED_COMPLEX_GENERA else ""),
                             t in config.RED_COMPLEX_GENERA) for t in to_label])
    style.save(fig, "08_volcano")


def _clr_boxes(axes, clr_values: pd.DataFrame, taxa: list[str], meta: pd.DataFrame,
               notes: dict[str, str]) -> None:
    for ax, taxon in zip(axes, taxa):
        values = {config.GROUP_LABELS[c]: clr_values.loc[taxon, meta.index[meta["group_code"] == c]]
                  .to_numpy() for c in ("H", "P")}
        flagged = meta.get("flagged_atypical", pd.Series(False, index=meta.index))
        hollow = {config.GROUP_LABELS[c]: flagged.loc[meta.index[meta["group_code"] == c]].to_numpy(bool)
                  for c in ("H", "P")}
        style.box_with_points(ax, values, hollow)
        ax.set_title(taxon.replace("_", " "), fontsize=9.5, fontstyle="italic")
        ax.set_xlabel(notes[taxon], fontsize=8, color=style.INK_SECONDARY)
    for ax in axes[len(taxa):]:
        ax.set_visible(False)


def plot_top_taxa(clr_values: pd.DataFrame, table: pd.DataFrame, meta: pd.DataFrame,
                  n_show: int = 8) -> None:
    """CLR boxplots for the top genera by p-value."""
    taxa = list(table.index[:n_show])
    notes = {t: f"q = {table.loc[t, 'q_bh']:.3f}, r = {table.loc[t, 'rank_biserial']:+.2f}"
             for t in taxa}
    fig, axes = plt.subplots(2, 4, figsize=(13, 7))
    _clr_boxes(axes.ravel(), clr_values, taxa, meta, notes)
    for ax in axes[:, 0]:
        ax.set_ylabel("CLR abundance")
    fig.suptitle(f"{style.title_prefix()}Top {n_show} genera by Mann-Whitney p "
                 f"(BH q and rank-biserial r shown; open circles = flagged samples)", fontweight="bold")
    fig.tight_layout()
    style.save(fig, "09_top_taxa_boxplots")


def plot_red_complex(species_clr: pd.DataFrame, targeted: pd.DataFrame, meta: pd.DataFrame) -> None:
    """CLR boxplots for the three pre-declared red-complex species."""
    taxa = [t for t in config.RED_COMPLEX_SPECIES if t in species_clr.index]
    notes = {t: (f"p = {targeted.loc[t, 'p']:.3f}, q = {targeted.loc[t, 'q_bh_across_targets']:.3f}"
                 f"\ndetected in {targeted.loc[t, 'prevalence_reference']:.0%} H / "
                 f"{targeted.loc[t, 'prevalence_disease']:.0%} P") for t in taxa}
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.4))
    _clr_boxes(axes, species_clr, taxa, meta, notes)
    axes[0].set_ylabel("CLR abundance (species table)")
    fig.suptitle(f"{style.title_prefix()}Pre-declared test: red-complex species "
                 f"(open circles = flagged samples)",
                 fontweight="bold")
    fig.tight_layout()
    style.save(fig, "10_red_complex_species")


def high_richness_check(pre: Preprocessed, table: pd.DataFrame, reference: str = "H",
                        disease: str = "P", min_asvs: int = 1000, n_top: int = 10) -> pd.DataFrame:
    """Where do the top hits come from? Prevalence split by a high-richness subset.

    The QC step showed a few samples with > ``min_asvs`` observed features
    (about 5x the median), mostly in the healthy group. This table shows, for
    the top genera, how often each is detected in those samples versus the
    rest, so a reader can see whether a "group difference" is carried by that
    subset.
    """
    richness = (pre.asv_rarefied > 0).sum(axis=0)
    labels = pre.metadata["group_code"]
    subsets = {
        f"{reference}_high_richness": labels.index[(labels == reference) & (richness > min_asvs)],
        f"{reference}_other": labels.index[(labels == reference) & (richness <= min_asvs)],
        f"{disease}_high_richness": labels.index[(labels == disease) & (richness > min_asvs)],
        f"{disease}_other": labels.index[(labels == disease) & (richness <= min_asvs)],
    }
    rows = {}
    for taxon in table.index[:n_top]:
        present = pre.genus_counts.loc[taxon] > 0
        rows[taxon] = {f"detected_{k} (n={len(v)})": f"{int(present[v].sum())}/{len(v)}"
                       for k, v in subsets.items()}
        rows[taxon]["direction"] = table.loc[taxon, "direction"]
    return pd.DataFrame.from_dict(rows, orient="index")


def run(pre: Preprocessed, disease: str = "P", reference: str = "H") -> dict:
    """Run M5 and return its headline numbers."""
    meta, labels = pre.metadata, pre.metadata["group_code"]
    rel_all = relative_abundance(pre.genus_counts)
    table = differential_table(pre.genus_clr, rel_all, labels, disease, reference)
    table.to_csv(config.RESULTS / "m5_differential_abundance_genus.csv")

    null = permutation_null(pre.genus_clr, labels, disease=disease, reference=reference)
    null.to_csv(config.RESULTS / "m5_null_permutations.csv", index=False)

    targeted, species_clr = targeted_species_test(pre.species_counts, labels, disease=disease,
                                     reference=reference)
    targeted.to_csv(config.RESULTS / "m5_red_complex_species.csv")
    red_genera = table.loc[[g for g in config.RED_COMPLEX_GENERA if g in table.index],
                           ["rank", "clr_mean_diff", "rank_biserial", "p", "q_bh"]]
    red_genera.to_csv(config.RESULTS / "m5_red_complex_genera.csv")

    high_richness_check(pre, table).to_csv(config.RESULTS / "m5_high_richness_check.csv")

    n = {c: int((labels == c).sum()) for c in (reference, disease)}
    n = {"H": n[reference], "P": n[disease]}
    plot_volcano(table, n, float(null["n_q_below_fdr"].mean()))
    plot_top_taxa(pre.genus_clr, table, meta)
    plot_red_complex(species_clr, targeted, meta)

    sig = table["q_bh"] < config.ALPHA
    summary = {
        "n_genera_tested": len(table),
        "n_significant_fdr05": int(sig.sum()),
        "n_higher_in_disease": int((sig & (table["clr_mean_diff"] > 0)).sum()),
        "n_higher_in_reference": int((sig & (table["clr_mean_diff"] < 0)).sum()),
        "n_raw_p_below_005": int((table["p"] < 0.05).sum()),
        "expected_raw_p_below_005_by_chance": round(0.05 * len(table), 1),
        "null_fdr_hits_mean": float(null["n_q_below_fdr"].mean()),
        "null_fdr_hits_95th_percentile": float(null["n_q_below_fdr"].quantile(0.95)),
        "null_fraction_perms_with_any_fdr_hit": float((null["n_q_below_fdr"] > 0).mean()),
        "null_raw_p_hits_mean": float(null["n_p_below_0.05"].mean()),
        "n_null_permutations": len(null),
        "top10": table.head(10)[["clr_mean_diff", "rank_biserial", "p", "q_bh"]]
        .round(5).to_dict(orient="index"),
        "red_complex_species": targeted[["clr_mean_diff", "rank_biserial", "p",
                                         "q_bh_across_targets", "prevalence_reference",
                                         "prevalence_disease", "species_rank_in_scan"]]
        .round(5).to_dict(orient="index"),
        "red_complex_genera": red_genera.round(5).to_dict(orient="index"),
    }
    (config.RESULTS / "m5_summary.json").write_text(json.dumps(summary, indent=2))
    return summary
