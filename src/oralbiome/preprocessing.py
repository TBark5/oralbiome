"""M1 - Preprocessing: aggregation, prevalence filtering, zeros, rarefaction, CLR.

Two versions of the data are produced for different jobs:
* rarefied counts (every sample subsampled to the same depth) for alpha
  diversity, where richness depends on sequencing effort;
* relative abundances / centred log-ratios (CLR) of prevalence-filtered
  genera for beta diversity, differential abundance, networks and the
  classifier.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import config
from .data import Dataset


def aggregate(counts: pd.DataFrame, taxonomy: pd.DataFrame, level: str) -> pd.DataFrame:
    """Sum ASV counts into taxa at a taxonomic rank (e.g. 'genus')."""
    return counts.groupby(taxonomy.loc[counts.index, level]).sum()


def prevalence_filter(
    counts: pd.DataFrame,
    min_prevalence: float = config.MIN_PREVALENCE,
    min_mean_rel: float = config.MIN_MEAN_REL_ABUND,
) -> pd.Series:
    """Boolean mask of taxa to keep.

    A taxon is kept if it is non-zero in at least ``min_prevalence`` of samples
    AND its mean relative abundance is at least ``min_mean_rel``.
    """
    prevalence = (counts > 0).mean(axis=1)
    mean_rel = relative_abundance(counts).mean(axis=1)
    return (prevalence >= min_prevalence) & (mean_rel >= min_mean_rel)


def relative_abundance(counts: pd.DataFrame) -> pd.DataFrame:
    """Divide each sample (column) by its total so columns sum to 1."""
    totals = counts.sum(axis=0).replace(0, np.nan)
    return (counts / totals).fillna(0.0)


def rarefy(counts: pd.DataFrame, depth: int, seed: int = config.SEED) -> pd.DataFrame:
    """Subsample every sample to exactly ``depth`` reads without replacement.

    Samples with fewer reads than ``depth`` are dropped (none are in this
    dataset because the minimum depth is used).
    """
    rng = np.random.default_rng(seed)
    kept = counts.loc[:, counts.sum(axis=0) >= depth]
    out = np.column_stack(
        [rng.multivariate_hypergeometric(kept[c].to_numpy(dtype=np.int64), depth) for c in kept]
    )
    return pd.DataFrame(out, index=kept.index, columns=kept.columns)


def clr(counts: pd.DataFrame | np.ndarray, pseudocount: float = config.PSEUDOCOUNT):
    """Centred log-ratio transform per sample (samples are columns).

    clr(x)_i = log(x_i + c) - mean_j log(x_j + c). Zeros are handled by the
    pseudocount ``c``; this is the simplest standard choice.
    """
    logged = np.log(np.asarray(counts, dtype=float) + pseudocount)
    result = logged - logged.mean(axis=0, keepdims=True)
    if isinstance(counts, pd.DataFrame):
        return pd.DataFrame(result, index=counts.index, columns=counts.columns)
    return result


def flag_atypical_samples(full: Dataset, max_unique_fraction: float = 0.5) -> pd.Series:
    """Flag samples where most reads come from ASVs found in no other sample.

    In saliva, most bacteria are shared across people. A sample whose reads
    are dominated by ASVs seen nowhere else in the study (often soil- or
    gut-associated lineages here) is more likely to be contaminated or
    low-biomass. Flagged samples stay in the primary analysis and are dropped
    in a sensitivity analysis.
    """
    unique = (full.counts > 0).sum(axis=1) == 1
    fraction = full.counts[unique].sum(axis=0) / full.counts.sum(axis=0)
    return fraction > max_unique_fraction


@dataclass
class Preprocessed:
    """All derived tables for one two-group comparison."""

    metadata: pd.DataFrame
    genus_counts: pd.DataFrame          # all genera, raw counts
    genus_filtered: pd.DataFrame        # prevalence-filtered genera, raw counts
    genus_rarefied: pd.DataFrame        # all genera, rarefied (alpha diversity)
    asv_rarefied: pd.DataFrame          # all ASVs, rarefied (sensitivity)
    genus_rel: pd.DataFrame             # filtered genera, relative abundance
    genus_clr: pd.DataFrame             # filtered genera, CLR
    phylum_rel: pd.DataFrame            # all phyla, relative abundance
    species_counts: pd.DataFrame        # all species labels, raw counts
    rarefy_depth: int = 0
    summary: dict = field(default_factory=dict)


def preprocess(ds: Dataset, tag: str = "primary", write: bool = True) -> Preprocessed:
    """Run M1 on one comparison; optionally write the processed tables + summary."""
    counts, tax = ds.counts, ds.taxonomy
    genus = aggregate(counts, tax, "genus")
    keep = prevalence_filter(genus)
    genus_f = genus.loc[keep]
    depth = config.RAREFY_DEPTH or int(counts.sum(axis=0).min())

    pre = Preprocessed(
        metadata=ds.metadata.copy(),
        genus_counts=genus,
        genus_filtered=genus_f,
        genus_rarefied=rarefy(genus, depth),
        asv_rarefied=rarefy(counts, depth),
        genus_rel=relative_abundance(genus_f),
        genus_clr=clr(genus_f),
        phylum_rel=relative_abundance(aggregate(counts, tax, "phylum")),
        species_counts=aggregate(counts, tax, "species"),
        rarefy_depth=depth,
    )
    retained = genus_f.sum(axis=0) / genus.sum(axis=0)
    pre.summary = {
        "comparison": tag,
        "n_samples": int(counts.shape[1]),
        "n_asvs_input": int(counts.shape[0]),
        "n_genera_input": int(genus.shape[0]),
        "n_genera_after_filter": int(keep.sum()),
        "min_prevalence": config.MIN_PREVALENCE,
        "min_mean_rel_abund": config.MIN_MEAN_REL_ABUND,
        "pseudocount": config.PSEUDOCOUNT,
        "rarefy_depth": depth,
        "reads_retained_after_filter_median": float(retained.median()),
        "reads_retained_after_filter_min": float(retained.min()),
        "zero_fraction_filtered_genus_table": float((genus_f.to_numpy() == 0).mean()),
    }
    if not write:
        return pre
    out = config.RESULTS / f"m1_preprocessing_summary_{tag}.json"
    out.write_text(json.dumps(pre.summary, indent=2))
    if tag == "primary":
        pre.genus_rel.to_csv(config.DATA_PROCESSED / "genus_relative_abundance.csv")
        pre.genus_clr.to_csv(config.DATA_PROCESSED / "genus_clr.csv")
        pre.genus_rarefied.to_csv(config.DATA_PROCESSED / "genus_rarefied_counts.csv")
        pre.phylum_rel.to_csv(config.DATA_PROCESSED / "phylum_relative_abundance.csv")
        pre.metadata.to_csv(config.DATA_PROCESSED / "primary_metadata.csv")
    return pre
