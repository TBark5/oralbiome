"""M5 - Differential abundance: per-taxon tests, BH-FDR, null control, pathogens.

Each taxon is compared between groups with a two-sided Mann-Whitney U test
on CLR-transformed abundances. The CLR puts every sample on a log-ratio scale
relative to its own geometric mean, which is the standard way to reduce the
compositional artefacts of comparing raw proportions.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from . import config
from .preprocessing import clr, prevalence_filter, relative_abundance
from .stats import benjamini_hochberg, mean_difference_ci, rank_biserial


def differential_table(clr_values: pd.DataFrame, rel: pd.DataFrame, labels: pd.Series,
                       disease: str = "P", reference: str = "H") -> pd.DataFrame:
    """Per-taxon test table, sorted by p-value.

    Columns: CLR mean difference (disease - reference) with 95% Welch CI,
    rank-biserial r, Mann-Whitney p, BH q, prevalence and mean relative
    abundance in each group.
    """
    lab = labels.loc[clr_values.columns]
    x = clr_values.loc[:, lab == disease].to_numpy()
    y = clr_values.loc[:, lab == reference].to_numpy()
    _, p = stats.mannwhitneyu(x, y, axis=1, alternative="two-sided")
    rows = []
    for i, taxon in enumerate(clr_values.index):
        diff, lo, hi = mean_difference_ci(x[i], y[i])
        rows.append({"taxon": taxon, "clr_mean_diff": diff, "clr_diff_ci_low": lo,
                     "clr_diff_ci_high": hi, "rank_biserial": rank_biserial(x[i], y[i])})
    table = pd.DataFrame(rows).set_index("taxon")
    table["p"] = p
    table["q_bh"] = benjamini_hochberg(p)
    for code, name in ((reference, "reference"), (disease, "disease")):
        ids = lab.index[lab == code]
        table[f"prevalence_{name}"] = (rel.loc[clr_values.index, ids] > 0).mean(axis=1)
        table[f"mean_rel_abund_{name}"] = rel.loc[clr_values.index, ids].mean(axis=1)
    table["direction"] = np.where(table["clr_mean_diff"] > 0, "higher in disease",
                                  "higher in reference")
    table = table.sort_values("p")
    table.insert(0, "rank", np.arange(1, len(table) + 1))
    return table


def permutation_null(clr_values: pd.DataFrame, labels: pd.Series, n_perm: int = 200,
                     disease: str = "P", reference: str = "H", fdr: float = config.ALPHA,
                     seed: int = config.SEED) -> pd.DataFrame:
    """How many taxa pass FDR (and raw p < 0.05) when group labels are shuffled?

    With no true signal, BH should give ~0 discoveries; this is the negative
    control for M5.
    """
    rng = np.random.default_rng(seed)
    lab = labels.loc[clr_values.columns].to_numpy()
    values = clr_values.to_numpy()
    rows = []
    for i in range(n_perm):
        perm = rng.permutation(lab)
        _, p = stats.mannwhitneyu(values[:, perm == disease], values[:, perm == reference],
                                  axis=1, alternative="two-sided")
        rows.append({"permutation": i, "n_q_below_fdr": int((benjamini_hochberg(p) < fdr).sum()),
                     "n_p_below_0.05": int((p < 0.05).sum())})
    return pd.DataFrame(rows)


def targeted_species_test(species_counts: pd.DataFrame, labels: pd.Series,
                          targets: tuple[str, ...] = config.RED_COMPLEX_SPECIES,
                          disease: str = "P", reference: str = "H") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pre-declared test of the red-complex species, BH-corrected across targets.

    CLR is computed on the species table after the standard prevalence filter,
    with the target species always kept so that they can be tested even if
    they are rare. Returns the result table and the species CLR matrix used.
    """
    keep = prevalence_filter(species_counts)
    for t in targets:
        if t in keep.index:
            keep.loc[t] = True
    species_clr = clr(species_counts.loc[keep])
    rel = relative_abundance(species_counts)
    table = differential_table(species_clr, rel, labels, disease, reference)
    present = [t for t in targets if t in table.index]
    out = table.loc[present].drop(columns=["rank", "q_bh"])
    out["q_bh_across_targets"] = benjamini_hochberg(out["p"].to_numpy())
    out["species_rank_in_scan"] = [int(table.loc[t, "rank"]) for t in present]
    out["n_species_tested"] = len(table)
    missing = [t for t in targets if t not in table.index]
    for t in missing:
        out.loc[t, "note"] = "not detected in any sample"
    return out, species_clr
