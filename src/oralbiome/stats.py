"""Shared statistics: Mann-Whitney with effect sizes, bootstrap CIs, BH-FDR."""
from __future__ import annotations

import numpy as np
from scipy import stats

from . import config


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values).

    Sort p ascending, compute p_(i) * m / i, then take the running minimum
    from the largest rank down so q-values are monotone, and cap at 1.
    NaN p-values are passed through as NaN and do not count towards m.
    """
    p = np.asarray(pvalues, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    m = ok.sum()
    if m == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def rank_biserial(x: np.ndarray, y: np.ndarray) -> float:
    """Rank-biserial correlation r = 2*U_x/(n_x*n_y) - 1, in [-1, 1].

    r = P(X > Y) - P(X < Y) for a random pair (ties count half). Positive
    means values in ``x`` tend to be larger. Identical to Cliff's delta.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    diff = x[:, None] - y[None, :]
    return float((np.sum(diff > 0) - np.sum(diff < 0)) / (len(x) * len(y)))


def bootstrap_ci(x: np.ndarray, y: np.ndarray, func, n_boot: int = config.N_BOOTSTRAP,
                 seed: int = config.SEED, level: float = 0.95) -> tuple[float, float]:
    """Percentile bootstrap CI of func(x, y), resampling within each group."""
    rng = np.random.default_rng(seed)
    x, y = np.asarray(x, float), np.asarray(y, float)
    values = np.array([
        func(rng.choice(x, len(x), replace=True), rng.choice(y, len(y), replace=True))
        for _ in range(n_boot)
    ])
    tail = (1 - level) / 2 * 100
    return float(np.percentile(values, tail)), float(np.percentile(values, 100 - tail))


def compare_two_groups(x: np.ndarray, y: np.ndarray, with_ci: bool = True) -> dict:
    """Two-sided Mann-Whitney U test of x vs y with effect size and CI.

    ``x`` is the disease group and ``y`` the reference group, so a positive
    rank-biserial r means higher values in disease.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    result = {
        "n_disease": len(x), "n_reference": len(y),
        "median_disease": float(np.median(x)), "median_reference": float(np.median(y)),
        "U": float(u), "p": float(p), "rank_biserial": rank_biserial(x, y),
    }
    if with_ci:
        lo, hi = bootstrap_ci(x, y, rank_biserial)
        result.update({"rank_biserial_ci_low": lo, "rank_biserial_ci_high": hi})
    return result


def mean_difference_ci(x: np.ndarray, y: np.ndarray, level: float = 0.95) -> tuple[float, float, float]:
    """Difference of means (x - y) with a Welch t confidence interval."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    diff = x.mean() - y.mean()
    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y)
    se = np.sqrt(vx + vy)
    if se == 0:
        return float(diff), float(diff), float(diff)
    dof = (vx + vy) ** 2 / (vx**2 / (len(x) - 1) + vy**2 / (len(y) - 1))
    t = stats.t.ppf(0.5 + level / 2, dof)
    return float(diff), float(diff - t * se), float(diff + t * se)
