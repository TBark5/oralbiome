"""Tests for the shared statistics helpers."""
import numpy as np
from scipy import stats as sps

from oralbiome.stats import (benjamini_hochberg, bootstrap_ci, compare_two_groups,
                             mean_difference_ci, rank_biserial)


def test_bh_matches_scipy():
    rng = np.random.default_rng(0)
    p = rng.uniform(size=200) ** 2
    assert np.allclose(benjamini_hochberg(p), sps.false_discovery_control(p, method="bh"))


def test_bh_known_values_and_nan():
    q = benjamini_hochberg(np.array([0.01, 0.04, 0.03, np.nan]))
    assert np.isnan(q[3])
    assert np.allclose(q[:3], [0.03, 0.04, 0.04])


def test_rank_biserial_extremes_and_symmetry():
    assert rank_biserial([5, 6, 7], [1, 2, 3]) == 1.0
    assert rank_biserial([1, 2, 3], [5, 6, 7]) == -1.0
    assert rank_biserial([1, 2], [1, 2]) == 0.0
    x, y = np.array([1.0, 4, 6, 8]), np.array([2.0, 3, 5])
    assert np.isclose(rank_biserial(x, y), -rank_biserial(y, x))


def test_rank_biserial_matches_mann_whitney_u():
    rng = np.random.default_rng(1)
    x, y = rng.normal(size=12), rng.normal(0.5, size=15)
    u = sps.mannwhitneyu(x, y).statistic
    assert np.isclose(rank_biserial(x, y), 2 * u / (len(x) * len(y)) - 1)


def test_bootstrap_ci_contains_estimate():
    rng = np.random.default_rng(2)
    x, y = rng.normal(1, size=20), rng.normal(size=20)
    lo, hi = bootstrap_ci(x, y, rank_biserial, n_boot=500)
    assert lo <= rank_biserial(x, y) <= hi


def test_compare_two_groups_direction():
    res = compare_two_groups(np.arange(10, 20), np.arange(0, 10), with_ci=False)
    assert res["rank_biserial"] == 1.0 and res["p"] < 0.001


def test_mean_difference_ci():
    diff, lo, hi = mean_difference_ci(np.array([3.0, 4, 5]), np.array([1.0, 2, 3]))
    assert diff == 2.0 and lo < 2.0 < hi
