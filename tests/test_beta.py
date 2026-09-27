"""Tests for M4 distances, PCoA, PERMANOVA and PERMDISP."""
import numpy as np
import pandas as pd
from scipy import stats as sps
from scipy.spatial.distance import pdist, squareform

from oralbiome import beta


def _dist(points: np.ndarray) -> pd.DataFrame:
    names = [f"s{i}" for i in range(len(points))]
    return pd.DataFrame(squareform(pdist(points)), index=names, columns=names)


def test_bray_curtis_hand_example():
    rel = pd.DataFrame({"a": [0.5, 0.5, 0.0], "b": [0.5, 0.0, 0.5]})
    # BC = 1 - 2 * (sum of shared minimums) / (total a + total b) = 1 - 2 * 0.5 / 2
    assert np.isclose(beta.bray_curtis(rel).loc["a", "b"], 0.5)


def test_jaccard_hand_example():
    counts = pd.DataFrame({"a": [3, 1, 0], "b": [2, 0, 7]})
    # 1 shared taxon out of 3 present in either sample -> distance 2/3
    assert np.isclose(beta.jaccard(counts).loc["a", "b"], 2 / 3)


def test_pcoa_recovers_euclidean_distances():
    rng = np.random.default_rng(4)
    dist = _dist(rng.normal(size=(12, 3)))
    coords, explained, _ = beta.pcoa(dist)
    assert np.allclose(squareform(pdist(coords.to_numpy())), dist.to_numpy(), atol=1e-8)
    assert np.isclose(explained.sum(), 1.0)


def test_permanova_equals_anova_for_one_dimension():
    rng = np.random.default_rng(5)
    values = np.r_[rng.normal(0, 1, 10), rng.normal(1, 1, 12)]
    labels = pd.Series(["A"] * 10 + ["B"] * 12, index=[f"s{i}" for i in range(22)])
    res = beta.permanova(_dist(values[:, None]), labels, n_perm=99)
    f_anova = sps.f_oneway(values[:10], values[10:]).statistic
    assert np.isclose(res["pseudo_F"], f_anova)


def test_permanova_detects_shift_and_not_noise():
    rng = np.random.default_rng(6)
    labels = pd.Series(["A"] * 15 + ["B"] * 15, index=[f"s{i}" for i in range(30)])
    shifted = np.r_[rng.normal(0, 1, (15, 2)), rng.normal(3, 1, (15, 2))]
    noise = rng.normal(0, 1, (30, 2))
    assert beta.permanova(_dist(shifted), labels, n_perm=199)["p"] < 0.01
    assert beta.permanova(_dist(noise), labels, n_perm=199)["p"] > 0.05


def test_permdisp_detects_spread_difference():
    rng = np.random.default_rng(7)
    labels = pd.Series(["A"] * 20 + ["B"] * 20, index=[f"s{i}" for i in range(40)])
    pts = np.r_[rng.normal(0, 0.3, (20, 2)), rng.normal(0, 3, (20, 2))]
    res = beta.permdisp(_dist(pts), labels, n_perm=199)
    assert res["p"] < 0.01
    assert res["mean_distance_to_centroid"]["B"] > res["mean_distance_to_centroid"]["A"]
