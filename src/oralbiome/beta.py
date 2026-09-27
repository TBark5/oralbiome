"""M4 - Beta diversity: distances, PCoA, NMDS, PERMANOVA and PERMDISP.

PERMANOVA and PERMDISP are implemented directly (Anderson 2001, 2006) so the
code matches the formulas in METHODS.md line by line.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS

from . import config


def bray_curtis(rel: pd.DataFrame) -> pd.DataFrame:
    """Bray-Curtis dissimilarity between samples (columns of ``rel``)."""
    d = squareform(pdist(rel.T.to_numpy(), metric="braycurtis"))
    return pd.DataFrame(d, index=rel.columns, columns=rel.columns)


def jaccard(counts: pd.DataFrame) -> pd.DataFrame:
    """Jaccard distance on presence/absence (1 - shared / union)."""
    d = squareform(pdist((counts.T.to_numpy() > 0), metric="jaccard"))
    return pd.DataFrame(d, index=counts.columns, columns=counts.columns)


def pcoa(dist: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Principal coordinates analysis (classical MDS).

    Returns coordinates on axes with positive eigenvalues, the fraction of
    positive-eigenvalue variance each axis explains, and all eigenvalues
    (negative ones occur for non-Euclidean distances such as Bray-Curtis).
    """
    d = dist.to_numpy()
    n = d.shape[0]
    a = -0.5 * d**2
    centre = np.eye(n) - np.ones((n, n)) / n
    b = centre @ a @ centre
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos = eigvals > 1e-10
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    explained = eigvals[pos] / eigvals[pos].sum()
    cols = [f"PCo{i + 1}" for i in range(pos.sum())]
    return pd.DataFrame(coords, index=dist.index, columns=cols), explained, eigvals


def nmds(dist: pd.DataFrame, seed: int = config.SEED, n_init: int = 8) -> tuple[pd.DataFrame, float]:
    """Non-metric MDS in 2 dimensions; returns coordinates and Kruskal stress-1."""
    model = MDS(n_components=2, metric_mds=False, metric="precomputed", n_init=n_init,
                init="random", max_iter=1000, random_state=seed, normalized_stress=True)
    coords = model.fit_transform(dist.to_numpy())
    return pd.DataFrame(coords, index=dist.index, columns=["NMDS1", "NMDS2"]), float(model.stress_)


def _pseudo_f(d2: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    """PERMANOVA pseudo-F and R^2 from squared distances and group labels."""
    n = len(labels)
    groups = np.unique(labels)
    ss_total = d2[np.triu_indices(n, 1)].sum() / n
    ss_within = 0.0
    for g in groups:
        idx = np.where(labels == g)[0]
        sub = d2[np.ix_(idx, idx)]
        ss_within += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
    ss_between = ss_total - ss_within
    a = len(groups)
    f = (ss_between / (a - 1)) / (ss_within / (n - a))
    return float(f), float(ss_between / ss_total)


def permanova(dist: pd.DataFrame, labels: pd.Series, n_perm: int = config.N_PERMUTATIONS,
              seed: int = config.SEED) -> dict:
    """One-way PERMANOVA: are group centroids different in distance space?

    The p-value is (1 + #{F_perm >= F_obs}) / (1 + n_perm), labels shuffled.
    """
    rng = np.random.default_rng(seed)
    d2 = dist.to_numpy() ** 2
    y = labels.loc[dist.index].to_numpy()
    f_obs, r2 = _pseudo_f(d2, y)
    f_perm = np.array([_pseudo_f(d2, rng.permutation(y))[0] for _ in range(n_perm)])
    p = (1 + np.sum(f_perm >= f_obs)) / (1 + n_perm)
    return {"pseudo_F": f_obs, "R2": r2, "p": float(p), "n_perm": n_perm, "n": len(y)}


def distances_to_centroid(dist: pd.DataFrame, labels: pd.Series) -> pd.Series:
    """Distance of each sample to its group centroid in full PCoA space.

    Negative eigenvalues are handled as in Anderson (2006): their axes
    contribute negatively to the squared distance.
    """
    d = dist.to_numpy()
    n = d.shape[0]
    centre = np.eye(n) - np.ones((n, n)) / n
    b = centre @ (-0.5 * d**2) @ centre
    eigvals, eigvecs = np.linalg.eigh(b)
    keep = np.abs(eigvals) > 1e-10
    vecs = eigvecs[:, keep] * np.sqrt(np.abs(eigvals[keep]))
    sign = np.sign(eigvals[keep])
    y = labels.loc[dist.index].to_numpy()
    out = np.empty(n)
    for g in np.unique(y):
        idx = y == g
        diff = vecs[idx] - vecs[idx].mean(axis=0)
        out[idx] = np.sqrt(np.maximum((diff**2 * sign).sum(axis=1), 0))
    return pd.Series(out, index=dist.index)


def _anova_f(values: np.ndarray, labels: np.ndarray) -> float:
    groups = [values[labels == g] for g in np.unique(labels)]
    grand = values.mean()
    ss_b = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
    ss_w = sum(((g - g.mean()) ** 2).sum() for g in groups)
    return float((ss_b / (len(groups) - 1)) / (ss_w / (len(values) - len(groups))))


def permdisp(dist: pd.DataFrame, labels: pd.Series, n_perm: int = config.N_PERMUTATIONS,
             seed: int = config.SEED) -> dict:
    """PERMDISP: do groups differ in spread (distance to own centroid)?

    ANOVA F on distances-to-centroid; p-value from permuting group labels of
    those distances.
    """
    rng = np.random.default_rng(seed)
    z = distances_to_centroid(dist, labels)
    y = labels.loc[z.index].to_numpy()
    f_obs = _anova_f(z.to_numpy(), y)
    f_perm = np.array([_anova_f(z.to_numpy(), rng.permutation(y)) for _ in range(n_perm)])
    p = (1 + np.sum(f_perm >= f_obs)) / (1 + n_perm)
    means = {g: float(z[y == g].mean()) for g in np.unique(y)}
    return {"F": f_obs, "p": float(p), "n_perm": n_perm, "mean_distance_to_centroid": means}
