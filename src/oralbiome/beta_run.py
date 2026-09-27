"""M4 runner: compute beta-diversity statistics and draw the ordination figures."""
from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import chi2

from . import beta, config, style
from .preprocessing import Preprocessed


def confidence_ellipse(ax: plt.Axes, x: np.ndarray, y: np.ndarray, color: str,
                       level: float = 0.95) -> None:
    """Draw the ellipse expected to contain ``level`` of a bivariate normal fit."""
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    scale = np.sqrt(chi2.ppf(level, 2))
    width, height = 2 * scale * np.sqrt(vals)
    ax.add_patch(Ellipse((x.mean(), y.mean()), width, height, angle=angle, facecolor=color,
                         alpha=0.10, edgecolor=color, linewidth=1.5))


def _scatter_groups(ax: plt.Axes, coords: pd.DataFrame, meta: pd.DataFrame,
                    xcol: str, ycol: str) -> None:
    flagged = meta.get("flagged_atypical", pd.Series(False, index=meta.index))
    for label in [g for g in config.GROUP_LABELS.values() if g in set(meta["group"])]:
        ids = meta.index[meta["group"] == label]
        color = style.GROUP_COLORS[label]
        c = coords.loc[ids]
        f = flagged.loc[ids].to_numpy(bool)
        confidence_ellipse(ax, c[xcol].to_numpy(), c[ycol].to_numpy(), color)
        ax.scatter(c.loc[~f, xcol], c.loc[~f, ycol], s=38, color=color, edgecolor="white",
                   linewidth=0.6, label=style.group_label_with_n(label, len(ids)), zorder=3)
        ax.scatter(c.loc[f, xcol], c.loc[f, ycol], s=40, facecolor="white", edgecolor=color,
                   linewidth=1.4, zorder=3)
    ax.axhline(0, color=style.GRID, lw=0.8)
    ax.axvline(0, color=style.GRID, lw=0.8)


def _stat_text(res: dict) -> str:
    pm, pd_ = res["permanova"], res["permdisp"]
    return (f"PERMANOVA R² = {pm['R2']:.3f}, F = {pm['pseudo_F']:.2f}, p = {pm['p']:.3f}\n"
            f"PERMDISP F = {pd_['F']:.2f}, p = {pd_['p']:.3f}")


def plot_ordinations(ords: dict, results: dict, meta: pd.DataFrame, kind: str) -> None:
    """Side-by-side Bray-Curtis / Jaccard ordination (PCoA or NMDS)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, metric in zip(axes, ("bray_curtis", "jaccard")):
        o = ords[metric]
        if kind == "pcoa":
            coords, explained = o
            xcol, ycol = "PCo1", "PCo2"
            ax.set_xlabel(f"PCo1 ({explained[0]:.1%} of variance)")
            ax.set_ylabel(f"PCo2 ({explained[1]:.1%} of variance)")
            sub = ""
        else:
            coords, stress = o
            xcol, ycol = "NMDS1", "NMDS2"
            ax.set_xlabel("NMDS1")
            ax.set_ylabel("NMDS2")
            sub = f" (stress = {stress:.3f})"
        _scatter_groups(ax, coords, meta, xcol, ycol)
        name = "Bray-Curtis" if metric == "bray_curtis" else "Jaccard"
        ax.set_title(f"{name}{sub}")
        ax.text(0.02, 0.02, _stat_text(results[metric]), transform=ax.transAxes, fontsize=8.5,
                color=style.INK_SECONDARY, va="bottom",
                bbox={"facecolor": "white", "edgecolor": style.GRID, "alpha": 0.9})
        ax.legend(loc="upper right")
    title = "PCoA" if kind == "pcoa" else "Non-metric MDS"
    fig.suptitle(f"{style.title_prefix()}{title} of genus-level community composition "
                 f"(95% ellipses; open circles = flagged samples)", fontweight="bold")
    fig.tight_layout()
    style.save(fig, "05_pcoa" if kind == "pcoa" else "06_nmds")


def plot_heatmap(dist: pd.DataFrame, meta: pd.DataFrame) -> None:
    """Bray-Curtis matrix, samples grouped by status and clustered within group."""
    order = []
    groups = [g for g in config.GROUP_LABELS.values() if g in set(meta["group"])]
    for g in groups:
        ids = list(meta.index[meta["group"] == g])
        sub = dist.loc[ids, ids].to_numpy()
        leaves = leaves_list(linkage(squareform(sub, checks=False), "average"))
        order += [ids[i] for i in leaves]
    d = dist.loc[order, order]
    fig, ax = plt.subplots(figsize=(8.2, 7.4))
    im = ax.imshow(d.to_numpy(), cmap=style.SEQUENTIAL_CMAP.reversed(), vmin=0, vmax=1)
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_yticklabels(order, fontsize=6)
    ax.grid(False)
    pos = 0
    for g in groups:
        n = int((meta["group"] == g).sum())
        color = style.GROUP_COLORS[g]
        ax.add_patch(plt.Rectangle((pos - 0.5, -2.2), n, 1.2, color=color, clip_on=False))
        ax.text(pos + n / 2 - 0.5, -2.6, style.group_label_with_n(g, n), ha="center",
                va="bottom", fontsize=9, color=color, fontweight="bold")
        pos += n
        if pos < len(order):
            ax.axhline(pos - 0.5, color="white", lw=2)
            ax.axvline(pos - 0.5, color="white", lw=2)
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Bray-Curtis dissimilarity (0 = identical, 1 = no shared taxa)")
    ax.set_title(f"{style.title_prefix()}Pairwise Bray-Curtis dissimilarity (genus level)",
                 pad=40)
    style.save(fig, "07_distance_heatmap")


def run(pre: Preprocessed) -> dict:
    """Compute distances, ordinations and tests; write tables and figures."""
    meta = pre.metadata
    labels = meta["group_code"]
    dists = {"bray_curtis": beta.bray_curtis(pre.genus_rel),
             "jaccard": beta.jaccard(pre.genus_filtered)}
    results, pcoas, nmdss = {}, {}, {}
    flagged = meta.get("flagged_atypical", pd.Series(False, index=meta.index)).astype(bool)
    keep = meta.index[~flagged]
    for metric, dist in dists.items():
        dist.to_csv(config.RESULTS / f"m4_distance_{metric}.csv")
        coords, explained, _ = beta.pcoa(dist)
        pcoas[metric] = (coords, explained)
        nmdss[metric] = beta.nmds(dist)
        coords.iloc[:, :3].join(nmdss[metric][0]).to_csv(
            config.RESULTS / f"m4_ordination_{metric}.csv")
        results[metric] = {
            "permanova": beta.permanova(dist, labels),
            "permdisp": beta.permdisp(dist, labels),
            "permanova_excluding_flagged": beta.permanova(dist.loc[keep, keep], labels),
            "pcoa_explained_axis1": float(explained[0]),
            "pcoa_explained_axis2": float(explained[1]),
            "nmds_stress": nmdss[metric][1],
        }
    (config.RESULTS / "m4_beta_tests.json").write_text(json.dumps(results, indent=2))
    plot_ordinations(pcoas, results, meta, "pcoa")
    plot_ordinations(nmdss, results, meta, "nmds")
    plot_heatmap(dists["bray_curtis"], meta)
    return results
