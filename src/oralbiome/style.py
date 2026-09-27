"""Shared matplotlib style so every figure looks like part of one set.

Palettes were checked with a colour-vision-deficiency validator (protan,
deutan and tritan simulations): the two group colours are the Okabe-Ito blue
and vermillion, and the eight taxon colours are a validated categorical order.
Anything beyond eight taxa is folded into a grey "Other" band rather than
given a generated ninth hue.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from . import config  # noqa: E402

GROUP_COLORS: dict[str, str] = {
    "Healthy": "#0072B2",
    "Periodontitis": "#D55E00",
    "Hypertension": "#56B4E9",
    "Hypertension + periodontitis": "#E69F00",
}
TAXON_COLORS: list[str] = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#4a3aa7", "#e34948",
]
OTHER_COLOR: str = "#c3c2b7"
NEUTRAL_COLOR: str = "#898781"
INK: str = "#0b0b0b"
INK_SECONDARY: str = "#52514e"
GRID: str = "#e1e0d9"

# Sequential (one hue, light -> dark) and diverging (two hues, grey midpoint).
SEQUENTIAL_CMAP = LinearSegmentedColormap.from_list(
    "oral_seq", ["#f4f8fd", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"]
)
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "oral_div", ["#1c5cab", "#86b6ef", "#f0efec", "#f2a27f", "#b3441a"]
)

DPI: int = 300


def apply_style() -> None:
    """Set global rcParams used by every figure in the project."""
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": DPI,
            "savefig.bbox": "tight",
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "axes.edgecolor": "#c3c2b7",
            "axes.labelcolor": INK,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.axisbelow": True,
            "xtick.color": INK_SECONDARY,
            "ytick.color": INK_SECONDARY,
            "legend.frameon": False,
            "legend.fontsize": 9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def group_label_with_n(label: str, n: int) -> str:
    """Return a legend/axis label such as 'Healthy (n=16)'."""
    return f"{label} (n={n})"


def title_prefix() -> str:
    """Return 'SYNTHETIC - ' when the pipeline is running on synthetic data."""
    marker = config.DATA_PROCESSED / "DATA_IS_SYNTHETIC"
    return f"{config.SYNTHETIC_TAG} - " if marker.exists() else ""


def save(fig: plt.Figure, name: str, out_dir: Path | None = None) -> Path:
    """Save a figure as a 300 dpi PNG and close it."""
    out_dir = out_dir or config.FIGURES
    out_dir.mkdir(parents=True, exist_ok=True)
    if title_prefix():
        name = f"{config.SYNTHETIC_TAG}_{name}"
    path = out_dir / f"{name}.png"
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def box_with_points(ax: plt.Axes, values: dict[str, "np.ndarray"],
                    hollow: dict[str, "np.ndarray"] | None = None) -> None:
    """Boxplot per group with every sample drawn as a jittered point.

    ``values`` maps group label -> values. ``hollow`` optionally maps group
    label -> boolean mask of points to draw as open circles (flagged samples).
    Tick labels include the group size.
    """
    import numpy as np

    rng = np.random.default_rng(config.SEED)
    labels = list(values)
    for i, label in enumerate(labels):
        v = np.asarray(values[label], dtype=float)
        color = GROUP_COLORS.get(label, NEUTRAL_COLOR)
        ax.boxplot(v, positions=[i], widths=0.55, showfliers=False, patch_artist=True,
                   boxprops={"facecolor": color + "22", "edgecolor": color},
                   medianprops={"color": color, "linewidth": 2},
                   whiskerprops={"color": color}, capprops={"color": color})
        jitter = i + rng.uniform(-0.13, 0.13, len(v))
        mask = np.zeros(len(v), bool) if hollow is None else np.asarray(hollow[label], bool)
        ax.scatter(jitter[~mask], v[~mask], s=24, color=color, edgecolor="white",
                   linewidth=0.5, zorder=3)
        ax.scatter(jitter[mask], v[mask], s=26, facecolor="white", edgecolor=color,
                   linewidth=1.2, zorder=3)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([f"{lab}\n(n={len(values[lab])})" for lab in labels])
    ax.grid(axis="x", visible=False)


def place_labels(ax: plt.Axes, points: list[tuple[float, float, str, bool]],
                 boxed: bool = False, spread: float = 1.0) -> None:
    """Annotate points while avoiding overlapping text boxes.

    ``points`` holds (x, y, text, bold). For each label, candidate offsets are
    tried in order and the first whose bounding box does not overlap an
    already placed label is kept. ``boxed`` draws a white background behind
    each label (for busy plots) and ``spread`` scales the offsets.
    """
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    placed = []
    candidates = [(6, 6), (6, -8), (-6, 6), (-6, -8), (6, 16), (-6, 16), (6, -18), (-6, -18),
                  (6, 26), (-6, 26), (6, -28), (-6, -28), (20, 40), (-20, 40), (20, -40),
                  (-20, -40)]
    candidates = [(dx * spread, dy * spread) for dx, dy in candidates]
    bbox = {"facecolor": "white", "edgecolor": GRID, "alpha": 0.92, "pad": 1.5,
            "boxstyle": "round,pad=0.25"} if boxed else None
    for x, y, _, _ in points:
        ax.scatter([x], [y], s=34, facecolor="none", edgecolor=INK, linewidth=0.9, zorder=4)
    for x, y, text, bold in points:
        for dx, dy in candidates:
            ann = ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points",
                              fontsize=8, ha="left" if dx > 0 else "right", va="center",
                              fontweight="bold" if bold else "normal", color=INK, bbox=bbox,
                              arrowprops={"arrowstyle": "-", "color": INK_SECONDARY, "lw": 0.6,
                                          "shrinkA": 0, "shrinkB": 3})
            box = ann.get_window_extent(renderer).expanded(1.05, 1.15)
            inside = ax.get_window_extent(renderer)
            fits = (box.x0 >= inside.x0 and box.x1 <= inside.x1
                    and box.y0 >= inside.y0 and box.y1 <= inside.y1)
            if fits and not any(box.overlaps(b) for b in placed):
                placed.append(box)
                break
            ann.remove()
        else:
            ann = ax.annotate(text, (x, y), xytext=(6, 6), textcoords="offset points",
                              fontsize=8, color=INK, fontweight="bold" if bold else "normal")
            placed.append(ann.get_window_extent(renderer))
