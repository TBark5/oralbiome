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
    "Hypertension": "#0072B2",
    "Hypertension + periodontitis": "#D55E00",
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
