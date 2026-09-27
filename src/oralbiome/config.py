"""Project-wide paths, random seed and analysis thresholds.

Every threshold used anywhere in the pipeline lives here so it can be cited
from DECISIONS.md and changed in one place.
"""
from __future__ import annotations

from pathlib import Path

ROOT: Path = Path(__file__).resolve().parents[2]
DATA_RAW: Path = ROOT / "data" / "raw"
DATA_PROCESSED: Path = ROOT / "data" / "processed"
RESULTS: Path = ROOT / "results"
FIGURES: Path = ROOT / "figures"

SEED: int = 42

# Source data (Guo et al. 2026, figshare, CC BY 4.0).
DATA_URL: str = "https://ndownloader.figshare.com/files/57145445"
DATA_DOI: str = "10.6084/m9.figshare.29897750.v1"
RAW_ZIP: Path = DATA_RAW / "guo2025_figshare_29897750.zip"
SYNTHETIC_TAG: str = "SYNTHETIC"

# Group codes in the source metadata and human-readable labels.
GROUP_LABELS: dict[str, str] = {
    "H": "Healthy",
    "P": "Periodontitis",
    "T": "Hypertension",
    "TP": "Hypertension + periodontitis",
}
PRIMARY_GROUPS: tuple[str, str] = ("H", "P")
REPLICATION_GROUPS: tuple[str, str] = ("T", "TP")

# M1 preprocessing thresholds (justified in DECISIONS.md).
MIN_PREVALENCE: float = 0.10       # taxon must be non-zero in >= 10% of samples
MIN_MEAN_REL_ABUND: float = 1e-4   # and have mean relative abundance >= 0.01%
PSEUDOCOUNT: float = 0.5           # added to counts before log-ratio transforms
RAREFY_DEPTH: int | None = None    # None -> use the minimum sample depth

# M2 core microbiome: taxon present in >= 90% of samples of a group.
CORE_PREVALENCE: float = 0.90
TOP_N_TAXA_PLOT: int = 12

# Statistics.
ALPHA: float = 0.05
N_PERMUTATIONS: int = 999
N_BOOTSTRAP: int = 2000

# M6 network thresholds.
NETWORK_MIN_ABS_RHO: float = 0.6
NETWORK_MAX_Q: float = 0.05

# M7 classifier.
CV_FOLDS: int = 5
CV_REPEATS: int = 10
N_LABEL_PERMUTATIONS: int = 100
RF_TREES: int = 300

RED_COMPLEX_SPECIES: tuple[str, ...] = (
    "Porphyromonas_gingivalis",
    "Tannerella_forsythia",
    "Treponema_denticola",
)
RED_COMPLEX_GENERA: tuple[str, ...] = ("Porphyromonas", "Tannerella", "Treponema")


def ensure_dirs() -> None:
    """Create the output directories if they do not exist."""
    for path in (DATA_RAW, DATA_PROCESSED, RESULTS, FIGURES):
        path.mkdir(parents=True, exist_ok=True)
