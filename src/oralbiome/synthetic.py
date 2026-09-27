"""Dirichlet-multinomial synthetic abundance tables with planted differences.

Used in two ways:
1. As the mandatory fallback when the real dataset cannot be downloaded.
2. As a positive control: the recovery test in ``tests/`` and the validation
   step in ``run_all.py`` check that the pipeline finds the taxa that were
   deliberately made more or less abundant in the "Periodontitis" group.

Model: each group g has a mean composition pi_g. Every sample draws its own
composition p ~ Dirichlet(theta * pi_g) (theta controls between-sample
overdispersion, smaller = noisier), then counts ~ Multinomial(depth, p).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from .data import RANKS, Dataset

_REAL_GENERA = [
    "Streptococcus", "Prevotella", "Veillonella", "Neisseria", "Haemophilus",
    "Rothia", "Porphyromonas", "Fusobacterium", "Leptotrichia", "Actinomyces",
    "Gemella", "Granulicatella", "Alloprevotella", "Treponema", "Tannerella",
    "Campylobacter", "Capnocytophaga", "Filifactor", "Parvimonas", "Selenomonas",
]
_PHYLA = ["Firmicutes", "Bacteroidota", "Proteobacteria", "Actinobacteriota",
          "Fusobacteriota", "Spirochaetota", "Campilobacterota", "Patescibacteria"]


def generate_dm(
    n_healthy: int = 16,
    n_disease: int = 18,
    n_taxa: int = 150,
    n_up: int = 10,
    n_down: int = 5,
    fold_change: float = 4.0,
    theta: float = 300.0,
    depth_range: tuple[int, int] = (30_000, 70_000),
    seed: int = config.SEED,
) -> tuple[Dataset, pd.DataFrame]:
    """Simulate a two-group count table and return it with the ground truth.

    The three red-complex genera are always among the planted "up" taxa so the
    pathogen cross-check in M5 has something real to find.
    """
    rng = np.random.default_rng(seed)
    names = _REAL_GENERA + [f"SynGenus_{i:03d}" for i in range(n_taxa - len(_REAL_GENERA))]
    names = names[:n_taxa]

    # Long-tailed baseline, like real communities: a few dominant taxa.
    base = np.sort(rng.lognormal(mean=0.0, sigma=2.0, size=n_taxa))[::-1]
    base[names.index("Streptococcus")] = base.max() * 1.5
    healthy = base / base.sum()

    effect = np.ones(n_taxa)
    candidates = [i for i in range(5, n_taxa) if i not in _red_idx(names)]
    chosen = rng.choice(candidates, size=n_up + n_down - 3, replace=False)
    up_idx = list(_red_idx(names)) + list(chosen[: n_up - 3])
    down_idx = list(chosen[n_up - 3:])
    effect[up_idx] = fold_change
    effect[down_idx] = 1.0 / fold_change
    disease = healthy * effect
    disease /= disease.sum()

    columns, groups, blocks = [], [], []
    for code, n, pi in (("H", n_healthy, healthy), ("P", n_disease, disease)):
        for i in range(n):
            p = rng.dirichlet(theta * pi + 1e-3)
            depth = int(rng.integers(*depth_range))
            blocks.append(rng.multinomial(depth, p))
            columns.append(f"{code}{i + 1}")
            groups.append(code)
    counts = pd.DataFrame(np.column_stack(blocks), index=names, columns=columns)

    taxonomy = pd.DataFrame(index=names, columns=RANKS)
    taxonomy["kingdom"] = "Bacteria"
    taxonomy["phylum"] = [_PHYLA[i % len(_PHYLA)] for i in range(n_taxa)]
    taxonomy["class"] = "Unclassified " + taxonomy["phylum"]
    taxonomy["order"] = taxonomy["class"]
    taxonomy["family"] = taxonomy["class"]
    taxonomy["genus"] = names
    species = {g: s for g, s in zip(config.RED_COMPLEX_GENERA, config.RED_COMPLEX_SPECIES)}
    taxonomy["species"] = [species.get(g, f"Unclassified {g}") for g in names]

    meta = pd.DataFrame({"group_code": groups}, index=pd.Index(columns, name="sample_id"))
    meta["group"] = meta["group_code"].map(config.GROUP_LABELS)
    meta["read_depth"] = counts.sum(axis=0).to_numpy()

    truth = pd.DataFrame({"taxon": names, "log2_fold_change": np.log2(effect)})
    truth["planted"] = np.where(effect > 1, "up", np.where(effect < 1, "down", "none"))
    return Dataset(counts, taxonomy, meta, is_synthetic=True), truth.set_index("taxon")


def _red_idx(names: list[str]) -> list[int]:
    return [names.index(g) for g in config.RED_COMPLEX_GENERA]


def generate_synthetic_dataset() -> Dataset:
    """Fallback dataset used when the real download fails (saves ground truth)."""
    ds, truth = generate_dm()
    config.ensure_dirs()
    truth.to_csv(config.DATA_PROCESSED / f"{config.SYNTHETIC_TAG}_ground_truth.csv")
    return ds
