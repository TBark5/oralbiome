"""Shared fixtures: a small synthetic dataset reused by several tests."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oralbiome import style  # noqa: E402
from oralbiome.synthetic import generate_dm  # noqa: E402

style.apply_style()


@pytest.fixture(scope="session")
def synthetic():
    """(Dataset, truth) with the same group sizes as the real primary comparison."""
    return generate_dm(seed=7)
