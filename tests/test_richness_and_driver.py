"""Tests for the Chao1 / ACE richness estimators and the driver-taxon report."""
import numpy as np
import pytest

from oralbiome import alpha, config

needs_results = pytest.mark.skipif(not (config.RESULTS / "m7_driver_taxon.json").exists(),
                                   reason="run_all.py has not been run")


def test_chao1_hand_example():
    # S_obs = 5, F1 = 2, F2 = 1 -> 5 + 2*1 / (2*2) = 5.5
    assert np.isclose(alpha.chao1(np.array([1, 1, 2, 3, 5])), 5.5)


def test_chao1_equals_observed_without_singletons():
    counts = np.array([2, 3, 4, 10])
    assert alpha.chao1(counts) == alpha.observed(counts)


def test_ace_hand_example():
    # rare = [1, 1, 2, 3]: S_rare 4, N_rare 7, F1 2, C = 5/7,
    # gamma^2 = 4/(5/7) * (2*1*1 + 3*2*1) / (7*6) - 1 = 0.0667; ACE = 1 + 5.6 + 2.8*0.0667
    expected = 1 + 4 / (5 / 7) + 2 / (5 / 7) * (4 / (5 / 7) * 8 / 42 - 1)
    assert np.isclose(alpha.ace(np.array([1, 1, 2, 3, 15])), expected)


def test_ace_without_singletons_equals_observed():
    counts = np.array([2, 5, 8, 40, 100])
    assert np.isclose(alpha.ace(counts), alpha.observed(counts))


def test_ace_falls_back_to_chao1_when_all_rare_reads_are_singletons():
    counts = np.array([1, 1, 1, 50])
    assert np.isclose(alpha.ace(counts), alpha.chao1(counts))


def test_richness_estimators_are_at_least_observed():
    rng = np.random.default_rng(9)
    for _ in range(20):
        counts = rng.negative_binomial(1, 0.2, size=60)
        if counts.sum() == 0:
            continue
        assert alpha.chao1(counts) >= alpha.observed(counts)
        assert alpha.ace(counts) >= alpha.observed(counts) - 1e-9


def test_alpha_table_has_six_metrics():
    import pandas as pd

    table = alpha.alpha_table(pd.DataFrame({"s1": [1, 2, 3, 0], "s2": [5, 5, 5, 5]}))
    assert list(table.columns) == list(alpha.METRICS)
    assert len(alpha.METRICS) == 6


@needs_results
def test_driver_report_matches_deposit():
    import json

    summary = json.loads((config.RESULTS / "m7_driver_taxon.json").read_text())
    assert summary["driver"] == "Unclassified Bacilli"
    assert summary["has_sequences"] is False
    assert not (config.RESULTS / "m7_driver_taxon.fasta").exists()
    assert sum(summary["detected"].values()) > 0
