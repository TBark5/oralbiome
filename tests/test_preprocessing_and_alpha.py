"""Tests for M1 transforms and M3 diversity formulas."""
import numpy as np
import pandas as pd

from oralbiome import alpha
from oralbiome.preprocessing import clr, prevalence_filter, rarefy, relative_abundance


def test_diversity_of_uniform_community():
    counts = np.full(8, 25.0)
    assert np.isclose(alpha.shannon(counts), np.log(8))
    assert np.isclose(alpha.simpson(counts), 1 - 1 / 8)
    assert alpha.observed(counts) == 8
    assert np.isclose(alpha.pielou(counts), 1.0)


def test_diversity_of_single_taxon():
    counts = np.array([100.0, 0, 0])
    assert alpha.shannon(counts) == 0 and alpha.simpson(counts) == 0
    assert alpha.observed(counts) == 1 and alpha.pielou(counts) == 0


def test_shannon_hand_calculation():
    p = np.array([0.5, 0.25, 0.25])
    assert np.isclose(alpha.shannon(p * 100), -(p * np.log(p)).sum())


def test_clr_centres_each_sample():
    df = pd.DataFrame({"s1": [10, 0, 5], "s2": [1, 2, 3]})
    assert np.allclose(clr(df).sum(axis=0), 0)


def test_clr_is_scale_invariant_without_zeros():
    x = pd.DataFrame({"s": [1.0, 2.0, 4.0]})
    assert np.allclose(clr(x, pseudocount=0), clr(x * 10, pseudocount=0))


def test_relative_abundance_columns_sum_to_one():
    df = pd.DataFrame({"a": [1, 3], "b": [0, 5]})
    assert np.allclose(relative_abundance(df).sum(axis=0), 1)


def test_rarefy_gives_exact_depth_and_never_adds_reads():
    rng = np.random.default_rng(3)
    df = pd.DataFrame(rng.integers(0, 50, size=(30, 4)), columns=list("abcd"))
    depth = int(df.sum().min())
    out = rarefy(df, depth)
    assert (out.sum(axis=0) == depth).all()
    assert (out.to_numpy() <= df.to_numpy()).all()


def test_prevalence_filter():
    df = pd.DataFrame({"s1": [10, 0, 1000], "s2": [10, 0, 1000], "s3": [10, 1, 1000]},
                      index=["common", "rare", "dominant"])
    keep = prevalence_filter(df, min_prevalence=0.5, min_mean_rel=1e-3)
    assert keep["common"] and keep["dominant"] and not keep["rare"]
