"""Tests for M7: in-fold preprocessing and a no-signal sanity check."""
import numpy as np

from oralbiome import classifier as clf


def test_prevalence_clr_learns_columns_from_training_data_only():
    train = np.array([[10, 0, 5], [12, 0, 4], [9, 0, 6]], dtype=float)
    test = np.array([[1, 50, 1]], dtype=float)
    step = clf.PrevalenceCLR(min_prevalence=0.5, min_mean_rel=0.0).fit(train)
    assert step.keep_.tolist() == [True, False, True]
    out = step.transform(test)
    assert out.shape == (1, 2) and np.isclose(out.sum(), 0)


def test_random_labels_give_chance_auc():
    rng = np.random.default_rng(8)
    X = rng.poisson(20, size=(40, 30)).astype(float)
    y = np.r_[np.zeros(20), np.ones(20)].astype(int)
    model = clf.make_models()["L1 logistic regression"]
    res = clf.repeated_cv(model, X, y, np.arange(30).astype(str), n_repeats=2, n_jobs=1,
                          want_importance=False)
    assert 0.2 < res["aucs"].mean() < 0.8


def test_bootstrap_auc_perfect_separation():
    y = np.r_[np.zeros(10), np.ones(10)].astype(int)
    res = clf.bootstrap_auc(y, np.r_[np.zeros(10), np.ones(10)].astype(float), n_boot=200)
    assert res["auc"] == 1.0 and res["ci_low"] == 1.0
