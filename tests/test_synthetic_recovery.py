"""Positive control: the pipeline must recover differences planted in synthetic data.

The Dirichlet-multinomial generator makes 10 genera 4x more abundant (always
including the three red-complex genera) and 5 genera 4x less abundant in the
"Periodontitis" group. If the pipeline cannot find these, its results on real
data would mean nothing.
"""
import numpy as np

from oralbiome import beta, config
from oralbiome import classifier as clf
from oralbiome.differential import differential_table, permutation_null
from oralbiome.preprocessing import preprocess, relative_abundance
from oralbiome.validation import recovery_metrics


def test_generator_shapes_and_truth(synthetic):
    ds, truth = synthetic
    assert ds.counts.shape == (150, 34)
    assert (ds.metadata["group_code"] == "H").sum() == 16
    assert (truth["planted"] == "up").sum() == 10 and (truth["planted"] == "down").sum() == 5
    assert all(truth.loc[g, "planted"] == "up" for g in config.RED_COMPLEX_GENERA)


def test_differential_abundance_recovers_planted_taxa(synthetic):
    ds, truth = synthetic
    pre = preprocess(ds, tag="test", write=False)
    labels = pre.metadata["group_code"]
    table = differential_table(pre.genus_clr, relative_abundance(pre.genus_counts), labels)
    table = table.join(truth)
    called = table["q_bh"] < config.ALPHA
    planted = table["planted"] != "none"
    sensitivity = (called & planted).sum() / planted.sum()
    fdp = (called & ~planted).sum() / max(called.sum(), 1)
    assert sensitivity >= 0.5
    assert fdp <= 0.25
    # every call on a planted taxon must point the right way
    right_way = np.sign(table["clr_mean_diff"]) == np.sign(table["log2_fold_change"])
    assert right_way[called & planted].all()


def test_null_control_finds_almost_nothing(synthetic):
    ds, _ = synthetic
    pre = preprocess(ds, tag="test", write=False)
    null = permutation_null(pre.genus_clr, pre.metadata["group_code"], n_perm=30)
    assert null["n_q_below_fdr"].mean() < 1.0


def test_permanova_detects_planted_shift(synthetic):
    ds, _ = synthetic
    pre = preprocess(ds, tag="test", write=False)
    res = beta.permanova(beta.bray_curtis(pre.genus_rel), pre.metadata["group_code"], n_perm=199)
    assert res["p"] < 0.05


def test_classifier_separates_planted_groups(synthetic):
    ds, _ = synthetic
    X = ds.counts.T.to_numpy(dtype=float)
    y = (ds.metadata["group_code"] == "P").to_numpy(dtype=int)
    model = clf.make_models()["L1 logistic regression"]
    res = clf.repeated_cv(model, X, y, ds.counts.index.to_numpy(), n_repeats=2, n_jobs=1,
                          want_importance=False)
    assert res["aucs"].mean() > 0.8


def test_recovery_metrics_default_seed():
    res = recovery_metrics(write=False)
    assert res["sensitivity"] >= 0.5 and res["false_discovery_proportion"] <= 0.25
    assert len(res["red_complex_recovered"]) >= 2
