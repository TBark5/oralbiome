"""Positive control: does the pipeline recover differences we planted ourselves?

A Dirichlet-multinomial dataset with the same group sizes as the real data
is simulated with 10 genera made 4x more abundant and 5 genera 4x less
abundant in the "Periodontitis" group. M1 + M4 + M5 are run on it and the
calls are compared with the known truth.
"""
from __future__ import annotations

import json

import numpy as np

from . import beta, config
from .differential import differential_table
from .preprocessing import preprocess, relative_abundance
from .synthetic import generate_dm


def recovery_metrics(seed: int = config.SEED, write: bool = True) -> dict:
    """Sensitivity and false discovery proportion of M5 on planted data."""
    ds, truth = generate_dm(seed=seed)
    pre = preprocess(ds, tag="synthetic_validation", write=False)
    labels = pre.metadata["group_code"]
    table = differential_table(pre.genus_clr, relative_abundance(pre.genus_counts), labels)
    table = table.join(truth, how="left")
    called = table["q_bh"] < config.ALPHA
    planted = table["planted"] != "none"
    direction_ok = np.sign(table["clr_mean_diff"]) == np.sign(table["log2_fold_change"])
    true_pos = int((called & planted & direction_ok).sum())
    false_pos = int((called & ~planted).sum())
    perm = beta.permanova(beta.bray_curtis(pre.genus_rel), labels, n_perm=499, seed=seed)
    result = {
        "n_taxa_tested": int(len(table)),
        "n_planted_tested": int(planted.sum()),
        "n_called_fdr05": int(called.sum()),
        "true_positives_correct_direction": true_pos,
        "false_positives": false_pos,
        "sensitivity": true_pos / max(int(planted.sum()), 1),
        "false_discovery_proportion": false_pos / max(int(called.sum()), 1),
        "red_complex_recovered": [g for g in config.RED_COMPLEX_GENERA
                                  if g in table.index and called.get(g, False)],
        "permanova_R2": perm["R2"],
        "permanova_p": perm["p"],
    }
    if write:
        (config.RESULTS / "synthetic_validation.json").write_text(json.dumps(result, indent=2))
    return result
