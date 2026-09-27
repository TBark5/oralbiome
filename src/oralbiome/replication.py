"""Secondary analysis: does the H vs P signal replicate in T vs TP?

The hypertensive participants without (T) and with (TP) periodontitis are a
second, independent set of people from the same study. The same pipeline is
re-run on them and effect sizes are compared side by side with the primary
comparison in a forest plot.
"""
from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import alpha, beta, classifier_run, config, style
from .data import Dataset
from .differential import targeted_species_test
from .preprocessing import Preprocessed, preprocess
from .stats import bootstrap_ci, rank_biserial

COHORTS = {"primary": ("P", "H"), "replication": ("TP", "T")}


def effect_rows(pre: Preprocessed, disease: str, reference: str, cohort: str) -> list[dict]:
    """Rank-biserial r (+ bootstrap CI) for alpha metrics and red-complex species."""
    labels = pre.metadata["group_code"]
    rows = []
    alpha_table = alpha.alpha_table(pre.genus_rarefied)
    for metric in ("shannon", "observed", "pielou"):
        x = alpha_table.loc[labels == disease, metric].to_numpy()
        y = alpha_table.loc[labels == reference, metric].to_numpy()
        lo, hi = bootstrap_ci(x, y, rank_biserial)
        rows.append({"cohort": cohort, "measure": f"alpha: {alpha.METRICS[metric]}",
                     "r": rank_biserial(x, y), "ci_low": lo, "ci_high": hi})
    targeted, species_clr = targeted_species_test(pre.species_counts, labels, disease=disease,
                                                  reference=reference)
    for sp in targeted.index:
        x = species_clr.loc[sp, labels.index[labels == disease]].to_numpy()
        y = species_clr.loc[sp, labels.index[labels == reference]].to_numpy()
        lo, hi = bootstrap_ci(x, y, rank_biserial)
        rows.append({"cohort": cohort, "measure": sp.replace("_", " "), "r": rank_biserial(x, y),
                     "ci_low": lo, "ci_high": hi, "p": targeted.loc[sp, "p"],
                     "q_bh_across_targets": targeted.loc[sp, "q_bh_across_targets"]})
    return rows


def plot_forest(table: pd.DataFrame, headline: dict) -> None:
    """Forest plot of effect sizes in both cohorts."""
    measures = list(dict.fromkeys(table["measure"]))
    fig, ax = plt.subplots(figsize=(9, 5.6))
    offsets = {"primary": -0.15, "replication": 0.15}
    colors = {"primary": "#4a3aa7", "replication": "#1baf7a"}
    for cohort, off in offsets.items():
        sub = table[table["cohort"] == cohort].set_index("measure").loc[measures]
        ypos = np.arange(len(measures)) + off
        ax.errorbar(sub["r"], ypos, xerr=[sub["r"] - sub["ci_low"], sub["ci_high"] - sub["r"]],
                    fmt="o", color=colors[cohort], ms=7, capsize=3, lw=1.5,
                    label=headline[cohort]["label"])
    ax.axvline(0, color=style.INK_SECONDARY, lw=1)
    ax.set_yticks(range(len(measures)))
    ax.set_yticklabels(measures)
    ax.invert_yaxis()
    ax.set_xlim(-1, 1)
    ax.set_xlabel("Rank-biserial r (disease vs reference) with bootstrap 95% CI\n"
                  "< 0: lower in periodontitis      > 0: higher in periodontitis")
    text = "\n".join(
        f"{h['label']}: Bray-Curtis PERMANOVA R2={h['permanova_R2']:.3f}, p={h['permanova_p']:.3f}; "
        f"L1-LR AUC={h['lasso_auc']:.2f} (perm. p={h['lasso_perm_p']:.3f})" for h in headline.values())
    ax.set_title(f"{style.title_prefix()}Does the periodontitis signal replicate in a second "
                 f"cohort?\n{text}", fontsize=9.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=8.5)
    fig.tight_layout()
    style.save(fig, "15_replication_forest")


def run(full: Dataset, primary_pre: Preprocessed, primary_m4: dict, primary_m7: dict) -> dict:
    """Run the replication cohort and build the comparison table and figure."""
    disease, reference = COHORTS["replication"]
    rep_ds = full.subset((reference, disease))
    rep = preprocess(rep_ds, "replication")
    labels = rep.metadata["group_code"]

    bc = beta.bray_curtis(rep.genus_rel)
    perm = beta.permanova(bc, labels)
    disp = beta.permdisp(bc, labels)
    alpha_tests = alpha.compare_alpha(alpha.alpha_table(rep.genus_rarefied), rep.metadata,
                                      disease, reference, "replication_genus")
    m7 = classifier_run.run(rep, disease=disease, reference=reference, figures=False,
                            tag="replication")

    rows = effect_rows(primary_pre, *COHORTS["primary"], "primary") + \
        effect_rows(rep, disease, reference, "replication")
    table = pd.DataFrame(rows)
    table.to_csv(config.RESULTS / "replication_effect_sizes.csv", index=False)
    alpha_tests.to_csv(config.RESULTS / "replication_alpha_tests.csv", index=False)

    n_rep = {c: int((labels == c).sum()) for c in (reference, disease)}
    n_pri = {c: int((primary_pre.metadata["group_code"] == c).sum()) for c in COHORTS["primary"]}
    headline = {
        "primary": {"label": f"Primary: H (n={n_pri['H']}) vs P (n={n_pri['P']})",
                    "permanova_R2": primary_m4["bray_curtis"]["permanova"]["R2"],
                    "permanova_p": primary_m4["bray_curtis"]["permanova"]["p"],
                    "lasso_auc": primary_m7["L1 logistic regression"]["auc_mean_over_repeats"],
                    "lasso_perm_p": primary_m7["L1 logistic regression"]["permutation_p"]},
        "replication": {"label": f"Replication: T (n={n_rep['T']}) vs TP (n={n_rep['TP']})",
                        "permanova_R2": perm["R2"], "permanova_p": perm["p"],
                        "lasso_auc": m7["L1 logistic regression"]["auc_mean_over_repeats"],
                        "lasso_perm_p": m7["L1 logistic regression"]["permutation_p"]},
    }
    plot_forest(table, headline)
    summary = {"n": n_rep, "permanova_bray_curtis": perm, "permdisp_bray_curtis": disp,
               "alpha": alpha_tests.set_index("metric")[["p", "rank_biserial"]].to_dict(orient="index"),
               "classifier": {k: m7[k] for k in ("L1 logistic regression", "Random forest")},
               "classifier_ablation": m7["ablation_without_top_genus"],
               "effect_sizes": table.round(4).to_dict(orient="records")}
    (config.RESULTS / "replication_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    return summary
