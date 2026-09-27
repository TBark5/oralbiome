"""Characterise the taxon that drives the M7 classifier.

The deposited data contain feature IDs and taxonomy strings but no
representative sequences, so the taxon cannot be identified by sequence
search from these files. This module documents everything that *can* be
measured: which features make up the label, where they are detected, how
abundant they are, whether they sit in the flagged or high-richness samples,
and what they correlate with in the M6 network.
"""
from __future__ import annotations

import json
import zipfile

import numpy as np
import pandas as pd
from scipy import stats

from . import config
from .data import Dataset

GROUPS = ("H", "P", "T", "TP")


def _has_sequences() -> bool:
    """True if the raw deposit contains any FASTA / representative-sequence file."""
    with zipfile.ZipFile(config.RAW_ZIP) as zf:
        names = [n.lower() for n in zf.namelist()]
    return any(n.endswith((".fa", ".fasta", ".fna", ".qza")) or "seq" in n for n in names)


def profile(full: Dataset, driver: str, high_richness: set[str]) -> dict:
    """Prevalence, abundance and feature make-up of ``driver`` across all four groups."""
    members = full.taxonomy.index[full.taxonomy["genus"] == driver]
    counts = full.counts.loc[members]
    rel = counts.sum(axis=0) / full.counts.sum(axis=0)
    meta = full.metadata
    flagged = meta.get("flagged_atypical", pd.Series(False, index=meta.index)).astype(bool)
    rows = []
    for code in GROUPS:
        ids = meta.index[meta["group_code"] == code]
        r = rel[ids]
        rows.append({"group": code, "n": len(ids), "detected": int((r > 0).sum()),
                     "median_rel_abund_pct": float(r.median() * 100),
                     "max_rel_abund_pct": float(r.max() * 100),
                     "total_reads": int(counts[ids].sum().sum())})
    by_group = pd.DataFrame(rows).set_index("group")
    top = counts.sum(axis=1).sort_values(ascending=False)
    features = pd.DataFrame({"total_reads": top})
    for code in GROUPS:
        ids = meta.index[meta["group_code"] == code]
        features[f"detected_{code}"] = (counts.loc[top.index, ids] > 0).sum(axis=1)

    special = {
        "flagged": list(meta.index[flagged]),
        "high_richness_healthy": sorted(high_richness),
    }
    detected_in = {k: [s for s in v if rel.get(s, 0) > 0] for k, v in special.items()}
    table = [[by_group.loc["P", "detected"], by_group.loc["P", "n"] - by_group.loc["P", "detected"]],
             [by_group.loc["H", "detected"], by_group.loc["H", "n"] - by_group.loc["H", "detected"]]]
    rep = [[by_group.loc["TP", "detected"], by_group.loc["TP", "n"] - by_group.loc["TP", "detected"]],
           [by_group.loc["T", "detected"], by_group.loc["T", "n"] - by_group.loc["T", "detected"]]]
    return {
        "by_group": by_group, "features": features, "detected_in": detected_in,
        "fisher_presence_P_vs_H_p": float(stats.fisher_exact(table)[1]),
        "fisher_presence_TP_vs_T_p": float(stats.fisher_exact(rep)[1]),
        "lineage": full.taxonomy.loc[members[0]].to_dict(),
    }


def network_partners(driver: str) -> pd.DataFrame:
    """Edges of ``driver`` in the M6 co-occurrence network."""
    edges = pd.read_csv(config.RESULTS / "m6_network_edges.csv")
    mine = edges[(edges["source"] == driver) | (edges["target"] == driver)].copy()
    mine["partner"] = np.where(mine["source"] == driver, mine["target"], mine["source"])
    return mine[["partner", "rho", "q_bh"]].sort_values("rho", key=np.abs, ascending=False)


def run(full: Dataset) -> dict:
    """Write results/m7_driver_taxon.md and the supporting CSV."""
    importance = pd.read_csv(config.RESULTS / "m7_importance_coef_primary.csv", index_col=0)
    driver = importance.index[0]
    alpha_ps = pd.read_csv(config.RESULTS / "m3_alpha_per_sample.csv", index_col=0)
    high = set(alpha_ps.index[(alpha_ps["asv_observed"] > 1000) & (alpha_ps["group"] == "Healthy")])
    prof = profile(full, driver, high)
    partners = network_partners(driver)
    has_seq = _has_sequences()
    prof["by_group"].to_csv(config.RESULTS / "m7_driver_taxon_by_group.csv")
    prof["features"].to_csv(config.RESULTS / "m7_driver_taxon_features.csv")

    bg, feats = prof["by_group"], prof["features"]
    top_id = feats.index[0]
    lines = [
        f"# Driver taxon of the M7 classifier: {driver}",
        "",
        f"Selected in {importance.iloc[0]['selection_frequency']:.0%} of L1 folds "
        f"(mean standardised coefficient {importance.iloc[0]['mean']:+.3f}). Lineage in the deposit: "
        + "; ".join(f"{k} = {v}" for k, v in prof["lineage"].items()) + ".",
        "",
        "## Can it be identified?",
        "",
        ("The deposit contains representative sequences." if has_seq else
         "No. The figshare deposit contains a feature table, taxonomy strings and sample "
         "groups but no representative sequences, so there is nothing to extract to FASTA "
         "and nothing to search against a 16S reference. `m7_driver_taxon.fasta` was "
         "therefore not created. Identification would require re-processing the raw reads "
         "(SRA BioProject PRJNA1304526) to recover the sequences of these features."),
        "",
        "## Prevalence and abundance by group",
        "",
        "| Group | n | Detected | Median rel. abundance (%) | Max rel. abundance (%) | Reads |",
        "|---|---|---|---|---|---|",
    ]
    for code, r in bg.iterrows():
        n, det = int(r["n"]), int(r["detected"])
        lines.append(f"| {code} | {n} | {det}/{n} | {r['median_rel_abund_pct']:.4f} "
                     f"| {r['max_rel_abund_pct']:.4f} | {int(r['total_reads'])} |")
    lines += [
        "",
        f"Presence/absence, Fisher's exact test (descriptive, not multiplicity-corrected): "
        f"P vs H p = {prof['fisher_presence_P_vs_H_p']:.2g}; TP vs T p = "
        f"{prof['fisher_presence_TP_vs_T_p']:.2g}.",
        "",
        "## Which features make up the label",
        "",
        f"{len(feats)} features carry this label. The most abundant, {top_id}, holds "
        f"{int(feats.iloc[0]['total_reads'])} of {int(feats['total_reads'].sum())} reads and is "
        f"detected in H {feats.iloc[0]['detected_H']}/16, P {feats.iloc[0]['detected_P']}/18, "
        f"T {feats.iloc[0]['detected_T']}/17, TP {feats.iloc[0]['detected_TP']}/16 samples. "
        f"The other {len(feats) - 1} features contribute the remaining reads; "
        f"{int((feats[[f'detected_{g}' for g in GROUPS]].sum(axis=1) == 1).sum())} "
        "of all features are seen in a single sample only. Full table: `m7_driver_taxon_features.csv`.",
        "",
        "## Overlap with flagged and high-richness samples",
        "",
        f"- Flagged samples where it is detected: {', '.join(prof['detected_in']['flagged']) or 'none'} "
        f"(of {', '.join(sorted(full.metadata.index[full.metadata['flagged_atypical'].astype(bool)]))}).",
        f"- High-richness healthy samples where it is detected: "
        f"{', '.join(prof['detected_in']['high_richness_healthy']) or 'none'} "
        f"(of {', '.join(sorted(high))}).",
        "",
        "## Network partners (M6, |rho| >= 0.6, q < 0.05)",
        "",
    ]
    if len(partners):
        lines += ["| Partner | rho | q |", "|---|---|---|"]
        lines += [f"| {p.partner} | {p.rho:+.3f} | {p.q_bh:.3g} |" for p in partners.itertuples()]
    else:
        lines.append("No edges passed the thresholds.")
    (config.RESULTS / "m7_driver_taxon.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    summary = {"driver": driver, "has_sequences": has_seq, "n_features": int(len(feats)),
               "top_feature": top_id, "top_feature_reads": int(feats.iloc[0]["total_reads"]),
               "total_reads": int(feats["total_reads"].sum()),
               "detected": bg["detected"].to_dict(),
               "fisher_presence_P_vs_H_p": prof["fisher_presence_P_vs_H_p"],
               "fisher_presence_TP_vs_T_p": prof["fisher_presence_TP_vs_T_p"],
               "network_partners": partners.round(4).to_dict(orient="records")}
    (config.RESULTS / "m7_driver_taxon.json").write_text(json.dumps(summary, indent=2, default=int))
    return summary
