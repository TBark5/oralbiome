"""ORALBIOME dashboard: one tab per analysis module (M1-M7) plus replication.

Run with:  streamlit run app.py
All numbers are read from /results and /data/processed, which run_all.py
writes; the app does not re-run the statistics.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from oralbiome import config  # noqa: E402
from oralbiome.preprocessing import aggregate, relative_abundance  # noqa: E402

RESULTS, FIGURES, PROCESSED = config.RESULTS, config.FIGURES, config.DATA_PROCESSED
GROUP_DOMAIN = ["Healthy", "Periodontitis"]
GROUP_RANGE = ["#0072B2", "#D55E00"]
TAXON_RANGE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
               "#4a3aa7", "#e34948", "#c3c2b7"]
LEVELS = ["phylum", "class", "order", "family", "genus"]

st.set_page_config(page_title="ORALBIOME", page_icon=":material/biotech:", layout="wide")


@st.cache_data(max_entries=4)
def load_json(name: str) -> dict:
    return json.loads((RESULTS / name).read_text())


@st.cache_data(max_entries=32)
def load_csv(name: str, folder: str = "results", index_col: int | None = None) -> pd.DataFrame:
    base = RESULTS if folder == "results" else PROCESSED
    return pd.read_csv(base / name, index_col=index_col)


@st.cache_data(max_entries=8)
def level_table(level: str) -> pd.DataFrame:
    """Relative abundance (taxa x primary samples) at any taxonomic rank."""
    meta = load_csv("primary_metadata.csv", "processed", index_col=0)
    counts = pd.read_csv(PROCESSED / "asv_counts.csv.gz", index_col=0)[meta.index]
    taxonomy = pd.read_csv(PROCESSED / "asv_taxonomy.csv.gz", index_col=0)
    return relative_abundance(aggregate(counts, taxonomy, level))


def group_color() -> alt.Color:
    return alt.Color("group:N", scale=alt.Scale(domain=GROUP_DOMAIN, range=GROUP_RANGE),
                     legend=alt.Legend(title="Group", orient="top"))


def tidy(df: pd.DataFrame) -> pd.DataFrame:
    """Round floats for display; keep p/q values to 3 significant digits."""
    out = df.copy()
    for col in out.select_dtypes("float").columns:
        if col == "p" or col.startswith(("p_", "q_")):
            out[col] = out[col].map(lambda v: float(f"{v:.3g}"))
        else:
            out[col] = out[col].round(4)
    return out


def figure(name: str, caption: str) -> None:
    path = FIGURES / name
    if path.exists():
        st.image(str(path), caption=caption)


if not (RESULTS / "summary.json").exists():
    st.error("No results found. Run `python run_all.py` first, then reload this page.",
             icon=":material/error:")
    st.stop()

summary = load_json("summary.json")
meta_all = load_csv("primary_metadata.csv", "processed", index_col=0)

# ---------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("Filters")
    groups = st.pills("Groups", GROUP_DOMAIN, selection_mode="multi", default=GROUP_DOMAIN,
                      key="groups")
    include_flagged = st.toggle("Include flagged atypical samples", value=True, key="flagged",
                                help="Samples where most reads come from features seen in no "
                                     "other sample (possible contamination).")
    level = st.selectbox("Taxonomic level (composition tab)", LEVELS, index=4, key="level")
    top_n = st.slider("Taxa shown in composition", 3, 8, 8, key="top_n")
    st.caption("Filters change what is plotted. Test statistics always come from the full "
               "primary comparison (healthy n=16, periodontitis n=18) saved in /results.")

groups = groups or GROUP_DOMAIN
meta = meta_all[meta_all["group"].isin(groups)]
if not include_flagged:
    meta = meta[~meta["flagged_atypical"].astype(bool)]
counts_n = meta["group"].value_counts().reindex(GROUP_DOMAIN).dropna().astype(int)
n_text = ", ".join(f"{g} n={n}" for g, n in counts_n.items())

st.title("ORALBIOME: oral microbiome dysbiosis in periodontitis")
st.caption(f"Saliva 16S rRNA data from the source study (figshare, CC BY 4.0; full citation in "
           f"DATA_SOURCE.md). Showing: {n_text}. "
           "All results are associations from a cross-sectional study, not causal findings "
           "and not a diagnostic test.")

tabs = st.tabs(["Overview", "M1 Preprocessing", "M2 Composition", "M3 Alpha diversity",
                "M4 Beta diversity", "M5 Differential abundance", "M6 Network",
                "M7 Classifier", "Replication"])

# ---------------------------------------------------------------- overview
with tabs[0]:
    m4, m5, m7 = summary["m4"]["bray_curtis"], summary["m5"], summary["m7"]
    with st.container(horizontal=True):
        st.metric("Samples (H / P)", f"{m7['n']['H']} / {m7['n']['P']}", border=True)
        st.metric("PERMANOVA R² (Bray-Curtis)", f"{m4['permanova']['R2']:.3f}",
                  f"p = {m4['permanova']['p']:.3f}", delta_color="off", delta_arrow="off", border=True)
        st.metric("Genera passing FDR 5%", f"{m5['n_significant_fdr05']} / {m5['n_genera_tested']}",
                  border=True)
        st.metric("L1-LR AUC (95% CI)",
                  f"{m7['L1 logistic regression']['auc_pooled_oof']:.2f}",
                  "[{:.2f}, {:.2f}]".format(*m7["L1 logistic regression"]["auc_bootstrap_ci95"]),
                  delta_color="off", delta_arrow="off", border=True)
    st.markdown("**Question.** Does salivary bacterial composition differ between periodontally "
                "healthy people and people with periodontitis, and can composition alone "
                "predict the group? See `HYPOTHESIS.md` and `RESULTS_DISCUSSION.md`.")
    figure("00_data_quality.png", "Read depth and rarefaction curves for all 67 samples.")

# ---------------------------------------------------------------- M1
with tabs[1]:
    m1 = summary["m1"]
    with st.container(horizontal=True):
        st.metric("Genera before filter", m1["n_genera_input"], border=True)
        st.metric("Genera after filter", m1["n_genera_after_filter"], border=True)
        st.metric("Median reads kept", f"{m1['reads_retained_after_filter_median']:.1%}",
                  border=True)
        st.metric("Rarefaction depth", f"{m1['rarefy_depth']:,}", border=True)
    st.markdown(f"A genus is kept if it is detected in at least {m1['min_prevalence']:.0%} of "
                f"samples and has mean relative abundance of at least {m1['min_mean_rel_abund']:.2%}. "
                f"Zeros get a pseudocount of {m1['pseudocount']} before the CLR transform.")
    st.dataframe(meta_all[["group", "read_depth", "flagged_atypical"]], height=300)

# ---------------------------------------------------------------- M2
with tabs[2]:
    rel = level_table(level)[meta.index]
    top = rel.mean(axis=1).sort_values(ascending=False).index[:top_n]
    shown = pd.concat([rel.loc[top], rel.drop(top).sum().to_frame("Other").T])
    long = shown.T.join(meta["group"]).reset_index(names="sample").melt(
        id_vars=["sample", "group"], var_name="taxon", value_name="abundance")
    order = list(top) + ["Other"]
    color = alt.Color("taxon:N", scale=alt.Scale(domain=order, range=TAXON_RANGE[:top_n] + ["#c3c2b7"]),
                      sort=order, legend=alt.Legend(title=level.capitalize(), orient="right"))
    with st.container(border=True):
        st.subheader(f"Per-sample composition ({level}; {n_text})")
        st.altair_chart(alt.Chart(long).mark_bar().encode(
            x=alt.X("sample:N", sort=list(meta.sort_values("group").index), title=None),
            y=alt.Y("abundance:Q", stack="normalize", axis=alt.Axis(format="%"), title="Relative abundance"),
            color=color, order=alt.Order("taxon_order:Q"),
            tooltip=["sample", "group", "taxon", alt.Tooltip("abundance:Q", format=".1%")],
        ).transform_calculate(taxon_order=f"indexof({json.dumps(order)}, datum.taxon)").properties(height=360))
    means = long.groupby(["group", "taxon"], as_index=False)["abundance"].mean()
    with st.container(border=True):
        st.subheader(f"Group means ({n_text})")
        st.altair_chart(alt.Chart(means).mark_bar().encode(
            x=alt.X("abundance:Q", stack="normalize", axis=alt.Axis(format="%"), title="Mean relative abundance"),
            y=alt.Y("group:N", title=None), color=color,
            tooltip=["group", "taxon", alt.Tooltip("abundance:Q", format=".1%")],
        ).properties(height=140))
    st.markdown("**Core genera** (detected in at least 90% of a group)")
    st.dataframe(load_csv("m2_core_microbiome.csv", index_col=0), height=280)

# ---------------------------------------------------------------- M3
with tabs[3]:
    alpha = load_csv("m3_alpha_per_sample.csv", index_col=0).loc[meta.index]
    metric = st.segmented_control("Metric", ["shannon", "simpson", "observed", "pielou", "chao1", "ace"],
                                  default="shannon", key="alpha_metric") or "shannon"
    st.subheader(f"Alpha diversity, genus level ({n_text})")
    data = alpha.reset_index(names="sample")
    base = alt.Chart(data).encode(x=alt.X("group:N", title=None, sort=GROUP_DOMAIN,
                                          axis=alt.Axis(labelAngle=0)))
    chart = base.mark_boxplot(size=60, opacity=0.35, extent="min-max").encode(
        y=alt.Y(f"genus_{metric}:Q", title=metric.capitalize(), scale=alt.Scale(zero=False)),
        color=group_color()) + base.mark_circle(size=60).encode(
        y=f"genus_{metric}:Q", color=group_color(),
        xOffset=alt.XOffset("jitter:Q", scale=alt.Scale(domain=[-2.5, 2.5])),
        tooltip=["sample", "group", alt.Tooltip(f"genus_{metric}:Q", format=".3f")]
    ).transform_calculate(jitter="random() - 0.5")
    st.altair_chart(chart.properties(height=380))
    tests = load_csv("m3_alpha_tests.csv")
    st.markdown("**Mann-Whitney U tests** (rank-biserial r > 0 means higher in periodontitis)")
    st.dataframe(tidy(tests), hide_index=True)

# ---------------------------------------------------------------- M4
with tabs[4]:
    dist = st.segmented_control("Distance", ["bray_curtis", "jaccard"], default="bray_curtis",
                                format_func=lambda d: d.replace("_", "-").title(),
                                key="distance") or "bray_curtis"
    coords = load_csv(f"m4_ordination_{dist}.csv", index_col=0).loc[meta.index]
    method = st.segmented_control("Ordination", ["PCoA", "NMDS"], default="PCoA",
                                  key="ordination") or "PCoA"
    xcol, ycol = ("PCo1", "PCo2") if method == "PCoA" else ("NMDS1", "NMDS2")
    st.subheader(f"{method} of {dist.replace('_', '-')} distances ({n_text})")
    pts = coords.join(meta[["group", "flagged_atypical"]]).reset_index(names="sample")
    st.altair_chart(alt.Chart(pts).mark_point(size=90, filled=True).encode(
        x=alt.X(f"{xcol}:Q"), y=alt.Y(f"{ycol}:Q"), color=group_color(),
        shape=alt.Shape("flagged_atypical:N", title="Flagged"),
        tooltip=["sample", "group", "flagged_atypical"]).properties(height=420))
    res = summary["m4"][dist]
    with st.container(horizontal=True):
        st.metric("PERMANOVA R²", f"{res['permanova']['R2']:.3f}", border=True)
        st.metric("PERMANOVA p", f"{res['permanova']['p']:.3f}", border=True)
        st.metric("PERMDISP p", f"{res['permdisp']['p']:.3f}", border=True)
        st.metric("PERMANOVA p without flagged", f"{res['permanova_excluding_flagged']['p']:.3f}",
                  border=True)
    figure("07_distance_heatmap.png", "Bray-Curtis dissimilarity between all primary samples.")

# ---------------------------------------------------------------- M5
with tabs[5]:
    da = load_csv("m5_differential_abundance_genus.csv", index_col=0).reset_index()
    st.subheader("Volcano plot, 354 genera (healthy n=16, periodontitis n=18)")
    da["neg_log10_p"] = -np.log10(da["p"])
    da["status"] = da.apply(lambda r: "q < 0.05" if r["q_bh"] < 0.05 else
                            ("p < 0.05" if r["p"] < 0.05 else "n.s."), axis=1)
    st.altair_chart(alt.Chart(da).mark_circle(size=45).encode(
        x=alt.X("clr_mean_diff:Q", title="CLR difference (periodontitis - healthy)"),
        y=alt.Y("neg_log10_p:Q", title="-log10 p"),
        color=alt.Color("status:N", scale=alt.Scale(domain=["q < 0.05", "p < 0.05", "n.s."],
                                                   range=["#D55E00", "#eda100", "#c3c2b7"])),
        tooltip=["taxon", alt.Tooltip("clr_mean_diff:Q", format=".2f"),
                 alt.Tooltip("p:Q", format=".4f"), alt.Tooltip("q_bh:Q", format=".3f")],
    ).properties(height=400))
    st.caption(f"Null control: with shuffled labels, on average {m5['null_fdr_hits_mean']:.2f} genera "
               f"pass FDR ({m5['n_null_permutations']} shuffles). Observed raw p < 0.05: "
               f"{m5['n_raw_p_below_005']} vs {m5['expected_raw_p_below_005_by_chance']} expected by chance.")
    st.markdown("**Pre-declared red-complex species test**")
    st.dataframe(tidy(load_csv("m5_red_complex_species.csv", index_col=0)))
    query = st.text_input("Filter ranked table by genus name", key="da_query")
    table = da if not query else da[da["taxon"].str.contains(query, case=False)]
    st.dataframe(tidy(table.drop(columns=["neg_log10_p", "status"])), hide_index=True, height=320)

# ---------------------------------------------------------------- M6
with tabs[6]:
    m6 = summary["m6"]
    with st.container(horizontal=True):
        st.metric("Genera in network", m6["n_nodes"], border=True)
        st.metric("Edges", m6["n_edges"], border=True)
        st.metric("Positive / negative", f"{m6['n_positive_edges']} / {m6['n_negative_edges']}",
                  border=True)
    figure("11_cooccurrence_network.png", "Spearman co-occurrence network on CLR abundances.")
    st.markdown("**Hub genera** (by degree) and red-complex pairwise correlations")
    st.dataframe(tidy(load_csv("m6_hub_taxa.csv", index_col=0).head(20)))
    st.dataframe(tidy(load_csv("m6_red_complex_pairs.csv")), hide_index=True)

# ---------------------------------------------------------------- M7
with tabs[7]:
    for name in ("L1 logistic regression", "Random forest"):
        r = m7[name]
        with st.container(horizontal=True):
            st.metric(f"{name}: AUC", f"{r['auc_pooled_oof']:.2f}",
                      "95% CI [{:.2f}, {:.2f}]".format(*r["auc_bootstrap_ci95"]), delta_color="off", delta_arrow="off",
                      border=True)
            st.metric("Permuted-label AUC (mean)", f"{r['null_auc_mean']:.2f}", border=True)
            st.metric("Permutation p", f"{r['permutation_p']:.4f}", border=True)
            st.metric("AUC without top genus", f"{m7['ablation_without_top_genus'][name]:.2f}",
                      border=True)
    left, right = st.columns(2)
    with left:
        figure("12_roc_curves.png", "Out-of-fold ROC curves with bootstrap 95% bands.")
    with right:
        figure("14_permutation_null.png", "AUCs from label-shuffled data (negative control).")
    figure("13_feature_importance.png", "Genera used by the fold-level models.")

# ---------------------------------------------------------------- replication
with tabs[8]:
    rep = summary["replication"]
    st.markdown("The same pipeline re-run on hypertensive participants without (T) and with (TP) "
                "periodontitis, a second, independent set of people from the same study.")
    with st.container(horizontal=True):
        st.metric("Samples (T / TP)", f"{rep['n']['T']} / {rep['n']['TP']}", border=True)
        st.metric("PERMANOVA p", f"{rep['permanova_bray_curtis']['p']:.3f}", border=True)
        st.metric("L1-LR AUC", f"{rep['classifier']['L1 logistic regression']['auc_pooled_oof']:.2f}",
                  border=True)
    figure("15_replication_forest.png", "Effect sizes in both cohorts.")
    st.dataframe(tidy(load_csv("replication_effect_sizes.csv")), hide_index=True)
