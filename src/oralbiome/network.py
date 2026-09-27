"""M6 - Co-occurrence network: Spearman correlation on CLR abundances.

Correlating raw proportions creates spurious negative correlations (if one
taxon goes up, the others must go down). Correlating CLR values reduces that
artefact; it is the approach used by tools such as SparCC and SPIEC-EASI in
spirit, but much simpler. Edges are kept only if |rho| >= 0.6 and the
BH-corrected p-value < 0.05.
"""
from __future__ import annotations

import json

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

from . import config, style
from .preprocessing import Preprocessed
from .stats import benjamini_hochberg

MIN_NETWORK_PREVALENCE: float = 0.30


def correlation_edges(clr_values: pd.DataFrame) -> pd.DataFrame:
    """All taxon pairs with Spearman rho, p and BH q (over all pairs)."""
    rho, p = stats.spearmanr(clr_values.T.to_numpy())
    rho, p = np.atleast_2d(rho), np.atleast_2d(p)
    i, j = np.triu_indices(len(clr_values), 1)
    edges = pd.DataFrame({"source": clr_values.index[i], "target": clr_values.index[j],
                          "rho": rho[i, j], "p": p[i, j]})
    edges["q_bh"] = benjamini_hochberg(edges["p"].to_numpy())
    return edges


def build_graph(edges: pd.DataFrame, min_abs_rho: float = config.NETWORK_MIN_ABS_RHO,
                max_q: float = config.NETWORK_MAX_Q) -> nx.Graph:
    """Graph of edges passing both thresholds (weight = rho)."""
    keep = edges[(edges["rho"].abs() >= min_abs_rho) & (edges["q_bh"] < max_q)]
    graph = nx.Graph()
    for row in keep.itertuples():
        graph.add_edge(row.source, row.target, rho=row.rho, q=row.q_bh)
    return graph


def hub_table(graph: nx.Graph, da: pd.DataFrame | None = None) -> pd.DataFrame:
    """Degree, positive/negative degree and betweenness per node, sorted by degree."""
    between = nx.betweenness_centrality(graph)
    rows = []
    for node in graph.nodes:
        signs = [graph.edges[node, nb]["rho"] > 0 for nb in graph.neighbors(node)]
        rows.append({"taxon": node, "degree": graph.degree[node],
                     "positive_edges": int(sum(signs)), "negative_edges": int(len(signs) - sum(signs)),
                     "betweenness": between[node]})
    table = pd.DataFrame(rows).set_index("taxon").sort_values(["degree", "betweenness"],
                                                              ascending=False)
    if da is not None:
        table["m5_direction"] = da["direction"].reindex(table.index)
        table["m5_p"] = da["p"].reindex(table.index)
    return table


def red_complex_pairs(clr_values: pd.DataFrame) -> pd.DataFrame:
    """Pre-declared P5 test: Spearman rho between the three red-complex genera."""
    genera = [g for g in config.RED_COMPLEX_GENERA if g in clr_values.index]
    rows = []
    for a in range(len(genera)):
        for b in range(a + 1, len(genera)):
            rho, p = stats.spearmanr(clr_values.loc[genera[a]], clr_values.loc[genera[b]])
            rows.append({"pair": f"{genera[a]} - {genera[b]}", "rho": float(rho), "p": float(p)})
    table = pd.DataFrame(rows)
    if len(table):
        table["q_bh_across_pairs"] = benjamini_hochberg(table["p"].to_numpy())
    return table


def _node_color(taxon: str, da: pd.DataFrame) -> str:
    if taxon not in da.index or da.loc[taxon, "p"] >= 0.05:
        return style.OTHER_COLOR
    return (style.GROUP_COLORS["Periodontitis"] if da.loc[taxon, "clr_mean_diff"] > 0
            else style.GROUP_COLORS["Healthy"])


def plot_network(graph: nx.Graph, hubs: pd.DataFrame, da: pd.DataFrame, rel: pd.DataFrame,
                 n: int) -> None:
    """Spring-layout network; node colour = M5 direction, size = mean abundance."""
    full = graph
    graph = graph.subgraph(max(nx.connected_components(graph), key=len)).copy()
    fig, ax = plt.subplots(figsize=(11, 9.5))
    pos = nx.kamada_kawai_layout(graph, weight=None)
    pos_edges = [e for e in graph.edges if graph.edges[e]["rho"] > 0]
    neg_edges = [e for e in graph.edges if graph.edges[e]["rho"] < 0]
    nx.draw_networkx_edges(graph, pos, edgelist=pos_edges, ax=ax, width=0.6, alpha=0.25,
                           edge_color="#52514e")
    nx.draw_networkx_edges(graph, pos, edgelist=neg_edges, ax=ax, width=0.6, alpha=0.45,
                           edge_color="#4a3aa7", style="dashed")
    sizes = [30 + 3000 * np.sqrt(rel.loc[t].mean()) if t in rel.index else 30 for t in graph]
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_size=sizes, linewidths=0.6,
                           node_color=[_node_color(t, da) for t in graph], edgecolors="white")
    label_nodes = list(hubs.index[:10]) + [g for g in config.RED_COMPLEX_GENERA
                                            if g in graph and g not in hubs.index[:10]]
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=style.GROUP_COLORS["Periodontitis"],
                   label="higher in periodontitis (M5 p < 0.05)"),
        plt.Line2D([], [], marker="o", ls="", color=style.GROUP_COLORS["Healthy"],
                   label="higher in healthy (M5 p < 0.05)"),
        plt.Line2D([], [], marker="o", ls="", color=style.OTHER_COLOR, label="no group difference"),
        plt.Line2D([], [], color="#52514e", label="positive correlation"),
        plt.Line2D([], [], color="#4a3aa7", ls="--", label="negative correlation"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
    ax.set_title(f"{style.title_prefix()}Genus co-occurrence network (Spearman on CLR, "
                 f"|rho| >= {config.NETWORK_MIN_ABS_RHO}, q < {config.NETWORK_MAX_Q}; "
                 f"all {n} samples pooled)\nLargest component: {graph.number_of_nodes()} of "
                 f"{full.number_of_nodes()} genera, {graph.number_of_edges()} of "
                 f"{full.number_of_edges()} edges; labels = top-10 hubs by degree and "
                 f"red-complex genera", fontsize=10.5)
    ax.axis("off")
    fig.tight_layout()
    # Numbered markers + a key: names do not fit inside the dense core.
    style.place_labels(ax, [(pos[t][0], pos[t][1], str(i + 1), t in config.RED_COMPLEX_GENERA)
                            for i, t in enumerate(label_nodes)], boxed=True, spread=1.6)
    key = "\n".join(f"{i + 1:>2}. {t}{' (red complex)' if t in config.RED_COMPLEX_GENERA else ''}"
                    f" - degree {graph.degree[t]}" for i, t in enumerate(label_nodes))
    ax.text(1.0, 0.72, "Labelled genera\n" + key, transform=ax.transAxes, va="top", ha="left",
            fontsize=8.5, color=style.INK, family="DejaVu Sans Mono", linespacing=1.4)
    style.save(fig, "11_cooccurrence_network")


def run(pre: Preprocessed, da: pd.DataFrame) -> dict:
    """Build the network, write edges/hubs, test red-complex co-occurrence."""
    prevalence = (pre.genus_filtered > 0).mean(axis=1)
    clr_values = pre.genus_clr.loc[prevalence >= MIN_NETWORK_PREVALENCE]
    edges = correlation_edges(clr_values)
    graph = build_graph(edges)
    hubs = hub_table(graph, da)
    edges[(edges["rho"].abs() >= config.NETWORK_MIN_ABS_RHO) & (edges["q_bh"] < config.NETWORK_MAX_Q)] \
        .to_csv(config.RESULTS / "m6_network_edges.csv", index=False)
    hubs.to_csv(config.RESULTS / "m6_hub_taxa.csv")
    pairs = red_complex_pairs(pre.genus_clr)
    pairs.to_csv(config.RESULTS / "m6_red_complex_pairs.csv", index=False)
    plot_network(graph, hubs, da, pre.genus_rel, pre.genus_clr.shape[1])

    summary = {
        "n_genera_input": int(len(clr_values)),
        "min_prevalence_for_network": MIN_NETWORK_PREVALENCE,
        "n_pairs_tested": int(len(edges)),
        "n_nodes": graph.number_of_nodes(),
        "n_edges": graph.number_of_edges(),
        "n_positive_edges": int(sum(graph.edges[e]["rho"] > 0 for e in graph.edges)),
        "n_negative_edges": int(sum(graph.edges[e]["rho"] < 0 for e in graph.edges)),
        "n_components": nx.number_connected_components(graph) if len(graph) else 0,
        "top_hubs": hubs.head(10)[["degree", "positive_edges", "negative_edges"]]
        .to_dict(orient="index"),
        "red_complex_pairs": pairs.round(5).to_dict(orient="records"),
    }
    (config.RESULTS / "m6_summary.json").write_text(json.dumps(summary, indent=2))
    return summary
