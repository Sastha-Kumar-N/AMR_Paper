#!/usr/bin/env python3
"""Build a stringent, readable functional association network."""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict

import networkx as nx
import pandas as pd

from analysis_common import ROOT, atomic_write_tsv


PHI_THRESHOLD = 0.50
Q_THRESHOLD = 0.05


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    files = [
        "AMP_ARG_associations.tsv",
        "AMP_BGC_associations.tsv",
        "ARG_BGC_associations.tsv",
        "AMR_gene_cooccurrence.tsv",
    ]
    colocalization = pd.read_csv(ROOT / "04_statistics/feature_colocalization.tsv", sep="\t", dtype=str, keep_default_na=False)
    coloc_lookup = {
        (row.association_source, row.feature_1, row.feature_2): row
        for row in colocalization.itertuples(index=False)
    }
    graph = nx.Graph()
    edge_rows = []
    node_metadata: dict[str, dict[str, object]] = {}
    for filename in files:
        frame = pd.read_csv(ROOT / f"04_statistics/{filename}", sep="\t")
        frame = frame[(frame["q_value"] < Q_THRESHOLD) & (frame["Phi_coefficient"].abs() >= PHI_THRESHOLD)]
        for row in frame.itertuples(index=False):
            for feature, category, annotation, prevalence, n_genomes in (
                (row.feature_1, row.category_1, row.annotation_1, row.prevalence_percent_1, row.n_genomes_1),
                (row.feature_2, row.category_2, row.annotation_2, row.prevalence_percent_2, row.n_genomes_2),
            ):
                node_metadata[str(feature)] = {
                    "functional_category": category,
                    "annotation": annotation,
                    "prevalence_percent": prevalence,
                    "n_genomes": n_genomes,
                }
            weight = abs(float(row.Phi_coefficient))
            graph.add_edge(
                str(row.feature_1),
                str(row.feature_2),
                weight=weight,
                signed_weight=float(row.Phi_coefficient),
                distance=1.0 / weight,
            )
            coloc = coloc_lookup.get((filename, str(row.feature_1), str(row.feature_2)))
            edge_rows.append(
                {
                    "source": row.feature_1,
                    "target": row.feature_2,
                    "source_category": row.category_1,
                    "target_category": row.category_2,
                    "association_family": filename.removesuffix(".tsv"),
                    "association_sign": "positive" if row.Phi_coefficient > 0 else "negative",
                    "Phi_coefficient": row.Phi_coefficient,
                    "absolute_Phi": weight,
                    "Jaccard_index": row.Jaccard_index,
                    "n_both": row.n_both,
                    "Fisher_odds_ratio": row.Fisher_odds_ratio,
                    "q_value": row.q_value,
                    "same_contig_genomes": getattr(coloc, "n_genomes_same_contig", "") if coloc else "",
                    "within_10kb_genomes": getattr(coloc, "n_genomes_distance_le_10kb", "") if coloc else "",
                    "physical_linkage_interpretation": getattr(coloc, "physical_linkage_interpretation", "") if coloc else "",
                    "network_filter": f"q<{Q_THRESHOLD} and |Phi|>={PHI_THRESHOLD}",
                }
            )

    if graph.number_of_nodes() == 0:
        raise ValueError("No edges survived the stringent network filter")
    communities = nx.community.louvain_communities(graph, weight="weight", seed=1350)
    ordered_communities = sorted(communities, key=lambda nodes: (-len(nodes), sorted(nodes)[0]))
    community_for_node = {
        node: f"COMMUNITY_{index:03d}"
        for index, nodes in enumerate(ordered_communities, start=1)
        for node in nodes
    }
    betweenness = nx.betweenness_centrality(graph, weight="distance", normalized=True)
    weighted_degree = dict(graph.degree(weight="weight"))
    node_rows = []
    for node in graph.nodes:
        metadata = node_metadata[node]
        node_rows.append(
            {
                "node": node,
                **metadata,
                "degree": graph.degree(node),
                "weighted_degree": weighted_degree[node],
                "betweenness_centrality": betweenness[node],
                "community": community_for_node[node],
                "network_filter": f"q<{Q_THRESHOLD} and |Phi|>={PHI_THRESHOLD}",
            }
        )
    node_rows.sort(key=lambda row: (-float(row["weighted_degree"]), -int(row["degree"]), str(row["node"])))
    edge_rows.sort(key=lambda row: (-float(row["absolute_Phi"]), float(row["q_value"]), str(row["source"]), str(row["target"])))
    atomic_write_tsv(ROOT / "06_tables/functional_network_nodes.tsv", node_rows, list(node_rows[0]))
    atomic_write_tsv(ROOT / "06_tables/functional_network_edges.tsv", edge_rows, list(edge_rows[0]))

    community_rows = []
    for index, nodes in enumerate(ordered_communities, start=1):
        community = f"COMMUNITY_{index:03d}"
        subgraph = graph.subgraph(nodes)
        categories = Counter(node_metadata[node]["functional_category"] for node in nodes)
        internal_weights = [abs(data["signed_weight"]) for _, _, data in subgraph.edges(data=True)]
        top_nodes = sorted(nodes, key=lambda node: (-weighted_degree[node], node))[:10]
        community_rows.append(
            {
                "community": community,
                "number_nodes": len(nodes),
                "number_internal_edges": subgraph.number_of_edges(),
                "functional_category_counts": ";".join(f"{key}:{value}" for key, value in sorted(categories.items())),
                "top_nodes_by_weighted_degree": ";".join(top_nodes),
                "mean_absolute_Phi_internal_edges": sum(internal_weights) / len(internal_weights) if internal_weights else "",
                "positive_internal_edges": sum(data["signed_weight"] > 0 for _, _, data in subgraph.edges(data=True)),
                "negative_internal_edges": sum(data["signed_weight"] < 0 for _, _, data in subgraph.edges(data=True)),
                "community_method": "NetworkX Louvain on absolute Phi edge weights; signs retained for interpretation",
            }
        )
    atomic_write_tsv(ROOT / "06_tables/functional_network_communities.tsv", community_rows, list(community_rows[0]))
    print(
        f"Functional network: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} strong edges, "
        f"{len(ordered_communities)} Louvain communities"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
