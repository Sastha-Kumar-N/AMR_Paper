"""Louvain, preselected Fisher tests, and typed NetworkX knowledge-graph export."""
from __future__ import annotations

import json
import logging
import random
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact, false_discovery_control
from tqdm import tqdm

LOG = logging.getLogger('panaroo_figures')


def candidate_pairs(ids: list[int], count: int, seed: int) -> list[tuple[int, int]]:
    """Uniformly sample unique unordered pairs without inspecting joint presence."""
    ids = sorted(ids)
    n = len(ids)
    population = n * (n - 1) // 2
    count = min(count, population)
    # Triangular indexing avoids materializing the potentially billion-pair universe.
    ranks = sorted(random.Random(seed).sample(range(population), count))
    starts = np.arange(n, dtype=np.int64)
    starts = starts * (2 * n - starts - 1) // 2
    pairs = []
    for rank in ranks:
        i = int(np.searchsorted(starts, rank, side='right') - 1)
        j = i + 1 + int(rank - starts[i])
        pairs.append((ids[i], ids[j]))
    return pairs


def contingency(a: int, b: int, total: int) -> tuple[int, int, int, int]:
    """Return both present, A only, B only, neither from genome-presence bitsets."""
    both = (a & b).bit_count()
    a_only = a.bit_count() - both
    b_only = b.bit_count() - both
    neither = total - both - a_only - b_only
    if neither < 0:
        raise ValueError('Membership exceeds cohort denominator')
    return both, a_only, b_only, neither


def presence_bitsets(context: Any, isolate_names: list[str]) -> dict[int, int]:
    """Persist binary CSR presence and its exact node/genome ordering."""
    ids = sorted(context.graph)
    members = sorted({m for n in ids for m in context.graph.nodes[n]['genome_members']})
    if len(members) > context.total:
        raise ValueError('Too many genome IDs')
    if isolate_names:
        if any(not m.isdigit() or not 0 <= int(m) < len(isolate_names) for m in members):
            raise ValueError('Genome member IDs cannot be mapped to isolateNames indices')
        genome_ids = [str(i) for i in range(context.total)]
        names = isolate_names
    else:
        if len(members) != context.total:
            raise ValueError('Provide graph isolateNames to map genomes absent from all clusters')
        genome_ids = members
        names = [''] * context.total
    mapping = {value: i for i, value in enumerate(genome_ids)}
    indptr = [0]
    indices = []
    bitsets = {}
    for n in ids:
        present = sorted(mapping[m] for m in context.graph.nodes[n]['genome_members'])
        if len(present) != context.graph.nodes[n]['size']:
            raise ValueError(f'Incomplete genome membership for node {n}')
        bitset = 0
        for index in present:
            bitset |= 1 << index
        bitsets[n] = bitset
        indices.extend(present)
        indptr.append(len(indices))
    matrix = sparse.csr_matrix((np.ones(len(indices), dtype=np.uint8),
                               np.asarray(indices, dtype=np.int32), np.asarray(indptr, dtype=np.int64)),
                              shape=(len(ids), context.total))
    sparse.save_npz(context.outdir / 'gene_presence_absence.npz', matrix)
    pd.DataFrame({'row_index': range(len(ids)), 'node_id': ids}).to_csv(context.outdir/'presence_node_order.tsv', sep='\t', index=False)
    pd.DataFrame({'column_index': range(context.total), 'genome_id': genome_ids, 'isolate_name': names}).to_csv(context.outdir/'presence_genome_order.tsv', sep='\t', index=False)
    return bitsets


def detect_communities(context: Any) -> tuple[list[set[int]], pd.DataFrame, float]:
    """Use full adjacency support as Louvain weight; canonicalize community IDs."""
    from pangenome_graph_figures import CATEGORIES
    config = context.config['statistics']
    LOG.info('Detecting weighted Louvain communities on the complete graph')
    if context.graph.number_of_edges():
        communities = nx.community.louvain_communities(context.graph, weight='size',
            resolution=config['louvain_resolution'], seed=context.config['seed'])
        modularity = nx.community.modularity(context.graph, communities, weight='size', resolution=config['louvain_resolution'])
    else:
        communities = [{n} for n in context.graph]
        modularity = 0.0
    communities = sorted(communities, key=lambda nodes: (-len(nodes), min(nodes)))
    assignment = {n: cid for cid, nodes in enumerate(communities) for n in nodes}
    context.table['community_id'] = [assignment[n] for n in context.table.node_id]
    rows = []
    for cid, nodes in enumerate(communities):
        frame = context.table.loc[sorted(nodes)]
        subgraph = context.graph.subgraph(nodes)
        rows.append({'community_id': cid, 'node_count': len(nodes), 'internal_edges': subgraph.number_of_edges(),
                     'internal_support_sum': sum(a['size'] for _, _, a in subgraph.edges(data=True)),
                     **{f'{c.replace("-", "_")}_count': int((frame.category == c).sum()) for c in CATEGORIES},
                     'evidence_type': 'graph_detected'})
    table = pd.DataFrame(rows)
    table.to_csv(context.outdir/'community_table.tsv', sep='\t', index=False)
    context.table[['node_id', 'community_id', 'component_id', 'evidence_type']].to_csv(context.outdir/'community_assignments.tsv', sep='\t', index=False)
    return communities, table, float(modularity)


def test_pairs(context: Any, bitsets: dict[int, int]) -> pd.DataFrame:
    """One Fisher test per prespecified pair; BH correction across the whole family."""
    config = context.config['statistics']
    eligible = [n for n in sorted(context.graph) if 0 < context.graph.nodes[n]['size'] < context.total]
    pairs = candidate_pairs(eligible, config['candidate_pairs'], context.config['seed'])
    if len(pairs) < config['candidate_pairs']:
        context.warn(f'Only {len(pairs)} variable-cluster pairs exist; testing all instead of {config["candidate_pairs"]}.')
    pd.DataFrame(pairs, columns=['source_node_id', 'target_node_id']).to_csv(context.outdir/'candidate_pairs.tsv', sep='\t', index=False)
    rows = []
    for u, v in tqdm(pairs, desc='Fisher exact tests', mininterval=2):
        both, a_only, b_only, neither = contingency(bitsets[u], bitsets[v], context.total)
        odds, p = fisher_exact([[both, a_only], [b_only, neither]], alternative=config['fisher_alternative'])
        expected = (both+a_only)*(both+b_only)/context.total
        rows.append({'source_node_id': u, 'target_node_id': v, 'both_present': both,
                     'a_only': a_only, 'b_only': b_only, 'neither_present': neither,
                     'expected_both_independent': expected, 'odds_ratio': odds, 'p_value': p,
                     'source_community_id': int(context.table.loc[u, 'community_id']),
                     'target_community_id': int(context.table.loc[v, 'community_id']),
                     'is_panaroo_adjacent': context.graph.has_edge(u, v), 'evidence_type': 'graph_detected'})
    columns = ['source_node_id','target_node_id','both_present','a_only','b_only','neither_present',
               'expected_both_independent','odds_ratio','p_value','source_community_id','target_community_id',
               'is_panaroo_adjacent','evidence_type']
    table = pd.DataFrame(rows, columns=columns)
    table['q_value_bh'] = false_discovery_control(table.p_value.to_numpy(), method='bh') if len(table) else np.array([])
    table['significant_positive'] = (table.q_value_bh <= config['fdr']) & (table.odds_ratio > 1)
    table.to_csv(context.outdir/'cooccurrence_all_pairs.tsv', sep='\t', index=False)
    table[table.significant_positive].sort_values(['q_value_bh','source_node_id','target_node_id']).to_csv(context.outdir/'cooccurrence_significant.tsv', sep='\t', index=False)
    return table


def export_knowledge_graph(context: Any, communities: pd.DataFrame, pairs: pd.DataFrame) -> dict:
    """Integrate attributed gene/community nodes and typed evidence edges in a MultiGraph."""
    graph = nx.MultiGraph(input_sha256=context.fingerprint, evidence_scope='graph and annotation; no phenotype validation',
                          fisher_alternative=context.config['statistics']['fisher_alternative'],
                          fdr_threshold=context.config['statistics']['fdr'])
    for row in context.table.to_dict('records'):
        attrs = {k: v for k, v in row.items() if v is not None and not (isinstance(v, float) and np.isnan(v))}
        graph.add_node(f'gene:{row["node_id"]}', entity_type='gene_cluster', **attrs)
    for row in communities.to_dict('records'):
        graph.add_node(f'community:{row["community_id"]}', entity_type='louvain_community', **row)
    edge_rows = []
    def add_edge(u: str, v: str, relationship: str, **attrs: Any) -> None:
        record = dict(relationship=relationship, evidence_type='graph_detected', **attrs)
        graph.add_edge(u, v, key=relationship, **record)
        edge_rows.append({'source': u, 'target': v, **record})
    for u, v, attrs in context.graph.edges(data=True):
        add_edge(f'gene:{u}', f'gene:{v}', 'panaroo_adjacency', size=attrs['size'], support=attrs['size']/context.total)
    for row in context.table.itertuples():
        add_edge(f'gene:{row.node_id}', f'community:{row.community_id}', 'member_of_community')
    for row in pairs[pairs.significant_positive].itertuples():
        add_edge(f'gene:{row.source_node_id}', f'gene:{row.target_node_id}', 'significant_cooccurrence',
                 odds_ratio=float(row.odds_ratio), p_value=float(row.p_value), q_value_bh=float(row.q_value_bh),
                 both_present=int(row.both_present))
    nx.write_graphml(graph, context.outdir/'pangenome_knowledge_graph.graphml.gz')
    pd.DataFrame([{'entity_id': n, **a} for n, a in graph.nodes(data=True)]).to_csv(context.outdir/'knowledge_graph_nodes.tsv', sep='\t', index=False)
    pd.DataFrame(edge_rows).to_csv(context.outdir/'knowledge_graph_edges.tsv', sep='\t', index=False)
    return {'nodes': graph.number_of_nodes(), 'edges': graph.number_of_edges(), 'type': 'NetworkX MultiGraph',
            'file': 'pangenome_knowledge_graph.graphml.gz'}


def plot_statistics(context: Any, communities: pd.DataFrame, pairs: pd.DataFrame) -> None:
    """Export community composition and co-occurrence evidence, with exact source rows."""
    from pangenome_graph_figures import CATEGORIES, COLORS, new_figure, save_figure, category_handles, export_source
    fig = new_figure(context, 95)
    ax = fig.add_axes([.10,.26,.85,.58])
    top = communities.head(12)
    base = np.zeros(len(top))
    for category in CATEGORIES:
        values = top[f'{category.replace("-", "_")}_count'].to_numpy()
        ax.bar(np.arange(len(top)), values, bottom=base, color=COLORS[category], width=.72)
        base += values
    ax.set_xticks(np.arange(len(top)), [str(v) for v in top.community_id])
    ax.set_xlabel('Louvain community ID (12 largest communities)')
    ax.set_ylabel('Gene clusters')
    fig.suptitle('Weighted Louvain community composition', fontsize=11, y=.97)
    fig.legend(handles=category_handles(context), loc='lower center', ncol=2, frameon=False, bbox_to_anchor=(.5,.03))
    export_source('panel_g_communities', context, list(context.graph),
                  node_metadata={n: {'community_id': int(context.table.loc[n,'community_id'])} for n in context.graph},
                  extras=[{'record_type':'community_summary', **row} for row in communities.to_dict('records')])
    save_figure(fig, 'fig_g_louvain_communities', context)
    fig = new_figure(context, 100)
    ax = fig.add_axes([.12,.23,.82,.62])
    for significant, color, label in [(False, '#AAAAAA', 'Not significant'), (True, '#0072B2', 'Positive; BH q ≤ 0.05')]:
        subset = pairs[pairs.significant_positive == significant]
        observed = np.log2((subset.both_present + .5)/(subset.expected_both_independent + .5))
        score = -np.log10(np.maximum(subset.q_value_bh.to_numpy(dtype=float), np.finfo(float).tiny))
        ax.scatter(observed, score, s=4, c=color, alpha=.55, linewidths=0, rasterized=True, label=label.replace('0.05',f'{context.config["statistics"]["fdr"]:g}'))
    ax.axhline(-np.log10(context.config['statistics']['fdr']), color='#555555', ls='--', lw=.7)
    ax.set_xlabel('log2[(observed + 0.5) / (expected + 0.5)]')
    ax.set_ylabel('-log10(BH-adjusted p-value)')
    ax.legend(loc='upper left', frameon=False)
    fig.suptitle(f'Gene co-occurrence · {len(pairs):,} tested pairs', fontsize=11, y=.97)
    fig.text(.5,.06,'One-sided Fisher tests; all sampled variable-cluster pairs enter BH correction.\nThe 0.5 offset is for visualization only. Associations are unadjusted for population structure.',ha='center',fontsize=7)
    source = pairs.copy()
    source['record_type'] = 'tested_pair'
    source['display_log2_observed_expected'] = np.log2((source.both_present+.5)/(source.expected_both_independent+.5))
    source['display_minus_log10_q'] = -np.log10(np.maximum(source.q_value_bh.to_numpy(dtype=float),np.finfo(float).tiny))
    source.to_csv(context.outdir/'figure_source_data'/'panel_h_cooccurrence.csv', index=False)
    save_figure(fig, 'fig_h_gene_cooccurrence', context)


def run_statistics(context: Any, isolate_names: list[str]) -> dict:
    """Run and export the specified analyses without using significance to select tests."""
    options = context.config['statistics']
    if options['candidate_pairs'] < 1 or not 0 < options['fdr'] < 1 or options['louvain_resolution'] <= 0:
        raise ValueError('Invalid statistical configuration')
    if options['fisher_alternative'] != 'greater':
        raise ValueError('Positive co-occurrence analysis requires fisher_alternative: greater')
    bitsets = presence_bitsets(context, isolate_names)
    partition, communities, modularity = detect_communities(context)
    pairs = test_pairs(context, bitsets)
    knowledge = export_knowledge_graph(context, communities, pairs)
    plot_statistics(context, communities, pairs)
    eligible = sum(0 < a['size'] < context.total for _, a in context.graph.nodes(data=True))
    summary = {'communities': len(partition), 'weighted_modularity': modularity,
               'louvain_weight': 'size', 'louvain_resolution': options['louvain_resolution'],
               'seed': context.config['seed'], 'candidate_selection': 'uniform sample without replacement of unordered pairs among variable gene clusters; no joint-presence screening',
               'eligible_gene_clusters': eligible, 'eligible_pair_universe': eligible*(eligible-1)//2,
               'requested_pairs': options['candidate_pairs'], 'tested_pairs': len(pairs),
               'fisher_alternative': 'greater', 'fdr_method': 'Benjamini-Hochberg', 'fdr_threshold': options['fdr'],
               'significant_pairs': int(pairs.significant_positive.sum()), 'knowledge_graph': knowledge,
               'population_structure_adjusted': False}
    (context.outdir/'statistics_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    context.warn('Fisher co-occurrence associations are unadjusted for population structure; BH significance is restricted to the sampled candidate family and does not establish functional interaction.')
    return summary


def write_statistical_methods(context: Any) -> None:
    """Append actual statistical design, results, and evidence limitations to Methods."""
    s = context.statistics
    text = (f"\n\nCommunity structure was detected using NetworkX's Louvain modularity optimization algorithm on the full Panaroo adjacency graph, "
            f"weighted by edge genome support (size), with resolution {s['louvain_resolution']} and seed {s['seed']}. "
            f"This yielded {s['communities']:,} communities and weighted modularity {s['weighted_modularity']:.6f}. "
            f"Binary genome presence was reconstructed from validated node membership lists. A uniform sample without replacement of {s['tested_pairs']:,} unordered "
            f"candidate pairs was selected from {s['eligible_gene_clusters']:,} variable clusters ({s['eligible_pair_universe']:,} possible pairs), using seed {s['seed']}, "
            "without screening joint-presence counts, adjacency or community membership. Invariant clusters were excluded. "
            "One-sided Fisher exact tests (alternative greater) used 2×2 counts of both present, A only, B only and neither present across all genomes. "
            f"Benjamini-Hochberg correction was applied jointly to all {s['tested_pairs']:,} tested p-values; {s['significant_pairs']:,} positive associations "
            f"met BH-adjusted p ≤ {s['fdr_threshold']}. These exploratory associations are unadjusted for population structure or phylogenetic relatedness, "
            "and the reported FDR threshold applies to the sampled test family, not all possible pairs. BH assumptions concerning dependence between tests were not established. "
            "Node, adjacency, community and significant co-occurrence attributes were integrated into a typed NetworkX MultiGraph, exported as compressed GraphML "
            "and node/edge tables. Community membership, Panaroo adjacency and statistical co-occurrence are distinct relationships. "
            "Names and descriptions remain annotation-inferred, and no phenotypic or causal functional claims are made.\n")
    with (context.outdir/'methods.md').open('a') as handle:
        handle.write(text)
    with (context.outdir/'figure_legends.md').open('a') as handle:
        handle.write('\n\n**g. Louvain communities.** Category composition of the 12 largest communities in the complete adjacency graph; weights are genome support. All assignments and community summaries are exported.\n\n**h. Gene co-occurrence.** Each point is one tested pair. Blue points meet the configured BH threshold and have odds ratio >1. All candidates enter the correction family. The horizontal coordinate uses a 0.5 visualization offset; Fisher tests use unmodified integer counts. Values of q=0 from numerical underflow are plotted at the smallest positive floating-point number. These associations are not adjusted for population structure.\n')
