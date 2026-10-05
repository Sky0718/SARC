from collections import defaultdict
from pathlib import Path

import networkx as nx

def collect_memberships(mutation_matrix, outlier_matrix):
    samples = mutation_matrix.columns.tolist()
    observed_mutations = mutation_matrix.isin([1])
    mutations = {
        sample: list(
            set(observed_mutations.index[observed_mutations[sample] == True].tolist())
        )
        for (sample) in (samples)
    }
    outliers = []
    for (sample) in (samples):
        for (gene) in (outlier_matrix.index.tolist()):
            if (
                sample in outlier_matrix.columns
                and outlier_matrix[sample][gene] == True
            ):
                outliers.append(str(sample + "_" + gene))
    return samples, mutations, outliers

def index_outliers(outliers):
    by_gene = defaultdict(list)
    for (position, outlier) in (enumerate(outliers)):
        by_gene[outlier.split("_")[1]].append(position)
    return by_gene

def build_graph(sample, mutations, mutation_sets, outliers, by_gene, ppi):
    graph = nx.Graph()
    graph.add_nodes_from(mutations[sample], bipartite = 0)
    graph.add_nodes_from(outliers, bipartite = 1)
    for (mutation) in (mutations[sample]):
        if (mutation in ppi):
            positions = sorted(
                position
                for (gene) in (ppi[mutation])
                if (gene in ppi)
                for (position) in (by_gene.get(gene, ()))
            )
            for (position) in (positions):
                outlier = outliers[position]
                other_sample, outlier_gene = outlier.split("_")[:2]
                if (mutation in mutation_sets[other_sample] and outlier_gene != mutation):
                    graph.add_edge(mutation, outlier)
    mutation_nodes = {
        node for ((node, data)) in (graph.nodes(data = True)) if (data["bipartite"] == 0)
    }
    outlier_nodes = {
        node for ((node, data)) in (graph.nodes(data = True)) if (data["bipartite"] == 1)
    }
    for (outlier) in (outlier_nodes):
        if (int(graph.degree(outlier)) == 0):
            graph.remove_node(outlier)
    for (mutation) in (mutation_nodes):
        if (int(graph.degree(mutation)) == 0):
            graph.remove_node(mutation)
    return graph

def construct_pbns(mutation_matrix, outlier_matrix, ppi, cancer, dataset, network):
    output = Path("graphs") / dataset / (cancer + "_" + network)
    output.mkdir(parents = True, exist_ok = False)
    samples, mutations, outliers = collect_memberships(mutation_matrix, outlier_matrix)
    mutation_sets = {sample: set(genes) for ((sample, genes)) in (mutations.items())}
    by_gene = index_outliers(outliers)
    for (sample) in (samples):
        graph = build_graph(sample, mutations, mutation_sets, outliers, by_gene, ppi)
        nx.write_gml(graph, output / ("BP_" + sample + ".gml"))
