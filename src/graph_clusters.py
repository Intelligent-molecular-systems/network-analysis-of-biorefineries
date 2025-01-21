import networkx as nx
import itertools

import numpy as np
from matplotlib import pyplot as plt
from networkx import k_core, core_number
from pyvis.network import Network

from src import utils

COLOURS = [
    "#ff8c00",
    "#ffff00",
    "#00ff00",
    "#008000",
    "#00ffff",
    "#ff00ff",
    "#2f4f4f",
    "#ff69b4",
    "#7f0000",
    "#00008b",
    "#1e90ff",
    "#ffdead",
    "gray",
]


def analyse(df):
    kcore_clustering(df)

    betwenness_clustering(df, 100)


def kcore_clustering(df):
    """
    Performs k-core clustering on the given DataFrame and saves visualisations for the top 4 k-cores.

    This function creates an undirected graph from the given DataFrame, calculates the core numbers for each node,
    and visualizes the top 4 k-cores based on the highest core numbers. The visualisations are saved as HTML files.

    :param df: DataFrame containing the data to create the graph.
    :return: None
    """
    # Create a network
    graph, nodes_dict = utils.create_weighted_undirected_graph(df, use_chemical_names=True, allow_multiple_edges=False)
    graph.remove_edges_from(nx.selfloop_edges(graph))

    number = core_number(graph)
    possible_cores = list(utils.count_by(1, list(number.items()))[0])
    possible_cores.sort()
    print("Molecules with these k numbers for kcore were found:", possible_cores)
    print("Visualising top 4 kcores.")
    for i in range(max(possible_cores) - 3, max(possible_cores) + 1):
        result = k_core(graph, k=i)
        # Visualise the network
        net = Network(directed=False, notebook=False)
        net.from_nx(result)

        html = net.generate_html()
        with open("../visualisations/graph_clusters_" + str(i) + "core.html", mode="w", encoding="utf-8") as fp:
            fp.write(html)
            print("Visualisation saved in the: visualisations/graph_clusters_" + str(i) + "core.html")
    print()


def betwenness_clustering(df, k=8, visualise=True):
    """
    Method to generate clusters using girvan-newman method based on betweenness centrality. Creates also
    file visualisation for it.
    Note, this method takes about 5 minutes to run since betweenness centrality has to be
    recalculated after every edge deletion. This method is able to visualise at most 13 clusters (limited by the
    number of colors)
    :param df: DataFrame with appropriate columns for graph creation.
    :param k: Number of graph splits before stopping.
    :param visualise: Whether to create a visualisation of the best clustering.
    :return:
    """
    # Extract the largest connected component
    graph_lcc, _ = utils.create_weighted_undirected_graph(df, use_chemical_names=True, allow_multiple_edges=False)
    largest_connected_component = graph_lcc.subgraph(max(nx.connected_components(graph_lcc), key=len)).copy()

    # Create a network using the largest connected component
    graph, nodes = utils.create_weighted_directed_graph(df, use_chemical_names=True, allow_multiple_edges=False)
    graph = graph.subgraph(largest_connected_component)

    net = Network(directed=True, notebook=False, select_menu=True, cdn_resources="remote")
    net.from_nx(graph)

    # Compute clusters based on girvan-neuman and betweenness centrality
    comp = nx.community.girvan_newman(graph)

    print("Calculating betweenness centrality (it can take a few minutes).")
    cluster_performance = []
    best_performance = 0
    best_partition = []
    for communities in itertools.islice(comp, k):
        components = sorted(tuple(sorted(c) for c in communities), key=lambda x: -len(x))
        performance = calculate_cluster_performance(components, graph)
        # print("Cluster rating for", len(components), "components:", performance)
        cluster_performance.append(performance)
        if performance > best_performance:
            best_performance = performance
            best_partition = components

    print("Best performance found:", best_performance, "for", len(best_partition), "clusters.")
    y_points = np.array(cluster_performance)
    x_points = np.arange(2, len(y_points) + 2)

    # Set font sizes for the paper
    plt.rc("axes", labelsize=16)
    plt.rc("legend", fontsize=13)
    plt.rc("xtick", labelsize=12)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=12)  # fontsize of the tick labels

    # plt.title("Clustering score by the number of clusters")
    plt.xlabel("Number of clusters")
    plt.ylabel("Cluster score")
    plt.plot(x_points, y_points)
    plt.show()

    if visualise:
        create_html_visualisation(best_partition, k, net)


def create_html_visualisation(best_partition, k, net):
    """
    Creates an HTML visualization of the best partition by coloring each component with a different color.

    This helper method colors the nodes of the best partition found by the clustering algorithm and saves the
    visualization as an HTML file.

    :param best_partition: List of components (each component is a list of nodes) representing the best partition.
    :param k: Number of clusters to visualize.
    :param net: Pyvis Network object used for visualization.
    :return: None
    """
    # Color the graph of the best partition
    component_number = 0
    for component in best_partition:
        for node_id in component:
            if component_number >= k:
                net.get_node(node_id)["color"] = COLOURS[-1]
            else:
                net.get_node(node_id)["color"] = COLOURS[component_number]
        component_number += 1
    best_num = len(best_partition)
    html = net.generate_html()
    with open("../visualisations/graph_clusters_betweenness" + str(best_num) + ".html", mode="w", encoding="utf-8") as fp:
        fp.write(html)
        print("Visualisation saved in the: visualisations/graph_clusters_betweenness" + str(best_num) + ".html")


def calculate_cluster_performance(communities, graph):
    """
    Calculates the clustering performance for the given communities in the graph.

    This function calculates the intra-community (in the same community) and inter-community
    (between different communities) edge counts to determine a clustering performance metric.
    If there is only one community, it returns -1.

    :param communities: List of communities (each community is a list of nodes).
    :param graph: NetworkX graph object representing the network.
    :return: A float value representing the clustering performance. Returns -1 if there is only one community.
    """
    if len(communities) == 1:
        return -1
    node_community = {}
    intra_edges = {}
    inter_edges = {}
    for i, community in enumerate(communities):
        intra_edges[i] = 0
        inter_edges[i] = 0
        for node in community:
            node_community[node] = i

    for e in graph.edges():
        community1 = node_community[e[0]]
        community2 = node_community[e[1]]
        if community1 == community2:
            intra_edges[community1] += 1
        else:
            inter_edges[community1] += 1
            inter_edges[community2] += 1

    result = 0
    for i, community in enumerate(communities):
        cluster_rating = 2.0 * intra_edges[i] / inter_edges[i]
        result += cluster_rating
    result /= len(communities)

    return result
