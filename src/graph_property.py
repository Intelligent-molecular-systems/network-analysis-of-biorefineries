from collections import Counter

import networkx as nx
from matplotlib import pyplot as plt
from numpy import inf

from src import utils
from networkx.algorithms.smallworld import omega

from src.visualise import plot_comparison_with_NOC_avg_shortest_path_length


def analyse(df):
    """
    Creates graph to analyse if it has small-world property. Calculates average shortest path length,
    clustering coefficient and omega metric. ASPL is calculated using weighted, directed graph, while CC and omega
    are calculated using the largest connected component in the weighted, undirected graph.
    :param df: DataFrame with appropriate columns
    :return: None
    """
    # Create a directed network to analyse molecules
    graph, nodes = utils.create_weighted_directed_graph(df, use_chemical_names=True, allow_multiple_edges=False)

    # Calculate average shortest path length
    x, y, avg = average_shortest_path(graph)
    print("Average shortest path length =", avg)

    # Plot average shortest path length using log scale plot
    plt.plot(x, y, "o-b")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylim(top=1)
    plt.xlim(right=100)
    plt.xlabel("Shortest path length, l")
    plt.ylabel("Probability, P(l)")
    plt.show()

    # Get largest connected component
    graph, nodes = utils.create_weighted_undirected_graph(df, use_chemical_names=True, allow_multiple_edges=False)
    largest_connected_component = graph.subgraph(max(nx.connected_components(graph), key=len)).copy()
    # Calculate clustering coefficient
    clustering_coefficient = nx.average_clustering(largest_connected_component)
    print("Average clustering coefficient is:", clustering_coefficient)
    # Calculate omega metric
    # With this settings (niter=3 and nrand=3) this part can run about 800 seconds. (Comment out if not needed
    # to speed up the calculations)
    graph_omega = omega(largest_connected_component, niter=3, nrand=3)
    print("Omega of this graph is:", graph_omega)


def average_shortest_path(graph):
    """
    Calculates the average shortest path length of the given graph and returns the distribution of path lengths.

    This method computes the shortest path lengths between all pairs of nodes using the Floyd-Warshall algorithm.
    It then calculates the average shortest path length and the distribution of path lengths,
    excluding infinite distances.

    :param graph: A NetworkX graph object for which to calculate the average shortest path length.
    :return: A tuple containing:
        - x: A list of path lengths
        - y: A list of normalized counts of the respective path lengths
        - avg: The average shortest path length
    """
    distance = nx.floyd_warshall_numpy(graph)
    distances = []
    max_length = 0
    for i in range(len(graph.nodes)):
        for j in range(len(graph.nodes)):
            if i != j and distance[i][j] != inf:
                distances.append(distance[i][j])
                if distance[i][j] > max_length:
                    max_length = int(distance[i][j])
    avg = 0
    counter = Counter(distances)
    x = []
    y = []
    num = 0
    for i in range(1, max_length):
        x.append(i)
        pom = counter[i]
        y.append(pom)
        num += pom
        avg += pom * i

    for i in range(1, max_length):
        y[i - 1] /= num

    avg /= num
    return x, y, avg


def compare_average_shortest_path_length_undirected_only_largest_component(df):
    """
    Method for plotting average shortest path length compared with the Network of Organic Chemistry (NOC).
    It is calculated on undirected graph using only the largest connected component.
    :param df: DataFrame with the network to plot.
    :return: None
    """
    # Create an undirected network to analyse molecules
    graph, nodes = utils.create_weighted_undirected_graph(df, use_chemical_names=True, allow_multiple_edges=False)

    largest_connected_component = graph.subgraph(max(nx.connected_components(graph), key=len)).copy()
    x, y, avg = average_shortest_path(largest_connected_component)

    plot_comparison_with_NOC_avg_shortest_path_length(x, y)
    plt.tight_layout()
    plt.show()


def compare_average_shortest_path_length_directed(df):
    """
    Method for plotting average shortest path length compared with the Network of Organic Chemistry (NOC).
    This function was used to produce the plot for the final paper.
    It is calculated on the full, directed graph.
    :param df: DataFrame with the network to plot.
    :return: None
    """
    # Set font sizes for the paper
    plt.rc("axes", labelsize=18)
    plt.rc("legend", fontsize=16)
    plt.rc("xtick", labelsize=16)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=16)  # fontsize of the tick labels

    # Create a directed network to analyse molecules
    graph, nodes = utils.create_weighted_directed_graph(df, use_chemical_names=True, allow_multiple_edges=False)

    x, y, avg = average_shortest_path(graph)
    print("Average shortest path length =", avg)

    plot_comparison_with_NOC_avg_shortest_path_length(x, y)
    plt.tight_layout()
    plt.show()


def compare_average_shortest_path_length_directed_only_largest_component(df):
    """
    Method for plotting average shortest path length compared with the Network of Organic Chemistry (NOC).
    It is calculated on directed graph using only the largest connected component.
    :param df: DataFrame with the network to plot.
    :return: None
    """
    # Get largest connected component
    graph, nodes = utils.create_weighted_undirected_graph(df, use_chemical_names=True, allow_multiple_edges=False)

    largest_connected_component = graph.subgraph(max(nx.connected_components(graph), key=len)).copy()

    # Create a directed network to analyse molecules
    graph, nodes = utils.create_weighted_directed_graph(df, use_chemical_names=True, allow_multiple_edges=False)
    graph = graph.subgraph(largest_connected_component)

    x, y, avg = average_shortest_path(graph)
    print("Average shortest path length =", avg)

    plot_comparison_with_NOC_avg_shortest_path_length(x, y)
    plt.tight_layout()
    plt.show()
