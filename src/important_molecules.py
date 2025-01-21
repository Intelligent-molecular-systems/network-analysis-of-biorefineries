import networkx as nx
import matplotlib.pyplot as plt
from networkx import set_edge_attributes
from pyvis.network import Network

import utils


def analyse(df, same_width_edges=False):
    """
    Creates a weighted, directed graph and performs analysis to display nodes
    with the highest degree and betweenness centrality.
    :param df: DataFrame with appropriate columns
    :param same_width_edges: boolean to set the edges to the same width. Defaults to False
    :return: None
    """
    # Create a directed network to analyse molecules
    graph, nodes = utils.create_weighted_directed_graph(df, use_chemical_names=True, allow_multiple_edges=False)
    print("There are", graph.number_of_edges(), "edges in this graph and", graph.number_of_nodes(), "nodes.")

    if same_width_edges:
        set_edge_attributes(graph, 1, "weight")

    k = 6
    num_to_plot = 5
    betweenness_centrality_analysis(graph, k, num_to_plot)

    inout_degree_analysis(graph, k)


def betweenness_centrality_analysis(graph, k, num_to_plot):
    """
    Performs betweenness centrality analysis on the graph and visualizes the top nodes.

    :param graph: NetworkX graph on which to perform the analysis.
    :param k: Number of top nodes to display in the analysis.
    :param num_to_plot: Number of top nodes to plot in the bar chart.
    :return: None
    """
    # Calculate centrality
    node_centrality = nx.betweenness_centrality(graph, weight="weight", normalized=True)
    node_centrality = sorted(node_centrality.items(), key=lambda node: node[1], reverse=True)

    # Visualise the network
    net = Network(directed=True, notebook=False, select_menu=True, filter_menu=True)
    net.from_nx(graph)

    print("\nTop", k, "chemicals according to their betweenness centrality:")
    most_central_chemicals = get_best_chemicals(
        scores=node_centrality,
        net=net,
        scale=10000,
        file_name="important_molecules_centrality",
        k=k,
        num_to_plot=num_to_plot,
    )
    # Create plots of most important molecules
    plt.figure(1)
    plt.title("Most important molecules according to betweenness centrality")
    plt.ylabel("Betweenness centrality score")
    plt.bar(most_central_chemicals[0], most_central_chemicals[1])
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


def inout_degree_analysis(graph, k):
    """
    Performs total degree analysis on the graph and visualizes the top nodes.

    :param graph: NetworkX graph on which to perform the analysis.
    :param k: Number of top nodes to display in the analysis.
    :return: None
    """
    # Restore the weight of the edges (because pyvis doesn't recognise if only 'width' is provided)
    for e in graph.edges(data=True):
        e[2]["weight"] = e[2]["width"]

    # Visualise the network
    net = Network(directed=True, notebook=False, select_menu=True, filter_menu=True)
    net.from_nx(graph)

    # Calculate degree
    degree_value = sorted(graph.degree, key=lambda node: node[1], reverse=True)
    print("\nTop", k, "chemicals according to their degree and their score:")
    highest_degree_chemicals = get_best_chemicals(
        scores=degree_value, net=net, file_name="important_molecules_degree", scale=1.5, k=k, num_to_plot=k
    )

    plt.figure(2)
    plt.title("Most important molecules according to degree")
    plt.ylabel("Total Degree")
    plt.bar(highest_degree_chemicals[0], highest_degree_chemicals[1])
    plt.xticks(rotation=60)
    plt.tight_layout()
    plt.show()


def get_best_chemicals(scores, net, scale, file_name, k, num_to_plot):
    """
    Identifies and visualizes the top k chemicals based on the provided scores.

    :param scores: List of tuples containing node names and their scores.
    :param net: Pyvis Network object for visualization.
    :param scale: Scaling factor for node size in the visualization.
    :param file_name: Name of the file to save the visualization.
    :param k: Number of top nodes to display in the analysis.
    :param num_to_plot: Number of top nodes to plot in the bar chart.
    :return: List containing two lists - (0) top k chemical names and (1) their scores.
    """
    best_chemicals = [[], []]
    for i in range(k):
        print("{0:80}\t{1}".format(scores[i][0], scores[i][1]))
        if i < num_to_plot:
            best_chemicals[0].append(scores[i][0])
            best_chemicals[1].append(scores[i][1])

        # Set options for visualisation
        net.get_node(scores[i][0])["shape"] = "star"
        net.get_node(scores[i][0])["size"] = scale * scores[i][1]
        net.get_node(scores[i][0])["color"] = {"background": "red", "border": "#648FC9"}

    # Save visualisation in the file
    html = net.generate_html()
    with open("../visualisations/" + file_name + ".html", mode="w", encoding="utf-8") as fp:
        fp.write(html)
    print("Visualisation saved in the: visualisations/" + file_name + ".html")

    return best_chemicals
