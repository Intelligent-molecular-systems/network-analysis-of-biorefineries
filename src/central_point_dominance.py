import networkx as nx

from src import utils


def analyse(df):
    """
    Create graph from the df and calculate central point dominance.
    :param df: DataFrame to calculate central point dominance for.
    :return: None
    """
    # Create graph
    graph, nodes = utils.create_weighted_directed_graph(df, use_chemical_names=True, allow_multiple_edges=False)
    # Calculate betweenness
    node_centrality = nx.betweenness_centrality(graph)
    node_centrality = sorted(node_centrality.items(), key=lambda node: node[1], reverse=True)
    # Extract the node with the highest centrality
    highest_centrality = node_centrality[0][1]

    # Calculate central point dominance
    result = 0
    for centrality in node_centrality:
        result += highest_centrality - centrality[1]
    result /= len(node_centrality) - 1
    print("Central point dominance =", result)
