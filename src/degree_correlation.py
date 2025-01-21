import networkx as nx

from src import utils
from src.visualise import plot_correlation


def analyse(df, degree_type="out"):
    """
    Creates weighted, directed graph. Calculates the average nearest neighbor degree of nodes with degree k. Plots it.
    :param df: DataFrame to plot the degree correlation for.
    :param degree_type: Whether to plot "in" or "out" degree. Defaults to "out".
    :return: None
    """
    graph, nodes = utils.create_weighted_directed_graph(df, allow_multiple_edges=False)
    avg_deg = nx.average_degree_connectivity(graph, source=degree_type, target=degree_type, weight=None)
    plot_correlation(avg_deg, degree_type)
