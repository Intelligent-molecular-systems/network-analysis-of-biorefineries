import networkx as nx
from matplotlib import pyplot as plt
from pyvis.network import Network

import utils
from src.visualise import plot_components, plot_components_size_distribution


def analyse(df, file_name="graph_fragmentation"):
    """
    Build the network for analysing the graph fragmentation.
    Calculate how many components are there. Sort them decreasingly by size and plot their distribution.
    Calculate the proportion of the largest component size to all nodes.
    :param df: DataFrame with appropriate columns
    :param file_name: File name to save the visualisation
    :return: nothing
    """
    # Create a network for analysing fragmentation
    graph, nodes_dict = utils.create_weighted_undirected_graph(df, use_chemical_names=True)

    print("There are", graph.number_of_edges(), "edges in this graph and", graph.number_of_nodes(), "nodes.")
    print("This graph is fully connected:", nx.is_connected(graph))
    print()
    print("The number of components:", nx.number_connected_components(graph))

    # Visualise the network
    net = Network(directed=True, notebook=False)
    net.from_nx(graph)

    html = net.generate_html()
    with open("../visualisations/" + file_name + ".html", mode="w", encoding="utf-8") as fp:
        fp.write(html)
        print("Visualisation saved in the: visualisations/" + file_name + ".html")

    # Get the sizes of the components in the decreasing order
    component_sizes = [len(c) for c in sorted(nx.connected_components(graph), key=len, reverse=True)]

    figure, axis = plt.subplots(2, 1)

    # Plot the components
    plot_components(component_sizes, axis, 0)

    # Plot the number of components with each size
    plot_components_size_distribution(component_sizes, axis, 1)

    figure.tight_layout()
    plt.show()

    percent_of_nodes_in_biggest_component = round(component_sizes[0] * 100 / graph.number_of_nodes(), 2)
    print("The largest components contains:", percent_of_nodes_in_biggest_component, "% of all the nodes")
