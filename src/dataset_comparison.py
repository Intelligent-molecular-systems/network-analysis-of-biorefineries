from pyvis.network import Network

from src import utils
import pandas


def analyse(df1, df2):
    """
    Analyse the commonality and uniqueness of molecules between two datasets
    and compute the Jaccard similarity coefficient. Further, create HTML visualisation of the results.

    This method creates weighted, undirected graphs from two DataFrames, identifies the common and unique molecules
    between the two sets, checks the presence of important molecules (THESE VALUES ARE HARDCODED) in the second set,
    and calculates the Jaccard similarity coefficient. It also calls a helper method to visualize the combined datasets.

    :param df1: First DataFrame with appropriate columns for graph creation.
    :param df2: Second DataFrame with appropriate columns for graph creation.
    :return: None
    """
    _, nodes_dict1 = utils.create_weighted_undirected_graph(df1, use_chemical_names=True)
    _, nodes_dict2 = utils.create_weighted_undirected_graph(df2, use_chemical_names=True)

    nodes1 = set(nodes_dict1.keys())
    nodes2 = set(nodes_dict2.keys())
    intersection = nodes1.intersection(nodes2)

    print("\nThere are", len(intersection), "common molecules,"
                                            " that is", len(intersection) * 100 / len(nodes2), "% of the 2nd network")
    print("Common molecules are:", intersection, "\n")

    print("Checking if important molecules are in the other set")
    important_molecules = [
        "2-methoxy-phenol",
        "formic acid",
        "methanol",
        "syringic aldehyde",
        "carbon dioxide",
        "1-(4-hydroxy-3,5-dimethoxyphenyl)-2-(2'-methoxyphenoxy)-1,3-propanediol",
        "5-hydroxymethyl-2-furfuraldehyde",
        "methanol",
        "levulinic acid",
        "furfural",
        "vanillin",
    ]
    for molecule in important_molecules:
        if molecule not in intersection:
            print(molecule, "doesn't appear in the other dataset")

    union = nodes1.union(nodes2)
    jaccard_similarity_coefficient = float(len(intersection)) / float(len(union))
    print("\nJaccard similarity coefficient is equal to", jaccard_similarity_coefficient)

    create_visualisation_of_combined_datasets(df1, df2, nodes_dict1, nodes_dict2)


def create_visualisation_of_combined_datasets(df1, df2, nodes_dict1, nodes_dict2):
    """
    Creates an HTML visualization of the combined datasets, highlighting nodes from each dataset and their intersection.

    This helper method combines the two DataFrames, creates a weighted, undirected graph, and visualizes it using the
    Pyvis library. Nodes from the first dataset are colored teal, nodes from the second dataset are colored pink,
    and the intersection nodes are colored lavender. The visualization is saved as an HTML file.

    :param df1: First DataFrame with appropriate columns for graph creation.
    :param df2: Second DataFrame with appropriate columns for graph creation.
    :param intersection: Set of nodes that are common between the two datasets.
    :param nodes1: Set of nodes from the first dataset.
    :param nodes2: Set of nodes from the second dataset.
    :return: None
    """
    nodes1 = set(nodes_dict1.keys())
    nodes2 = set(nodes_dict2.keys())
    intersection = nodes1.intersection(nodes2)

    df = pandas.concat([df1, df2], axis="rows")
    graph, _ = utils.create_weighted_undirected_graph(df, use_chemical_names=True)
    net = Network(directed=True, notebook=False)
    net.from_nx(graph)
    options = {"physics": {"barnesHut": {"springLength": 10, "springConstant": 0.015}, "minVelocity": 0.75}}
    net.options = options
    # Color nodes from 1st dataset as TEAL (#1abc9c) / border #648FC9
    for node in nodes1:
        net.get_node(node)["color"] = {"background": "#1abc9c", "border": "#1abc9c"}
    # Color nodes from 2nd dataset as PINK (#ff6b81)
    for node in nodes2:
        net.get_node(node)["color"] = {"background": "#ff6b81", "border": "#ff6b81"}
    # Color the intersection of the datasets as LAVENDER (#a29bfe)
    for node in intersection:
        net.get_node(node)["color"] = {"background": "#a29bfe", "border": "#a29bfe"}
    html = net.generate_html()
    with open("../visualisations/graph_combined.html", mode="w", encoding="utf-8") as fp:
        fp.write(html)
    print("Visualisation saved in the:", "visualisations/graph_combined.html")
