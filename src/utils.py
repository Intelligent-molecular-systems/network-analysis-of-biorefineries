import networkx as nx
import numpy as np


def create_weighted_undirected_graph(df, use_chemical_names=False, allow_multiple_edges=False):
    """
    Create a weighted, undirected, network out of the data frame with Reactants and Products.
    Weights on edges are assigned according to the "Number of Reaction Steps" field for each reaction.
    If DataFrame doesn't have the "Number of Reaction Steps" field, it is going to be set to 1 by default.
    :param df: data frame with products and reactants
    :param use_chemical_names: boolean whether to use chemical names or numbers. Defaults to False
    :param allow_multiple_edges: boolean whether to allow for multiple edges between two nodes. Defaults to False
    :return: networkX graph and dictionary (node name -> index)
    """

    # Extract chemicals as nodes from Reactant and Product columns
    nodes = get_distinct_chemicals(df)
    print("There are " + str(len(nodes)) + " chemicals in this graph.")

    # Create edges from Reactants to Products
    edges = get_edges(df, nodes, use_chemical_names, allow_multiple_edges)

    # Build a final graph
    graph = nx.Graph()
    if not use_chemical_names:
        graph.add_nodes_from(nodes.values())
    else:
        graph.add_nodes_from(nodes.keys())
    graph.add_weighted_edges_from(edges)
    return graph, nodes


def create_weighted_directed_graph(df, use_chemical_names=False, allow_multiple_edges=False):
    """
    Create a weighted, directed network out of the dataframe with Reactants and Products.
    Weights on edges are assigned according to the "Number of Reaction Steps" field for each reaction.
    :param df: data frame with reactants and products
    :param use_chemical_names: boolean whether to use chemical names or numbers. Defaults to False
    :param allow_multiple_edges: boolean whether to allow for multiple edges between two nodes. Defaults to False
    :return: networkX directed multi graph and dictionary (node name -> index)
    """
    # Extract chemicals as nodes
    nodes = get_distinct_chemicals(df)
    print("There are " + str(len(nodes)) + " chemicals in this graph.")

    # Create a set of edges
    edges = get_edges(df, nodes, use_chemical_names, allow_multiple_edges)
    # Build a final graph
    graph = nx.MultiDiGraph()
    if not use_chemical_names:
        graph.add_nodes_from(nodes.values())
    else:
        graph.add_nodes_from(nodes.keys())
    graph.add_weighted_edges_from(edges)
    return graph, nodes


def get_distinct_chemicals(df):
    """
    Extract distinct chemicals from Reactants and Products fields, based on the chemical name.
    :param df: dataframe with products and reactants
    :return: Dictionary of pairs (chemical name, index)
    """
    nodes = {}
    idx = 1
    list_of_chemicals = df["Reactant"].tolist() + df["Product"].tolist()
    for chem_string in list_of_chemicals:
        chemicals = chem_string.split(sep="; ")
        for chemical in chemicals:
            if chemical.strip() not in nodes:
                nodes[chemical.strip()] = idx
                idx += 1
    return nodes


def get_edges(df, nodes, use_chemical_names=False, allow_multiple_edges=True):
    """
    Extract all the edges from the data frame as tuples.
    The tuple have the following format (reactant_index, product_index, number_of_steps).
    :param df: dataframe with reactants and products
    :param nodes: dictionary of (chemical name, index)
    :param use_chemical_names: boolean whether to use chemical names or numbers. Defaults to False
    :param allow_multiple_edges: boolean whether to allow for multiple edges between two nodes. Defaults to True
    :return: Set of edges from products to reactants
    """
    if "Number of Reaction Steps" not in df.columns:
        print("Number of Reaction Steps was not given. Feeling in with ones.")
        ones = np.ones(len(df), dtype=int).tolist()
        df["Number of Reaction Steps"] = ones

    added_edges = {}
    edges = set()
    for x, y, z in zip(df["Reactant"], df["Product"], df["Number of Reaction Steps"]):
        if not use_chemical_names:
            reacts = [nodes[temp] for temp in x.split(sep="; ")]
            products = [nodes[temp] for temp in y.split(sep="; ")]
        else:
            reacts = [temp for temp in x.split(sep="; ")]
            products = [temp for temp in y.split(sep="; ")]
        for react in reacts:
            for product in products:
                if allow_multiple_edges:
                    edges.add((react, product, z))
                    added_edges[(react, product)] = z
                else:
                    if (react, product) in added_edges and added_edges[(react, product)] > z:
                        added_edges[(react, product)] = z
                    elif (react, product) not in added_edges:
                        added_edges[(react, product)] = z
    if not allow_multiple_edges:
        edges = set()
        for key1, key2 in added_edges:
            edges.add((key1, key2, added_edges[(key1, key2)]))
    return edges


def count_by(atr_number, tuple_list):
    """
    Groups the elements in the list by specified attribute.
    Counts how many times each element has occurred
    :param atr_number: Attribute number to group by with
    :param tuple_list: List of analysed tuples
    :return: list of distinct elements for the specified field
    :return: list of counts the element from the corresponding place in the other list appeared
    """
    only_attributes = [element[atr_number] for element in tuple_list]
    counter = {}
    for x in only_attributes:
        if x in counter:
            counter[x] = counter[x] + 1
        else:
            counter[x] = 1
    x = counter.keys()
    y = counter.values()
    return x, y
