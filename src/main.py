import graph_fragmentation
import important_molecules
import degree_distribution
import graph_property
import pandas
from preprocess import preprocess
from src import graph_clusters, degree_correlation, central_point_dominance, dataset_comparison
from src.read_from_file import read_columns

weighted_data = [
    "important_molecules",
    "degree_distribution",
    "graph_property",
    "graph_clusters",
    "degree_correlation",
    "central_point_dominance",
]


def create_columns(columns, option):
    """
    Method for returning, which columns from the DataFrame should be extracted.
    :param columns: Columns to extract if already specified.
    :param option: Type of analysis.
    :return: List with indexes of columns to extract.
    """
    # Specific columns are specified. Return them
    if columns is not None:
        return columns
    # The analysis uses weighted data. Extract columns corresponding to Reactants, Products and Number of Reaction Steps
    elif option in weighted_data:
        return ["Reactant", "Product", "Number of Reaction Steps"]
    # The analysis uses unweighted data. xtract columns corresponding to Reactants and Products
    else:
        return ["Reactant", "Product"]


def analyse(
    option,
    columns=None,
    my_file="../resources/reaction_data.tsv",
    file_merge="../resources/lignocellulosic_reactions.tsv",
    file_compare="../resources/lignocellulosic_reactions.tsv",
):
    """
    Method providing access to various analysis methods.
    :param option: String specifying what kind of analysis to perform. Possible options are:
    - "merge_and_graph_fragmentation" - merges data from `my_file` with the `file_merge` data. Then, performs
        "graph_fragmentation" analysis on such combined data.
    - "graph_fragmentation" - analyses the network fragmentation. Provides info about connected components like
        their size distribution. Creates HTML visualisation of the whole network.
    - "important_molecules" - analyses important molecules in the network according to the highest degree and
        betweenness centrality. Plots the results. Creates the HTML visualisation where the most important molecules
        are highlighted.
    - "degree_distribution" - analyses the degree distribution of the network. Creates plot specified as the `plot_type`
        and finds the exponent for the power law best fitting the data.
    - "graph_property" - calculates network metrics (average shortest path length, clustering coefficient, omega metric)
        and plots the average shortest path length as probability, compared to the Network of Organic Chemistry.
    - "graph_clusters" - performs k-core clustering and girvan-neuman clustering. Finds the best girvan-neuman
        clustering according to created metric and creates its HTML visualisation.
    - "degree_correlation" - calculates and plots the degree correlation for both in and out-degree.
    - "central_point_dominance" - calculates central point dominance
    - "dataset_comparison" - compares the dataset with the data from `file_compare` by calculating Jaccard similarity
        and checking whether important molecules from Biorefinery Network appear in the other dataset.
    :param columns: List of columns to extract for specific analysis. When None it is set to [6, 7, 14] for analysis
        with weighted data and [6, 7] for unweighted. Column numbers are consistent with REAXYS output.
    :param my_file: .tsv file with data to analyse.
    :param file_merge: .tsv file with data that should be merged with `my_file` data for `merge_and_graph_fragmentation`
        analysis.
    :param file_compare: .tsv file with data that should be compared to `my_file` data.
    :return:
    """

    columns = create_columns(columns, option)
    df = read_columns(my_file, columns)
    df = preprocess(df)

    if option == "merge_and_graph_fragmentation":
        df_add = read_columns(file_merge, ["Reactant", "Product"])
        df_add = preprocess(df_add)
        df = pandas.concat([df, df_add], axis="rows")
        graph_fragmentation.analyse(df, "graph_fragmentation_merged_data")
    if option == "graph_fragmentation":
        graph_fragmentation.analyse(df, "graph_fragmentation")
    if option == "important_molecules":
        important_molecules.analyse(df)
    if option == "degree_distribution":
        degree_distribution.analyse(df, "distribution_comparison_plot")
    if option == "graph_property":
        # Calculate metrics (average shortest path length, clustering coefficient, omega) to analyse network properties
        graph_property.analyse(df)
        # Plot the comparison of Biorefinery Network average shortest path length and the Network of Organic Chemistry
        graph_property.compare_average_shortest_path_length_directed(df)
    if option == "graph_clusters":
        graph_clusters.analyse(df)
    if option == "degree_correlation":
        degree_correlation.analyse(df, degree_type="out")
        degree_correlation.analyse(df, degree_type="in")
    if option == "central_point_dominance":
        central_point_dominance.analyse(df)
    if option == "dataset_comparison":
        df_to_compare_to = read_columns(file_compare, ["Reactant", "Product"])
        df_to_compare_to = preprocess(df_to_compare_to)
        dataset_comparison.analyse(df, df_to_compare_to)
    return 0


if __name__ == "__main__":
    # Change option param to run different analysis.
    option = "degree_distribution"

    # Below parameters are set to default. Look for `analyse` docstring for their descriptions.
    columns = None
    my_file = "../resources/reaction_data.tsv"
    file_merge = "../resources/lignocellulosic_reactions.tsv"
    file_compare = "../resources/lignocellulosic_reactions.tsv"

    analyse(option, columns, my_file, file_merge, file_compare)
