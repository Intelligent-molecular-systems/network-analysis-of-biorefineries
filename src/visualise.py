import itertools

import pandas as pd
from matplotlib import pyplot as plt


def plot_components(component_sizes, axis, a):
    """
    Helper method. Used in `graph_fragmentation.py`. Scatter plots the size of each connected component
    on the given axis and subplot.
    """
    axis[a].plot(component_sizes, "o", markersize=3)
    axis[a].set_xlabel("Index of the component")
    axis[a].set_ylabel("Size of the x-th largest component")
    axis[a].set_title("Sizes of the connected components")


def plot_components_size_distribution(component_sizes, axis, a):
    """
    Helper method. Used in `graph_fragmentation.py`. Plots the distribution of the connected component sizes
    on the given axis and subplot.
    """
    x = [(k, len(list(g))) for k, g in itertools.groupby(component_sizes)]
    component_size = []
    number_of_components = []
    for size, num in x:
        component_size.append(str(size))
        number_of_components.append(num)

    axis[a].bar(component_size, number_of_components)
    axis[a].set_title("Distribution of component sizes")
    axis[a].set_xlabel("Size of the component")
    axis[a].set_ylabel("The number of components")
    axis[a].set_ylim(0, max(number_of_components) * 1.2)
    axis[a].yaxis.grid(True)


def plot_comparison_with_NOC_avg_shortest_path_length(x, y, data_file="../resources/NOC data/average_path_length.csv"):
    """
    Helper method. Used in `graph_property.py`. Plots the average shortest path length using given x and y.
    Plots the average shortest path length using data from the data_file.
    :param x: List with X coordinates for the data to plot.
    :param y: List with Y coordinates for the data to plot.
    :param data_file: Path to the file with the data for the compared plot.
    :return: None
    """
    df = pd.read_csv(data_file, header=None)
    plt.loglog(df[0], df[1], "o-b")
    plt.loglog(x, y, "o-r")
    plt.ylim(top=1)
    plt.xlim(right=100)
    plt.xlabel("Shortest path length, l")
    plt.ylabel("Probability, P(l)")
    plt.legend(["Network of organic chemistry", "Biorefinery network"])


def plot_correlation(avg_deg, degree_type):
    x = avg_deg.keys()
    y = avg_deg.values()
    fig, ax = plt.subplots()
    ax.set_title("Degree correlation")
    ax.set_ylabel("Average target " + degree_type + "-degree")
    ax.set_xlabel("Source " + degree_type + "-degree")
    ax.scatter(x, y)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim(right=100)
    ax.set_ylim(top=100)
    ax.axline((0, 0), (5, 5))
    plt.show()
