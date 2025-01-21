import numpy as np
import pandas as pd
import powerlaw
from matplotlib import pyplot as plt
from src import utils


def analyse(df, plot_type="distribution_comparison_plot"):
    """
    Analyzes the degree distribution of a graph created from the DataFrame `df`.

    Depending on the `plot_type`, different plots or calculations are performed.

    :param df: DataFrame with appropriate columns to create a graph.
    :param plot_type: Type of plot or calculation to perform. Defaults to "distribution_comparison_plot".
        Options include:
        - "bar" - creates bar plot, representing degree distribution from df.
        - "scatter" - creates scatter plot for in, out and total-degree distribution.
        - "exponent" - finds the best exponent and kmin using the powerlaw package to best fit the distribution
            approximation.
        - "scatter_comparison_plot" - creates scatter plot comparing Biorefinery Network and Network of Organic
            Chemistry. It's done for both in and out degree.
        - "distribution_comparison_plot" - creates plots comparing degree distribution in the Biorefinery Network and
            Network of Organic Chemistry.
    :return: None
    """
    # Create a directed network to analyse molecules
    graph, nodes = utils.create_weighted_directed_graph(df, allow_multiple_edges=False)
    if plot_type == "bar":
        create_bar_plot(graph)
    elif plot_type == "scatter":
        create_scatter_plot(graph)
    elif plot_type == "exponent":
        print("\nAnalysing exponent for the total-degree...")
        find_exponent([int(num[1]) for num in graph.degree()], "total")
        print("\nAnalysing exponent for the in-degree...")
        find_exponent([int(num[1]) for num in graph.in_degree()], "in")
        print("\nAnalysing exponent for the out-degree...")
        find_exponent([int(num[1]) for num in graph.out_degree()], "out")
    elif plot_type == "scatter_comparison_plot":
        create_comparison_scatter_plot_option1(graph, "out")
        create_comparison_scatter_plot_option1(graph, "in")
    elif plot_type == "distribution_comparison_plot":
        create_pdf_comparison(graph, "out", "Comparison of degree distribution in biorefinery network and NOC")
        create_pdf_comparison(graph, "in", "Comparison of degree distribution in biorefinery network and NOC")
        create_pdf_comparison_in_out(graph)


def create_scatter_plot(graph):
    """
    Creates scatter plots for in-degree, out-degree, and total-degree distributions of a graph.

    :param graph: Directed graph to analyze.
    :return: None
    """
    figure, axis = plt.subplots(3, 1)
    figure.set_figwidth(4)
    figure.suptitle("Degree distribution")
    plt.tight_layout(h_pad=2, pad=2)
    for i in range(3):
        axis[i].set_yscale("log")
        axis[i].set_xscale("log")

    # Plot the in-degree
    x, y = utils.count_by(1, graph.in_degree())
    axis[0].set_xlabel("In-degree, $k_{in}$", labelpad=-2)
    axis[0].set_ylabel("Count")
    axis[0].set_xlim(right=100)
    axis[0].set_ylim(bottom=0.1, top=1000)
    axis[0].scatter(x, y, marker=".")

    # Plot the out-degree
    x, y = utils.count_by(1, graph.out_degree())
    axis[1].set_xlabel("Out-degree, $k_{out}$", labelpad=-2)
    axis[1].set_ylabel("Count")
    axis[1].set_xlim(right=100)
    axis[1].set_ylim(bottom=0.1, top=1000)
    axis[1].scatter(x, y, marker=".")

    # Plot the in&out degree
    x, y = utils.count_by(1, graph.degree())
    axis[2].set_xlabel("Total-degree, $k_{total}$", labelpad=-2)
    axis[2].set_ylabel("Count")
    axis[2].set_xlim(right=100)
    axis[2].set_ylim(bottom=0.1, top=1000)
    axis[2].scatter(x, y, marker=".")

    plt.subplots_adjust(left=0.2)
    plt.show()


def create_bar_plot(graph):
    """
    Creates bar plots for in-degree, out-degree, and total-degree distributions of a graph.

    :param graph: Directed graph to analyze.
    :return: None
    """
    figure, axis = plt.subplots(3, 2)
    figure.suptitle("Degree distribution")
    plt.tight_layout(h_pad=2)

    # Plot the in-degree
    x, y = utils.count_by(1, graph.in_degree())
    maximum_in_degree = max(x)
    axis[0, 0].set_title("In-degree")
    axis[0, 0].set_xlim(0, 20)
    axis[0, 0].bar(x, y, align="edge")

    axis[0, 1].set_title("In-degree")
    axis[0, 1].set_xlim(10, maximum_in_degree)
    axis[0, 1].set_ylim(0, 10)
    axis[0, 1].bar(x, y, align="edge")

    # Plot the out-degree
    x, y = utils.count_by(1, graph.out_degree())
    maximum_out_degree = max(x)
    axis[1, 0].set_title("Out-degree")
    axis[1, 0].set_xlim(0, 20)
    axis[1, 0].bar(x, y, align="edge")

    axis[1, 1].set_title("Out-degree")
    axis[1, 1].set_xlim(10, maximum_out_degree)
    axis[1, 1].set_ylim(0, 10)
    axis[1, 1].bar(x, y, align="edge")

    # Plot the in&out degree
    x, y = utils.count_by(1, graph.degree())
    maximum_inout_degree = max(x)
    axis[2, 0].set_title("In&Out-degree")
    axis[2, 0].set_xlim(0, 20)
    axis[2, 0].bar(x, y, align="edge")

    axis[2, 1].set_title("In&Out-degree")
    axis[2, 1].set_xlim(10, maximum_inout_degree)
    axis[2, 1].set_ylim(0, 15)
    axis[2, 1].bar(x, y, align="edge")

    plt.show()


def find_exponent(data_init, text):
    """
    Finds and displays the exponent of the power-law distribution for the given degree data.
    Prints comparison of how well different distributions approximate the given data.

    :param data_init: List of degree values.
    :param text: Type of degree (total, in, out) for labeling purposes.
    :return: None
    """
    data = []
    for x in data_init:
        if x != 0:
            data.append(x)
    data.sort()
    results = powerlaw.Fit(data, discrete=True, estimate_discrete=False, xmin=2)
    print("Exponent =", results.power_law.alpha)
    print("K_min =", results.power_law.xmin)

    plt.xlim(1, 100)
    fig = results.plot_pdf(original_data=True, color="orange", linewidth=2)
    results.power_law.plot_pdf(color="r", linestyle="--", ax=fig)

    fig.set_title("$discrete, y=" + str(np.round(results.alpha, 4)) + ", k_{min}=" + str(results.xmin) + "$")
    fig.set_xlabel("$k_{" + text + "}$")
    fig.set_ylabel("$P(k_{" + text + "})$")
    fig.legend(["pdf discrete", "pdf empirical"])
    plt.show()

    print("\nCandidate 1 vs Candidate 2 \t\t\t\t\t\tR\t\t\t\t\tp")
    print("Power law vs exponential:\t\t\t\t", results.distribution_compare("power_law", "exponential"))
    print("Power law vs lognormal:\t\t\t\t\t", results.distribution_compare("power_law", "lognormal"))
    print("Power law vs truncated power law:\t\t", results.distribution_compare("power_law", "truncated_power_law"))
    print("Power law vs stretched exponential:\t\t", results.distribution_compare("power_law", "stretched_exponential"))
    print("Power law vs lognormal positive:\t\t", results.distribution_compare("power_law", "lognormal_positive"))
    print()
    print("Truncated power law vs lognormal:\t\t", results.distribution_compare("truncated_power_law", "lognormal"))
    print("Truncated power law vs exponential:\t\t", results.distribution_compare("truncated_power_law", "exponential"))
    print(
        "Truncated power law vs stretched exponential:\t\t",
        results.distribution_compare("truncated_power_law", "stretched_exponential"),
    )
    print(
        "Truncated power law vs lognormal positive:\t\t",
        results.distribution_compare("truncated_power_law", "lognormal_positive"),
    )
    print()
    print(
        "Lognormal vs stretched exponential:\t\t\t\t",
        results.distribution_compare("lognormal", "stretched_exponential"),
    )


def create_comparison_scatter_plot_option1(graph, option):
    """
    Creates a scatter plot comparing the degree distribution of a graph with pre-existing data.
    Chosen as the used default.

    :param graph: Directed graph to analyze.
    :param option: Degree type to plot ("in" or "out").
    :return: None
    """
    figure, axis = plt.subplots(1, 1)
    figure.suptitle("Degree distribution")
    axis.set_ylabel("Count")
    axis.set_yscale("log")
    axis.set_xscale("log")

    if option == "out":
        data_file = "../resources/NOC data/out_degree_distribution.csv"
        x, y = utils.count_by(1, graph.out_degree())
        axis.set_xlabel("$Out-degree, k_{out}$")
    elif option == "in":
        data_file = "../resources/NOC data/in_degree_distribution.csv"
        x, y = utils.count_by(1, graph.in_degree())
        axis.set_xlabel("$In-degree, k_{in}$")
    else:
        raise Exception("Invalid option provided")

    df = pd.read_csv(data_file, header=None)
    axis.scatter(df[0], df[1], marker=".")
    axis.scatter(x, y, marker=".", c="red")
    plt.show()


def create_comparison_scatter_plot_option2(graph, option):
    """
    Creates a scatter plot comparing the degree distribution of a graph with pre-existing data.
    Alternative to the above method. Not used in the final version.

    :param graph: Directed graph to analyze.
    :param option: Degree type to plot ("in" or "out").
    :return: None
    """
    figure, axis = plt.subplots(1, 2)
    figure.suptitle("Degree distribution")
    figure.set_figwidth(9.6)
    for i in range(2):
        axis[i].set_yscale("log")
        axis[i].set_xscale("log")
        axis[i].set_ylabel("Count")

    if option == "out":
        data_file = "../resources/NOC data/out_degree_distribution.csv"
        x, y = utils.count_by(1, graph.out_degree())
        for i in range(2):
            axis[i].set_xlabel("$Out-degree, k_{out}$")
    elif option == "in":
        data_file = "../resources/NOC data/in_degree_distribution.csv"
        x, y = utils.count_by(1, graph.in_degree())
        for i in range(2):
            axis[i].set_xlabel("$In-degree, k_{in}$")
    else:
        raise Exception("Invalid option provided")

    df = pd.read_csv(data_file, header=None)
    axis[0].scatter(df[0], df[1], marker=".")
    axis[1].scatter(x, y, marker=".")

    plt.show()


def create_pdf_comparison(graph, option, title):
    """
    Creates a PDF comparison plot for the given degree distribution type (in or out).

    :param graph: Directed graph to analyze.
    :param option: Degree type to plot ("in" or "out").
    :param title: Super title to be used for the plot.
    :return: None
    """
    fig, ax = plt.subplots(1, 2, sharey="all", sharex="all")
    fig.set_figwidth(9.6)
    fig.supylabel("$P(k_{" + option + "})$")
    fig.supxlabel("$k_{" + option + "}$")
    if option == "out":
        data_init = [int(num[1]) for num in graph.out_degree()]
        data_file = "../resources/NOC data/out_pdf_empirical.csv"
        data_file_exp = "../resources/NOC data/out_pdf_discrete.csv"
        ax[1].set_title("$discrete, y=2.0935, k_{min}=1.0$")
        xmin = 2
    elif option == "in":
        data_init = [int(num[1]) for num in graph.in_degree()]
        data_file = "../resources/NOC data/in_pdf_empirical.csv"
        data_file_exp = "../resources/NOC data/in_pdf_discrete.csv"
        ax[1].set_title("$discrete, y=3.0436, k_{min}=3.0$")
        xmin = None
    else:
        raise Exception("Invalid option provided")

    df = pd.read_csv(data_file, header=None)
    df2 = pd.read_csv(data_file_exp, header=None)

    data = []
    for x in data_init:
        if x != 0:
            data.append(x)
    data.sort()

    results = powerlaw.Fit(data, discrete=True, estimate_discrete=False, xmin=xmin)
    ax[0].set_title(
        "$discrete, y=" + str(np.round(results.power_law.alpha, 4)) + ", k_{min}=" + str(results.xmin) + "$"
    )

    results.plot_pdf(original_data=True, color="orange", linewidth=2, ax=ax[0])
    results.power_law.plot_pdf(color="r", linestyle="--", ax=ax[0])
    ax[1].loglog(df[0], df[1], color="blue")
    ax[1].loglog(df2[0], df2[1], color="green", linestyle="--")

    ax[0].legend(["pdf empirical", "pdf discrete"])
    ax[1].legend(["pdf empirical", "pdf discrete"])
    fig.suptitle(title)
    plt.show()


def create_pdf_comparison_in_out(graph):
    """
    Creates a PDF comparison plot for the in and out-degree distribution displayed side by side.

    :param graph: Directed graph to analyze.
    :return: None
    """
    # Set fonts for the paper plots
    plt.rc("axes", labelsize=20)
    plt.rc("legend", fontsize=16)
    plt.rc("xtick", labelsize=16)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=16)  # fontsize of the tick labels

    fig, ax = plt.subplots(1, 2, sharey="all", sharex="all")
    fig.set_figheight(6.0)
    fig.set_figwidth(10.6)

    data_init_out = [int(num[1]) for num in graph.out_degree()]
    data_file_out = "../resources/NOC data/out_pdf_empirical.csv"
    data_file_exp_out = "../resources/NOC data/out_pdf_discrete.csv"
    xmin = 2

    df = pd.read_csv(data_file_out, header=None)
    df2 = pd.read_csv(data_file_exp_out, header=None)

    data = []
    for x in data_init_out:
        if x != 0:
            data.append(x)
    data.sort()

    results = powerlaw.Fit(data, discrete=True, estimate_discrete=False, xmin=xmin)
    ax[0].set_ylabel("$P(k_{out})$")
    ax[0].set_xlabel("$k_{out}$")
    ax[0].loglog(df[0], df[1], color="blue")
    ax[0].loglog(df2[0], df2[1], color="green", linestyle="--")
    results.plot_pdf(original_data=True, color="orange", linewidth=2, ax=ax[0])
    results.power_law.plot_pdf(color="r", linestyle="-.", ax=ax[0])

    ax[0].legend(["pdf empirical NOC", "pdf discrete NOC", "pdf empirical biorefinery", "pdf discrete biorefinery"])

    data_init_in = [int(num[1]) for num in graph.in_degree()]
    data_file_in = "../resources/NOC data/in_pdf_empirical.csv"
    data_file_exp_in = "../resources/NOC data/in_pdf_discrete.csv"
    xmin = None

    df = pd.read_csv(data_file_in, header=None)
    df2 = pd.read_csv(data_file_exp_in, header=None)

    data = []
    for x in data_init_in:
        if x != 0:
            data.append(x)
    data.sort()

    results = powerlaw.Fit(data, discrete=True, estimate_discrete=False, xmin=xmin)
    ax[1].set_ylabel("$P(k_{in})$")
    ax[1].set_xlabel("$k_{in}$")
    ax[1].loglog(df[0], df[1], color="blue")
    ax[1].loglog(df2[0], df2[1], color="green", linestyle="--")
    results.plot_pdf(original_data=True, color="orange", linewidth=2, ax=ax[1])
    results.power_law.plot_pdf(color="r", linestyle="-.", ax=ax[1])

    ax[1].legend(["pdf empirical NOC", "pdf discrete NOC", "pdf empirical biorefinery", "pdf discrete biorefinery"])

    plt.margins = 0
    plt.tight_layout()
    plt.show()
