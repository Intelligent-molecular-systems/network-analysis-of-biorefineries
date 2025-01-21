import pandas as pd
import pytest
from networkx import MultiDiGraph, is_weighted

from src import utils


@pytest.fixture
def df_1():
    data = [
        ["calcium; oxygen; nitrogen", "oxygen; hydrogen"],
        ["carbon-dioxide", "carbon-dioxide"],
        ["nitrogen", "nitrogen"],
        ["oxygen", "oxygen"],
        ["oxygen", "oxygen"],
    ]
    return pd.DataFrame(data=data, columns=["Reactant", "Product"])


def get_expected_edges_df_1(nodes):
    return {
        (nodes["calcium"], nodes["oxygen"], 1),
        (nodes["oxygen"], nodes["oxygen"], 1),
        (nodes["nitrogen"], nodes["oxygen"], 1),
        (nodes["calcium"], nodes["hydrogen"], 1),
        (nodes["oxygen"], nodes["hydrogen"], 1),
        (nodes["nitrogen"], nodes["hydrogen"], 1),
        (nodes["carbon-dioxide"], nodes["carbon-dioxide"], 1),
        (nodes["nitrogen"], nodes["nitrogen"], 1),
    }


def get_expected_edges_df_with_steps(nodes):
    return {
        (nodes["calcium"], nodes["oxygen"], 2),
        (nodes["oxygen"], nodes["oxygen"], 2),
        (nodes["nitrogen"], nodes["oxygen"], 2),
        (nodes["calcium"], nodes["hydrogen"], 2),
        (nodes["oxygen"], nodes["hydrogen"], 2),
        (nodes["nitrogen"], nodes["hydrogen"], 2),
        (nodes["carbon-dioxide"], nodes["carbon-dioxide"], 1),
        (nodes["nitrogen"], nodes["nitrogen"], 3),
        (nodes["oxygen"], nodes["oxygen"], 7),
        (nodes["oxygen"], nodes["oxygen"], 4),
    }


def get_expected_edges_without_multiples(nodes):
    return {
        (nodes["calcium"], nodes["oxygen"], 2),
        (nodes["oxygen"], nodes["oxygen"], 2),
        (nodes["nitrogen"], nodes["oxygen"], 2),
        (nodes["calcium"], nodes["hydrogen"], 2),
        (nodes["oxygen"], nodes["hydrogen"], 2),
        (nodes["nitrogen"], nodes["hydrogen"], 2),
        (nodes["carbon-dioxide"], nodes["carbon-dioxide"], 1),
        (nodes["nitrogen"], nodes["nitrogen"], 3),
    }


@pytest.fixture
def expected_edges_with_chemical_names():
    return {
        ("calcium", "oxygen", 2),
        ("oxygen", "oxygen", 2),
        ("nitrogen", "oxygen", 2),
        ("calcium", "hydrogen", 2),
        ("oxygen", "hydrogen", 2),
        ("nitrogen", "hydrogen", 2),
        ("carbon-dioxide", "carbon-dioxide", 1),
        ("nitrogen", "nitrogen", 3),
        ("oxygen", "oxygen", 7),
        ("oxygen", "oxygen", 4),
    }


@pytest.fixture
def df_with_steps():
    data = [
        ["calcium; oxygen; nitrogen", "oxygen; hydrogen", 2],
        ["carbon-dioxide", "carbon-dioxide", 1],
        ["nitrogen", "nitrogen", 3],
        ["oxygen", "oxygen", 7],
        ["oxygen", "oxygen", 4],
    ]
    return pd.DataFrame(data=data, columns=["Reactant", "Product", "Number of Reaction Steps"])


@pytest.fixture
def list_of_tuples():
    data = [
        (1, 0),
        (2, 0),
        (3, 0),
        (1, 1),
        (2, 2),
    ]
    return data


def test_get_distinct_chemicals_1(df_1):
    actual_chemicals = utils.get_distinct_chemicals(df_1)
    expected_chemicals = ["oxygen", "hydrogen", "calcium", "nitrogen", "carbon-dioxide"]
    assert len(actual_chemicals.keys()) == len(expected_chemicals)
    assert set(actual_chemicals.keys()) == set(expected_chemicals)
    assert len(actual_chemicals.values()) == len(expected_chemicals)


def test_get_edges_without_number_of_reaction_steps(df_1):
    nodes = utils.get_distinct_chemicals(df_1)
    expected_edges = get_expected_edges_df_1(nodes)
    actual_edges = utils.get_edges(df_1, nodes)
    assert expected_edges == actual_edges


def test_get_edges_with_number_of_reaction_steps(df_with_steps):
    nodes = utils.get_distinct_chemicals(df_with_steps)
    expected_edges = get_expected_edges_df_with_steps(nodes)
    actual_edges = utils.get_edges(df_with_steps, nodes)
    assert expected_edges == actual_edges


def test_get_edges_with_chemical_names(df_with_steps, expected_edges_with_chemical_names):
    nodes = utils.get_distinct_chemicals(df_with_steps)
    expected_edges = expected_edges_with_chemical_names
    actual_edges = utils.get_edges(df_with_steps, nodes, use_chemical_names=True)
    assert expected_edges == actual_edges


def test_get_edges_not_allow_multiples(df_with_steps):
    nodes = utils.get_distinct_chemicals(df_with_steps)
    expected_edges = get_expected_edges_without_multiples(nodes)
    actual_edges = utils.get_edges(df_with_steps, nodes, allow_multiple_edges=False)
    assert expected_edges == actual_edges


def test_create_weighted_directed_graph(df_with_steps):
    graph, nodes = utils.create_weighted_directed_graph(df_with_steps)
    assert isinstance(graph, MultiDiGraph)
    assert is_weighted(graph)
    assert set(graph.edges.data("weight")) == set(get_expected_edges_df_with_steps(nodes))


def test_group_by_first_element(list_of_tuples):
    x, y = utils.count_by(0, list_of_tuples)
    expected = {1: 2, 2: 2, 3: 1}
    expected_x = expected.keys()
    expected_y = expected.values()
    assert x == expected_x
    assert sorted(y) == sorted(expected_y)


def test_group_by_second_element(list_of_tuples):
    x, y = utils.count_by(1, list_of_tuples)
    expected = {0: 3, 1: 1, 2: 1}
    expected_x = expected.keys()
    expected_y = expected.values()
    assert x == expected_x
    assert sorted(y) == sorted(expected_y)
