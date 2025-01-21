import pandas as pd


def read_columns(my_file, columns):
    """
    Read the data from csv file and drop NA fields
    :param my_file: File from which to read data
    :param columns: List of columns to be read
    :return: DataFrame with only columns requested
    """
    # Read file
    df = pd.read_csv(my_file, sep="\t", usecols=columns)
    df = df.dropna()
    return df
