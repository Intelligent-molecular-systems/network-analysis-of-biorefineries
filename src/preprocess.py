def preprocess(df):
    """
    Method for specifying the preprocessing pipeline.
    :param df: DataFrame to preprocess
    :return: Preprocessed DataFrame
    """
    step1_lowercase(df)
    return df


def step1_lowercase(df):
    """
    Lowercase all names in Reactant and Product fields.
    :param df: DataFrame to preprocess
    :return: Preprocessed DataFrame
    """
    df["Reactant"] = df["Reactant"].apply(lambda x: x.lower() if isinstance(x, str) else x)
    df["Product"] = df["Product"].apply(lambda x: x.lower() if isinstance(x, str) else x)
