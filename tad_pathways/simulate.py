"""
Gregory Way 2019
tad_pathways/simulate.py

Description:
Functions to simulate TAD_Pathways analysis

Usage: Import only

    from tad_pathways import simulate
"""

import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.stats import zscore

tad_column_musts = [
    "TAD_id",
    "TAD_end",
    "TAD_start",
    "gene_type",
    "gene_name",
    "chromosome",
]
tad_df_err_msg = ", ".join(x for x in tad_column_musts)
tad_df_err_msg = "TAD file must include certain columns: [{}]".format(tad_df_err_msg)


def define_tad_similarity(
    tad_index_file="None", tad_df="None", gene_type_list="all", metric="euclidean"
):
    """
    Given a TAD Index File or pandas DataFrame, build a TAD similarity matrix. The
    similarity is defined based on TAD length and TAD content (number of genes).

    Arguments:
    tad_index_file - a string indicating the file path to the tad index file
    tad_df - a pandas DataFrame storing genic information about all specified TADs
    gene_type_list - a list of genetypes to consider when building the similarity
                     matrix. Defaults to "all", which considers all available genetypes.
    metric - similarity metric as input into scipy.spatial.distance.cdist

    Output:
    a ranked similarity matrix where each column and row represents each TAD. The
    entries represent the relative similarity rank of the TAD in the column.
    """
    # Read in data
    try:
        tad_df = pd.read_table(tad_index_file, index_col=0)
    except FileNotFoundError:
        if type(tad_df) == str:
            raise ValueError("User must specify a correct path tad_file or tad_df")

    # Assert that specific columns are present in the input tad DataFrame
    assert all(x in tad_df.columns for x in tad_column_musts), tad_df_err_msg

    if gene_type_list != "all":
        # Check genetypes
        available_gene_types = sorted(tad_df.gene_type.unique())
        genetype_err_msg = ", ".join(x for x in available_gene_types)
        genetype_err_msg = "Must use valid genetypes: [{}]".format(genetype_err_msg)
        assert all(x in available_gene_types for x in gene_type_list), genetype_err_msg

        # Subset to the gene types of interest
        tad_df = tad_df.query("gene_type in @gene_type_list")

    # Define TAD length
    tad_df = tad_df.assign(TAD_length=tad_df.TAD_end - tad_df.TAD_start)

    # Count genetypes
    tad_df = (
        tad_df.groupby(["TAD_id", "TAD_length", "gene_type"])
        .agg("count")
        .reset_index()
        .sort_values(by="TAD_length", ascending=False)
        .reset_index(drop=True)
        .loc[:, ["TAD_id", "TAD_length", "gene_type", "chromosome"]]
        .rename({"chromosome": "num_genes"}, axis="columns")
    )

    # Remove excess info
    tad_df.index = tad_df.TAD_id
    tad_df = tad_df.drop("TAD_id", axis="columns").drop("Boundary", axis="rows")

    # Split out the gene types
    gene_type_df = (
        tad_df.reset_index()
        .pivot_table(values="num_genes", index="TAD_id", columns="gene_type")
        .fillna(0)
        .astype(int)
    )

    # Merge the two dataframes
    tad_df = (
        tad_df.loc[:, ["TAD_length"]]
        .drop_duplicates()
        .merge(gene_type_df, how="outer", left_index=True, right_index=True)
    )

    # Scale data
    tad_df = pd.DataFrame(
        zscore(tad_df, axis=0), index=tad_df.index, columns=tad_df.columns
    ).fillna(0)

    # Get the distance
    tad_distance_df = distance.cdist(tad_df, tad_df, metric=metric)

    # Build new dataframe
    tad_distance_df = pd.DataFrame(
        tad_distance_df, index=tad_df.index, columns=tad_df.index
    )

    # Get similarity rank
    tad_similarity_df = tad_distance_df.rank(method="first", axis="rows")
    tad_similarity_df = tad_similarity_df - 1  # Self similarity equal to 0
    tad_similarity_df = tad_similarity_df.astype(int)

    return tad_similarity_df


def draw_random_tad(similarity_matrix_df, tad_id, p=0.2):
    """
    Randomly draw a similarity matched TAD

    Arguments:
    similarity_matrix_df - a symmetrical matrix indicating similarity ranked TADs
                           output from function `define_tad_similarity()`
    tad_id - an identifier in the similarity_matrix_df that indicates the TAD to match
    p - parameter in the np.random.negative_binomial distribution. On the scale [0,1].
        p represents the probability of success. A higher number will sample more
        closely matched TADs, where a lower number will tend to select more distantly
        similar TADs.

    Output:
    A similarity matched TAD id
    """
    tad_rank_array = similarity_matrix_df.loc[:, tad_id].sort_values()
    random_rank = np.random.negative_binomial(n=1, p=p, size=1) + 1
    random_tad = tad_rank_array[random_rank]
    return random_tad.index[0]


def simulate_tad_list(tad_list, similarity_matrix_df, p=0.2, with_replacement=False,
                      restrict_list="None"):
    """
    Randomly sample similarity matched TADs. Makes repeated calls to `draw_random_tad()`

    Arguments:
    tad_list - a list of TAD ids to match with random TADs
    similarity_matrix_df - a symmetrical matrix indicating similarity ranked TADs
                           output from function `define_tad_similarity()`
    tad_id - an identifier in the similarity_matrix_df that indicates the TAD to match
    p - parameter in the np.random.negative_binomial distribution. On the scale [0,1].
        p represents the probability of success. A higher number will sample more
        closely matched TADs, where a lower number will tend to select more distantly
        similar TADs.
    with_replacement - a boolean indicating if TADs should be sampled with or without
                       replacement. Sampling with replacement means TADs can be
                       sampled more than once.
    restrict_list - a list of TADs that cannot be selected to sample

    Output:
    A list of similarity matched TAD ids
    """
    simulated_tad_list = []
    pass_restrict = True

    for tad in tad_list:

        while pass_restrict:

            random_tad_id = draw_random_tad(
                similarity_matrix_df=similarity_matrix_df, tad_id=tad, p=p
            )

            if random_tad_id not in restrict_list:
                pass_restrict = False

        simulated_tad_list.append(random_tad_id)
        if not with_replacement:
            try:
                similarity_matrix_df = similarity_matrix_df.drop(
                    random_tad_id, axis="rows"
                ).copy()
            except KeyError:
                continue

        pass_restrict = True

    return simulated_tad_list


def simulate_tad_genelist(
    tad_df,
    tad_list,
    similarity_matrix_df,
    p=0.2,
    with_replacement=False,
    return_full_dataframe=False,
    restrict_list="None"
):
    """
    Randomly sample similarity matched TADs. Wrapper for `simulate_tad_list()`. Will
    return the genelist of interest, or, the full TAD index DataFrame.

    Arguments:
    tad_df - the TAD index dataframe
    tad_list - a list of TAD ids to match with random TADs
    similarity_matrix_df - a symmetrical matrix indicating similarity ranked TADs
                           output from function `define_tad_similarity()`
    tad_id - an identifier in the similarity_matrix_df that indicates the TAD to match
    p - parameter in the np.random.negative_binomial distribution. On the scale [0,1].
        p represents the probability of success. A higher number will sample more
        closely matched TADs, where a lower number will tend to select more distantly
        similar TADs.
    with_replacement - a boolean indicating if TADs should be sampled with or without
                       replacement. Sampling with replacement means TADs can be
                       sampled more than once.
    return_full_dataframe - a boolean if the full TAD index data should be returned
    restrict_list - a list of TADs that cannot be selected to sample

    Output:
    Either a genelist or a subsetted TAD index of similarity matched random TADs.
    """
    # Assert that specific columns are present in the input tad DataFrame
    assert all(x in tad_df.columns for x in tad_column_musts), tad_df_err_msg

    # Simulate a list of tads
    simulated_tad_list = simulate_tad_list(
        tad_list=tad_list,
        similarity_matrix_df=similarity_matrix_df,
        p=p,
        with_replacement=with_replacement,
        restrict_list=restrict_list
    )

    tad_df = tad_df.query("TAD_id in @simulated_tad_list")

    if return_full_dataframe:
        return tad_df
    else:
        return tad_df.gene_name.tolist()
