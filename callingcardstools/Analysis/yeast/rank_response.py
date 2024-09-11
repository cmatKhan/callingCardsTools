"""
The critical function is `create_rank_response_table`. This takes an argparse object
with the appropriate arguments and returns a tuple that includes the random expectation
dataframe, the merged binding/expression data with bins and responsive labels according
to the threshold settings, and the summarized dataframe. The `main` method currently
only returns the summarized dataframe.
"""

import argparse
import json
import logging
import os
import re
from typing import Union

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from scipy.stats._result_classes import BinomTestResult

from .read_in_data import combine_data, read_in_data

logger = logging.getLogger(__name__)

__all__ = [
    "create_partitions",
    "bin_by_binding_rank",
    "parse_binomtest_results",
    "compute_rank_response",
    "label_responsive_genes",
    "create_rank_response_table",
    "set_none_str_to_none",
    "validate_config",
    "parse_args",
    "main",
]


def create_partitions(vector_length, equal_parts=100):
    """
    Splits a vector of a specified length into nearly equal partitions.

    This function creates a partition vector where each partition is of equal
    size, except the last partition which may be smaller depending on the
    vector length and the number of equal parts specified. Each element in
    the partition vector represents the partition number.

    Args:
        vector_length (int): The total length of the vector to be partitioned.
        equal_parts (int, optional): The number of equal parts to divide the
            vector. Defaults to 100.

    Returns:
        numpy.ndarray: An array where each element represents the partition
            number for each element in the original vector.

    Examples:
        >>> create_partitions(10, 3)
        array([1, 1, 1, 2, 2, 2, 3, 3, 3, 3])
    """
    quotient, remainder = divmod(vector_length, equal_parts)
    return np.concatenate(
        [
            np.repeat(np.arange(1, quotient + 1), equal_parts),
            np.repeat(quotient + 1, remainder),
        ]
    )


def bin_by_binding_rank(
    df: pd.DataFrame, bin_size: int, rank_by_binding_effect: bool = False
):
    """
    Assigns a rank bin to each row in a DataFrame based on binding signal.

    This function divides the DataFrame into partitions based on the specified
    bin size, assigns a rank to each row within these partitions, and then
    sorts the DataFrame based on the 'effect' and 'binding_pvalue' columns. The
    ranking is assigned such that rows within each bin get the same rank, and
    the rank value is determined by the bin size.

    Args:
        df (pd.DataFrame): The DataFrame to be ranked and sorted.
            It must contain 'effect' and 'binding_pvalue' columns.
        bin_size (int): The size of each bin for partitioning the DataFrame
            for ranking.
        rank_by_binding_effect (bool, optional): If True, the DataFrame is sorted by
            abs('effect') in descending order first with ties broken by pvalue.
            If False, sort by pvalue first with ties broken by effect size.
            Defaults to False

    Returns:
        pd.DataFrame: The input DataFrame with an added 'rank' column, sorted
            by 'effect' in descending order or 'binding_pvalue' in
            ascending order depending on `rank_by_binding_effect`.

    Example:
        >>> df = pd.DataFrame({'effect': [1.2, 0.5, 0.8],
        ...                    'binding_pvalue': [5, 3, 4]})
        >>> bin_by_binding_rank(df, 2)
        # Returns a DataFrame with added 'rank' column and sorted as per
        # the specified criteria.
    """
    if "binding_pvalue" not in df.columns:
        raise KeyError("Column 'binding_pvalue' is not in the data")
    if "binding_effect" not in df.columns:
        raise KeyError("Column 'binding_effect' is not in the data")

    parts = min(len(df), bin_size)
    df_abs = df.assign(abs_binding_effect=df["binding_effect"].abs())

    df_sorted = df_abs.sort_values(
        by=(
            ["abs_binding_effect", "binding_pvalue"]
            if rank_by_binding_effect
            else ["binding_pvalue", "abs_binding_effect"]
        ),
        ascending=[False, True] if rank_by_binding_effect else [True, False],
    )

    return (
        df_sorted.drop(columns=["abs_binding_effect"])
        .reset_index(drop=True)
        .assign(rank_bin=create_partitions(len(df_sorted), parts) * parts)
    )


def parse_binomtest_results(binomtest_obj: BinomTestResult, **kwargs):
    """
    Parses the results of a binomtest into a tuple of floats.

    This function takes the results of a binomtest and returns a tuple of
    floats containing the response ratio, p-value, and confidence interval
    bounds.

    Args:
        binomtest_obj (scipy.stats.BinomTestResult): The results of a binomtest
            for a single rank bin.
        Additional keyword arguments: Additional keyword arguments are passed
            to the proportional_ci method of the binomtest object.

    Returns:
        tuple: A tuple of floats containing the response ratio, p-value, and
            confidence interval bounds.

    Example:
        >>> parse_binomtest_results(binomtest(1, 2, 0.5, alternative='greater')
        (0.5, 0.75, 0.2, 0.8)
    """
    return (
        binomtest_obj.statistic,
        binomtest_obj.pvalue,
        binomtest_obj.proportion_ci(
            confidence_level=kwargs.get("confidence_level", 0.95),
            method=kwargs.get("method", "exact"),
        ).low,
        binomtest_obj.proportion_ci(
            confidence_level=kwargs.get("confidence_level", 0.95),
            method=kwargs.get("method", "exact"),
        ).high,
    )


def compute_rank_response(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """
    Computes rank-based statistics and binomial test results for a DataFrame.

    This function groups the DataFrame by 'rank_bin' and aggregates it to
    calculate the number of responsive items in each rank bin, as well as
    various statistics related to a binomial test.  It calculates the
    cumulative number of successes, response ratio, p-value, and confidence
    intervals for each rank bin.

    Args:
        df (pd.DataFrame): DataFrame containing the columns 'rank_bin',
            'responsive', and 'random'. 'rank_bin' is an integer representing
            the rank bin, 'responsive' is a boolean indicating responsiveness,
            and 'random' is a float representing the random expectation.
        Additional keyword arguments: Additional keyword arguments are passed
            to the binomtest function, including arguments to the
            proportional_ci method of the BinomTestResults object (see scipy
            documentation for details)

    Returns:
        pd.DataFrame: A DataFrame indexed by 'rank_bin' with columns for the
            number of responsive items in each bin ('n_responsive_in_rank'),
            cumulative number of successes ('n_successes'), response ratio
            ('response_ratio'), p-value ('p_value'), and confidence interval
            bounds ('ci_lower' and 'ci_upper').

    Example:
        >>> df = pd.DataFrame({'rank_bin': [1, 1, 2],
        ...                    'responsive': [True, False, True],
        ...                    'random': [0.5, 0.5, 0.5]})
        >>> compute_rank_response(df)
        # Returns a DataFrame with rank-based statistics and binomial
        # test results.
    """
    rank_response_df = (
        df.groupby("rank_bin")
        .agg(
            n_responsive_in_rank=pd.NamedAgg(column="responsive", aggfunc="sum"),
            random=pd.NamedAgg(column="random", aggfunc="first"),
        )
        .reset_index()
    )

    rank_response_df["n_successes"] = rank_response_df["n_responsive_in_rank"].cumsum()

    # Binomial Test and Confidence Interval
    rank_response_df[["response_ratio", "pvalue", "ci_lower", "ci_upper"]] = (
        rank_response_df.apply(
            lambda row: parse_binomtest_results(
                binomtest(
                    int(row["n_successes"]),
                    int(row.rank_bin),
                    float(row["random"]),
                    alternative=kwargs.get("alternative", "two-sided"),
                ),
                **kwargs,
            ),
            axis=1,
            result_type="expand",
        )
    )

    return rank_response_df


def label_responsive_genes(
    df,
    abs_expression_effect_threshold,
    expression_pvalue_threshold,
    normalization_cutoff: int = -1,
):
    """
    Labels genes in a DataFrame as responsive or not based on thresholds for
    expression effect and p-value. Note that the comparisons on the thresholds
    are strictly greater than for the abs_expression_effect_threshold and
    strictly less than for the expression_pvalue_threshold.

    The function adds a new boolean column 'responsive' to the DataFrame, where
    each gene is labeled as responsive if its absolute effect expression is
    strictly greater than a threshold and its p-value is strictly less than
    a specified threshold. If normalization is enabled, only the top genes
    meeting the criteria up to the minimum number found in the normalized
    subset are labeled as responsive.

    Args:
        df (pd.DataFrame): DataFrame containing gene data. Must include
            'expression_effect' and 'expression_pvalue' columns.
        abs_expression_effect_threshold (float): Absolute value threshold
            for the absolute value of the expression effect. Values strictly
            greater than this threshold are considered responsive if the pvalue
            threshold passes.
        expression_pvalue_threshold (float): Threshold for the expression
            p-value. Values strictly less than this threshold are considered
            responsive if the effect threshold passes.
        normalization_cutoff (int, optional): The maximum number of responsive
            genes to consider prior to labelling. This serves to normalize
            rank response across expression data sets. Defaults to -1, which
            disables normalization.

    Returns:
        pd.DataFrame: The input DataFrame with an added 'responsive' column.

    Raises:
        KeyError: If 'expression_effect' or 'expression_pvalue' are not in
        the DataFrame.

    Examples:
        >>> df = pd.DataFrame({'effect_expression': [0.5, 0.7, 1.2],
                               'p_expression': [0.01, 0.05, 0.2]})
        >>> label_responsive_genes(df, 0.6, 0.05).responsive
        [False, True, False]
    """
    if "expression_effect" not in df.columns:
        raise KeyError("Column 'effect_expression' is not in the data")
    if "expression_pvalue" not in df.columns:
        raise KeyError("Column 'effect_pvalue' is not in the data")

    expression_effect_rank_cutoff = (
        normalization_cutoff if normalization_cutoff > 0 else len(df) + 1
    )

    df_abs = df.assign(abs_expression_effect=df["expression_effect"].abs())

    # if either the effect or p-value threshold is `None`, then set
    # the threshold to the appropiate boundary to prevent filtering on that
    # column
    abs_expression_effect_threshold = (
        abs_expression_effect_threshold
        if abs_expression_effect_threshold is not None
        else min(df_abs["abs_expression_effect"]) - 1
    )

    expression_pvalue_threshold = (
        expression_pvalue_threshold
        if expression_pvalue_threshold is not None
        else max(df_abs["expression_pvalue"]) + 1
    )

    df_ranked = (
        df_abs.sort_values(
            by=["abs_expression_effect", "expression_pvalue"], ascending=[False, True]
        ).reset_index(drop=True)
        # Add 1 to start ranking from 1 instead of 0
        .assign(rank=lambda x: x.index + 1)
    )

    df_ranked["responsive"] = (
        (df_ranked["abs_expression_effect"] > abs_expression_effect_threshold)  # noqa
        & (df_ranked["expression_pvalue"] < expression_pvalue_threshold)
        & (df_ranked["rank"] <= expression_effect_rank_cutoff)
    )

    return df_ranked.drop(columns=["rank", "abs_expression_effect"])


def create_rank_response_table(
    config_dict: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Create a rank repsonse table from a dictionary which contains the
    configuration parameters. See docs at
    https://cmatkhan.github.io/callingCardsTools/file_format_specs/yeast_rank_response/ # noqa
    for details

    Args:
        config_dict (dict): A dictionary containing the configuration
            parameters

    Returns:
        tuple: A tuple containing three DataFrames
            (see rank_response_summarize):
               1. A dataframe summarized where hte responsive_ratio is summarized by
                    rank_bin
               2. The merged, labelled and sorted dataframe with both binding and
                    expression data
               3. The random expectation dataframe

    Raises:
        KeyError: if the configuration dictionary is missing any of the
            required keys
        FileExistsError: if the data files do not exist
        AttributeError: if there are NA values in the effect or pvalue columns
        ValueError: if there are incomplete cases in the data
    """
    # validate the configuration key/value pairs
    args = validate_config(config_dict)

    try:
        if len(args["binding_data_path"]) > 1:
            binding_data = combine_data(
                data_paths=args["binding_data_path"],
                identifier_col=args["binding_identifier_col"],
                effect_col=args["binding_effect_col"],
                pval_col=args["binding_pvalue_col"],
                source=args["binding_source"],
                data_type="binding",
            )
        else:
            binding_data = read_in_data(
                data_path=args["binding_data_path"][0],
                identifier_col=args["binding_identifier_col"],
                effect_col=args["binding_effect_col"],
                pval_col=args["binding_pvalue_col"],
                source=args["binding_source"],
                data_type="binding",
            )
    except (KeyError, FileExistsError, AttributeError) as exc:
        logger.error("Error reading in binding data: %s", exc)
        raise

    try:
        if len(args["expression_data_path"]) > 1:
            expression_data = combine_data(
                data_paths=args["expression_data_path"],
                identifier_col=args["expression_identifier_col"],
                effect_col=args["expression_effect_col"],
                pval_col=args["expression_pvalue_col"],
                source=args["expression_source"],
                data_type="expression",
            )
        else:
            expression_data = read_in_data(
                data_path=args["expression_data_path"][0],
                identifier_col=args["expression_identifier_col"],
                effect_col=args["expression_effect_col"],
                pval_col=args["expression_pvalue_col"],
                source=args["expression_source"],
                data_type="expression",
            )
    except (KeyError, FileExistsError, AttributeError) as exc:
        logger.error("Error reading in expression data: %s", exc)
        raise

    labeled_expression_data = label_responsive_genes(
        expression_data,
        args["expression_effect_thres"],
        args["expression_pvalue_thres"],
        args["normalization_cutoff"],
    )

    # Calculate counts for responsive and unresponsive
    responsive_unresponsive_counts = labeled_expression_data[
        "responsive"
    ].value_counts()

    # Create the DataFrame
    random_expectation_df = pd.DataFrame(
        {
            "unresponsive": [responsive_unresponsive_counts.get(False, 0)],
            "responsive": [responsive_unresponsive_counts.get(True, 0)],
        }
    )

    # Calculate the 'random' column
    total_expression_genes = random_expectation_df.sum(axis=1)
    random_expectation_df["random"] = (
        random_expectation_df["responsive"] / total_expression_genes
    )

    df = labeled_expression_data.merge(
        binding_data[["binding_effect", "binding_pvalue", "binding_source", "feature"]],
        how="inner",
        on="feature",
    )
    # test that there no incomplete cases. raise an error if there are
    if df.isnull().values.any():
        raise ValueError("There are incomplete cases in the data")

    logger.debug(
        "There are %s genes in the data after merging "
        "the %s binding data and "
        " %s expression data",
        str(df.shape[0]),
        args["binding_source"],
        args["expression_source"],
    )

    df_expression_labeled_binding_ranked = bin_by_binding_rank(
        df, args["rank_bin_size"], args["rank_by_binding_effect"]
    )

    df_expression_labeled_binding_ranked["random"] = random_expectation_df[
        "random"
    ].iloc[0]

    rank_response_df = compute_rank_response(df_expression_labeled_binding_ranked)

    return rank_response_df, df_expression_labeled_binding_ranked, random_expectation_df


def set_none_str_to_none(
    value: Union[str, None, int, float]
) -> Union[str, None, int, float]:
    """
    Test whether a string matches 'none' in any case. Return None if it does.
    Otherwise, return the original value

    Args:
        value (str, type(None), int, float): the string to test, or a value
            already set to None.

    Returns:
        str, type(None), int, float: the original value if it does not
            match 'none' in any case, or None if it does.

    Raises:
        TypeError: if the value is not a string, None or base python numeric.
    """
    if not isinstance(value, (str, type(None), int, float)):
        raise TypeError("value must be a string, None or base python numeric")

    none_pattern = r"(?i)^none$"

    # if the value is a string
    if isinstance(value, str):
        # test whether it matches 'none' in any case and return None if it does
        if bool(re.match(none_pattern, value)):
            return None
        # else, return the original value
        else:
            return value
    # if the value is not a string, it must be None, so return it
    return value


def validate_config(config: dict) -> dict:
    """
    Validate the yeast rank_response input configuration file.

    Args:
        config (dict): the configuration dictionary.

    Returns:
        dict: the validated configuration dictionary.

    Raises:
        KeyError: if the configuration is invalid due to either a missing
            key or an invalid value.
        TypeError: if the configuration is invalid due to an invalid type.
        FileNotFoundError: if the configuration is invalid due to a missing
    """
    # set default values if they are not in the config file
    config.setdefault("rank_by_binding_effect", False)
    config.setdefault("rank_bin_size", 5)
    config.setdefault("normalization_cutoff", -1)

    # this is used to check if a column is set to 'none' and replace it
    # with None. currently set for the expression_{effect/pvalue}_{col/thres}

    try:
        if not isinstance(config["binding_data_path"], list):
            raise TypeError("binding_data_path must be a list of strings")
        for path in config["binding_data_path"]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Binding data file {path} does not exist")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_data_path' in config") from exc

    try:
        config["binding_source"] = str(config["binding_source"])
    except KeyError as exc:
        raise KeyError("Missing key 'binding_source' in config") from exc

    try:
        if not isinstance(config["binding_identifier_col"], str):
            raise TypeError("binding_identifier_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_identifier_col' in config") from exc

    try:
        if not isinstance(config["binding_effect_col"], str):
            raise TypeError("binding_effect_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_effect_col' in config") from exc

    try:
        if not isinstance(config["binding_pvalue_col"], str):
            raise TypeError("binding_pvalue_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_pvalue_col' in config") from exc

    try:
        if not isinstance(config["rank_by_binding_effect"], bool):
            raise TypeError("rank_by_binding_effect must be a boolean")
    except KeyError as exc:
        raise KeyError("Missing key 'rank_by_binding_effect' in config") from exc

    try:
        if not isinstance(config["expression_data_path"], list):
            raise TypeError("expression_data_path must be a list of strings")
        for path in config["expression_data_path"]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Expression data file {path} does not exist")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_data_path' in config") from exc

    try:
        config["expression_source"] = str(config["expression_source"])
    except KeyError as exc:
        raise KeyError("Missing key 'expression_source' in config") from exc

    try:
        if not isinstance(config["expression_identifier_col"], str):
            raise TypeError("expression_identifier_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_identifier_col' in config") from exc

    try:
        if not isinstance(config["expression_effect_col"], (str, type(None))):
            raise TypeError("expression_effect_col must be a string or None")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_effect_col' in config") from exc

    try:
        if not isinstance(config["expression_effect_thres"], (int, float, type(None))):
            raise TypeError("expression_effect_thres must be numeric or None")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_effect_thres' in config") from exc

    for key in ["expression_effect_col", "expression_effect_thres"]:
        config[key] = set_none_str_to_none(config[key])

    if (
        config["expression_effect_col"] is None
        and config["expression_effect_thres"] is not None
    ) or (
        config["expression_effect_col"] is not None
        and config["expression_effect_thres"] is None
    ):
        raise KeyError(
            "expression_effect_thres must be None if " "expression_effect_col is None"
        )

    try:
        if not isinstance(config["expression_pvalue_col"], (str, type(None))):
            raise TypeError("expression_pvalue_col must be a string or None")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_pvalue_col' in config") from exc

    try:
        if not isinstance(config["expression_pvalue_thres"], (int, float, type(None))):
            raise TypeError("expression_pvalue_thres must be numeric or None")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_pvalue_thres' in config") from exc

    for key in ["expression_pvalue_col", "expression_pvalue_thres"]:
        config[key] = set_none_str_to_none(config[key])

    if (
        config["expression_pvalue_col"] is None
        and config["expression_pvalue_thres"] is not None
    ) or (
        config["expression_pvalue_col"] is not None
        and config["expression_pvalue_thres"] is None
    ):
        raise KeyError(
            "expression_pvalue_thres must be None if " "expression_pvalue_col is None"
        )

    if (
        config["expression_pvalue_col"] is None
        and config["expression_effect_col"] is None
    ):
        raise KeyError(
            "expression_pvalue_col and expression_effect_col " "cannot both be None"
        )

    try:
        if not isinstance(config["rank_bin_size"], int):
            raise TypeError("rank_bin_size must be an integer")
    except KeyError as exc:
        raise KeyError("Missing key 'rank_bin_size' in config") from exc

    try:
        if not isinstance(config["normalization_cutoff"], int):
            raise TypeError(
                "`normalization_cutoff` must be an integer >= -1. "
                + "Default is -1, which disables normalization"
            )
    except KeyError as exc:
        raise KeyError("Missing key 'normalization_cutoff' in config") from exc

    return config


def parse_args(
    subparser: argparse.ArgumentParser,
    script_desc: str,
    common_args: argparse.ArgumentParser,
) -> argparse.ArgumentParser:
    """
    Parse the command line arguments.

    :param subparser: the subparser object.
    :type subparser: argparse.ArgumentParser
    :param script_desc: the description of the script.
    :type script_desc: str
    :param common_args: the common arguments.
    :type common_args: argparse.ArgumentParser
    :return: the parser.
    :rtype: argparse.ArgumentParser
    """

    parser = subparser.add_parser(
        "yeast_rank_response",
        help=script_desc,
        prog="yeast_rank_response",
        parents=[common_args],
    )

    parser.set_defaults(func=main)

    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to the configuration json file. "
        "For details, see "
        "https://cmatkhan.github.io/callingCardsTools/file_format_specs/yeast_rank_response/",  # noqa
    )
    parser.add_argument(
        "--output_file",
        type=str,
        default="rank_response.csv",
        help="Path to the output file. Default is rank_response.csv",
    )
    parser.add_argument(
        "--compress",
        action="store_true",
        help="Compress the output file using gzip",
    )

    return subparser


def main(args: argparse.Namespace):
    # Load the JSON configuration file
    with open(args.config, "r", encoding="utf-8") as config_file:
        config_dict = json.load(config_file)

    try:
        config_dict = validate_config(config_dict)
    except (KeyError, TypeError, FileNotFoundError) as exc:
        logger.error("Error in configuration file: %s", exc)
        raise

    if os.path.exists(args.output_file):
        # warn that the file will be overwritten
        logger.warning(
            "The output file %s already exists and will be overwritten",
            args.output_file,
        )

    output_path = (
        args.output_file + ".gz"
        if args.compress and not args.output_file.endswith(".gz")
        else args.output_file
    )

    rank_response_df, labelled_binding_expression_df, random_expectation_df = (
        create_rank_response_table(config_dict)
    )

    compress = "gzip" if args.compress else None
    rank_response_df.to_csv(output_path, compression=compress, index=False)
