import argparse
import logging
import os
import re

from callingcardstools.Analysis.yeast.rank_response import (
    label_responsive_genes,
)
from callingcardstools.Analysis.yeast.read_in_data import read_in_data

logger = logging.getLogger(__name__)

__all__ = ["parse_args", "main"]


def find_min_responsive(
    data_path_list: list,
    identifier_col_list: list,
    effect_col_list: list,
    effect_thres_list: list,
    pval_col_list: list,
    pval_thres_list: list,
) -> int:
    """
    Finds the minimum number of responsive genes in a list of DataFrames.

    This function takes a list of DataFrames and finds the minimum number of
    responsive genes in any of the DataFrames. This is used to normalize the
    rank response across expression data sets.

    Args:
        data_path_list (list): A list of paths to expression dataframes
        identifier_col_list (list): A list of column names for the feature
            identifier in each DataFrame
        effect_col_list (list): A list of column names for the effect in each
            DataFrame. If there is no effect column in the dataframe at the
            same index, enter `None`
        effect_thres_list (list): A list of effect thresholds in each DataFrame
        pval_col_list (list): A list of column names for the p-value in each
            DataFrame. If no threshold is to be applied to the dataframe at
            the same index, enter `None`
        pval_thres_list (list): A list of p-value thresholds in each DataFrame.
            If no threshold is to be applied to the dataframe at the same
            index, enter `None`

    Returns:
        int: The minimum number of responsive genes in any of the DataFrames.

    Raises:
        TypeError: if data_path_list, identifier_col_list, effect_col_list, or
            pval_col_list is not a list. Also raised if there is an error
            in the `min()` statement after tallying number of responsive
            genes in each dataframe.
        ValueError: if the length of data_path_list, identifier_col_list,
            effect_col_list, or pval_col_list is not equal.
    """
    if not isinstance(data_path_list, list):
        raise TypeError("data_path_list must be a list")
    if not isinstance(identifier_col_list, list):
        raise TypeError("identifier_col_list must be a list")
    if not isinstance(effect_col_list, list):
        raise TypeError("effect_col_list must be a list")
    if not isinstance(pval_col_list, list):
        raise TypeError("pval_col_list must be a list")

    if (
        len(data_path_list)
        != len(identifier_col_list)
        != len(data_path_list)
        != len(effect_col_list)
        != len(data_path_list)
        != len(pval_col_list)
        != len(effect_thres_list)
        != len(pval_thres_list)
    ):
        raise ValueError(
            "Length of data_path_list, identifier_col_list, "
            "effect_col_list, and pval_col_list must be equal"
        )

    for data_path in data_path_list:
        if not os.path.exists(data_path):
            raise FileNotFoundError(f"{data_path} does not exist")

    for identifier_col in identifier_col_list:
        if not isinstance(identifier_col, str):
            raise TypeError("identifier_col_list must contain only strings")

    for effect_col in effect_col_list:
        if not isinstance(effect_col, (str, type(None))):
            raise TypeError("effect_col_list must contain only strings or " "`None`")

    for pval_col in pval_col_list:
        if not isinstance(pval_col, (str, type(None))):
            raise TypeError("pval_col_list must contain only strings or " "`None`")

    for effect_thres in effect_thres_list:
        if not isinstance(effect_thres, (int, float, type(None))):
            raise TypeError("effect_thres_list must contain only numbers or " "`None`")

    for pval_thres in pval_thres_list:
        if not isinstance(pval_thres, (int, float, type(None))):
            raise TypeError("pval_thres_list must contain only numbers or " "`None`")

    df_list = [
        label_responsive_genes(
            read_in_data(
                data_path, identifier_col, effect_col, pval_col, "source", "expression"
            ),
            effect_thres,
            pval_thres,
        )
        for data_path, identifier_col, effect_col, effect_thres, pval_col, pval_thres in zip(
            data_path_list,
            identifier_col_list,
            effect_col_list,
            effect_thres_list,
            pval_col_list,
            pval_thres_list,
        )
    ]

    try:
        min_responsive = min([df["responsive"].sum() for df in df_list])
    except TypeError as exc:
        logger.error("Error in `find_min_responsive()`: %s", exc)
        raise

    return min_responsive


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
        "yeast_min_responsive",
        help=script_desc,
        prog="yeast_min_responsive",
        parents=[common_args],
    )

    parser.set_defaults(func=main)

    parser.add_argument(
        "--data_path_list",
        nargs="+",
        type=list,
        help="A list of paths to expression dataframes",
        required=True,
    )
    parser.add_argument(
        "--identifier_col_list",
        nargs="+",
        type=list,
        help="A list of column names for the feature identifier " "in each DataFrame",
        required=True,
    )
    parser.add_argument(
        "--effect_col_list",
        nargs="+",
        type=list,
        help="A list of column names for the effect in each DataFrame",
        required=True,
    )
    parser.add_argument(
        "--effect_thres_list",
        nargs="+",
        type=list,
        help="A list of effect thresholds in each DataFrame. "
        "Enter `None` for no threshold on the effect for the dataframe "
        "at the same index",
        required=True,
    )
    parser.add_argument(
        "--pval_col_list",
        nargs="+",
        type=list,
        help="A list of column names for the p-value in each DataFrame",
        required=True,
    )
    parser.add_argument(
        "--pval_thres_list",
        nargs="+",
        type=list,
        help="A list of p-value thresholds in each DataFrame. "
        "Enter `None` for no threshold on the pvalue for the dataframe "
        "at the same index",
        required=True,
    )

    return subparser


def main(args: argparse.Namespace) -> None:
    """
    Find the minimum number of responsive genes in a list of DataFrames.

    This function takes a list of DataFrames and finds the minimum number of
    responsive genes in any of the DataFrames. This is used to normalize the
    rank response across expression data sets.

    Args:
        args (argparse.Namespace): The parsed command line arguments.

    Returns:
        None. prints the minimum number of responsive genes in any of the
        DataFrames to stdout.
    """
    none_pattern = r"(?i)^none$"

    try:
        effect_col_list = [
            None if bool(re.match(none_pattern, x)) else str(x)
            for x in args.effect_col_list
        ]
    except ValueError as exc:
        logger.error("Error in `find_min_responsive_main()`: %s", exc)
        raise

    try:
        effect_thres_list = [
            None if bool(re.match(none_pattern, x)) else float(x)
            for x in args.effect_thres_list
        ]
    except ValueError as exc:
        logger.error("Error in `find_min_responsive_main()`: %s", exc)
        raise

    try:
        pval_col_list = [
            None if bool(re.match(none_pattern, x)) else str(x)
            for x in args.pval_col_list
        ]
    except ValueError as exc:
        logger.error("Error in `find_min_responsive_main()`: %s", exc)
        raise

    try:
        pval_thres_list = [
            None if bool(re.match(none_pattern, x)) else float(x)
            for x in args.pval_thres_list
        ]
    except ValueError as exc:
        logger.error("Error in `find_min_responsive_main()`: %s", exc)
        raise

    try:
        min_responsive = find_min_responsive(
            args.data_path_list,
            args.identifier_col_list,
            effect_col_list,
            effect_thres_list,
            pval_col_list,
            pval_thres_list,
        )

        print(min_responsive)
    except (KeyError, TypeError, FileNotFoundError) as exc:
        logger.error("Error in `find_min_responsive_main()`: %s", exc)
        raise
