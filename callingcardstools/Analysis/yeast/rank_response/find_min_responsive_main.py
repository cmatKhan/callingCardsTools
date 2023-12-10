import logging
import argparse
import re
from .find_min_responsive import find_min_responsive

logger = logging.getLogger(__name__)

__all__ = ['parse_args', 'main']


def parse_args(
        subparser: argparse.ArgumentParser,
        script_desc: str,
        common_args: argparse.ArgumentParser) -> argparse.ArgumentParser:
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
        'yeast_min_responsive',
        help=script_desc,
        prog='yeast_min_responsive',
        parents=[common_args]
    )

    parser.set_defaults(func=main)

    parser.add_argument(
        '--data_path_list',
        nargs='+',
        type=list,
        help='A list of paths to expression dataframes',
        required=True
    )
    parser.add_argument(
        '--identifier_col_list',
        nargs='+',
        type=list,
        help='A list of column names for the feature identifier '
        'in each DataFrame',
        required=True
    )
    parser.add_argument(
        '--effect_col_list',
        nargs='+',
        type=list,
        help='A list of column names for the effect in each DataFrame',
        required=True
    )
    parser.add_argument(
        '--effect_thres_list',
        nargs='+',
        type=list,
        help='A list of effect thresholds in each DataFrame. '
        'Enter `None` for no threshold on the effect for the dataframe '
        'at the same index',
        required=True
    )
    parser.add_argument(
        '--pval_col_list',
        nargs='+',
        type=list,
        help='A list of column names for the p-value in each DataFrame',
        required=True
    )
    parser.add_argument(
        '--pval_thres_list',
        nargs='+',
        type=list,
        help='A list of p-value thresholds in each DataFrame. '
        'Enter `None` for no threshold on the pvalue for the dataframe '
        'at the same index',
        required=True
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
    none_pattern = r'(?i)^none$'

    try:
        effect_col_list = \
            [None if bool(re.match(none_pattern, x)) else str(x)
             for x in args.effect_col_list]
    except ValueError as exc:
        logger.error("Error in `find_min_responsive_main()`: %s", exc)
        raise

    try:
        effect_thres_list = \
            [None if bool(re.match(none_pattern, x)) else float(x)
             for x in args.effect_thres_list]
    except ValueError as exc:
        logger.error("Error in `find_min_responsive_main()`: %s", exc)
        raise

    try:
        pval_col_list = \
            [None if bool(re.match(none_pattern, x)) else str(x)
             for x in args.pval_col_list]
    except ValueError as exc:
        logger.error("Error in `find_min_responsive_main()`: %s", exc)
        raise

    try:
        pval_thres_list = \
            [None if bool(re.match(none_pattern, x)) else float(x)
             for x in args.pval_thres_list]
    except ValueError as exc:
        logger.error("Error in `find_min_responsive_main()`: %s", exc)
        raise

    try:
        min_responsive = find_min_responsive(args.data_path_list,
                                             args.identifier_col_list,
                                             effect_col_list,
                                             effect_thres_list,
                                             pval_col_list,
                                             pval_thres_list)

        print(min_responsive)
    except (KeyError, TypeError, FileNotFoundError) as exc:
        logger.error("Error in `find_min_responsive_main()`: %s", exc)
        raise
