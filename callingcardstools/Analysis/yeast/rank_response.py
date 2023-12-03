import logging
import warnings
from typing import Literal
import os
import argparse
import json
from math import inf
from scipy.stats import binomtest
from scipy.stats._result_classes import BinomTestResult
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def read_in_data(
        data_path: str,
        identifier_col: str,
        effect_col: str,
        pval_col: str,
        source: str,
        data_type: Literal['binding', 'expression']) -> pd.DataFrame:
    """
    Read in data from a file and return a dataframe with
    the following columns: gene_id, {binding/expression}_effect,
    {binding/expression}_pvalue, source

    Args:
        data_path (str): path to the data file
        identifier_col (str): name of the feature identifier column in the data
        effect_col (str): name of the effect column in the data
        pval_col (str): name of the pvalue column in the data
        source (str): source of the data
        data_type (str): type of data, either 'binding' or 'expression'

    Returns:
        pd.DataFrame: dataframe with the following columns:
            feature, {binding/expression}_effect, {binding/expression}_pvalue,
            source

    Raises:
        FileExistsError: if data_path does not exist
        KeyError: if identifier_col, effect_col, or pval_col is not in the
            data, or if the `identifier_col` is something other than `feature`
            and the column `feature` also exists in the data
    """

    if not os.path.exists(data_path):
        raise FileExistsError(f"{data_path} does not exist")

    compressed = data_path.endswith('.gz')
    logger.debug("data compressed: %s", compressed)

    sep = '\t' if any(data_path.endswith(ext)
                      for ext in ['.tsv', '.txt', '.tsv.gz', 'txt.gz']) \
        else ','
    logger.debug("data separator: %s", sep)

    df = pd.read_csv(
        data_path,
        sep=sep,
        compression='gzip' if compressed else None)

    if identifier_col not in df.columns:
        raise KeyError(f"Column {identifier_col} is not in {data_path}")
    if 'feature' in df.columns and identifier_col != 'feature':
        raise KeyError(f"Column `feature` exists in the data, but is not the "
                       f"`identifier_col` {identifier_col}. Please rename the "
                       f"current `feature` column to avoid confusion.")

    try:
        effect_colname = data_type + '_effect'
        df[effect_colname] = inf \
            if effect_col == 'none' \
            else df[effect_col]
    except KeyError as exc:
        raise KeyError(f"Column {effect_col} is not `none` and "
                       "does not exist in {data_path}") from exc

    try:
        pval_colname = data_type + '_pvalue'
        df[pval_colname] = 0 \
            if pval_col == 'none' \
            else df[pval_col]
    except KeyError as exc:
        raise KeyError(f"Column {pval_col} is not `none` and "
                       "does not exist in {data_path}") from exc

    source_colname = data_type + '_source'
    df[source_colname] = source

    df.rename(columns={identifier_col: 'feature'}, inplace=True)

    return df[['feature', effect_colname, pval_colname, source_colname]]


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
    return np.concatenate([np.repeat(np.arange(1, quotient + 1), equal_parts),
                           np.repeat(quotient + 1, remainder)])


def label_responsive_genes(df,
                           expression_effect_threshold,
                           expression_pvalue_threshold,
                           normalization_cutoff: int = -1):
    """
    Labels genes in a DataFrame as responsive or not based on thresholds for
    expression effect and p-value.

    The function adds a new boolean column 'responsive' to the DataFrame, where
    each gene is labeled as responsive if its absolute effect expression is
    greater than a threshold and its p-value is less than a specified
    threshold. If normalization is enabled, only the top genes meeting the
    criteria up to the minimum number found in the normalized subset are
    labeled as responsive.

    Args:
        df (pd.DataFrame): DataFrame containing gene data. Must include
            'expression_effect' and 'expression_pvalue' columns.
        expression_effect_threshold (float): Threshold for the absolute value
            of the expression effect.
        expression_pvalue_threshold (float): Threshold for the expression
            p-value.
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
    if 'expression_effect' not in df.columns:
        raise KeyError("Column 'effect_expression' is not in the data")
    if 'expression_pvalue' not in df.columns:
        raise KeyError("Column 'effect_pvalue' is not in the data")

    expression_effect_rank_cutoff = normalization_cutoff \
        if normalization_cutoff > 0 else len(df)+1

    df_abs = df.assign(abs_expression_effect=df['expression_effect'].abs())

    df_ranked = (df_abs.sort_values(by=['abs_expression_effect',
                                    'expression_pvalue'],
                                    ascending=[False, True])
                 .reset_index(drop=True)
                 # Add 1 to start ranking from 1 instead of 0
                 .assign(rank=lambda x: x.index + 1))

    df_ranked['responsive'] = \
        ((df_ranked['abs_expression_effect'] >= expression_effect_threshold)
         & (df_ranked['expression_pvalue'] <= expression_pvalue_threshold)
         & (df_ranked['rank'] <= expression_effect_rank_cutoff))

    return df_ranked.drop(columns=['rank', 'abs_expression_effect'])


def calculate_random_expectation(df):
    """
    Calculates the random expectation of responsiveness in a DataFrame.

    This function takes a DataFrame that contains a 'responsive' column with
    boolean values.  It calculates the proportion of rows marked as responsive
    and unresponsive, and then computes the expected random proportion of
    responsiveness.

    Args:
        df (pd.DataFrame): A DataFrame containing at least a 'responsive'
        column with boolean values.

    Returns:
        pd.DataFrame: A DataFrame with columns 'unresponsive', 'responsive',
        and 'random', where 'unresponsive' and 'responsive' are counts of each
        category, and 'random' is the proportion of responsive rows over the
        total number of rows.

    Example:
        >>> df = pd.DataFrame({'responsive': [True, False, True, False]})
        >>> calculate_random_expectation(df)
        # Returns a DataFrame with the counts of responsive and unresponsive
        # rows and the proportion
        # of responsive rows.
    """
    # Calculate counts for responsive and unresponsive
    counts = df['responsive'].value_counts()

    # Create the DataFrame
    random_expectation_df = pd.DataFrame({
        'unresponsive': [counts.get(False, 0)],
        'responsive': [counts.get(True, 0)]
    })

    # Calculate the 'random' column
    total = random_expectation_df.sum(axis=1)
    random_expectation_df['random'] = \
        random_expectation_df['responsive'] / total

    return random_expectation_df


def bin_by_binding_rank(df: pd.DataFrame,
                        bin_size: int,
                        order_by_effect: bool = False):
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
        order_by_effect (bool, optional): If True, the DataFrame is sorted by
            abs('effect') in descending order first with ties broken by pvalue.
            If False, sort by pvalue first with ties broken by effect size.
            Defaults to False

    Returns:
        pd.DataFrame: The input DataFrame with an added 'rank' column, sorted
            by 'effect' in descending order and 'binding_pvalue' in
            ascending order.

    Example:
        >>> df = pd.DataFrame({'effect': [1.2, 0.5, 0.8], 
        ...                    'binding_pvalue': [5, 3, 4]})
        >>> bin_by_binding_rank(df, 2)
        # Returns a DataFrame with added 'rank' column and sorted as per
        # the specified criteria.
    """
    if 'binding_pvalue' not in df.columns:
        raise KeyError("Column 'binding_pvalue' is not in the data")
    if 'binding_effect' not in df.columns:
        raise KeyError("Column 'binding_effect' is not in the data")

    parts = min(len(df), bin_size)
    df_abs = df.assign(abs_binding_effect=df['binding_effect'].abs())

    df_sorted = df_abs.sort_values(
        by=['abs_binding_effect', 'binding_pvalue']
        if order_by_effect
        else ['binding_pvalue', 'abs_binding_effect'],
        ascending=[False, True]
        if order_by_effect
        else [True, False])

    return df_sorted\
        .drop(columns=['abs_binding_effect'])\
        .reset_index(drop=True)\
        .assign(rank_bin=create_partitions(len(df_sorted), parts) * parts)


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
    return (binomtest_obj.statistic,
            binomtest_obj.pvalue,
            binomtest_obj.proportion_ci(
                confidence_level=kwargs.get('confidence_level', 0.95),
                method=kwargs.get('method', 'exact')).low,
            binomtest_obj.proportion_ci(
                confidence_level=kwargs.get('confidence_level', 0.95),
                method=kwargs.get('method', 'exact')).high)


def compute_rank_response(df: pd.DataFrame, **kwargs):
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
    rank_response_df = df\
        .groupby('rank_bin')\
        .agg(
            n_responsive_in_rank=pd.NamedAgg(
                column='responsive', aggfunc='sum'),
            random=pd.NamedAgg(column='random', aggfunc='first'))\
        .reset_index()

    rank_response_df['n_successes'] = \
        rank_response_df['n_responsive_in_rank'].cumsum()

    # Binomial Test and Confidence Interval
    rank_response_df[['response_ratio', 'pvalue', 'ci_lower', 'ci_upper']] = \
        rank_response_df\
        .apply(lambda row: parse_binomtest_results(binomtest(
            int(row['n_successes']),
            int(row.rank_bin),
            float(row['random']),
            alternative=kwargs.get('alternative', 'two-sided')),
            **kwargs),
            axis=1, result_type='expand')

    return rank_response_df


def rank_response_ratio_summarize(
        df: pd.DataFrame,
        effect_expression_thres: float = 0,
        p_expression_thres: float = 0.05,
        normalize: bool = False,
        bin_size: int = 5) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    """
    Processes a DataFrame to compute and summarize rank response ratios.

    This function applies several processing steps on the input DataFrame,
    including labeling responsive genes, calculating random expectations,
    binning by binding rank, and computing rank responses. It returns three
    DataFrames containing various processed results.

    Args:
        df (pd.DataFrame): DataFrame to process.
        effect_expression_thres (float, optional): Threshold for effect
            expression. Defaults to 0.
        p_expression_thres (float, optional): Threshold for expression p-value.
            Defaults to 0.05.
        normalize (bool, optional): Whether to normalize the data. Defaults to
            False.
        bin_size (int, optional): Size of each bin for binding rank. Defaults
            to 5.

    Returns:
        tuple: A tuple containing three DataFrames:
               1. The input DataFrame with additional processing,
               2. A DataFrame of random expectations,
               3. A DataFrame of rank response calculations.

    Example:
        >>> test_df = pd.DataFrame({'gene_id': ['gene1', 'gene2', 'gene3'],
                                    'effect_expression': [0.5, -0.7, 1.2],
                                    'p_expression': [0.04, 0.07, 0.01],
                                    'binding_signal': [10, 20, 30]})
        >>> df, random_expectation_df, rank_response_df = \
        ...                  rank_response_ratio_summarize(test_df)
        >>> df.shape
        (3, x)  # x depends on the processing steps
        >>> random_expectation_df.shape
        (y, z)  # y and z depend on the structure of random expectations
        >>> rank_response_df.shape
        (a, b)  # a and b depend on the structure of rank response calculations
    """
    df_expression_labeled = label_responsive_genes(
        df,
        effect_expression_thres,
        p_expression_thres, normalize)

    random_expectation_df = calculate_random_expectation(df_expression_labeled)

    df_expression_labeled_binding_ranked = \
        bin_by_binding_rank(df_expression_labeled, bin_size)

    df_expression_labeled_binding_ranked_with_random = \
        df_expression_labeled_binding_ranked\
        .assign(random=float(random_expectation_df['random']))

    rank_response_df = compute_rank_response(
        df_expression_labeled_binding_ranked_with_random)

    return (df_expression_labeled_binding_ranked_with_random,
            random_expectation_df,
            rank_response_df)


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
        'yeast_rank_response',
        help=script_desc,
        prog='yeast_rank_response',
        parents=[common_args]
    )

    parser.set_defaults(func=main)

    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help="Path to the configuration json file. "
        "For details, see "
        "https://cmatkhan.github.io/callingCardsTools/file_format_specs/yeast_rank_response/"  # noqa
    )

    # parser.add_argument(
    #     '--binding_data_path',
    #     type=str,
    #     help='path to the binding data file. The `binding_effect_col`, '
    #     '`binding_pval_col`, and \'gene_id\' are required.',
    #     required=True
    # )
    # parser.add_argument(
    #     '--binding_identifier_col',
    #     type=str,
    #     help='name of the feature identifier column in the binding data',
    #     required=True
    # )
    # parser.add_argument(
    #     '--binding_effect_col',
    #     type=str,
    #     help='name of the effect column in the binding data. Set to '
    #     '`none` if an effect column does not exist',
    #     required=True
    # )
    # parser.add_argument(
    #     '--binding_pvalue_col',
    #     type=str,
    #     help='name of the pvalue column in the binding data. Set to `none` '
    #     'if a pvalue column does not exist'
    # )
    # parser.add_argument(
    #     '--rank_by_effect',
    #     action='store_true',
    #     help='The default is to rank by binding pvalue. A smaller pvalue '
    #     'indicates less probably that the effect value occurs by chance. '
    #     'If `--rank_by_effect` is set, then rank by the binding '
    #     'effect instead of binding pvalue. binding effect is '
    #     'expected to be a ratio of case/control where the minimum is 0 and '
    #     'larger values are more enriched.'
    # )
    # parser.add_argument(
    #     '--expression_data_path',
    #     type=str,
    #     help='path to the effect data file. The `expression_effect_col`, '
    #     '`expression_pval_col`, and \'gene_id\' are required.`',
    #     required=True
    # )
    # parser.add_argument(
    #     '--expression_identifier_col',
    #     type=str,
    #     help='name of the feature identifier column in the binding data',
    #     required=True
    # )
    # parser.add_argument(
    #     '--expression_effect_col',
    #     type=str,
    #     help='name of the effect column in the gene expression data. '
    #     'Set to `none` if an effect column does not exist',
    #     required=True
    # )
    # parser.add_argument(
    #     '--expression_effect_thres',
    #     type=float,
    #     help='threshold for effect expression',
    #     default=0
    # )
    # parser.add_argument(
    #     '--expression_pvalue_col',
    #     type=str,
    #     help='name of the pvalue column in the gene expression data. Set to '
    #     '`none` if a pvalue column does not exist'
    # )
    # parser.add_argument(
    #     '--expression_pvalue_thres',
    #     type=float,
    #     help='threshold for pvalue of effect expression',
    #     default=inf
    # )
    # parser.add_argument(
    #     '--rank_bin_size',
    #     type=int,
    #     help='bin size for rank response',
    #     default=5
    # )
    # parser.add_argument(
    #     '--normalize',
    #     action='store_true',
    #     help='normalize the number of responsive genes in each rank. '
    #     'This is not currently implemented -- it is a placeholder for '
    #     'future development when list input of binding/effect data is '
    #     'supported',
    #     default=False,
    #     hidden=True
    # )
    # parser.add_argument(
    #     '--output_file',
    #     type=str,
    #     default="rank_response.tsv",
    #     help='Path to the output file.'
    # )
    # parser.add_argument(
    #     '--compress',
    #     action='store_true',
    #     help='Set this flag to gzip the output file.'
    # )
    return subparser


def validate_config(config: dict) -> None:
    """
    Validate the yeast rank_response input configuration file.

    Args:
        config (dict): the configuration dictionary.

    Returns:
        None

    Raises:
        KeyError: if the configuration is invalid due to either a missing
            key or an invalid value.
        TypeError: if the configuration is invalid due to an invalid type.
        FileNotFoundError: if the configuration is invalid due to a missing
    """
    try:
        if not os.path.exists(config['binding_data_path']):
            raise FileNotFoundError(f"Binding data file "
                                    f"{config['binding_data_path']} "
                                    "does not exist")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_data_path' in config") from exc

    try:
        if not isinstance(config['binding_identifier_col'], str):
            raise TypeError("binding_identifier_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_identifier_col' in config") \
            from exc

    try:
        if not isinstance(config['binding_effect_col'], str):
            raise TypeError("binding_effect_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_effect_col' in config") \
            from exc

    try:
        if not isinstance(config['binding_pvalue_col'], str):
            raise TypeError("binding_pvalue_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_pvalue_col' in config") \
            from exc

    try:
        if not isinstance(config['rank_by_effect'], bool):
            raise TypeError("rank_by_effect must be a boolean")
    except KeyError as exc:
        raise KeyError("Missing key 'rank_by_effect' in config") \
            from exc

    try:
        if not os.path.exists(config['expression_data_path']):
            raise FileNotFoundError(f"Expression data file "
                                    f"{config['expression_data_path']} "
                                    "does not exist")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_data_path' in config") \
            from exc

    try:
        if not isinstance(config['expression_identifier_col'], str):
            raise TypeError("expression_identifier_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_identifier_col' in config") \
            from exc

    try:
        if not isinstance(config['expression_effect_col'], str):
            raise TypeError("expression_effect_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_effect_col' in config") \
            from exc

    try:
        if not isinstance(config['expression_effect_thres'], (int, float)):
            raise TypeError("expression_effect_thres must be numeric")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_effect_thres' in config") \
            from exc

    try:
        if not isinstance(config['expression_pvalue_col'], str):
            raise TypeError("expression_pvalue_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_pvalue_col' in config") \
            from exc

    try:
        if not isinstance(config['expression_pvalue_thres'], (int, float)):
            raise TypeError("expression_pvalue_thres must be numeric")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_pvalue_thres' in config") \
            from exc

    try:
        if not isinstance(config['rank_bin_size'], int):
            raise TypeError("rank_bin_size must be an integer")
    except KeyError as exc:
        raise KeyError("Missing key 'rank_bin_size' in config") \
            from exc

    try:
        if not isinstance(config['normalize'], bool):
            raise TypeError("normalize must be a boolean")
    except KeyError as exc:
        raise KeyError("Missing key 'normalize' in config") \
            from exc

    try:
        if not isinstance(config['output_file'], str):
            raise TypeError("output_file must be a string")
        if os.path.exists(config['output_file']):
            warnings.warn("The output file already exists. "
                          "It will be overwritten.")
        if not config['output_file'].endswith('.csv'):
            warnings.warn("We suggest that the output file end with .csv")
    except KeyError as exc:
        raise KeyError("Missing key 'output_file' in config") \
            from exc

    try:
        if not isinstance(config['compress'], bool):
            raise TypeError("compress must be a boolean")
    except KeyError as exc:
        raise KeyError("Missing key 'compress' in config") \
            from exc


def main(args: argparse.Namespace):
    # Load the JSON configuration file
    with open(args.config, 'r', encoding='utf-8') as config_file:
        config = json.load(config_file)

    # set default values if they are not in the config file
    config.setdefault('rank_by_effect', False)
    config.setdefault('rank_bin_size', 5)
    config.setdefault('normalize', False)
    config.setdefault('output_file', 'rank_response.csv')
    config.setdefault('compress', False)

    try:
        validate_config(config)
    except (KeyError, TypeError, FileNotFoundError) as exc:
        logger.error("Error in configuration file: %s", exc)
        raise

    # TODO: offer normalization method
    # if normalize:
    #     min_responsive = \
    #         df[(df['effect_expression'].abs() > effect_expression_thres)
    #            & (df['p_expression'] < p_expression_thres)].shape[0]
    # else:
    #     min_responsive = int('inf')

    # try:
    #     binding_data = read_in_data(
    #         args.binding_data_path,
    #         args.identifier_col,
    #         args.binding_effect_col,
    #         args.binding_pval_col,
    #         'binding',
    #         args.binding_source)
    # except (KeyError, FileExistsError) as exc:
    #     logger.error("Error reading in binding data: %s", exc)
    #     raise

    # try:
    #     expression_data = read_in_data(
    #         args.expression_data_path,
    #         args.expression_effect_col,
    #         args.expression_pval_col,
    #         'expression',
    #         args.expression_source)
    # except (KeyError, FileExistsError) as exc:
    #     logger.error("Error reading in expression data: %s", exc)
    #     raise

    # df = expression_data.merge(binding_data[['binding_pvalue', 'gene_id']],
    #                            how='left',
    #                            left_on='gene_id',
    #                            right_on='gene_id')

    # df, random_expectation_df, rank_response_df = \
    #     rank_response_ratio_summarize(
    #         df,
    #         effect_expression_thres=args.effect_expression_thres,
    #         p_expression_thres=args.p_expression_thres,
    #         normalize=args.normalize,
    #         bin_size=args.bin_size)

    # df.to_csv(args.output_prefix + '_rank_response.tsv', sep='\t', index=False)
    # random_expectation_df.to_csv(args.output_prefix +
    #                              '_random_expectation.tsv',
    #                              sep='\t',
    #                              index=False)
    # rank_response_df.to_csv(args.output_prefix +
    #                         '_rank_response.tsv',
    #                         sep='\t',
    #                         compression='gzip' if args.compress else None,
    #                         index=True)
