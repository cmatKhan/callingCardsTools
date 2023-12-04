import logging
from math import inf
from scipy.stats import binomtest
import pandas as pd
from .parse_binomtest_results import parse_binomtest_results

logger = logging.getLogger(__name__)


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
