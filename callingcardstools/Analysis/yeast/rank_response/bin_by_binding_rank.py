import logging
import pandas as pd
from .create_partitions import create_partitions

logger = logging.getLogger(__name__)


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