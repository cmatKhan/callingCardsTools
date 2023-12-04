import logging
import pandas as pd
from .label_responsive_genes import label_responsive_genes
from .calculate_random_expectation import calculate_random_expectation
from .bin_by_binding_rank import bin_by_binding_rank
from .compute_rank_response import compute_rank_response

logger = logging.getLogger(__name__)


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
