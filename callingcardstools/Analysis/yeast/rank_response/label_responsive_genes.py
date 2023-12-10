import logging

logger = logging.getLogger(__name__)


def label_responsive_genes(df,
                           abs_expression_effect_threshold,
                           expression_pvalue_threshold,
                           normalization_cutoff: int = -1):
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
    if 'expression_effect' not in df.columns:
        raise KeyError("Column 'effect_expression' is not in the data")
    if 'expression_pvalue' not in df.columns:
        raise KeyError("Column 'effect_pvalue' is not in the data")

    expression_effect_rank_cutoff = normalization_cutoff \
        if normalization_cutoff > 0 else len(df)+1

    df_abs = df.assign(abs_expression_effect=df['expression_effect'].abs())

    # if either the effect or p-value threshold is `None`, then set
    # the threshold to the appropiate boundary to prevent filtering on that
    # column
    abs_expression_effect_threshold = abs_expression_effect_threshold \
        if abs_expression_effect_threshold is not None \
        else min(df_abs['abs_expression_effect'])-1
    
    expression_pvalue_threshold = expression_pvalue_threshold \
        if expression_pvalue_threshold is not None \
        else max(df_abs['expression_pvalue'])+1

    df_ranked = (df_abs.sort_values(by=['abs_expression_effect',
                                    'expression_pvalue'],
                                    ascending=[False, True])
                 .reset_index(drop=True)
                 # Add 1 to start ranking from 1 instead of 0
                 .assign(rank=lambda x: x.index + 1))

    df_ranked['responsive'] = \
        ((df_ranked['abs_expression_effect'] > abs_expression_effect_threshold) # noqa
         & (df_ranked['expression_pvalue']
            < expression_pvalue_threshold)
         & (df_ranked['rank'] <= expression_effect_rank_cutoff))

    return df_ranked.drop(columns=['rank', 'abs_expression_effect'])
