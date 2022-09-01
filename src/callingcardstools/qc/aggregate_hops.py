
def aggregate_hops(df, coordinate_grouping_fields, qbed_col_order):
    """aggregate a qbed-style bed format dataframe. Note that this expects that
    the column with number of reads at a give spot be called 'reads'

    :param df: a qbed-style bed format dataframe
    :type df: pandas DataFrame
    :param coordinate_grouping_fields: list of fields by which to group
    :type coordinate_grouping_fields: list
    :param qbed_col_order: order of columns in return dataframe
    :type qbed_col_order: list

    :returns: an aggregated dataframe in expected qbed col order
    :rtype: pandas DataFrame
    """

    if not sum([True if x in df.columns else False \
        for x in coordinate_grouping_fields]) == len(coordinate_grouping_fields):
        ValueError('coordinate_grouping_fields not subset of df.columns')
    elif not sum([True if x in coordinate_grouping_fields else False \
        for x in qbed_col_order]) == len(qbed_col_order):
            ValueError('qbed_col_order not in coordinate_grouping_fields')

    agg_df = df\
        .groupby(coordinate_grouping_fields)['reads']\
        .agg(['sum'])\
        .reset_index()\
        .rename(columns={'sum':'reads'})\
        [qbed_COL_ORDER]

    return agg_df
