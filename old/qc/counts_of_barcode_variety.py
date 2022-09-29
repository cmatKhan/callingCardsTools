
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

def barcode_qc(bed_df,
               barcode_details,
               fltr_bed_output_name):
    """parse the calling cards bed file name column into components and perform
    This expects a dataframe with AT LEAST the columns (there may be more):
    ['chrom','chromStart','chromEnd','barcode', 'strand', 'insert_seq']

    :param bed_df:
    :param barcode_details:
    :param fltr_bed_output_name:

    """

    # Using the keys of the barcode_components dict, make a dictionary
    # with structure {component1: [], component2:[], ...} where the keys are
    # barcode components, and the values are empty lists
    barcode_vector_dict = { key:[] for key in barcode_details['components'].keys()}

    # iterate over the barcodes in the 'name' column of the cc bed file
    # and parse the barcode into the barcode components
    for barcode in bed_df.loc[:,'barcode']:
        [barcode_vector_dict[k].append(barcode[v[0]:v[1]]) \
        for k,v in barcode_details['indicies'].items()]

    # TODO rename the function below to reflect the fact that both barcodes
    # and insertion sites may be handled

    # create position probability matricies of the barcode components
    barcode_components_ppm(barcode_vector_dict)
    barcode_components_ppm({'insert_seq': bed_df.loc[:,'insert_seq']})

    # transform barcode_vector_dict to a dataframe
    barcode_component_df = pd.DataFrame(data=barcode_vector_dict)

    # TODO rename the function below to reflect the fact that both barcodes
    # and insertion sites may be handled

    # create tables which describe the variety of barcode components (eg,
    # group by the srt column and count how many of each unique srt seq
    # there is)
    counts_of_barcode_variety(barcode_component_df)
    counts_of_barcode_variety(pd.DataFrame({'insert_seq': bed_df.loc[:,'insert_seq']}))

    # In the following section, create a vector to filter out
    # barcodes which do not meet barcode expectation OR!! the insert sequence
    # expectation

    # initialize a boolean vector with length equal to the number of rows
    # in the bed dataframe
    barcode_fltr_vector = [True]*len(bed_df)
    # iterate over the parsed barcode vector dataframe
    for index,row in barcode_component_df.iterrows():
        # initialize variable keep_index to True
        keep_index = True
        # iterate over the barcode components
        for barcode_component in barcode_vector_dict.keys():
            # if the current barcode component is not in the list of expected
            # values for that component
            if row[barcode_component] not in \
                barcode_details['components'][barcode_component]:
                keep_index = False
            # if there is a specific insert seq on which to filter, also
            # exclude reads that don't match it
            try:
                if not bed_df.loc[index,'insert_seq'] in barcode_details['insert_seq']:
                    keep_index = False
            except KeyError:
                pass
        # if keep_index if false,
        if not keep_index:
            barcode_fltr_vector[index] = False

    # reduce the cc bed file to only those rows which have expected
    # barcodes, and write out
    col_list = list(bed_df.columns)

    # filter the bed_df for those reads which did not pass the barcode_fltr
    fltr_bed_df = bed_df[barcode_fltr_vector]
    try:
        # try to read the field tf_map into a dataframe. if there is a keyerror,
        # handle in the except block
        barcode_map = pd.DataFrame(barcode_details['tf_map'])

        # group by TF. If there is more than one field to group on, add to the
        # list
        barcode_grouping_fields = ["TF"]

        # row filter, join, transform mismatches to 'other'
        barcode_to_tf_df = pd.merge(barcode_component_df,barcode_map,
def counts_of_barcode_variety(barcode_component_df):
    """Group and tally by each of the barcode components, and write the output
    to file

    :param barcode_component_df: a dataframe which has columns corresponding to
    the expected components of the barcode in the name column of the calling
    cards bed file
    """
    # iterate over columns in the input df, group by that column, tally
    # and then write to file
    for barcode_component in barcode_component_df.columns:
        tally_series = barcode_component_df\
                            .groupby(by=barcode_component)\
                            .size()\
                            .sort_values(ascending=False)
        tally_series.name = "tally"
        tally_series.to_csv(barcode_component+"_tally.tsv", sep = "\t")