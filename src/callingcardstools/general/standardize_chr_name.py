"""Translate between chromosome naming conventions"""

# third party
import pandas as pd

def standardize_chr_names(df, chr_map_df, standard_chr_format, chr_colname = "chrom"):
    """Translate between chromosome naming conventions. Input frame must have 
    the `df[chr_colname]` levels contained entirely within one column of 
    chr_map_df.

    Args:
        df (Pandas DataFrame): dataframe with at least column chrom. all 
        entries in chrom field must be in one of the columns of chr_map_df
        chr_map_df (Pandas DataFrame): fields correspond to chromosome naming conventions, eg
        maybe refseq,ucsc,ensembl
        standard_chr_format (Str): the naming convention to which to translate 
        `df[chr_colname]`
        chr_colname (Str): name of the chromosome column in `df`. Default 'chrom'

    Raises:
        AttributeError: checks that all factor levels of the chrom col are 
        present in the factor levels of the 

    Returns:
        Pandas DataFrame: the input `df` with the `df[chrom_col]` translated 
        to `chr_map_df[standard_chr_format]`
    """

    # instantiate sentinel
    curr_chrom_format = -1
    # loop over colnames in chr_map_df. HALT if the current naming convention
    # is discovered
    for naming_convention in chr_map_df.columns:
        # check if all levels in the current chr_map_df[naming_convention]
        # contain the `df[chrom_colname]`.unique() levels
        if sum([True if x in chr_map_df[naming_convention].unique() else False \
            for x in df[chr_colname].unique()]) == len(df[chr_colname].unique()):
            curr_chrom_format = naming_convention
            break
    # if the current chromosome naming convention is not intuited above,
    # raise an error
    if curr_chrom_format == -1:
        raise AttributeError("Chromosome names are not "+\
            "recognized in a recognized format. Unique chr names which cause "+\
                " error are: %s." %",".join(df['chrom'].unique()))

    # if the current names are already in the requested format, return
    if curr_chrom_format == standard_chr_format:
        return df
    # else, use the chr_map_df to map to the standardized_chr_format convention
    else:
        return pd.merge(df,
                        chr_map_df.reindex(columns=[curr_chrom_format,
                                                    standard_chr_format]),
                        how = 'left',
                        left_on = 'chrom',
                        right_on = curr_chrom_format)\
                        .drop(['chrom', curr_chrom_format], axis=1)\
                        .rename(columns = {standard_chr_format:'chrom'})\
                        [df.columns]
