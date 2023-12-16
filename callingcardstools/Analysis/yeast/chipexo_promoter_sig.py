import argparse
import logging
import os

import pandas as pd

from callingcardstools.PeakCalling.yeast import (read_in_chrmap,
                                                 read_in_promoter_data,
                                                 relabel_chr_column)

logger = logging.getLogger(__name__)


def read_in_chipexo_data(
        chipexo_data_path: str,
        curr_chr_convention: str,
        new_chr_convention: str,
        chrmap_df: pd.DataFrame) -> pd.DataFrame:
    """
    Read in the data from the chipexo file. This is data parsed from
    yeastepigenome.org. see yeastepigenome.org and
    https://github.com/cmatKhan/parsing_yeast_database_data 

    Args:
            chipexo_allevents_data_path (str): path to the chipexo data
            chipexo_orig_chr_convention (str): chromosome convention of the
                chipexo allevents file
            unified_chr_convention (str): chromosome convention to convert to

    Returns:
            pandas.DataFrame: A pandas DataFrame containing the chipexo
                allevents data

    Raises:
            AttributeError: If the chipexo table does not contain at least the
                following columns: `chr`, `start`, `end`, `YPD_log2Fold`,
                `YPD_log2P`. Note that the `start` column is the original
                `coord` column from the yeastepigenome.org data and `end` is
                simply `coord` + 1. It is in this format to make it somewhat
                easier to input to other processes that accept bed-like files.
    """
    df = pd.read_csv(chipexo_data_path,
                     header=0,
                     index_col=False)

    if not {'chr', 'start', 'end',
            'YPD_log2Fold', 'YPD_log2P'}.issubset(df.columns):
        raise AttributeError('The chipexo table must contain at least the '
                             'following columns: `chr`, `start`, `end`, '
                             '`YPD_log2Fold`, `YPD_log2P`. Note that the '
                             '`start` column is the original `coord` column '
                             'from the yeastepigenome.org data and `end` '
                             'is simply `coord` + 1. It is in this format '
                             'to make it somewhat easier to input to other '
                             'processes that accept bed-like files.')

    df.rename(columns={'start': 'chipexo_start',
                       'end': 'chipexo_end'},
              inplace=True)

    return relabel_chr_column(df,
                              chrmap_df,
                              curr_chr_convention,
                              new_chr_convention)


def chipexo_promoter_sig(chipexo_data_path: str,
                         chipexo_orig_chr_convention: str,
                         promoter_data_path: str,
                         promoter_orig_chr_convention: str,
                         chrmap_data_path: str,
                         unified_chr_convention: str) -> pd.DataFrame:
    """
    Find the promoter signature of the chipexo data. This is calculated as
    the most significant peak in each promoter region.

    Args:
        chipexo_data_path (str): path to the chipexo allevents file.
        chipexo_orig_chr_convention (str): chromosome convention of the
            chipexo allevents file.
        promoter_data_path (str): path to the promoter data file.
        promoter_orig_chr_convention (str): chromosome convention of the
            promoter data file.
        chrmap_data_path (str): path to the chromosome map file.
        unified_chr_convention (str): chromosome convention to convert to.

    Returns:
        pandas.DataFrame: A pandas DataFrame containing the promoter
            signature of the chipexo data.

    Example:
        >>> import pandas as pd
        >>> import tempfile
        >>> # Create temporary chipexo data file
        >>> with tempfile.NamedTemporaryFile(mode='w+',
        ...                                  suffix='.tsv') as chipexo_file:
        ...     _ = chipexo_file.write('chr\\tcoord\\tYPD_log2Fold\\t'
        ...                        ' YPD_log2P\\nchr1\\t150\\t2.0\\t0.05\\n')
        >>> # Create temporary promoter data file
        >>> with tempfile.NamedTemporaryFile(mode='w+',
        ...                                  suffix='.tsv') as promoter_file:
        ...     _ = promoter_file.write('chr\\tstart\\tend\\t'
        ...                         'associated_feature\\nchr1\\t100\\t'
        ...                         '200\\tpromoter1\\n')
        >>> # Create temporary chromosome map file
        >>> with tempfile.NamedTemporaryFile(mode='w+',
        ...                                  suffix='.tsv') as chrmap_file:
        ...     - = chrmap_file.write('chr\\tucsc\\nchr1\\tchr1\\n')
        >>> # Call the function
        >>> result = chipexo_promoter_sig(chipexo_file.name, 'chr',
        ...                               promoter_file.name, 'chr',
        ...                               chrmap_file.name, 'ucsc')
        >>> isinstance(result, pd.DataFrame)
        True
    """
    # read in chrmap data
    chrmap_df = read_in_chrmap(chrmap_data_path,
                               {chipexo_orig_chr_convention,
                                promoter_orig_chr_convention,
                                unified_chr_convention})
    # read in promoter data
    promoter_df = read_in_promoter_data(promoter_data_path,
                                        promoter_orig_chr_convention,
                                        unified_chr_convention,
                                        chrmap_df)
    # read in chipexo data
    chipexo_df = read_in_chipexo_data(chipexo_data_path,
                                      chipexo_orig_chr_convention,
                                      unified_chr_convention,
                                      chrmap_df)

    # Step 1: Inner Join
    return pd.merge(promoter_df, chipexo_df,
                    on='chr',
                    how='inner')\
        .query('start <= chipexo_start <= end')\
        .groupby(['chr', 'start', 'end', 'name', 'strand'])\
        .agg(
        n_sig_peaks=pd.NamedAgg(column='chr',
                                aggfunc='count'),
        max_fc=pd.NamedAgg(column='YPD_log2Fold',
                           aggfunc='max'),
        min_pval=pd.NamedAgg(column='YPD_log2P',
                             aggfunc='min'))\
        .reset_index()


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
        'chipexo_promoter_sig',
        help=script_desc,
        prog='chipexo_promoter_sig',
        parents=[common_args]
    )

    parser.set_defaults(func=main)

    parser.add_argument(
        '--chipexo_data_path',
        help='Path to the chipexo data file.',
        required=True
    )
    parser.add_argument(
        '--chipexo_orig_chr_convention',
        help='Chromosome convention of the chipexo data file.',
        required=True
    )
    parser.add_argument(
        '--promoter_data_path',
        help='Path to the promoter data file.',
        required=True
    )
    parser.add_argument(
        '--promoter_orig_chr_convention',
        help='Chromosome convention of the promoter data file.',
        required=True
    )
    parser.add_argument(
        '--chrmap_data_path',
        help='Path to the chromosome map file.',
        required=True
    )
    parser.add_argument(
        '--unified_chr_convention',
        help='Chromosome convention to convert to.',
        required=True
    )
    parser.add_argument(
        '--output_file',
        default="chipexo_promoter_sig.csv",
        help='Path to the output file.'
    )
    parser.add_argument(
        '--compress',
        action='store_true',
        help='Set this flag to gzip the output file.'
    )

    return subparser


def main(args: argparse.Namespace) -> None:
    """
    Given the allevents file from yeastepigenome.org, which is the (parsed) 
    output from chexmix, find the most significant peak in each interval,
    if a peak exists, in a set of promoter regions

    Args:
        args (argparse.Namespace): arguments from command line
    """
    check_files = [args.chipexo_data_path,
                   args.promoter_data_path,
                   args.chrmap_data_path]
    for file in check_files:
        if not os.path.isfile(file):
            raise FileNotFoundError(f'{file} not found')

    result = chipexo_promoter_sig(args.chipexo_data_path,
                                  args.chipexo_orig_chr_convention,
                                  args.promoter_data_path,
                                  args.promoter_orig_chr_convention,
                                  args.chrmap_data_path,
                                  args.unified_chr_convention)

    result.to_csv(args.output_file,
                  compression='gzip' if args.compress else None,
                  index=False)
