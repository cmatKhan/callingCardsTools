"""
.. module:: callingcards_with_metrics
   :synopsis: Module for calling cards quantification functions.

This module contains functions for calculating various statistical values
related to the analysis of Calling Cards data. It includes functions for
computing Calling Cards effect (enrichment), Poisson p-value, and
hypergeometric p-value, as well as a function for processing and aggregating
data from multiple sources to obtain these values.

Functions
---------
- callingcards_with_metrics
- add_metrics
- parse_args
- main

.. author:: Chase Mateusiak
.. date:: 2023-11-23
"""
import argparse
import logging
import os
import time

import pandas as pd

from callingcardstools.PeakCalling.yeast import (read_in_background_data,
                                                 read_in_chrmap,
                                                 read_in_experiment_data,
                                                 read_in_promoter_data)
from callingcardstools.PeakCalling.yeast.enrichment_vectorized import \
    enrichment_vectorized
from callingcardstools.PeakCalling.yeast.hypergeom_pval_vectorized import \
    hypergeom_pval_vectorized
from callingcardstools.PeakCalling.yeast.poisson_pval_vectorized import \
    poisson_pval_vectorized

logger = logging.getLogger(__name__)


def count_hops(promoter_df: pd.DataFrame,
               qbed_df: pd.DataFrame,
               hop_colname: str,
               consider_strand: bool) -> pd.DataFrame:
    """
    Count the number of hops in the qbed_df for each promoter in the
    promoter_df.

    :param promoter_df: a pandas DataFrame of promoter regions.
    :type promoter_df: DataFrame
    :param qbed_df: a pandas DataFrame of qbed data.
    :type qbed_df: DataFrame
    :param hop_colname: the name of the column in the output DataFrame
        containing the number of hops.
    :type hop_colname: str
    :param consider_strand: whether to consider strand when counting hops.
    :type consider_strand: bool
    :return: a pandas DataFrame of promoter regions with a column containing
        the number of hops in the qbed_df for each promoter.
    :rtype: DataFrame

    :Example:

    >>> import pandas as pd
    >>> promoter_df = pd.DataFrame({
    ...     'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
    ...     'start': [100, 200, 300, 400, 500],
    ...     'end': [200, 300, 400, 500, 600],
    ...     'strand': ['+', '-', '+', '-', '+']
    ... })
    >>> qbed_df = pd.DataFrame({
    ...     'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
    ...     'start': [150, 250, 350, 450, 550],
    ...     'end': [200, 300, 400, 500, 600],
    ...     'depth': [1, 1, 1, 1, 1],
    ...     'strand': ['+', '-', '+', '-', '+']
    ... })
    >>> count_hops(promoter_df, qbed_df, 'hops', True)
         chr  start  end strand  hops
    0   chr1    100  200      +     1
    1   chr1    200  300      -     1
    2   chr1    300  400      +     1
    3   chr1    400  500      -     1
    4   chr1    500  600      +     1
    """
    if consider_strand:
        query_str = '(start <= qbed_start <= end) and strand == qbed_strand'
    else:
        # if consider_strand is false, then combine rows with the same
        # coordinates but different strand values and sum the depth. Set the
        # strand to "*" for all rows
        qbed_df = qbed_df\
            .groupby(['chr', 'start', 'end'])\
            .agg({'depth': 'sum'})\
            .reset_index()\
            .assign(strand='*')
        query_str = 'start <= qbed_start <= end'

    return promoter_df\
        .merge(qbed_df.rename(columns={'start': 'qbed_start',
                                       'end': 'qbed_end',
                                       'depth': 'qbed_hops',
                                       'strand': 'qbed_strand'}),
               on=['chr'], how='inner')\
        .query(query_str)\
        .drop(columns=['qbed_start',
                       'qbed_end',
                       'qbed_strand'])\
        .groupby(['chr', 'start', 'end', 'name', 'strand'])\
        .agg({'qbed_hops': 'count'})\
        .reset_index()\
        .rename(columns={'qbed_hops': hop_colname})


def call_peaks(
        experiment_data_path: str,
        experiment_orig_chr_convention: str,
        promoter_data_path: str,
        promoter_orig_chr_convention: str,
        background_data_path: str,
        background_orig_chr_convention: str,
        chrmap_data_path: str,
        consider_strand: bool,
        unified_chr_convention: str = 'ucsc') -> pd.DataFrame:
    """
    Call peaks for the given Calling Cards data.

    :param experiment_data_path: path to the experiment data file.
    :type experiment_data_path: str
    :param experiment_orig_chr_convention: the chromosome naming convention
        used in the experiment data file.
    :type experiment_orig_chr_convention: str
    :param promoter_data_path: path to the promoter data file.
    :type promoter_data_path: str
    :param promoter_orig_chr_convention: the chromosome naming convention
        used in the promoter data file.
    :type promoter_orig_chr_convention: str
    :param background_data_path: path to the background data file.
    :type background_data_path: str
    :param background_orig_chr_convention: the chromosome naming convention
        used in the background data file.
    :type background_orig_chr_convention: str
    :param chrmap_data_path: path to the chromosome map file.
    :type chrmap_data_path: str
    :param consider_strand: whether to consider strand when counting hops.
    :type consider_strand: bool
    :param unified_chr_convention: the chromosome naming convention
        to use in the output DataFrame.
    :type unified_chr_convention: str
    :return: a pandas DataFrame of promoter regions with Calling Cards
        metrics.
    :rtype: DataFrame
    """
    # read in the chr map
    chrmap_df = read_in_chrmap(
        chrmap_data_path,
        {experiment_orig_chr_convention,
         promoter_orig_chr_convention,
         background_orig_chr_convention,
         unified_chr_convention})

    # read in the experiment, promoter and background data
    experiment_df, experiment_total_hops = read_in_experiment_data(
        experiment_data_path,
        experiment_orig_chr_convention,
        unified_chr_convention,
        chrmap_df)
    promoter_df = read_in_promoter_data(promoter_data_path,
                                        promoter_orig_chr_convention,
                                        unified_chr_convention,
                                        chrmap_df)
    background_df, background_total_hops = read_in_background_data(
        background_data_path,
        background_orig_chr_convention,
        unified_chr_convention,
        chrmap_df)

    background_hops_df = count_hops(promoter_df,
                                    background_df,
                                    'background_hops',
                                    consider_strand)
    experiment_hops_df = count_hops(promoter_df,
                                    experiment_df,
                                    'experiment_hops',
                                    consider_strand)

    promoter_hops_df = experiment_hops_df\
        .merge(background_hops_df, on=['chr',
                                       'start',
                                       'end',
                                       'strand',
                                       'name'], how='left')\
        .fillna(0)\
        .assign(background_total_hops=background_total_hops,
                experiment_total_hops=experiment_total_hops)

    promoter_hops_df['background_hops'] = \
        promoter_hops_df['background_hops'].astype('int64')

    start_time = time.time()
    result_df = add_metrics(promoter_hops_df)
    logger.info("Time taken to process %s promoters: %s seconds",
                len(promoter_hops_df), time.time() - start_time)

    return result_df


def add_metrics(dataframe: pd.DataFrame,
                pseudocount: float = 0.2) -> pd.DataFrame:
    """
    Add Calling Cards metrics to the given DataFrame.

    :param dataframe: a pandas DataFrame of promoter regions.
    :type dataframe: DataFrame
    :param pseudocount: pseudocount to use when calculating Calling Cards
        metrics.
    :type pseudocount: float
    :return: a pandas DataFrame of promoter regions with Calling Cards
        metrics.
    :rtype: DataFrame
    """
    dataframe['callingcards_enrichment'] = enrichment_vectorized(
        dataframe['background_total_hops'],
        dataframe['experiment_total_hops'],
        dataframe['background_hops'],
        dataframe['experiment_hops'],
        pseudocount
    )

    dataframe['poisson_pval'] = poisson_pval_vectorized(
        dataframe['background_total_hops'],
        dataframe['experiment_total_hops'],
        dataframe['background_hops'],
        dataframe['experiment_hops'],
        pseudocount
    )

    dataframe['hypergeometric_pval'] = hypergeom_pval_vectorized(
        dataframe['background_total_hops'],
        dataframe['experiment_total_hops'],
        dataframe['background_hops'],
        dataframe['experiment_hops']
    )

    return dataframe


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
        'yeast_call_peaks',
        help=script_desc,
        prog='yeast_call_peaks',
        parents=[common_args]
    )

    parser.set_defaults(func=main)

    parser.add_argument(
        '--experiment_data_path',
        type=str,
        help='path to the experiment data file.',
        required=True
    )
    parser.add_argument(
        '--experiment_orig_chr_convention',
        type=str,
        help='the chromosome naming convention used in the experiment data '
             'file.',
        required=True
    )
    parser.add_argument(
        '--promoter_data_path',
        type=str,
        help='path to the promoter data file.',
        required=True
    )
    parser.add_argument(
        '--promoter_orig_chr_convention',
        type=str,
        help='the chromosome naming convention used in the promoter data '
             'file.',
        required=True
    )
    parser.add_argument(
        '--background_data_path',
        type=str,
        help='path to the background data file.',
        required=True
    )
    parser.add_argument(
        '--background_orig_chr_convention',
        type=str,
        help='the chromosome naming convention used in the background data '
             'file.',
        required=True
    )
    parser.add_argument(
        '--chrmap_data_path',
        type=str,
        help="path to the chromosome map file. this must include the data "
        "files' current naming conventions, the desired naming, and a column "
        "`type` that indicates whether the chromosome is 'genomic' or "
        "something else, eg 'mitochondrial' or 'plasmid'.",
        required=True
    )
    parser.add_argument(
        '--consider_strand',
        action='store_true',
        help='whether to consider strand when counting hops.'
    )
    parser.add_argument(
        '--unified_chr_convention',
        type=str,
        help='the chromosome naming convention to use in the output '
             'DataFrame.',
        required=False,
        default='ucsc'
    )
    parser.add_argument(
        '--output_path',
        default='sig_results.csv',
        type=str,
        help='path to the output file.'
    )
    parser.add_argument(
        '--pseudocount',
        type=float,
        help='pseudocount to use when calculating Calling Cards metrics.',
        required=False,
        default=0.2
    )
    parser.add_argument(
        '--compress_output',
        action='store_true',
        help='set this flag to gzip the output csv file.'
    )

    return subparser


def main(args: argparse.Namespace) -> None:
    """
    Call peaks for the given Calling Cards data.

    :param args: the command line arguments.
    :type args: Namespace
    """
    check_files = [args.experiment_data_path,
                   args.promoter_data_path,
                   args.background_data_path,
                   args.chrmap_data_path]
    for file in check_files:
        if not os.path.isfile(file):
            raise FileNotFoundError('The following path '
                                    f'does not exist: {file}')

    result_df = call_peaks(
        args.experiment_data_path,
        args.experiment_orig_chr_convention,
        args.promoter_data_path,
        args.promoter_orig_chr_convention,
        args.background_data_path,
        args.background_orig_chr_convention,
        args.chrmap_data_path,
        args.consider_strand,
        args.unified_chr_convention
    )

    result_df.to_csv(args.output_path,
                     compression='gzip' if args.compress_output else None,
                     index=False)
