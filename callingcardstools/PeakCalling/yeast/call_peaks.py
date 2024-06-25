"""
.. module:: call_peaks
   :synopsis: Module for calling cards quantification functions.

This module contains functions for calculating various statistical values
related to the analysis of Calling Cards data. It includes functions for
computing Calling Cards effect (enrichment), Poisson p-value, and
hypergeometric p-value, as well as a function for processing and aggregating
data from multiple sources to obtain these values.

Functions
---------
- count_hops
- call_peaks
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
import pyranges as pr

from callingcardstools.PeakCalling.yeast import (
    read_in_background_data,
    read_in_chrmap,
    read_in_experiment_data,
    read_in_promoter_data,
)
from callingcardstools.PeakCalling.yeast.enrichment_vectorized import (
    enrichment_vectorized,
)
from callingcardstools.PeakCalling.yeast.hypergeom_pval_vectorized import (
    hypergeom_pval_vectorized,
)
from callingcardstools.PeakCalling.yeast.poisson_pval_vectorized import (
    poisson_pval_vectorized,
)

# from memory_profiler import profile


logger = logging.getLogger(__name__)


def count_hops(
    promoters_pr: pr.PyRanges,
    qbed_pr: pr.PyRanges,
    hops_colname: str,
    **kwargs,
) -> pd.DataFrame:
    """
    Use pyranges to join the promoter regions with the qbed data and count the
        number of qbed records that overlap with each promoter.

    additional keyword arguments are passed to the join method of the
      PyRanges object. Currently, the following are configured:
      - slack: which defaults to 0
      - suffix: which defaults to "_b"
      - strandedness: which defaults to False

    :param promoter_pr: a PyRanges of promoter regions.
    :type promoter_df: pr.PyRanges
    :param qbed_pr: a pandas DataFrame of qbed data from the
        experiment.
    :type qbed_pr: pr.PyRanges
    :param hops_colname: the name of the column in the qbed_df that
        contains the number of hops.

    :return: a pandas DataFrame of promoter regions with a column containing
        the number of hops in the qbed_df for each promoter.
    :rtype: DataFrame
    """
    overlaps = promoters_pr.join(
        qbed_pr,
        how="left",
        slack=kwargs.get("slack", 0),
        suffix=kwargs.get("suffix", "_b"),
        strandedness=kwargs.get("strandedness", False),
    )

    # Group by 'name' and count the number of records in each group
    # `observed` set to true b/c grouping is over categorical variable. This is default
    # in pandas 2.0. Without this set, memory usage skyrockets.
    # Setting "Start_b >= 0" to remove rows where there is no overlap, which are
    # represented by -1 in the _b columns by pyranges.
    overlap_counts = (
        overlaps.df.query("Start_b >= 0")
        .groupby("name", observed=True)
        .size()
        .reset_index(name="Count")
        .rename(columns={"Count": hops_colname})
    )

    return overlap_counts


def promoter_pyranges(
    promoter_df: pd.DataFrame,
    pyranges_rename_dict: dict = {
        "chr": "Chromosome",
        "start": "Start",
        "end": "End",
        "strand": "Strand",
    },
) -> pr.PyRanges:
    """
    Create a PyRanges object from the given promoter DataFrame.

    :param promoter_df: a pandas DataFrame of promoter regions.
    :type promoter_df: DataFrame
    :param pyranges_rename_dict: a dictionary that maps the column names in the
        promoter data to the column names in the PyRanges object. This is used
        to rename the columns in the PyRanges object after the promoter data
        is read in. The default is {"chr": "Chromosome", "start": "Start",
        "end": "End", "strand": "Strand"}.
    :return: a PyRanges object of promoter regions.
    :rtype: pr.PyRanges
    """
    promoters_pr = pr.PyRanges(
        promoter_df.rename(
            pyranges_rename_dict,
            axis=1,
        )
    )
    # extend the End by 1 bp to entries that start on the endpoint to be counted
    return promoters_pr.apply(lambda df: df.assign(End=df.End + 1))


# @profile
def call_peaks(
    experiment_data_path: str,
    experiment_orig_chr_convention: str,
    promoter_data_path: str,
    promoter_orig_chr_convention: str,
    background_data_path: str,
    background_orig_chr_convention: str,
    chrmap_data_path: str,
    unified_chr_convention: str = "ucsc",
    deduplicate_experiment: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """
    Call peaks for the given Calling Cards data.

    The kwargs parameter is used to pass additional arguments into underlying
    functions. Currently, the following are configured:
    - pranges_rename_dict: a dictionary that maps the column names in the
        promoter data to the column names in the PyRanges object. This is used
        to rename the columns in the PyRanges object after the promoter data
        is read in. The default is {"chr": "Chromosome", "start": "Start",
        "end": "End", "strand": "Strand"}.
    - join_validate: the validation method to use when joining the promoter
        data with the experiment and background data. The default is
        "one_to_one".
    - background_total_hops: the total number of hops in the background data.
        The default is the number of hops in the background data, calculated from
        the input background data file
    - experiment_total_hops: the total number of hops in the experiment data.
        The default is the number of hops in the experiment data, calculated from
        the input experiment data file
    - genomic_only: set this flag to include only genomic chromosomes in the
        experiment and background. See read_in_<experiment/background>_data for
        more details. Passed in kwargs so that if none is passed to add_metrics,
        the default in enrichment_vectorized() and poisson_pval_vectorized() is used.
    - pseudocount: pseudocount to use when calculating enrichment and poisson
        pvalue. See either function for more documentation. Passed through
        kwargs so that if none is passed to add_metrics, the default in
        enrichment_vectorized() and poisson_pval_vectorized() is used.

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
    :param deduplicate_experiment: If this is true, the experiment data will be
        deduplicated based on `chr`, `start` and `end` such that if an insertion
        is found at the same coordinate on different strands, only one of those records
        will be retained. see `read_in_experiment_data` for more details.
    :type deduplicate_experiment: bool
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
        {
            experiment_orig_chr_convention,
            promoter_orig_chr_convention,
            background_orig_chr_convention,
            unified_chr_convention,
        },
    )

    # read in the experiment, promoter and background data
    promoter_df = read_in_promoter_data(
        promoter_data_path,
        promoter_orig_chr_convention,
        unified_chr_convention,
        chrmap_df,
    )

    read_in_data_kwargs = {}
    if "genomic_only" in kwargs:
        read_in_data_kwargs["genomic_only"] = kwargs["genomic_only"]

    experiment_pr, experiment_total_hops = read_in_experiment_data(
        experiment_data_path,
        experiment_orig_chr_convention,
        unified_chr_convention,
        chrmap_df,
        deduplicate_experiment,
        **read_in_data_kwargs,
    )
    background_pr, background_total_hops = read_in_background_data(
        background_data_path,
        background_orig_chr_convention,
        unified_chr_convention,
        chrmap_df,
        **read_in_data_kwargs,
    )

    pyranges_rename_dict = kwargs.get(
        "pranges_rename_dict",
        {"chr": "Chromosome", "start": "Start", "end": "End", "strand": "Strand"},
    )

    promoters_pr = promoter_pyranges(promoter_df, pyranges_rename_dict)

    experiment_hops_df = count_hops(
        promoters_pr, experiment_pr, "experiment_hops"
    ).set_index("name", drop=True)

    background_hops_df = count_hops(
        promoters_pr, background_pr, "background_hops"
    ).set_index("name", drop=True)

    promoter_hops_df = (
        promoter_df.drop("score", axis=1)
        .set_index("name")
        .join(
            [experiment_hops_df, background_hops_df],
            how="left",
            validate=kwargs.get("join_validate", "one_to_one"),
        )
        .fillna(0)
        .assign(
            background_total_hops=kwargs.get(
                "background_total_hops", background_total_hops
            ),
            experiment_total_hops=kwargs.get(
                "experiment_total_hops", experiment_total_hops
            ),
        )
        .astype(
            {
                "background_hops": "int64",
                "experiment_hops": "int64",
                "background_total_hops": "int64",
                "experiment_total_hops": "int64",
            }
        )
        .reset_index()
    )

    # extract the add_metric kwargs if they are present
    add_metric_kwargs = {}
    if "pseudocount" in kwargs:
        add_metric_kwargs["pseudocount"] = kwargs["pseudocount"]

    start_time = time.time()
    result_df = add_metrics(promoter_hops_df, **add_metric_kwargs)
    logger.debug(
        "Time taken to process %s promoters: %s seconds",
        len(promoter_hops_df),
        time.time() - start_time,
    )

    return result_df


def add_metrics(dataframe: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """
    Add Calling Cards metrics to the given DataFrame.

    The kwargs parameter is used to pass additional arguments into underlying
    functions. Currently, the following are configured:
    - pseudocount: pseudocount to use when calculating both the enrichment and
        poisson pvalue. See either function for more documentation. Passed through
        kwargs so that if none is passed to add_metrics, the default in
        enrichment_vectorized() and poisson_pval_vectorized() is used.

    :param dataframe: a pandas DataFrame of promoter regions.
    :type dataframe: DataFrame

    :return: a pandas DataFrame of promoter regions with Calling Cards
        metrics.
    :rtype: DataFrame
    """
    dataframe["callingcards_enrichment"] = enrichment_vectorized(
        dataframe["background_total_hops"],
        dataframe["experiment_total_hops"],
        dataframe["background_hops"],
        dataframe["experiment_hops"],
        **kwargs,
    )

    dataframe["poisson_pval"] = poisson_pval_vectorized(
        dataframe["background_total_hops"],
        dataframe["experiment_total_hops"],
        dataframe["background_hops"],
        dataframe["experiment_hops"],
        **kwargs,
    )

    dataframe["hypergeometric_pval"] = hypergeom_pval_vectorized(
        dataframe["background_total_hops"],
        dataframe["experiment_total_hops"],
        dataframe["background_hops"],
        dataframe["experiment_hops"],
    )

    return dataframe


def parse_args(
    subparser: argparse.ArgumentParser,
    script_desc: str,
    common_args: argparse.ArgumentParser,
) -> argparse.ArgumentParser:
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
        "yeast_call_peaks",
        help=script_desc,
        prog="yeast_call_peaks",
        parents=[common_args],
    )

    parser.set_defaults(func=main)

    parser.add_argument(
        "--experiment_data_path",
        type=str,
        help="path to the experiment data file.",
        required=True,
    )
    parser.add_argument(
        "--experiment_orig_chr_convention",
        type=str,
        help="the chromosome naming convention used in the experiment data " "file.",
        required=True,
    )
    parser.add_argument(
        "--promoter_data_path",
        type=str,
        help="path to the promoter data file.",
        required=True,
    )
    parser.add_argument(
        "--promoter_orig_chr_convention",
        type=str,
        help="the chromosome naming convention used in the promoter data " "file.",
        required=True,
    )
    parser.add_argument(
        "--background_data_path",
        type=str,
        help="path to the background data file.",
        required=True,
    )
    parser.add_argument(
        "--background_orig_chr_convention",
        type=str,
        help="the chromosome naming convention used in the background data " "file.",
        required=True,
    )
    parser.add_argument(
        "--chrmap_data_path",
        type=str,
        help="path to the chromosome map file. this must include the data "
        "files' current naming conventions, the desired naming, and a column "
        "`type` that indicates whether the chromosome is 'genomic' or "
        "something else, eg 'mitochondrial' or 'plasmid'.",
        required=True,
    )
    parser.add_argument(
        "--unified_chr_convention",
        type=str,
        help="the chromosome naming convention to use in the output " "DataFrame.",
        required=False,
        default="ucsc",
    )
    parser.add_argument(
        "--deduplicate_experiment",
        action="store_true",
        help="set this flag to deduplicate the experiment data based on `chr`, "
        "`start` and `end` such that if an insertion is found at the same "
        "coordinate on different strands, only one of those records will be "
        "retained.",
    )
    parser.add_argument(
        "--genomic_only",
        action="store_true",
        help="set this flag to include only genomic chromosomes in the "
        "experiment and background.",
    )
    parser.add_argument(
        "--output_path",
        default="sig_results.csv",
        type=str,
        help="path to the output file.",
    )
    parser.add_argument(
        "--pseudocount",
        type=float,
        help=(
            "pseudocount to use when calculating poisson pvalue. Note that ",
            "this is used only when the background hops are 0 for a given promoter.",
        ),
        required=False,
        default=0.1,
    )
    parser.add_argument(
        "--compress_output",
        action="store_true",
        help="set this flag to gzip the output csv file.",
    )

    return subparser


def main(args: argparse.Namespace) -> None:
    """
    Call peaks for the given Calling Cards data.

    :param args: the command line arguments.
    :type args: Namespace
    """
    check_files = [
        args.experiment_data_path,
        args.promoter_data_path,
        args.background_data_path,
        args.chrmap_data_path,
    ]
    for file in check_files:
        if not os.path.isfile(file):
            raise FileNotFoundError("The following path " f"does not exist: {file}")

    # If the peak calling fails, catch the error, log it, and print an error
    # message to the console
    try:
        result_df = call_peaks(
            args.experiment_data_path,
            args.experiment_orig_chr_convention,
            args.promoter_data_path,
            args.promoter_orig_chr_convention,
            args.background_data_path,
            args.background_orig_chr_convention,
            args.chrmap_data_path,
            args.unified_chr_convention,
            args.deduplicate_experiment,
            genomic_only=args.genomic_only,
            pseudocount=args.pseudocount,
        )

        result_df.to_csv(
            args.output_path,
            compression="gzip" if args.compress_output else None,
            index=False,
        )
    except Exception as e:
        logger.error(
            "Error processing experiment file: %s. Error: %s",
            args.experiment_data_path,
            e,
        )
        raise
