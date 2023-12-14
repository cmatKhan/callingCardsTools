import argparse
import logging
import sys
from importlib.metadata import version
from logging.config import dictConfig
from typing import Callable

from .Alignment.mammals import combine_qc as mammals_combine_qc
from .Alignment.mammals import process_alignments as process_mammals_bam
from .Alignment.yeast import legacy_makeccf
from .Alignment.yeast import process_alignments as process_yeast_bam
from .Analysis.yeast import chipexo_promoter_sig as yeast_chipexo_promoter_sig
from .Analysis.yeast.rank_response import rank_response_parse_args
from .BarcodeParser.yeast import barcode_table_to_json
from .BarcodeParser.yeast import combine_qc as yeast_combine_qc
from .PeakCalling.yeast import call_peaks as yeast_call_peaks
from .Reads import legacy_split_fastq, split_fastq


def parse_args() -> Callable[[list], argparse.Namespace]:
    """Create a cmd line argument parser for callingcardstools

    Returns:
            Callable[[list],Namespace]: This function returns the main
             argparser. If the subparser set_defaults method is set, it can
             make the correct decision as to which function to call to execute
             a given tool. See the main method below for usage
    """
    # create the top-level parser

    script_descriptions = {

        'barcode_to_json': 'parse a legacy mitra pipeline barcode table tsv '
        'to a barcode json file',

        'split_fastq': 'parse a (possibly multiplexed) batch of reads into '
        'expected barcode file(s), and a set of undetermined reads',

        'legacy_split_fastq': 'yeast cc_tools 3.0 version of parse_fastq',

        'legacy_makeccf': 'This function make .ccf files from mapped '
        '.bam files. ccf files have the following columns: '
        '[chr,start,end,reads,strand,barcode] but only the first 4 '
        'columns are required. The genome coordinates are 1-indexed',

        'yeast_combine_qc': 'Combine BarcodeQcCounter objects which may '
        'result from splitting the fastq files prior to demultiplexing',

        'process_yeast_bam': 'Iterate over yeast alignments to produce a '
        'qBed format file of the passing reads, and a qc file which '
        'allows finer exploration of the barcode and alignment metrics',

        'process_mammals_bam': 'Iterate over the reads in an alignment file '
        '(bam) and separate reads into passing.bam and failing.bam, a '
        'qBed format file of the passing reads, and a qc file which '
        'allows finer exploration of the barcode and alignment metrics',

        'mammals_combine_qc': 'Combine qbed and barcodeQC objects which may '
        'result from splitting the fastq and processing chunks in parallel',

        'yeast_call_peaks': 'Call peaks on yeast data',

        'yeast_chipexo_sig_promoter': 'call significant promoters from '
        'chipexo data from yeastepigenome.org',

        'yeast_find_min_responsive': 'Given a set of yeast expression data '
        'and thresholds on the effects and/or pvalues, find the minimum '
        'number of responsive genes in the data set given',

        'yeast_rank_response': 'Rank response analysis on yeast data'
    }

    # common options -- these can be applied to all scripts via the 'parent'---
    # argument -- see subparsers.add_parser for parse_bam below ---------------
    common_args = argparse.ArgumentParser(
        prog="callingcardstools", add_help=False)
    common_args_group = common_args.add_argument_group('General')
    common_args_group.add_argument(
        "-l",
        "--log_level",
        choices=("critical", "error", "warning", "info", "debug"),
        default="warning")

    # Create a top level parser -----------------------------------------------
    parser = argparse.ArgumentParser(
        prog='callingcardstools',
        description=f"callingcardstools: {version('callingcardstools')}")
    # add argument to get version
    parser.add_argument(
        "-v",
        "--version",
        action='version',
        version='%(prog)s '+f'{version("callingcardstools")}')

    # parse_bam subparser -----------------------------------------------------
    subparsers = parser.add_subparsers(
        help="Available Tools")

    subparsers = barcode_table_to_json.parse_args(
        subparsers,
        script_descriptions['barcode_to_json'],
        common_args
    )

    subparsers = legacy_split_fastq.parse_args(
        subparsers,
        script_descriptions['legacy_split_fastq'],
        common_args)

    subparsers = legacy_makeccf.parse_args(
        subparsers,
        script_descriptions['legacy_makeccf'],
        common_args)

    subparsers = split_fastq.parse_args(
        subparsers,
        script_descriptions['split_fastq'],
        common_args)

    subparsers = yeast_combine_qc.parse_args(
        subparsers,
        script_descriptions['yeast_combine_qc'],
        common_args)

    subparsers = process_yeast_bam.parse_args(
        subparsers,
        script_descriptions['process_yeast_bam'],
        common_args
    )

    subparsers = process_mammals_bam.parse_args(
        subparsers,
        script_descriptions['process_mammals_bam'],
        common_args
    )

    subparsers = mammals_combine_qc.parse_args(
        subparsers,
        script_descriptions['mammals_combine_qc'],
        common_args)

    subparsers = yeast_call_peaks.parse_args(
        subparsers,
        script_descriptions['yeast_call_peaks'],
        common_args)
    
    subparsers = yeast_chipexo_promoter_sig.parse_args(
        subparsers,
        script_descriptions['yeast_chipexo_sig_promoter'],
        common_args)

    subparsers = rank_response_parse_args(
        subparsers,
        script_descriptions['yeast_rank_response'],
        common_args)

    # return the top level parser to be used in the main method below
    return parser


def main(args=None) -> None:
    """Entry point to callingcardstools"""

    # parse the cmd line arguments
    arg_parser = parse_args()

    args = arg_parser.parse_args(args)

    # this is a default setting -- if it is not set, it means
    # that nothing was passed on the cmd line. Instead, print the
    # help message
    try:
        log_level = args.log_level.upper()
    except AttributeError:
        sys.exit(arg_parser.print_help())

    # set the logging details
    log_config = {
        "version": 1,
        "handlers": {
            "console": {
                "formatter": "std_out",
                "class": "logging.StreamHandler"
            }
        },
        "formatters": {
            "std_out": {
                "format": "%(asctime)s : %(module)s : "
                "%(funcName)s : line: %(lineno)d\n" +
                "\tprocess details : %(process)d, %(processName)s\n" +
                "\tthread details : %(thread)d, %(threadName)s\n" +
                "\t%(levelname)s : %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S"
            }
        },
        "loggers": {
            "": {
                "level": f"{log_level}",
                "handlers": ["console"]
            }
        }
    }
    dictConfig(log_config)
    # log the cmd line arguments at the debug level
    logging.debug(sys.argv)
    logging.debug(str(args))

    # note that this works b/c the subparser set_defaults function attribute
    # is set.
    # see https://docs.python.org/3/library/argparse.html#parser-defaults
    # scroll up from that point to see a usage example
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
