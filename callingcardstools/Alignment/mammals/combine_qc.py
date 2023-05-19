# pylint:disable=W1203
import os
import argparse
import logging

from callingcardstools.BarcodeParser.mammals.BarcodeQcCounter import BarcodeQcCounter  # noqa
from .Qbed import Qbed

__all__ = ['main']

logger = logging.getLogger(__name__)


def parse_args(
    subparser,
    script_description,
    common_args) -> argparse.ArgumentParser:  # noqa
    """This is intended to be used as a subparser for a parent parser passed
    from __main__.py. It adds the arguments required to combine mammals 
    Qbed and BarcodeQcCounter data that may result from splitting the reads 
    into multiple parts for parallel processing.

    Args:
        subparser (argparse.ArgumentParser): See __main__.py -- this is the
        subparser for the parent parser in __main__.py
        script_desc (str): Description of this script, which is set in
        __main__.py. The description is set in __main__.py so that all of
        the script descriptions are together in one spot and it is easier to
        write a unified cmd line interface
        common_args (argparse.ArgumentParser): These are the common arguments
        for all scripts in callingCardsTools, for instance logging level

    Returns:
        argparse.ArgumentParser: The subparser with the this additional
        cmd line tool added to it -- intended to be gathered in __main__.py
        to create a unified cmd line interface for the package
    """
    parser = subparser.add_parser(
        'mammals_combine_qc',
        help=script_description,
        prog='mammals_combine_qc',
        parents=[common_args]
    )

    # set the function to call when this subparser is used
    parser.set_defaults(func=main)

    input_group = parser.add_argument_group('input')
    input_group.add_argument(
        "-q",
        "--qbed_list",
        help="A list of pickled Qbed objects",
        nargs='+',
        required=True)
    input_group.add_argument(
        "-b",
        "--barcodeQcCounter_list",
        help="A list of pickled BarcodeQcCounter objects",
        nargs='+',
        required=True)

    parse_bam_output = parser.add_argument_group('output')
    parse_bam_output.add_argument(
        "-f",
        "--filename",
        help="base filename (no extension) for output files.",
        default="combined",
        required=False)
    parse_bam_output.add_argument(
        "-s",
        "--suffix",
        help="suffix to add to to the base filename. ",
        default="",
        required=False)
    parse_bam_output.add_argument(
        "-p",
        "--pickle",
        help="Set this flag to save the qbed and qc data as pickle files. "
        "this is useful when processing split files in parallel and then "
        "combining later. Defaults to False, which saves as qbed/tsv",
        action="store_true")

    return subparser


def main(args: argparse.Namespace) -> None:
    """Combine the Qbed and BarcodeQcCounter objects from multiple files
    into a single object and write to file"""

    for qbed_file in args.qbed_list:
        if not os.path.exists(qbed_file):
            raise FileNotFoundError(f"{qbed_file} does not exist")
        if os.path.splitext(qbed_file)[1] not in (".pkl"  ".pickle"):
            raise ValueError(f"{qbed_file} is not a pickled Qbed object")

    for barcodeqccounter_file in args.barcodeQcCounter_list:
        if not os.path.exists(barcodeqccounter_file):
            raise FileNotFoundError(f"{barcodeqccounter_file} does not exist")
        if os.path.splitext(barcodeqccounter_file)[1] not in (".pkl"  ".pickle"):
            raise ValueError(f"{barcodeqccounter_file} is not a "
                             "pickled BarcodeQcCounter object")

    qbed = Qbed()
    qbed_obj_list = [Qbed(qbed_file) for qbed_file in args.qbed_list]
    barcode_qc = BarcodeQcCounter()
    barcodeQcCounter_obj_list = [BarcodeQcCounter(barcodeQcCounter_file) for
                                 barcodeQcCounter_file in
                                 args.barcodeQcCounter_list]

    qbed = qbed.combine(qbed_obj_list)

    barcode_qc = barcode_qc.combine(barcodeQcCounter_obj_list)

    qbed.write(args.filename, args.suffix, args.pickle)
    barcode_qc.write(args.filename, args.suffix, args.pickle)

    return None
