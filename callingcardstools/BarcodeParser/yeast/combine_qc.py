# pylint:disable=C0206,W1514
import logging
import argparse

from .BarcodeQcCounter import BarcodeQcCounter
from ..BarcodeParser import BarcodeParser

__all__ = ['parse_args', 'combine_qc']

logger = logging.getLogger(__name__)


def parse_args(
        subparser: argparse.ArgumentParser,
        script_desc: str,
        common_args: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """This is intended to be used as a subparser for a parent parser passed
    from __main__.py. It adds the arguments required to combine 
    BarcodeQcCounter objects which may result from splitting the fastq files 
    prior to demultiplexing.

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
        'yeast_combine_qc',
        help=script_desc,
        prog='split_fastq',
        parents=[common_args]
    )

    parser.set_defaults(func=combine_qc)

    parser.add_argument('-i',
                        '--input',
                        help='a list of paths to BarcodeQcCounter object pickle files',
                        nargs='+',
                        required=True)

    parser.add_argument('-b',
                        '--barcode_details',
                        help='barcode filename (full path)',
                        required=True)

    parser.add_argument('-o',
                        '--output_dirpath',
                        help='a path to a directory where the output files '
                        'will be output. Defaults to the current directory',
                        default=".")

    parser.add_argument('-p',
                        '--prefix',
                        help='filename prefix for output files. Defaults to '
                        'barcode_qc',
                        default="barcode_qc")

    return subparser


def combine_qc(args):
    """
    Combine BarcodeQcCounter objects which may result from splitting the fastq
    files prior to demultiplexing.

    Args:
        args (argparse.Namespace): The cmd line arguments passed to this
        script

    Returns:
        None
    """
    bp = BarcodeParser(args.barcode_details)

    # TODO this is repeated code from split_fastq. Add to BarcodeParser?
    # construct the input to the BarcodeQcCounter summarize method
    component_dict = {k: [] for k in ['tf', 'r1_primer', 'r2_transposon']}
    r1_primer_start = bp.barcode_dict['r1']['primer']['index'][0]
    r1_primer_end = bp.barcode_dict['r1']['primer']['index'][1]
    r2_transposon_start = r1_primer_end + \
        bp.barcode_dict['r2']['transposon']['index'][0]
    r2_transposon_end = (r2_transposon_start +
                         bp.barcode_dict['r2']['transposon']['index'][1] -
                         bp.barcode_dict['r2']['transposon']['index'][0])
    for k, v in bp.barcode_dict['components']['tf']['map'].items():
        r1_primer_seq = k[r1_primer_start:r1_primer_end]
        r2_transposon_seq = k[r2_transposon_start:r2_transposon_end]
        component_dict['tf'].append(v)
        component_dict['r1_primer'].append(r1_primer_seq)
        component_dict['r2_transposon'].append(r2_transposon_seq)

    r1_bc_obj_list = [BarcodeQcCounter(x) for x in args.input]

    bc_combined = BarcodeQcCounter.combine(r1_bc_obj_list)

    bc_combined.write(component_dict=component_dict,
                      output_dirpath=args.output_dirpath,
                      filename=args.prefix)
