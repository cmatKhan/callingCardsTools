# pylint:disable=W1203
import os
import argparse
import logging
import tempfile

import pysam

from callingcardstools.Alignment.AlignmentTagger import AlignmentTagger
from callingcardstools.QC.create_status_coder import create_status_coder
from callingcardstools.BarcodeParser.mammals.BarcodeQcCounter \
    import BarcodeQcCounter
from .Qbed import Qbed

__all__ = ['main']

logger = logging.getLogger(__name__)


def parse_args(
    subparser,
    script_description,
    common_args) -> argparse.ArgumentParser:  # noqa
    """This is intended to be used as a subparser for a parent parser passed
    from __main__.py. It adds the arguments required to iterate over mammals
    alignments, set tags and create both a summary and the qbed (quantified
    hops file)

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
        'process_mammals_bam',
        help=script_description,
        prog='process_mammals_bam',
        parents=[common_args]
    )

    # set the function to call when this subparser is used
    parser.set_defaults(func=main)

    input_group = parser.add_argument_group('input')
    input_group.add_argument(
        "-i",
        "--input",
        help="path to bam file. Note that this must be "
        "sorted, and that an index .bai file must "
        "exist in the same directory",
        required=True)
    input_group.add_argument(
        "-b",
        "--barcode_details",
        help="path to the barcode details json",
        required=True)
    input_group.add_argument(
        "-g",
        "--genome",
        help="path to genome fasta file. "
        "Note that a index .fai must exist "
        "in the same directory",
        required=True)
    input_group.add_argument(
        "-q",
        "--mapq_threshold",
        help="Reads less than or equal to mapq_threshold "
        "will be marked as failed",
        type=int,
        default=10)

    parse_bam_output = parser.add_argument_group('output')
    parse_bam_output.add_argument(
        "-f",
        "--filename",
        help="Filename minus optional suffix "
        " and extension. Default is the input file basename "
        "minus the extension",
        default="",
        required=False)
    parse_bam_output.add_argument(
        "-s",
        "--suffix",
        help="suffix to add to output files.",
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


def process_chunk(bam_in: pysam.AlignmentFile,
                  barcode_details_path: str,
                  genome_path: str,
                  mapq_threshold: int,
                  **kwargs) -> dict:
    """This function is called when the subparser for this script is used.
    It parses the bam file, sets tags, and creates a summary and qbed file

    Args:
        bam_in (pysam.AlignmentFile): The bam file to parse
        contig (str): The contig to parse
        barcode_details_path (str): The path to the barcode details json
        genome_path (str): The path to the genome fasta file
        mapq_threshold (int): Reads less than or equal to mapq_threshold
        will be marked as failed
        kwargs (dict): Additional arguments passed to qbed.update

    Returns:
        dict: A dictionary with the following keys:
            'passing': A list of passing reads
            'failing': A list of failing reads
            'qbed': A Qbed object
            'barcode_qc': A BarcodeQcCounter object
    """
    output_dict = {
        'passing': [],
        'failing': [],
        'qbed': Qbed(),
        'barcode_qc': BarcodeQcCounter()
    }
    # create an AlignmentTagger object
    at = AlignmentTagger(barcode_details_path,
                         genome_path)  # pylint:disable=C0103
    status_coder = create_status_coder(at.insert_seq, mapq_threshold)

    # start iterating over the bam_chunk
    logger.debug("iterating over bam chunk...")
    for read in bam_in.fetch(until_eof=True):
        # parse the barcode, tag the read
        tagged_read = at.tag_read(read)
        # eval the read based on quality expectations, get the status
        status = status_coder(tagged_read)
        # add the data to the qbed and qc records
        output_dict['qbed'].update(tagged_read,
                                   status,
                                   **kwargs)
        if status == 0:
            # add the read to the passing_read list
            output_dict['passing'].append(tagged_read['read'])
        else:
            # add the read to the failing_read list
            output_dict['failing'].append(tagged_read['read'])
        # record barcode QC
        bc_counter_dict = {k: (v['query'], v['dist'] == 0) for
                           k, v in (tagged_read['barcode_details']
                                    ['details'].items())}
        bc_status = all([v[1] for v in bc_counter_dict.values()])
        output_dict['barcode_qc'].update(
            bc_counter_dict['r1_pb'][0],
            bc_counter_dict['r1_ltr1'][0],
            bc_counter_dict['r1_ltr2'][0],
            bc_counter_dict['r1_srt'][0],
            bc_status)

    return output_dict


def main(args: argparse.Namespace) -> None:

    logger.debug(args)

    # Check input paths
    logger.info('checking input...')
    input_path_list = [args.input,
                       args.barcode_details,
                       args.genome]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            error_msg = f"Input file DNE: {input_path}"
            logger.debug(error_msg)
            raise FileNotFoundError(error_msg)

    output_basename = args.filename if args.filename \
        else os.path.splitext(os.path.basename(args.input))[0]

    # open the bam file
    logger.info(f'beginning to parse {args.input}')
    bam_in = pysam.AlignmentFile(args.input)  # pylint:disable=E1101

    result_dict = process_chunk(bam_in,
                                args.barcode_details,
                                args.genome,
                                args.mapq_threshold)

    # write out
    with tempfile.TemporaryDirectory() as tmpdir:
        if len(result_dict['passing']) > 0:
            logger.info("writing unsorted passing bam...")
            # write the passing reads to a bam file
            tmp_passing_bam_output = os.path.join(tmpdir, 'passing.bam')
            with pysam.AlignmentFile(tmp_passing_bam_output,
                                     'wb',
                                     header=bam_in.header) as tmp_passing_bam:
                for read in result_dict['passing']:
                    tmp_passing_bam.write(read)
            # sort the passing and failing bam files with pysam
            logger.info("sorting passing bam...")

            try:
                passing_output_filename = \
                    output_basename + '_' + args.suffix + '_passing.bam' \
                    if args.suffix else output_basename + '_passing.bam'
                pysam.sort("-o", passing_output_filename,
                           tmp_passing_bam_output)
            except pysam.SamtoolsError as exc:
                raise pysam.SamtoolsError('Error sorting passing bam file') \
                    from exc

        else:
            logger.info(f"No passing reads found in {args.input}")

        if len(result_dict['failing']) > 0:

            logger.info("writing unsorted failing bam...")
            tmp_failing_bam_output = os.path.join(tmpdir, 'failing.bam')
            with pysam.AlignmentFile(tmp_failing_bam_output,
                                     'wb',
                                     header=bam_in.header) as tmp_failing_bam:
                for read in result_dict['failing']:
                    tmp_failing_bam.write(read)

            try:
                failing_output_filename = \
                    output_basename + '_' + args.suffix + '_failing.bam' \
                    if args.suffix else output_basename + '_failing.bam'
                logger.info("sorting failing bam...")
                pysam.sort("-o", failing_output_filename,
                           tmp_failing_bam_output)
            except pysam.SamtoolsError as exc:
                raise pysam.SamtoolsError('Error sorting failing bam file') \
                    from exc
        else:
            logger.info(f"No failing reads found in {args.input}")

    # write out
    result_dict['qbed'].write(output_basename, args.suffix, args.pickle)
    result_dict['barcode_qc'].write(output_basename, args.suffix, args.pickle)

    logger.info(f"file: {args.input} complete!")
