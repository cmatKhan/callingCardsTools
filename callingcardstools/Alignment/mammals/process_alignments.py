# pylint:disable=W1203
import os
import argparse
import logging
import tempfile

import pysam

from callingcardstools.Alignment.AlignmentTagger import AlignmentTagger
from callingcardstools.QC.StatusFlags import StatusFlags
from callingcardstools.QC.create_status_coder import create_status_coder # noqa
from .ReadRecords import ReadRecords

__all__ = ['parse_bam']

logging.getLogger(__name__).addHandler(logging.NullHandler())


def parse_args(
    subparser, 
    script_description, 
    common_args) -> argparse.ArgumentParser: # noqa
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
    parser.set_defaults(func=parse_bam)

    parser = parser.add_argument_group('input')
    parser.add_argument(
        "-i",
        "--input",
        help="path to bam file. Note that this must be " +
        "sorted, and that an index .bai file must " +
        "exist in the same directory",
        required=True)
    parser.add_argument(
        "-b",
        "--barcode_details",
        help="path to the barcode details json",
        required=True)
    parser.add_argument(
        "-g",
        "--genome",
        help="path to genome fasta file. " +
        "Note that a index .fai must exist " +
        "in the same directory",
        required=True)

    parse_bam_output = parser.add_argument_group('output')
    parse_bam_output.add_argument(
        "-o",
        "--output_prefix",
        help="path to output directory. if not provided, " +
        "output to current working directory",
        default="",
        required=False)
    parse_bam_output.add_argument(
        "-qc",
        "--record_barcode_qc",
        help="set to False to avoid recording barcode qc " +
        "details while processing the bam. Defaults to True, but this " +
        "can take an enormous amount of memory. Set to false if there " +
        "is a memory problem.",
        default=True,
        type=bool,
        required=False)
    parse_bam_output.add_argument(
        "-t",
        "--tally",
        help="tally srt and multi srt sites. this creates two additional " +
        "qc summaries",
        default=True,
        type=bool,
        required=False)

    parse_bam_settings = parser.add_argument_group('settings')
    parse_bam_settings.add_argument(
        "-q",
        "--mapq_threshold",
        help="Reads less than or equal to mapq_threshold " +
        "will be marked as failed",
        type=int,
        default=10)

    return subparser


def parse_bam(args: argparse.Namespace) -> None:
    """This function is called when the subparser for this script is used.
    It parses the bam file, sets tags, and creates a summary and qbed file
    
    Args:
        args (argparse.Namespace): The arguments passed from the cmd line
        
        Raises:
            FileNotFoundError: If the input file does not exist
            
        Returns:
            None
    """
    # Check input paths
    logging.info('checking input...')
    input_path_list = [args.input,
                       args.barcode_details,
                       args.genome]
    if args.output_prefix:
        input_path_list.append(args.output_prefix)
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            error_msg = f"Input file DNE: {input_path}"
            logging.debug(error_msg)
            raise FileNotFoundError(error_msg)

    logging.info(f'beginning to parse {args.input}')
    output_bampath_dict = \
        {k: os.path.join(
            args.output_prefix,
            os.path.splitext(os.path.basename(args.input))[0]+'_'+k+'.bam')
            for k in ['passing', 'failing']}

    logging.info(args)

    qbed_output = os.path.join(
        args.output_prefix,
        os.path.splitext(os.path.basename(args.input))[0] + '.qbed')
    qc_output = os.path.join(
        args.output_prefix,
        os.path.splitext(os.path.basename(args.input))[0] + '_qc.tsv')

    logging.info(f"qbed output: {qbed_output}")
    logging.info(f"qc_output: {qc_output}")

    # open the bam file
    bam_in = pysam.AlignmentFile(args.input)  # pylint:disable=E1101

    # open output files for tagged bam reads
    output_bam_dict = \
        {k: pysam.AlignmentFile(v, "wb", header=bam_in.header)  
            for k, v in output_bampath_dict.items()}

    # create an AlignmentTagger object
    at = AlignmentTagger(args.barcode_details,
                         args.genome)  # pylint:disable=C0103

    # create a temp directory which will be destroyed when this context block
    # exits
    with tempfile.TemporaryDirectory() as tmp_dir:
        # create tmp files and use them to instantiate the ReadRecords object
        qbed_tmp_file = os.path.join(tmp_dir, "tmp_qbed.tsv")
        qc_tmp_file = os.path.join(tmp_dir, "tmp_qc.tsv")
        # instantiate ReadRecords object to handle creating qbed, qc files
        read_records = ReadRecords(qbed_tmp_file, qc_tmp_file)
        # create a status coder object to handle adding status codes
        # to the reads
        status_coder = create_status_coder(at.insert_seq, args.mapq_threshold)
        logging.info("iterating over bam file...")
        for read in bam_in.fetch():
            # only look at mapped primary alignments
            if not read.is_secondary and \
                not read.is_supplementary and \
                    not read.is_unmapped:
                # parse the barcode, tag the read
                tagged_read = at.tag_read(read)
                # eval the read based on quality expectations, get the status
                status = status_coder(tagged_read)
                # add the data to the qbed and qc records
                read_records\
                    .add_read_info(
                        tagged_read,
                        status,
                        insert_offset=at.insert_length,
                        annotation_tags=at.annotation_tags,
                        record_barcode_qc=args.record_barcode_qc)
                # write the read out to the appropriate tagged bamfile
                if status == 0:
                    output_bam_dict['passing'].write(tagged_read['read'])
                else:
                    output_bam_dict['failing'].write(tagged_read['read'])

        logging.info("writing qBed...")
        read_records.to_qbed(qbed_output, args.tally)
        if args.record_barcode_qc:
            logging.info("writing qc...")
            read_records.summarize_qc(qc_output)

    logging.info(f"file: {args.input} complete!")
