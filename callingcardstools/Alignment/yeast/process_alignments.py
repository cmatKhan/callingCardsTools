# pylint:disable=W1203
import logging
import tempfile
import os
import argparse

# outside dependencies
import pysam
import pandas as pd

# from memory_profiler import profile
# local dependencies
from callingcardstools.Alignment.AlignmentTagger import AlignmentTagger
from callingcardstools.QC.create_status_coder import create_status_coder  # noqa
from callingcardstools.Alignment.SummaryParser import SummaryParser
from callingcardstools.QC.StatusFlags import StatusFlags

__all__ = ['parse_args', 'process_alignments']

logger = logging.getLogger(__name__)


def parse_args(
        subparser: argparse.ArgumentParser,
        script_desc: str,
        common_args: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """This is intended to be used as a subparser for a parent parser passed 
    from __main__.py. It adds the arguments required to iterate over yeast 
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
        'process_yeast_bam',
        help=script_desc,
        prog='process_yeast_bam',
        parents=[common_args]
    )

    parser.set_defaults(func=process_alignments)

    parser.add_argument("-i",
                        "--bampath",
                        help="path to the input bam file",
                        required=True)

    parser.add_argument("-g",
                        "--genome",
                        help=" ".join(["Path to a genome .fasta file.",
                                       "Note that an index .fai file must exist in the same path"]),
                        required=True)

    parser.add_argument("-j",
                        "--barcode_details",
                        help="Path to the barcode details json file",
                        required=True)

    parser.add_argument("-q",
                        "--mapq_threshold",
                        help="",
                        default=10,
                        type=int)

    parser.add_argument("-o",
                        "--output_dir",
                        help="path to the output directory"
                        " (default: current directory)",
                        default='.')

    parser.add_argument("-v",
                        "--verbose_qc",
                        help="save complete alignment summary",
                        action="store_true")

    return subparser


def process_alignments(args: argparse.Namespace) -> dict:
    """Iterate over a bam file, set tags and output updated bam with
     read groups added to the header, tags added to the reads.
     Also output a summary of the reads

    Args:
        args (argparse.Namespace): A argparse namespace object with 
        the following attributes:
        bampath (str): path to the alignment file (bam)
        insertion_length (int): Expected length of the insertion sequence
        barcode_details_json (str): Path to the barcode details json file
        genome(str): Path to a fasta file
        mapq_threshold (int, optional): mapq threshold below which to label a
         read as failing. Defaults to 10.
        out_suffix (str, optional): suffix to append to the augmented bam file
         output. Defaults to "_tagged.bam".
        nthreads (int): Number of threads which pysam.AlignmentFile may use to
         decompress lines
        output_dir (str): Path to the directory where the output files will be
            written


    Raises:
        FileNotFoundError: If one of the input files does not exist

    Returns:
        dict: A dictionary with keys summary and qbed and values the 
        corresponding dataframes
    """
    # Check inputs
    try:
        input_path_list = [args.bampath,
                           args.bampath + '.bai',
                           args.genome,
                           args.genome + '.fai',
                           args.barcode_details]
    except AttributeError as exc:
        raise AttributeError("Input file paths not specified") from exc

    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" % input_path)

    # set defaults
    try:
        mapq_threshold = args.mapq_threshold
    except AttributeError:
        mapq_threshold = 10

    try:
        out_suffix = args.out_suffix
    except AttributeError:
        out_suffix = '_tagged.bam'

    try:
        nthreads = args.nthreads
    except AttributeError:
        nthreads = 1

    try:
        output_dir = args.output_dir
    except AttributeError:
        output_dir = os.getcwd()

    output_basename = os.path.splitext(os.path.basename(args.bampath))[0]

    logger.info("tagging reads...")
    # temp_dir is automatically cleaned when context ends
    with tempfile.TemporaryDirectory() as temp_dir:
        # open files
        logger.info("opening input bam file: %s", args.bampath)
        input_bamfile = pysam.AlignmentFile(  # pylint:disable=E1101
            args.bampath, "rb",
            require_index=True,
            threads=nthreads)

        # create tmp output bam files
        logger.debug('creating tmp passing and failing bam in '
                      'temp dir: %s', temp_dir)
        tmp_output_bampath_dict = \
            {x: os.path.join(
                temp_dir,
                'tmp_' + output_basename + '_' + x + out_suffix)
                for x in ['passing', 'failing']}

        tmp_tagged_bam_dict = \
            {k: pysam.AlignmentFile(v, "wb", header=input_bamfile.header)
             for k, v in tmp_output_bampath_dict.items()}

        at = AlignmentTagger(args.barcode_details, args.genome)

        status_coder = create_status_coder(
            mapq_threshold=mapq_threshold,
            check_5_prime_clip=True,
            check_passing=False)

        read_group_set = set()
        read_summary = []
        # until_eof will include unmapped reads, also
        for read in input_bamfile.fetch(until_eof=True):
            tagged_read = at.tag_read(read, decompose_barcode=False)

            status_code = status_coder(tagged_read)

            summary_record = {"id": tagged_read.get('read').query_name,
                              "status": status_code,
                              "mapq": tagged_read.get('read').mapping_quality,
                              "flag": tagged_read.get('read').flag,
                              "chr": tagged_read.get('read').reference_name,
                              "strand": "*" if tagged_read.get('read').is_unmapped  # noqa
                              else "-" if tagged_read.get('read').is_reverse else '+',  # noqa
                              "five_prime": tagged_read.get('read').get_tag("XS"),  # noqa
                              "insert_start": tagged_read.get('read').get_tag("XI"),  # noqa
                              "insert_stop": tagged_read.get('read').get_tag("XE"),  # noqa
                              "insert_seq": tagged_read.get('read').get_tag("XZ")}  # noqa

            # add the additional tagged elements, defined in
            # the barcode_details json
            for k, v in at.tagged_components.items():
                summary_record[k] = tagged_read.get('read').get_tag(v)

            read_summary.append(summary_record)

            if status_code == 0:
                tmp_tagged_bam_dict['passing'].write(tagged_read.get('read'))
            else:
                tmp_tagged_bam_dict['failing'].write(tagged_read.get('read'))

        for k, v in tmp_tagged_bam_dict.items():
            v.close()
            pysam.index(tmp_output_bampath_dict[k])

        # copy alignments from the tmp file to the actual output so that
        # we can include the RG headers. It is frustrating that this
        # seems like the only way to do this in pysam.
        # TODO find a way to just add the header rather than having to
        # iterate over
        # the reads
        new_header = input_bamfile.header.to_dict()
        # Create new read group header. Note: this is used below in
        # the tagged_bam
        new_header['RG'] = [{'ID': rg} for rg in read_group_set]

        # open the tmp_tagged_bam for reading
        tmp_tagged_bam_dict = {k: pysam.AlignmentFile(v, "rb")
                               for k, v in tmp_output_bampath_dict.items()}

        output_bampath_dict = \
            {x: os.path.join(
                args.output_dir,
                output_basename + '_' + x + out_suffix)
                for x in ['passing', 'failing']}

        tagged_bam_output_dict = \
            {k: pysam.AlignmentFile(v, "wb", header=new_header)
                for k, v in output_bampath_dict.items()}

        # iterate over the reads to re-write
        logger.info("re-writing bam with updated header...")
        for k, tmp_tagged_bam in tmp_tagged_bam_dict.items():
            # until_eof will include unmapped reads, also
            for read in tmp_tagged_bam.fetch(until_eof=True):
                tagged_bam_output_dict[k].write(read)
            # close the temp bampath. Note that the whole temp directory
            # will be deleted when we leave the with
            # TempDirectory as ... clause
            tmp_tagged_bam.close()

    for k, v in tagged_bam_output_dict.items():
        v.close()
        pysam.index(output_bampath_dict[k])

    # Close input
    input_bamfile.close()

    logger.info('summarizing to summary_df and qbed...')
    aln_summary_df = pd.DataFrame(read_summary)
    sp = SummaryParser(aln_summary_df)
    qbed_df = sp.to_qbed()

    qbed_df.to_csv(os.path.join(output_dir, output_basename + '.qbed'),
                   sep='\t',
                   index=False)

    aln_summary_df = aln_summary_df\
        .assign(status_decomp=StatusFlags.decompose(aln_summary_df['status'],
                                                    as_str=True))
    # Convert the lists in the status_decomp column to strings
    # by joining the elements with a comma
    aln_summary_df['status_decomp'] = aln_summary_df['status_decomp']\
        .apply(lambda x: ', '.join(x))

    qc_summary = aln_summary_df\
        .groupby(['status_decomp'])\
        .agg({'id': 'count'})\
        .reset_index()\
        .rename(columns={'id': 'count'})

    qc_summary.to_csv(
        os.path.join(output_dir, output_basename + '_summary.tsv'),
        sep='\t',
        index=False)

    if args.verbose_qc:
        aln_summary_df.to_csv(
            os.path.join(output_dir, output_basename + '_aln_info.tsv'),
            sep='\t',
            index=False)

    logger.info('complete!')
    # return {'summary': aln_summary_df, 'qbed': qbed_df}
