# pylint:disable=C0206,W1514
import logging
import os
import argparse

from Bio import SeqIO

from callingcardstools.Reads.ReadParser import ReadParser
from callingcardstools.BarcodeParser.yeast import BarcodeQcCounter

__all__ = ['parse_args', 'split_fastq']

logger = logging.getLogger(__name__)


def parse_args(
        subparser: argparse.ArgumentParser,
        script_desc: str,
        common_args: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """This is intended to be used as a subparser for a parent parser passed
    from __main__.py. It adds the arguments required to iterate over yeast
    reads and demultiplex the fastq into separate files based on the TFs
    in the barcode details file.

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
        'split_fastq',
        help=script_desc,
        prog='split_fastq',
        parents=[common_args]
    )

    parser.set_defaults(func=split_fastq)

    parser.add_argument('-r1',
                        '--read1',
                        help='Read 1 filename (full path)',
                        required=True)
    parser.add_argument('-r2',
                        '--read2',
                        help='Read2 filename (full path)',
                        required=True)
    parser.add_argument('-b',
                        '--barcode_details',
                        help='barcode filename (full path)',
                        required=True)
    parser.add_argument('-s',
                        '--split_key',
                        help="Either a name of a key in " +
                        "barcode_details['components'], or just a string. "
                        "This will be used to create the passing "
                        "output fastq filenames. Defaults to 'tf' which is "
                        "appropriate for yeast data",
                        default="tf")
    parser.add_argument('-n',
                        '--split_suffix',
                        help='append this after the tf name and before _R1.fq '
                        'in the output fastq files. Defaults to "split"',
                        default="split")
    parser.add_argument('-o',
                        '--output_dirpath',
                        help='a path to a directory where the output files '
                        'will be output. Defaults to the current directory',
                        default=".")
    parser.add_argument('-v',
                        '--verbose_qc',
                        help='set this flag to output a file which contains '
                        'the barcode components for each read ID in the '
                        'fastq files associated with its barcode components',
                        action='store_true')
    parser.add_argument('-p',
                        '--pickle_qc',
                        help='set this flag to output a pickle file which '
                        'containing the BarcodeQcCounter object. This is '
                        'useful when splitting the fastq files prior to '
                        'demultiplexing',
                        action='store_true')

    return subparser


def split_fastq(args: argparse.Namespace):

    # Check inputs
    logger.info('checking input...')
    input_path_list = [args.read1,
                       args.read2,
                       args.barcode_details,
                       args.output_dirpath]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input file DNE: {input_path}")

    # create the BarcodeQcCounter object
    logger.info('creating BarcodeQcCounter object...')
    bc_counter = BarcodeQcCounter()

    # create the read parser object
    rp = ReadParser(args.barcode_details, args.read1, args.read2)
    logger.info('opening fq files')
    rp.open()
    # if the split_key isn't in the barcode_details, then every passing
    # read goes into a file with that name
    if args.split_key not in rp.barcode_dict['components']:
        msg = f"{args.split_key} not found in barcode_dict['components']. " \
            f"all output is directed to " \
            f"{args.split_key}_{args.split_suffix}_R1,2.fq"
        logger.info(msg)

        determined_out = {
            'r1': open(f"{args.split_key}_{args.split_suffix}_R1.fq", "w"),  # pylint:disable=W1514 # noqa
            'r2': open(f"{args.split_key}_{args.split_suffix}_R2.fq", "w")  # pylint:disable=W1514 # noqa
        }
    # else the split_key is in barcode_details, create/open a fq output file
    # for each of the keys in barcode[components][split_key]
    else:
        determined_out = {
            'r1': {tf: open(os.path.join(
                args.output_dirpath,
                f"{tf}_{args.split_suffix}_R1.fq"), "w") for tf in
                   rp.barcode_dict['components'][args.split_key]['map'].values()},  # noqa
            'r2': {tf: open(os.path.join(
                args.output_dirpath,
                f"{tf}_{args.split_suffix}_R2.fq"), "w") for tf in
                   rp.barcode_dict['components'][args.split_key]['map'].values()}  # noqa
        }
    # create/open undetermined read output -- these are reads which do not
    # match barcode expectations
    undetermined_out = {
        'r1': open(os.path.join(
            args.output_dirpath,
            f"undetermined_{args.split_suffix}_R1.fq"), "w"),
        'r2': open(os.path.join(
            args.output_dirpath,
            f"undetermined_{args.split_suffix}_R2.fq"), "w")
    }

    # if verbose_qc is true, for read paired read, record a line which
    # associates the fastq read ID with the barcode components
    if args.verbose_qc:
        logger.info('opening id to barcode map...')
        additional_components = ['tf', 'restriction_enzyme']
        id_bc_map = open(os.path.join(
            args.output_dirpath, "id_bc_map.tsv"), "w")
        # write header
        id_bc_map.write(
            "\t".join(['id'] + list(rp.components) + additional_components))
        id_bc_map.write("\n")

    logger.info('parsing fastq files...')
    # iterate over reads, split reads whose barcode components
    # match expectation into the appropriate file, and reads which don't
    # fulfill barcode expectations into undetermined.fq
    while True:
        try:
            rp.next()
        except StopIteration:
            break
        read_dict = rp.parse()
        # if the verbose_qc is on, record the line in the id_bc_map file
        if args.verbose_qc:
            tf = "_".join(
                [read_dict['status']['details'].get('tf', {}).get('name', "*"),
                str(read_dict['status']['details'].get('tf', {}).get('dist', ""))])  # noqa
            restriction_enzyme = \
                read_dict['status']['details']\
                .get('r2_restriction', {})\
                .get('name', "*")
            id_bc_line = \
                [read_dict['r1'].id, ] + \
                [read_dict['components'][comp] for comp in rp.components] + \
                [tf, restriction_enzyme]
            id_bc_map.write("\t".join(id_bc_line))
            id_bc_map.write("\n")
        # Determine to which fastq file the read should be written
        if read_dict['status']['passing'] is True:
            # check that a TF was actually found -- if the TF barcode had
            # a mismatch, then _3 for instance means that the closest match
            # had an edit distance of 3
            for read_end in ['r1', 'r2']:
                output_handle = determined_out[read_end]\
                    .get(
                        read_dict['status']['details'][args.split_key]['name'],  # noqa
                        determined_out[read_end])
                SeqIO.write(
                    read_dict[read_end],
                    output_handle,
                    'fastq')
        else:
            for read_end in ['r1', 'r2']:
                SeqIO.write(
                    read_dict[read_end],
                    undetermined_out[read_end],
                    'fastq')

        # extract component information to summarize as count metrics
        r1_primer_seq = (read_dict['status']
                         ['details']
                         ['r1_primer']
                         ['query'])
        r1_transposon_seq = (read_dict['status']
                             ['details']
                             ['r1_transposon']
                             ['query'])
        r1_transposon_dist = (read_dict['status']
                              ['details']
                              ['r1_transposon']
                              ['dist'])
        r2_transposon_seq = (read_dict['status']
                             ['details']
                             ['r2_transposon']
                             ['query'])
        bc_counter.update(
            (r1_primer_seq, r1_transposon_seq, r2_transposon_seq),
            r1_transposon_dist,
            read_dict['status']['details']['r2_restriction']['name'])

    # close the files
    for read_end in determined_out:
        # close the undetermined files
        undetermined_out[read_end].close()
        # close all tf files
        for write_handle in determined_out[read_end].values():
            write_handle.close()

    # construct the input to the BarcodeQcCounter summarize method
    component_dict = {k: [] for k in ['tf', 'r1_primer', 'r2_transposon']}
    r1_primer_start = rp.barcode_dict['r1']['primer']['index'][0]
    r1_primer_end = rp.barcode_dict['r1']['primer']['index'][1]
    r2_transposon_start = r1_primer_end + \
        rp.barcode_dict['r2']['transposon']['index'][0]
    r2_transposon_end = (r2_transposon_start +
                         rp.barcode_dict['r2']['transposon']['index'][1] -
                         rp.barcode_dict['r2']['transposon']['index'][0])
    for k, v in rp.barcode_dict['components']['tf']['map'].items():
        r1_primer_seq = k[r1_primer_start:r1_primer_end]
        r2_transposon_seq = k[r2_transposon_start:r2_transposon_end]
        component_dict['tf'].append(v)
        component_dict['r1_primer'].append(r1_primer_seq)
        component_dict['r2_transposon'].append(r2_transposon_seq)
        
    # summarize the barcode metrics
    if args.pickle_qc:
        bc_counter.write(raw=args.pickle_qc,
                         output_dirpath=args.output_dirpath,
                         suffix=args.split_suffix)
    else:
        bc_counter.write(component_dict=component_dict,
                         output_dirpath=args.output_dirpath,
                         suffix=args.split_suffix)

    logger.info('Done parsing the fastqs!')
