# pylint:disable=C0206,W1514
import logging
import os
import argparse

from Bio import SeqIO

from callingcardstools.Reads.ReadParser import ReadParser

__all__ = ['parse_args', 'split_fastq']

logging.getLogger(__name__).addHandler(logging.NullHandler())


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
                        "output fastq filenames",
                        required=True)
    parser.add_argument('-n',
                        '--split_suffix',
                        help='append this after the tf name and before _R1.fq '
                        'in the output fastq files',
                        default="split")
    parser.add_argument('-o',
                        '--output_prefix',
                        help='a path to a directory where the output files '
                        'will be output',
                        default=".")

    return subparser


def split_fastq(args: argparse.Namespace):

    # Check inputs
    logging.info('checking input...')
    input_path_list = [args.read1,
                       args.read2,
                       args.barcode_details,
                       args.output_prefix]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" % input_path)

    # create the read parser object
    rp = ReadParser(args.barcode_details, args.read1, args.read2)
    logging.info('opening fq files')
    rp.open()
    # if the split_key isn't in the barcode_details, then every passing
    # read goes into a file with that name
    if args.split_key not in rp.barcode_dict['components']:
        msg = f"{args.split_key} not found in barcode_dict['components']. " \
            f"all output is directed to " \
            f"{args.split_key}_{args.split_suffix}_R1,2.fq"
        logging.info(msg)

        determined_out = {
            'r1': open(f"{args.split_key}_{args.split_suffix}_R1.fq", "w"),  # pylint:disable=W1514 # noqa
            'r2': open(f"{args.split_key}_{args.split_suffix}_R2.fq", "w")  # pylint:disable=W1514 # noqa
        }
    # else the split_key is in barcode_details, create/open a fq output file
    # for each of the keys in barcode[components][split_key]
    else:
        determined_out = {
            'r1': {tf: open(os.path.join(
                args.output_prefix,
                f"{tf}_{args.split_suffix}_R1.fq"), "w") for tf in
                   rp.barcode_dict['components'][args.split_key]['map'].values()},  # noqa
            'r2': {tf: open(os.path.join(
                args.output_prefix,
                f"{tf}_{args.split_suffix}_R2.fq"), "w") for tf in
                   rp.barcode_dict['components'][args.split_key]['map'].values()}  # noqa
        }
    # create/open undetermined read output -- these are reads which do not
    # match barcode expectations
    undetermined_out = {
        'r1': open(os.path.join(
            args.output_prefix,
            f"undetermined_{args.split_suffix}_R1.fq"), "w"),
        'r2': open(os.path.join(
                args.output_prefix,
                f"undetermined_{args.split_suffix}_R2.fq"), "w")
    }

    # iterate over reads, split reads whose barcode components
    # match expectation into the appropriate file, and reads which don't
    # fulfill barcode expecations into undetermined.fq
    # also record each read and barcode details into the id_to_bc.csv file.
    # note that this will be pretty big (measured in GBs, not as big as R1,
    # but close)
    logging.info('opening id to barcode map...')
    additional_components = ['tf', 'restriction_enzyme']
    with open(os.path.join(
                args.output_prefix,
                "id_bc_map.tsv"), "w") as id_bc_map:  # pylint:disable=W1514
        id_bc_map.write(
            "\t".join(['id'] + list(rp.components) + additional_components))
        id_bc_map.write("\n")

        logging.info('parsing fastq files...')
        while True:
            try:
                rp.next()
            except StopIteration:
                break
            read_dict = rp.parse()
            # check that the barcode edit dist is 0 for each component
            if read_dict['status']['passing'] is True:
                # check that a TF was actually found -- if the TF barcode had
                # a mismatch, then _3 for instance means that the closest match
                # had an edit distance of 3
                for read_end in ['r1', 'r2']:
                    output_handle = determined_out[read_end]\
                        .get(
                            read_dict['status']['details'][args.split_key]['name'],
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

            # write line to id to bc map
            tf = "_".join(
                [read_dict['status']['details'].get('tf', {}).get('name', "*"),
                 str(read_dict['status']['details'].get('tf', {}).get('dist', ""))])
            restriction_enzyme = read_dict['status']['details']\
                .get('r2_restriction', {}).get('name', "*")
            id_bc_line = \
                [read_dict['r1'].id, ] + \
                [read_dict['components'][comp] for comp in rp.components] + \
                [tf, restriction_enzyme]
            id_bc_map.write("\t".join(id_bc_line))
            id_bc_map.write("\n")

    # close the files
    for read_end in determined_out:
        # close the undetermined files
        undetermined_out[read_end].close()
        # close all tf files
        for write_handle in determined_out[read_end].values():
            write_handle.close()

    logging.info('Done parsing the fastqs!')
