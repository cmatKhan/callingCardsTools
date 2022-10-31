#pylint:disable=C0206
import logging
from logging.config import dictConfig
import sys
import os
import argparse

from callingcardstools.ReadParser import ReadParser

from Bio import SeqIO
 
def parse_args(args=None):
    """parse command line arguments for tagBam

    Args:
        args (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    Description = "Demultiplex fastq into individual TF fastq files, "+\
        "and unidentifiable reads based on barcode"
    Epilog = ""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument('-r1', 
	                    '--read1', 
						help = 'Read 1 filename (full path)', 
						required=True)
    parser.add_argument('-r2',
	                    '--read2', 
						help = 'Read2 filename (full path)', 
						required=True)
    parser.add_argument('-b',
	                    '--barcode_details',
						help = 'barcode filename (full path)',
						required=True)
    parser.add_argument('-s',
	                    '--split_key',
						help = "Either a name of a key in "+\
                            "barcode_details['components'], or just a string. "+\
                                "This will be used to create the passing "+\
                                    "output fastq filenames",
						required=True)
    parser.add_argument("-l",
	                    "--log_level",
                        help="the desired log level (default warning).",
                        choices=("critical", "error", "warning", "info", "debug"),
                        default="warning")

    return parser.parse_args(args)

def main(args = None):

    args = parse_args(args)

    log_config = {
        "version": 1,
        "root": {
            "handlers": ["console"],
            "level": f"{args.log_level.upper()}"
        },
        "handlers": {
            "console": {
                "formatter": "std_out",
                "class": "logging.StreamHandler"
            }
        },
        "formatters": {
            "std_out": {
                "format": "%(asctime)s : %(module)s : %(funcName)s : line: %(lineno)d\n"+\
                    "process details : %(process)d, %(processName)s\n"+\
                        "thread details : %(thread)d, %(threadName)s\n"+\
                                "%(levelname)s : %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S"
            }
        }
    }
    dictConfig(log_config)

    # Check inputs
    logging.info('checking input...')
    input_path_list = [args.read1,
                       args.read2,
                       args.barcode_details]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)
    # create the read parser object
    rp = ReadParser(args.barcode_details, args.read1, args.read2)
    logging.info('opening fq files')
    rp.open()    
    # if the split_key isn't in the barcode_details, then every passing 
    # read goes into a file with that name
    if args.split_key not in rp.barcode_dict['components']:
        msg = f"{args.split_key} not found in barcode_dict['components']. "\
            f"all output is directed to {args.split_key}_R1,2.fq"
        logging.info(msg)

        determined_out = {
            'r1': open(f"{args.split_key}_R1.fq", "w"), #pylint:disable=W1514
            'r2': open(f"{args.split_key}_R2.fq", "w")  #pylint:disable=W1514
        }
    # elif the split_key is in barcode_details, create/open a fq output file for
    # each of the keys in barcode[components][split_key]
    else:
        determined_out = {
            'r1': {tf:open(f"{tf}_R1.fq", "w") for tf in\
                rp.barcode_dict['components'][args.split_key]['map'].values()}, #pylint:disable=W1514
            'r2': {tf:open(f"{tf}_R2.fq", "w") for tf in\
                rp.barcode_dict['components'][args.split_key]['map'].values()}  #pylint:disable=W1514
        }
    # create/open undetermined read output -- these are reads which do not 
    # match barcode expectations
    undetermined_out = {
        'r1': open("undetermined_R1.fq", "w"), #pylint:disable=W1514
        'r2': open("undetermined_R2.fq", "w")  #pylint:disable=W1514
    }

    # iterate over reads, split reads whose barcode components match expectation 
    # into the appropraite file, and reads which don't fulfill barcode expecations 
    # into undetermined.fq
    # also record each read and barcode details into the id_to_bc.csv file. 
    # note that this will be pretty big (measured in GBs, not as big as R1, but close)
    logging.info('opening id to barcode map...')
    additional_components = ['tf','restriction_enzyme']
    with open("id_bc_map.tsv", "w") as id_bc_map:  #pylint:disable=W1514
        id_bc_map.write("\t".join(['id'] + list(rp.components) + additional_components))
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
                for read_end in ['r1','r2']:
                    output_handle = determined_out[read_end]\
                        .get(
                            read_dict['status']['details'][args.split_key]['name'], 
                            determined_out[read_end])
                    SeqIO.write(
                        read_dict[read_end], 
                        output_handle, 
                        'fastq')
            else:
                for read_end in ['r1','r2']:
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
                [read_dict['r1'].id,] + \
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

if __name__ == '__main__':
	sys.exit(main())
