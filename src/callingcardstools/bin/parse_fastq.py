"""
cc_filter_reads.py
written 9/1/16 by RDM

Modified 4/3/16 to remove /n line from split files.

Modified 12/17/18 to update to include restriction enzyme QC
Modified 12/21/18 to add quality control path information and raw folder

Read #1
The first 5 base pairs will be the primer barcode (index 0 to 4), followed by
17bp of transposon sequence (index 5 to 21):  AATTCACTACGTCAACA, after which
is genomic sequence

Read #2
The first 8bp will be the transposon barcode, followed by the restriction 
enzyme bank and then genomic DNA.  So it will look like
NNNNNNNN TCGA GCGC CCGG where the genomic sequence could start at any 
of the three restriction sites.

The barcode file should be of the following format: 
#expt name \t primer barcode \t transposon barcode 

usage 
python cc_filter_reads.py -r1 <read1 file> -r2 <read2 file> 
 -b<barcode file> -o <output path> -t <path to temp folder>
-rp <path to raw folder> 
--hammp <hamming distance for primer barcode>
--hammt <hamming distance for transposon barcode>


    Required
    -r1 {read 1 filename (full path)}
    -r2 {read 2 filename (full path)}

    Not Required
    -b {barcode file = ../scripts/barcodes.txt}
    -t {path to temp folder = ../temp/}
    -rp {path to raw folder = ../raw/} 
    -hp {hamming distance for primer bc =0}
    -tp {hamming distance for transposon bc = 0}
    -o {output path} -t {path to temp folder} -rp {path to raw folder}  
    --hammp {hamming distance for primer barcode}
    --hammt {hamming distance for transposon barcode}
    
    Description:

	1. Reads barcodes and corresponding experiments into a dictionary.
	2. Opens the read 1 and checks for transposon sequence.
	3. If the tranposon sequence is present, it checks to see if the primer 
	barcode matches the transposon barcode.
	4. If both filters are passed, it prints the reads to a file of the format: 
	exptname_primerbc_transposonbc_R1.fasta (or R2 or I2) in the raw directory.  
	Undetermined reads are outputted to the temp directory.  
	5. It also prints a master file of all of the reads file in the temp
	directory.
	6. The program then outputs a brief QC that lists the fraction of reads 
	that have a transposon match, the fraction of reads that have matching
	primer and transposon barcodes and the number of reads for each experiment 
	and the total number of reads analyzed.  It also computes, for each 
	experiment,the number of reads generated for each restriction enzyme.  
	This qc file is outputted to the temp directory, and is used by 
	perform_QC.py.
"""
import logging
from logging import config
import sys
import os
import argparse
import re

from callingcardstools.ReadParser import ReadParser

from Bio import SeqIO

log_config = {
    "version":1,
    "root":{
        "handlers" : ["console"],
        "level": "INFO"
    },
    "handlers":{
        "console":{
            "formatter": "std_out",
            "class": "logging.StreamHandler",
            "level": "INFO"
        }
    },
    "formatters":{
        "std_out": {
            "format": "%(asctime)s : %(levelname)s : %(module)s : %(funcName)s : %(lineno)d : (Process Details : (%(process)d, %(processName)s), Thread Details : (%(thread)d, %(threadName)s))\nLog : %(message)s",
            "datefmt":"%d-%m-%Y %I:%M:%S"
        }
    },
}


config.dictConfig(log_config)
  
def parse_args(args=None):
    """parse command line arguments for tagBam

    Args:
        args (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    Description = "Demultiplex fastq into individual TF fastq files, and unidentifiable reads based on barcode"
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
    parser.add_argument('-o',
	                    '--output_id',
						help = "Either a name of a key in "+\
                            "barcode_details['components'], or just a string. "+\
                                "This will be used to create the passing "+\
                                    "output fastq filenames",
						required=True)

    return parser.parse_args(args)



def main(args = None):

    args = parse_args(args)

    # Check inputs
    logging.info('checking input...')
    input_path_list = [args.read1,
                       args.read2,
                       args.barcode_details]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)
    
    rp = ReadParser(args.barcode_details, args.read1, args.read2)
    rp.open()
    
    logging.info('opening fq files')

    if args.output_id not in rp.barcode_dict['components']:
        msg = f"{args.output_id} not found in barcode_dict['components']. "\
            f"all output is directed to {args.output_id}_R1,2.fq"
        logging.info(msg)

        determined_out = {
            'r1': open(f"{args.output_id}_R1.fq", "w"), #pylint:disable=W1514
            'r2': open(f"{args.output_id}_R2.fq", "w")  #pylint:disable=W1514
        }
    else:
        determined_out = {
            'r1': {tf:open(f"{tf}_R1.fq", "w") for tf in\
                rp.barcode_dict['components'][args.output_id]['map'].values()}, #pylint:disable=W1514
            'r2': {tf:open(f"{tf}_R2.fq", "w") for tf in\
                rp.barcode_dict['components'][args.output_id]['map'].values()}  #pylint:disable=W1514
        }

    undetermined_out = {
        'r1': open("undetermined_R1.fq", "w"), #pylint:disable=W1514
        'r2': open("undetermined_R2.fq", "w")  #pylint:disable=W1514
    }



    # create the id to barcode components (etc) map
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
                            read_dict['status']['details'][args.output_id]['name'], 
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
