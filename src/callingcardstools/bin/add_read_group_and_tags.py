#!/usr/bin/env python

"""
written by: chase mateusiak, chasem@wustl.edu

This script takes as input a bam file which has a barcode, probably added by
UMITools, of length n added to the end of each QNAME. The barcode is extracted
based on length (eg, if the barcode is 13 characters, then 13 characters are
extracted from the end of the QNAME string of each read). Currently, this
script loops over the bam file twice -- the first time, extracting the barcode
from each read, during which a unique set of barcodes are created. Then, a
write-able bamfile is opened, this unique set of barcodes is added as a RG
(read group) header, which makes parsing, eg splitting the bam into parts by
RG possible. Then, a second loop is performed during which the tags RG, XZ and
XI are added to each read line. RG is the read group (barcode), XZ is the
coordinate at which the transposon inserted on a given chromosome, and XI is the
sequence x bases upstream of the insertion site. If the read is unmapped,
XI and XZ are set to "*".
"""
# standard library
import os
import sys
import argparse

# local import
from callingcardstools.general import add_read_group_and_tags

def parse_args(args=None):
    Description = "Extract the barcode, added by UMItools, from each read id, \
                   Add the Read Group (RG) header, and add the appropriate \
                   RG tag to each alignment record, and add the XI and XZ tags \
                   which describe the insertion location and sequence at that \
                   location"
    Epilog = "Example usage: python add_read_group_and_tags.py \
        <input.bam> <output.bam> <genome.fasta> <genome.fasta.fai> \
        <id_length> <optional: insertion_length (default is 1)>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("bampath_in",
                         help="path to the input bam file")
    parser.add_argument("bampath_out",
                         help="path to the output bam file")
    parser.add_argument("genome_path",
                         help = "Path to the .fasta genome used in the \
                                 alignment")
    parser.add_argument("genome_index_path",
                         help = "path to the .fai produced by samtools \
                                 faidx from the genome .fasta")
    parser.add_argument("id_length",
                        help="length of the barcode in the bam ID line")
    parser.add_argument("insertion_length",
                        help = "The length of the transposase insertion. \
                                Default is 1 since pysam uses 0 based half \
                                    open intervals. [0,1) returns 1 base",
                        default = '1')

    return parser.parse_args(args)


def main(args=None):
    """_summary_

    Args:
        args (_type_, optional): _description_. Defaults to None.

    Raises:
        FileNotFoundError: _description_
    """
    args = parse_args(args)

    # Check inputs
    input_path_list = [args.bampath_in,
                       args.genome_path,
                       args.genome_index_path]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    # loop over the reads in the bam file and add the read group (header and tag)
    # and the XI and XZ tags
    add_read_group_and_tags(args.bampath_in,
                            args.bampath_out,
                            args.genome_path,
                            args.genome_index_path,
                            int(args.id_length),
                            int(args.insertion_length))

    sys.exit(0)


if __name__ == "__main__":
    sys.exit(main())
