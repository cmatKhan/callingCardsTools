#!/usr/bin/env python


"""
written by: chase mateusiak, chasem@wustl.edu

Loop over alignments. Only read1 are recorded.

Reads are considering passing if they are primary alignments which are not
unmapped, secondary, failing sequencer/aligner QC or supplementary alignment.
If paired end, then the reads must be correctly paired.

Reads may also be filtered by a mapq threshold.

If require_exactly_length is set to True, then only reads which
align without soft clipping are retained.

Output are two bed files in a format specified below. One file has suffix
_passing.bed, the other _failing.bed, and represent non overlapping paritions
of read1.

Bed format v1.0
https://github.com/samtools/hts-specs/blob/master/BEDv1.pdf

Historic Note: This is a 'modified' qbed file which more closely follows the
bed specifications. Each entry is a hop -- no aggregation at this point.

The output follows a Bed6+3 format with the following fields:

chrom, chromStart, chromEnd, name(barcode),
score(aln_mapq), strand, reads*, insert_seq, aln_flag

* reads are set to 1 -- these are not aggregated at this stage, so every row
represents a single hop. This field exists in order to accomodate previous
scripts
"""

# standard library
import os
import sys
import argparse
import logging

# outside dependencies
import pysam
import pandas as pd

# local dependencies
from callingcardstools.general import count_hops

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)


BED_6_3_COLNAMES = ['chrom', 'chromStart', 'chromEnd',
                    'barcode', 'aln_mapq', 'strand',
                    'reads', 'insert_seq', 'aln_flag']

def parse_args(args=None):
    Description = "Extract reads which potentially describe transpositions " +\
    "from a bam file, transform to bed6+3 format with column headers " +\
         " ".join(BED_6_3_COLNAMES)
    Epilog = "Example usage: python add_read_group_and_tags.py \
        <input.bam> <output.bam> <genome.fasta> <genome.fasta.fai> \
        <id_length> <optional: insertion_length (default is 1)>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("bampath",
                         help="path to the input bam file")
    parser.add_argument("require_exact_length",
                        help="True to filter out any soft-clipped reads, False "+\
                        "otherwise. Default is False",
                        default='False')
    parser.add_argument("mapq_filter",
                         help = "minimum value above which to accept alignments")

    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    def parse_bool(val):
        # translate from possible input to boolean
        switcher = {
            "0": False,
            "1": True,
            "true": True,
            "false": False
        }
        # throw error if cannot cast val to boolean based on switch statement
        if not isinstance(switcher.get(val.lower()),bool):
            raise ValueError("Invalid choice for argument single_end: %s. \
                Must be one of 0,1, true/True/TRUE, false/False/FALSE" %single_end)

        return switcher.get(val.lower())

    require_exact_length = parse_bool(args.require_exact_length)

    # Check inputs
    input_path_list = [args.bampath]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    # loop over the reads in the bam file and add the read group (header and tag)
    # and the XI and XZ tags
    count_hops(args.bampath,
               require_exact_length,
               int(args.mapq_filter))

    sys.exit(0)


if __name__ == "__main__":
    sys.exit(main())
