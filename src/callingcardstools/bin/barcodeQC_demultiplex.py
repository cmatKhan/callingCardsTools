#!/usr/bin/env python

# TODO better docstring, in particular, describe the 6x ? bed file format

"""
Given a 6x? bedfile detailing calling cards read insertions, split the
    name column into barcode components, tally in various ways (QC) and output
    a filtered bedfile with only those reads which pass barcode match standards.

    outputs the following files in the $PWD:
      - <input_filename>_<barcode_component>_tally.tsv
          for each barcode component, a tally of the unique elements in that
          component
      - <input_filename>_<barcode_component>_ppm.tsv
          a position probability matrix describing
          the percent probability of a given base in a given position given the
          multi-sequence alignment given a barcode component
      - <input_filename>_bc_fltr.bed
          a bed file in the same format as the input bedfile, but with a row
          filter applied such that only those records which conform to barcode
          AND insert_sequence specifications remain
"""
# standard library
import os
import sys
import argparse
import json
# outside dependencies
import pandas as pd

BED_6_3_COLNAMES = ['chrom', 'chromStart', 'chromEnd',
                    'barcode', 'aln_mapq', 'strand',
                    'reads', 'insert_seq', 'aln_flag']

# fields by which to group and count transpositions. Not that in the case
# of pooled barcodes, this occurs after splitting into individual barcode
# tables
COORDINATE_GROUPING_FIELDS = ['chrom','chromStart','chromEnd','strand']
# qbed fields and order
qbed_COL_ORDER = ['chrom', 'chromStart', 'chromEnd', 'reads', 'strand']

# NOTE that there is some hard coding with regards to the fields/keys of
# the barcode to TF mapping fields (this is for yeast right now)

def parse_args(args=None):
    Description = "Examine hop barcodes, output some QC metrics, and a bed file \
        which matches the columns in the input bed, but is row filtered for only \
            those reads which meet barcode expectations"
    Epilog = "Example usage: python mammal_barcode_qc.py \
        <cc_format_bedfile.bed> <barcode_details.json>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("bed_path",
                         help="path to the input bed file")
    parser.add_argument("barcode_details",
                         help="A json which describes components of the barcode \
                            which occurs in the name column of the calling cards \
                                bed file format")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    # Check inputs
    for input_path in [args.bed_path,
                       args.barcode_details]:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    bed_df = pd.read_csv(args.bed_path, sep = "\t", names=BED_6_3_COLNAMES)

    with open(args.barcode_details) as f1:
        barcode_details = json.load(f1)

    if not list(barcode_details['components'].keys()).sort() == \
        list(barcode_details['indicies'].keys()).sort():
        raise ValueError('The keys in the barcode_component json and the ' +
            'barcode_component_indicies json are not the same')

    barcode_components_length_check_dict = { k:v[1]-v[0] for k,v in \
                                              barcode_details['indicies'].items()}

    for k,v in barcode_details['components'].items():
        for comp in v:
            if not len(comp) == barcode_components_length_check_dict[k]:
                raise ValueError("There exists a component in " +
                "barcode_components which is not the length described by the " +
                "barcode_details['indicies']: %s, %s"
                %(comp, barcode_components_length_check_dict[k]))

    fltr_bed_output_name = os.path.splitext(os.path.basename(args.bed_path))[0]

    # parse out the name column of the bed file into its components and
    # perform some QC/filtering
    barcode_qc(bed_df,
               barcode_details,
               fltr_bed_output_name)

    sys.exit(0)


if __name__ == "__main__":
    sys.exit(main())
