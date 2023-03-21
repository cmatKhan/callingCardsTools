"""
    ---------------------
    make_ccffile.py
    ---------------------
    
    written 12/7/18 by RDM
    modified 10/21/22 by CAM -- runs on a single bam

    Description:
    This function make .ccf files from mapped .bam files.
    ccf files have the following columns: [chr,start,end,reads,strand,barcode]
    but only the first 4 columns are required. The genome coordinates are 
    1-indexed

    usage:
    make_ccffile -b <barcode filename>  -o <output path>

    required
    none

    not required
    -b {barcode filename}
    -p {path}
    -t,{temppath,default='../temp' path for temporary files}
"""

import argparse
import sys
import pysam
import pandas as pd
import csv
import re
import os
import pdb


def remove_suffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[:-len(suffix)]
    return input_string


def sort_ccf_file(ccffilename):
    sorter = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
              "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrM"]
    ccf_frame = pd.read_csv(ccffilename, delimiter='\t', header=None, names=[
                            'chr', 'start', 'end', 'reads', 'strand', 'barcode'])
    ccf_frame.chr = ccf_frame.chr.astype("category")
    ccf_frame.chr.cat.set_categories(sorter, inplace=True)
    ccf_frame = ccf_frame.sort_values(['chr', 'start'])
    ccf_frame.to_csv(ccffilename, sep='\t', header=False, index=False)
    return [len(ccf_frame), ccf_frame.reads.sum()]


def make_ccffile(sbamFilename, out_dir):
    # make chromosome list
    # mitrochondrial mappings break ccf viewer.  Uncomment this when we switch over to the
    # wustl browser.
    chr_list = {"NC_001133", "NC_001134", "NC_001135", "NC_001136",
                "NC_001137", "NC_001138", "NC_001139", "NC_001140",
                "NC_001141", "NC_001142", "NC_001143", "NC_001144",
                "NC_001145", "NC_001146", "NC_001147", "NC_001148", "NC_001224"}
    chr_dict = {"NC_001133": "I", "NC_001134": "II", "NC_001135": "III", "NC_001136": "IV",
                "NC_001137": "V", "NC_001138": "VI", "NC_001139": "VII", "NC_001140": "VIII",
                "NC_001141": "IX", "NC_001142": "X", "NC_001143": "XI", "NC_001144": "XII",
                "NC_001145": "XIII", "NC_001146": "XIV", "NC_001147": "XV", "NC_001148": "XVI", "NC_001224": "M"}

    # initialize quality control dictionary

    qc_dict = {}

    # inialize ccf dictionary
    for_ccf_dict = {}
    rev_ccf_dict = {}
    # make AlignmentFile object
    current_bamfile = pysam.AlignmentFile(sbamFilename, "rb")

    ccf_file = os.path.join(
        out_dir,
        remove_suffix(os.path.basename(sbamFilename), ".bam")+'.ccf')

    with open(ccf_file, 'w') as output_handle:
        # loop through the chromosomes and pileup start sites
        for chr in chr_list:
            # print chr
            aligned_reads_group = current_bamfile.fetch(chr)
            # now loop through each read and pile up start sites
            for aread in aligned_reads_group:
                if not (aread.is_read2):  # only count read1's
                    # is the read a reverse read?
                    if aread.is_reverse:
                        # read has to start right after primer
                        if aread.query_alignment_end == (aread.query_length):
                            pos = aread.get_reference_positions()[-1]+1
                        if (chr, pos) in rev_ccf_dict:
                            rev_ccf_dict[(chr, pos)] += 1
                        else:
                            rev_ccf_dict[(chr, pos)] = 1
                    elif aread.query_alignment_start == 0:  # forward read has to start right after primer
                        pos = aread.get_reference_positions()[0]+1
                        if (chr, pos) in for_ccf_dict:
                            for_ccf_dict[(chr, pos)] += 1
                        else:
                            for_ccf_dict[(chr, pos)] = 1
        # output dictionary to ccf file
        for key in for_ccf_dict:
            output_handle.write("%s\t%s\t%s\t%s\t%s\n" % (
                "chr"+chr_dict[key[0]], key[1], str(int(key[1])+1), for_ccf_dict[key], "+"))
        for key in rev_ccf_dict:
            output_handle.write("%s\t%s\t%s\t%s\t%s\n" % (
                "chr"+chr_dict[key[0]], key[1], str(int(key[1])+1), rev_ccf_dict[key], "-"))

    # OPEN GNASHY FILE AND SORT BY CHR THEN POS
    qc_list = sort_ccf_file(ccf_file)
    # after all experiments have been analyzed, print out qc
    qc_output = os.path.join(
        out_dir,
        remove_suffix(os.path.basename(ccf_file), ".ccf")+"_ccfQC.txt")
    with open(qc_output, 'w') as qc_handle:
        qc_handle.write("%s\t%s\t%s\n" % (remove_suffix(
            os.path.basename(ccf_file), '.ccf'), qc_list[0], qc_list[1]))


def parse_args(
        subparser: argparse.ArgumentParser,
        script_desc: str,
        common_args: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """This is intended to be used as a subparser for a parent parser passed 
    from __main__.py. It adds the arguments required to run the makeccf 
    script from the original calling cards pipeline

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
        'legacy_makeccf',
        help=script_desc,
        prog='legacy_makeccf',
        parents=[common_args]
    )

    parser.set_defaults(func=main)
    
    parser.add_argument('-s', '--sampath', help='path to sam/bam')
    parser.add_argument('-o', '--outputpath', help='output path')

    return subparser


def main(args=None):
    if not args.outputpath[-1] == "/":
        args.outputpath = args.outputpath+"/"
    make_ccffile(args.sampath, args.outputpath)


if __name__ == '__main__':
    sys.exit(main())
