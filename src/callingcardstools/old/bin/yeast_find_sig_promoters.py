#!/usr/bin/env python

"""Calculate poisson and hypergeomtric p-values for a calling cards experiment.
Given the experiment qbed file, background qbed file, promoter definitions,
a csv which maps various chromosome naming conventions to one another (eg
refseq to ucsc), and a choice of which chromosome naming convention to use,
aggregate hops over the promoter regions and calculate the poisson and
hypergeometric pvalues. Manipulations of the data frames are largely performed
within sqlite, which may be done either in memory or saved to disc. Choose to
save the sqlite database to disc for easy re-definition of promoter regions.

Note: currently no 'lax' score (lax just 'relaxes' the frame in which hops are
  counted by some amount. In the yeast 3.0 scripts, it is 200bp) is not
  calculated.

Note: the strand of the hop is not considered for the background or
  experimental data -- if a hop occurs on either strand in a given promoter
  region, it is counted.

written by: chase mateusiak, chasem@wustl.edu, 202207

:raises AttributeError: Various custom AttributeErrors are raised
:raises FileExistsError: All cmd line param filepaths are checked for
    existence
:raises KeyError: If a key in one of the dictionaries describing attributes
    related to the promoter, background or experimental table DNE

:output: a statistics data frame which will have poisson_pval
    and hypergeom_pval columns with filename
    <experimental_qbed_basename>_stats.csv. Optionally, the sqlite database
    may be saved to disc. A cmd line argument directs the location.
"""

# standard lib
import sys
import argparse
from os.path import exists,basename,splitext
import logging
# third party
import pandas as pd
# local
from callingcardstools.yeast_stats import create_hop_db, yeast_stats_frame
from callingcardstools.general import standardize_chr_names

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

def parse_args(args=None):
    Description = "create a sqlite database with background and expression data " +\
    "and calculate experimental significance. Optionally save the sqlite database."
    Epilog = "Example usage: yeast_find_sig_promoters.py \
        <input.bam> <output.bam> <genome.fasta> <genome.fasta.fai> \
        <id_length> <optional: insertion_length (default is 1)>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--promoter_bed_path", "-p",
                        type=str,
                        help="path to to a bed file format promoter file")
    parser.add_argument("--background_qbed_path", "-b",
                        type=str,
                        help="path` to the background data table")
    parser.add_argument("--experimental_qbed_path", "-e",
                        type=str,
                        help="path to the experimental data table")
    parser.add_argument("--chr_map", "-c",
                        type=str,
                        help="a csv which maps between various "+\
                            "chromosome naming conventions, eg ucsc, ignomes " +\
                                "and refseq")
    parser.add_argument("--standard_chr_format", "-s",
                        type=str,
                        help="this must be one of the column names in "+\
                            "the chr_map. All chromosome identifiers will be "+\
                                "translated to this naming convention. "+\
                                    "default is 'refseq'",
                        default="refseq")
    parser.add_argument("--sqlite_db_out", "-d",
                        type=str,
                        help="path to sqlite database output. "+\
                            "Default is ':memory:' for an in memory DB. " +\
                                "For instance, 'my_db.sqlite' would write to "+\
                                    "$PWD/my_db.sqlite.",
                        default=":memory:")
    parser.add_argument("--poisson_pseudocount", "-x",
                        type=float,
                        help="pseudocount to add to the poisson p-value "+\
                            "calculation. Default = 0.2",
                        default=0.2)

    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    # check cmd line input
    if args.sqlite_db_out != ":memory:" and\
         splitext(args.sqlite_db_out)[1] != '.sqlite':
        raise ValueError("sqlite must either be in memory -- ':memory:' -- or "+\
            "have the extension '.sqlite'.")
    if(float(args.poisson_pseudocount) < 0):
        raise ValueError('poisson_pseudocount must be 0 or greater')     
    if not exists(args.chr_map):
        raise FileExistsError(f"File Not Found: {args.chr_map}")
    # read in chr_map_df to check to make sure the standard_chr_format exists
    # in the colnames
    chr_map_df = pd.read_csv(args.chr_map)
    # check stadand_chr_format is in chr_map_df colnames
    if not args.standard_chr_format in chr_map_df.columns:
        raise AttributeError('Standard chr format not in chr_map columns')
    if not exists(args.promoter_bed_path):
        raise FileExistsError(f"File Not Found: {args.promoter_bed_path}")
    if not exists(args.background_qbed_path):
        raise FileExistsError(f"File Not Found: {args.background_qbed_path}")
    if not exists(args.experimental_qbed_path):
        raise FileExistsError(f"File Not Found: {args.experimental_qbed_path}")

    # setup colnames, etc for db tables
    PROMOTER_BED_COLNAMES = ["chrom", "chromStart", "chromEnd",
                             "name", "score", "strand"]
    QBED_COLNAMES = ['chrom', 'chromStart', 'chromEnd', 'reads', 'strand']
    INDEX_COL_STRING = '("chrom", "chromStart" ASC, "chromEnd" ASC, "strand");'

    promoter_dict = {
        'tablename': "promoters",
        'index_name': "promoter_index",
        'index_col_string': INDEX_COL_STRING
    }

    background_dict = {
        'tablename': "background",
        'hop_view': "background_hops",
        'index_name': "background_index",
        'index_col_string': INDEX_COL_STRING
    }

    experimental_dict = {
        'tablename': "experiment",
        'hop_view': "experiment_hops",
        'index_name': "experiment_index",
        'index_col_string': INDEX_COL_STRING
    }

    # read in promoter definitions
    promoter_df = pd.read_csv(args.promoter_bed_path,
                              sep = "\t",
                              names = PROMOTER_BED_COLNAMES)
    promoter_df = standardize_chr_names(promoter_df,
                                        chr_map_df,
                                        args.standard_chr_format)
    # read in background hop data
    background_df = pd.read_csv(args.background_qbed_path,
                                sep = "\t",
                                names = QBED_COLNAMES)
    background_df = standardize_chr_names(background_df,
                                          chr_map_df,
                                          args.standard_chr_format)
    # read in experimental hop data
    experimental_df = pd.read_csv(args.experimental_qbed_path,
                                  sep = "\t",
                                  names = QBED_COLNAMES)

    experimental_df = standardize_chr_names(experimental_df,
                                            chr_map_df,
                                            args.standard_chr_format)

    poisson_pseudocount = float(args.poisson_pseudocount)

    # create database
    create_hop_db(args.sqlite_db_out, 
                  args.promoter_bed_path, 
                  PROMOTER_BED_COLNAMES, 
                  args.background_qbed_path, 
                  args.experimental_qbed_path, 
                  QBED_COLNAMES,
                  promoter_dict,
                  background_dict,
                  experimental_dict,
                  chr_map_df,
                  args.standard_chr_format)

    # extract stats dataframe from db
    quant_df = yeast_stats_frame(args.sqlite_db_out, 
                                 background_dict['tablename'], 
                                 background_dict['hop_view'], 
                                 experimental_dict['tablename'],
                                 experimental_dict['hop_view'], 
                                 poisson_pseudocount)

    # write to CWD
    quant_output_path = splitext(basename(args.experimental_qbed_path))[0] + \
                            "_stats.csv"
    quant_df.to_csv(quant_output_path, index=False)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
