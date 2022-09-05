import pandas as pd
from callingcardstools.yeast_stats import create_hop_db, yeast_stats_frame
from callingcardstools.general import standardize_chr_names

def test_create_hop_db(tmp_path):

    db_out_path = "/home/oguzkhan/Desktop/cc_testing.sqlite"
    chr_map_df = pd.read_csv("tests/test_data/yeast_chr_map.csv")
    standard_chr_format = "ucsc"
    promoter_bed_path = "/home/oguzkhan/code/nf-core-callingcards/assets/yeast/notOrf_sacCer3_features.bed"
    background_qbed_path = "/home/oguzkhan/code/nf-core-callingcards/assets/yeast/S288C_dSir4_Background.qbed"
    experimental_qbed_path = "/mnt/scratch/calling_cards/nf_pipeline/results/run_5986_T1/quantification/passing/qbed/run_5986_T1_tagged_passing_ERT1.qbed"

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
    promoter_df = pd.read_csv(promoter_bed_path,
                              sep = "\t",
                              names = PROMOTER_BED_COLNAMES)
    promoter_df = standardize_chr_names(promoter_df,
                                        chr_map_df,
                                        standard_chr_format)
    # read in background hop data
    background_df = pd.read_csv(background_qbed_path,
                                sep = "\t",
                                names = QBED_COLNAMES)
    background_df = standardize_chr_names(background_df,
                                          chr_map_df,
                                          standard_chr_format)
    # read in experimental hop data
    experimental_df = pd.read_csv(experimental_qbed_path,
                                  sep = "\t",
                                  names = QBED_COLNAMES)

    experimental_df = standardize_chr_names(experimental_df,
                                            chr_map_df,
                                            standard_chr_format)

    # create database
    x = create_hop_db(db_out_path, 
                  promoter_bed_path, 
                  PROMOTER_BED_COLNAMES, 
                  background_qbed_path, 
                  experimental_qbed_path, 
                  QBED_COLNAMES,
                  promoter_dict,
                  background_dict,
                  experimental_dict,
                  chr_map_df,
                  standard_chr_format)

    assert x == 0

def test_yeast_stats_frame(tmp_path):
    db_path = "/home/oguzkhan/Desktop/cc_testing.sqlite"
    bg_table = "background"
    expr_table = "experiment"
    bg_viewname = 'background_hops'
    expr_viewname = 'experiment_hops'
    poisson_pseudocount = .2

    quant_df = yeast_stats_frame(db_path, bg_table, bg_viewname, expr_table, 
                                 expr_viewname, poisson_pseudocount)

    assert 1==1
