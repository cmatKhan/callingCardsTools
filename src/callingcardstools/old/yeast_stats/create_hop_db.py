import sqlite3

import pandas as pd

from callingcardstools.general import standardize_chr_names
from callingcardstools.yeast_stats import create_db_table,create_hop_view

def create_hop_db(db_path,
                  promoter_bed_path,
                  promoter_bed_colnames,
                  background_qbed_path, 
                  experimental_qbed_path, 
                  qbed_colnames,
                  promoter_dict,
                  background_dict, 
                  experimental_dict, 
                  chr_map_df,
                  standard_chr_format):

    con = sqlite3.connect(db_path)

    # read in promoter definitions
    promoter_df = pd.read_csv(promoter_bed_path,
                              sep = "\t",
                              names = promoter_bed_colnames)
    promoter_df = standardize_chr_names(promoter_df,
                                        chr_map_df,
                                        standard_chr_format)
    # read in background hop data
    background_df = pd.read_csv(background_qbed_path,
                                sep = "\t",
                                names = qbed_colnames)
    background_df = standardize_chr_names(background_df,
                                          chr_map_df,
                                          standard_chr_format)
    # read in experimental hop data
    experimental_df = pd.read_csv(experimental_qbed_path,
                                  sep = "\t",
                                  names = qbed_colnames)

    experimental_df = standardize_chr_names(experimental_df,
                                            chr_map_df,
                                            standard_chr_format)

    create_db_table(con, promoter_df, promoter_dict)

    create_db_table(con, background_df, background_dict)

    create_db_table(con, experimental_df, experimental_dict)

    # create background hop view
    create_hop_view(con,
                    background_dict['tablename'],
                    background_dict['hop_view'])
    # create experimental hop view
    create_hop_view(con,
                    experimental_dict['tablename'],
                    experimental_dict['hop_view'])

    # close db connection
    con.close()

    return 0