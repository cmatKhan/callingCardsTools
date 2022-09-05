import sqlite3
import pandas as pd
from scipy.stats import poisson as pois
from scipy.stats import hypergeom

def yeast_stats_frame(db_path, bg_table, bg_hops, expr_table, expr_hops, poisson_pseudocount):

    # connect to the database
    con = sqlite3.connect(db_path)
    # get some hop numbers
    total_background_hops = pd.read_sql_query(f"SELECT COUNT(*) as total FROM {bg_table}", con).total
    total_expr_hops       = pd.read_sql_query(f"SELECT COUNT(*) as total FROM {expr_table}", con).total
    hop_ratio             = float(total_expr_hops) / float(total_background_hops)
    bg_plus_expr_hops     = total_background_hops + total_expr_hops
    # get the tables
    promoter_bg_df   = pd.read_sql_query(f"SELECT * FROM {bg_hops}", con)
    promoter_expr_df = pd.read_sql_query(f"SELECT * FROM {expr_hops}", con)
    # close the connection
    con.close()

    # join background and experimental hop tables
    quant_df = pd.merge(promoter_bg_df, promoter_expr_df,
                        on = ['chrom', 'chromStart', 'chromEnd'])\
                 .rename(columns={'hops_x': 'bg_hops',
                                  'hops_y': 'expr_hops'})\
                 .dropna(axis = 'rows')

    # Calculate Statistics
    #usage: scistat.poisson.cdf(x,mu)
    # where X are the experimental hops,
    # and mu is the background hops * (total_expr_hops/total_background_hops)
    quant_df['poisson_pval'] = \
        [1-pois.cdf(quant_df.loc[index,'expr_hops'] + poisson_pseudocount,
                    (quant_df.loc[index, 'bg_hops'] * hop_ratio)+\
                        poisson_pseudocount)
        for index in range(len(quant_df))]

    # usage:
    # scistat.hypergeom.cdf(x,M,n,N)
    # where x is observed number of type I events
    # (white balls in draw) (experiment hops at locus)
    # M is total number of balls (total number of hops)
    # n is total number of white balls (total number of expeirment hops)
    # N is the number of balls drawn (total hops at a locus)
    quant_df['hypergeom_pval'] = \
        [1-hypergeom.cdf(quant_df.loc[index,'expr_hops']-1,
                         bg_plus_expr_hops,
                         total_expr_hops,
                         quant_df.loc[index,'expr_hops'] + \
                            quant_df.loc[index, 'bg_hops'])
        for index in range(len(quant_df))]
        
    return quant_df