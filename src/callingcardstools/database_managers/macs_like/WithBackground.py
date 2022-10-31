import os
import sys
import re
import logging 
import logging.config 
import yaml
from math import inf,floor,ceil

from callingcardstools import DatabaseApi

import pysqlite3 as sqlite3
import pandas as pd
import numpy as np
import scipy.stats as scistat

with open('logging.config.yaml', 'r') as f: #pylint:disable=W1514
    config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

logger = logging.getLogger(__name__)


# scale background
#    lambda_bg = \
#        ((num_bg_hops*(float(total_experiment_hops)/\
#            total_background_hops))/max(num_TTAAs,1))
#
#    pvalue = 1-scistat.poisson.cdf(
#        (num_exp_hops+pseudocounts),
#        lambda_f*max(num_TTAAs,1)+pseudocounts)

# yeast -- same as scale background
#    1-scistat.poisson.cdf((experiment_hops + pseudocount),
#                            (background_hops * (float(exp_hops)/float(bg_hops)) + pseudocount))

class HopsDb(DatabaseApi):
    """This extends DatabaseApi with functions to call significant peaks on mammals data using a MACS-like background free algorithm"""
     
    def region_significance_calculator(self, background_tbl, experiment_tbl, 
                                       ttaa_tbl, poisson_pseudocount):     
        # instantiate pvalue function from the factory function
        total_bg_hops = self.get_total_hops(background_tbl)
        total_expr_hops = self.get_total_hops(experiment_tbl)
        
        poisson_pval = self.poisson_pval_factory(total_bg_hops,
                                                 total_expr_hops,
                                                 poisson_pseudocount)
        # set up pvalue calculation options. Note that if ttaa_tbl is set to 
        # a tbl (rather than None), that field will be calculated/added in the 
        # loop
        pval_kwargs = {
            "scale_bg": True if total_bg_hops >= total_expr_hops else False,
            "scale_expr": True if total_bg_hops < total_expr_hops else False
        }
        
        sql_dict = {
            'aggregate_region': ' '.join(["SELECT COUNT(*) as hops"
                                   "FROM %s"
                                   "WHERE chr = %s AND start BETWEEN %s AND %s"])
        }
        def calculator(cur, chrom, region_start,region_stop):
            # extract hops in region
            expr_hops = cur.execute(sql_dict['aggregate_region'] \
                %(chrom,experiment_tbl,region_start, region_stop)).fetchall()
            bg_hops = cur.execute(sql_dict['aggregate_region'] \
                %(chrom,background_tbl,region_start, region_stop)).fetchall()
            # note that if this is set to None, num_ttaa is never set. if it 
            # is set to a valid tablename, then num_ttaa will be updated in 
            # each iteration
            if ttaa_tbl:
                # return the number of TTAA sites, or 1, whichever is greater
                pval_kwargs['num_ttaa'] = \
                    max(cur.execute(sql_dict['aggregate_region'] \
                %(chrom,ttaa_tbl,region_start,region_stop)).fetchall(),1)
            # calculate significance
            return poisson_pval(bg_hops, expr_hops, **pval_kwargs)              #pylint: disable=E1102
        return calculator

    # UNTESTED, UNREVISED, KONWN MISTAKES IN REGION SIG SETTINGS
    def range_score_macslike(self,background_tbl, experiment_tbl, ttaa_tbl, 
                             tile_width, sig_threshold, poisson_pseudocount = .2):      
        # check input to function
        for tbl in [background_tbl, experiment_tbl]:
            if not re.match(tbl, "^experiment_|^background_"):
                IOError(f"{tbl} must start with experiment_ or background_ "\
                        f"per HopsDB conventions")
        if tbl not in list(self.list_tables(self.con)):
            AttributeError(f"No table named {tbl}")
        # create significance table and table inserter
        sig_region_tablename = "_".join(["tileWidth", str(tile_width), 
                                         background_tbl,experiment_tbl])
        sig_region_col_dict = {'chr': 'TEXT', 'start':'INTEGER', 
                              'end':'INTEGER', 'poisson_pval':'NUMERIC'}        
        sig_table_inserter = self.new_table(
            sig_region_tablename, 
            sig_region_col_dict)
        # create poisson pvalue calculator
        region_sig_calculator = self.region_significance_calculator(
            background_tbl,
            experiment_tbl,
            ttaa_tbl,
            poisson_pseudocount)
        # iterate over chromosomes. Note that both chrom and seqlength are 
        # extracted from chr_map
        cur = self.con.cursor()
        seqinfo_sql = f"SELECT {self.standard_chr_format}, seqlength "\
                                   f"FROM {self.chr_map_table}"
        seqinfo = self.db_execute(cur,seqinfo_sql)
        for chrom, seqlength in seqinfo.items():
            
            min_expr_hop_position = cur.execute(f"SELECT MIN(end) as min_end "\
                                                f"FROM {experiment_tbl}").min_end
            max_expr_hop_position = cur.execute(f"SELECT MAX(end) as max_end "\
                                                f"FROM {experiment_tbl}").max_end

            iteration_start = min(0,min_expr_hop_position-tile_width)
            iteration_end = max(seqlength, max_expr_hop_position+tile_width)
            
            # instantiate variable to track significant regions
            sig_region_start = inf
            # iterate over tiles of width tile_width between iteration_start 
            # and iteration_end
            for window_start in range(iteration_start,iteration_end,tile_width):
                # last tile may be shorter than tile_width
                region_stop = min(window_start+tile_width, seqlength)
                region_pval = region_sig_calculator(                            #pylint: disable=E1102
                    cur,chrom,window_start,region_stop)
                # if the region pval is considered significant, update the 
                # sig region window (see initiation of this inner loop)
                if region_pval < sig_threshold:
                    sig_region_start = min(sig_region_start, window_start)
                # if the current region is not significant, but the previous 
                # region(s) were, add that (possibly combined) region, the hops 
                # and the significance to the database
                elif sig_region_start < window_start:
                    sig_region_pval = region_sig_calculator(                    #pylint: disable=E1102
                        cur,chrom,sig_region_start,window_start)
                    strand = "*"
                    values = [chrom, sig_region_start, window_start,strand,
                              sig_region_pval]
                    sig_table_inserter(values)                                  #pylint: disable=E1102
                # if the current region is not significant, and it is equal to 
                # or greater than the sig_region_start, then the last region(s) 
                # were also not significant and we just want to keep iterating 
                # over regions until we find a significant one. Nothing is done 
                # at this iteration in this event, hence the pass
                else:
                    pass
            # in the event that we reach the end of the seqlength without having 
            # entered the last (range of) significant region(s), add the region 
            # before moving on to the next seq/seqlength (outer loop)
            if sig_region_start < inf:
                sig_region_pval = region_sig_calculator(                        #pylint: disable=E1102
                    cur,chrom,sig_region_start,window_start)
                strand = "*"
                values = [chrom, sig_region_start, window_start,strand,
                            sig_region_pval]
                sig_table_inserter(values)                                      #pylint: disable=E1102
    
        return True