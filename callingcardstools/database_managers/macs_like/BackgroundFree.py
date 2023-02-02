from typing import Callable
import re
import logging 
import logging.config 
import yaml
from math import inf,floor,ceil

from callingcardstools.DatabaseApi import DatabaseApi

import pysqlite3 as sqlite3
import pandas as pd
import numpy as np
import scipy.stats as scistat

with open('logging.config.yaml', 'r') as f: #pylint:disable=W1514
    config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

logger = logging.getLogger(__name__)

# TODO move all sql statements into db_execute or wrap in the logging error handling

# discovered the 'row_factory' attribute in writing this. The connection 
# row factory is set to sqlite3.Row. This is another interesting one:
# see https://docs.python.org/3/library/sqlite3.html#connection-objects
# self.con.row_factory = lambda cursor, row: row[0]

class HopsDb(DatabaseApi):
    
    def significance_calculator(self,ttaa_tbl:str,experiment_tbl:str,
                                          background_region_size:int, 
                                          poisson_pseudocount:float) -> Callable[[str, int, int], float]:
        """A factory function to create a poisson pvalue calculator for a given 
        calling cards database.

        Args:
            ttaa_tbl (str): Name of the TTAA table
            experiment_tbl (str): Name of the experiment table for a given TF
            background_region_size (int): Width of the region to consider as background
            poisson_pseudocount (float): psuedocount to add when calculating the poisson pvalue

        Returns:
            Callable[[str, int, int], float]: A function which takes the chromosome, and a window_start and window_end 
            value and returns the poisson pvalue using constants set/calculated from the input to significance_calculator
        """
        aggregate_region = " ".join(["SELECT IFNULL(COUNT(*),0) as hops",
                                   "FROM %s",
                                   "WHERE chr = '%s' AND start BETWEEN %s AND %s"])
        def calculator(cur,chrom,window_start,window_stop):
            # NOTE there was a -1 on this in the original code, but I think 
            # that was to adjust for a 1 indexed TTAA (which the current is?) TTAA
            # may be. However, all qBed files should be 0 indexed, and that is 
            # the assumption here, hence no -1
            expected_region_start = window_start - (ceil(background_region_size/2))
            expected_region_stop  = window_stop + (ceil(background_region_size/2))
            sql_dict = {
                'background':{
                    'ttaa': aggregate_region \
                        %(ttaa_tbl,chrom, 
                        expected_region_start, expected_region_stop),
                    'hops': aggregate_region \
                        %(experiment_tbl, chrom, 
                        expected_region_start, expected_region_stop)
                },
                'window':{
                    'ttaa': aggregate_region \
                        %(ttaa_tbl,chrom,window_start, window_stop),
                    'hops': aggregate_region \
                        %(experiment_tbl,chrom,window_start, window_stop)

                }
            }
            # get ttaa in the expectation background
            ttaa_in_background = max(
                1,
                self.db_execute(
                    cur,
                    sql_dict['background']['ttaa']).fetchone()['hops'])
            # get hops in expectation background 
            hops_in_background = self.db_execute(
                cur,
                sql_dict['background']['hops']).fetchone()['hops'] 
            ttaa_in_window = max(
                1,
                self.db_execute(cur,sql_dict['window']['ttaa']).fetchone()['hops'])
            # extract hops in region
            hops_in_window = self.db_execute(
                cur,sql_dict['window']['hops']).fetchone()['hops']
            # set poisson parameters
            x = hops_in_window + poisson_pseudocount
            mu = ((hops_in_background / ttaa_in_background) * \
                ttaa_in_window) + poisson_pseudocount
            # return pvalue
            return 1-scistat.poisson.cdf(x, mu)
        return calculator

    def call_peaks(self, experiment_tbl, ttaa_tbl,
                                significance_threshold,
                                poisson_pseudocount = 0.2,
                                window_width = 1000, 
                                background_window_size = 100000, 
                                step_size = 500):
        # check input to function
        for tbl in [experiment_tbl, ttaa_tbl]:
            if not re.match("^experiment_|^background_|^ttaa_",tbl):
                IOError(f"{tbl} must start with experiment_ or background_ "\
                        f"per HopsDB conventions")
        if tbl not in list(self.list_tables(self.con)):
            AttributeError(f"No table named {tbl}")
        # create significance table and table inserter
        sig_region_tablename = "_".join(["tileWidth", str(window_width), 
                                         ttaa_tbl,experiment_tbl])
        sig_region_col_dict = {'chr': 'TEXT', 'start':'INTEGER', 
                              'end':'INTEGER', 'poisson_pval':'NUMERIC'}        
        sig_table_inserter = self.new_table(
            sig_region_tablename, 
            sig_region_col_dict)
        # create poisson pvalue calculator
        region_sig_calculator = self.significance_calculator(
            ttaa_tbl,
            experiment_tbl,
            background_window_size,
            poisson_pseudocount)
        # iterate over chromosomes. Note that both chrom and seqlength are 
        # extracted from chr_map
        cur = self.con.cursor()
        seqinfo_sql = f"SELECT {self.standard_chr_format}, seqlength "\
                                   f"FROM {self.chr_map_table}"
        seqinfo = self.db_execute(cur,seqinfo_sql).fetchall()
        for chrom, seqlength in seqinfo:
            
            min_max_sql = ' '.join(["SELECT MIN(end) AS min_end,",
                                    "max(END) as max_end",
                                    f"FROM {experiment_tbl}"])
            min_max_res = self.db_execute(cur,min_max_sql).fetchall()[0]
            # extract values
            min_expr_hop_position = min_max_res['min_end']
            max_expr_hop_position = min_max_res['max_end']
  
            iteration_start = max(
                0,
                floor(min_expr_hop_position-background_window_size/2))
            iteration_end = min(
                seqlength, 
                ceil(max_expr_hop_position+background_window_size/2))
            
            # instantiate variable to track significant regions
            sig_region_start = inf
            sig_region_end = 0
            # iterate over tiles of width tile_width between iteration_start 
            # and iteration_end
            print(f"working on {chrom}...")
            for window_start in range(iteration_start,iteration_end,step_size):
                window_end = window_start + window_width
                region_pval = region_sig_calculator(
                    cur, chrom, window_start,window_end)
                if region_pval <= significance_threshold:
                    # update if sig_region_start is infinity. otherwise, 
                    # leave it alone
                    if sig_region_start is inf:
                        sig_region_start = window_start
                    # update the sig window end
                    sig_region_end = window_end
                elif sig_region_start < window_start:
                    if sig_region_end < sig_region_start:
                        raise ValueError("Region end cannot be less than region start")
                    sig_region_pval = region_sig_calculator(
                        cur,chrom,sig_region_start,sig_region_end
                    )
                    values = ['"'+chrom+'"', sig_region_start, sig_region_end, sig_region_pval]
                    sig_table_inserter(values)
                    # reset sig_region boundaries
                    sig_region_start = inf
                    sig_region_end = 0
            # check to see if there is a sig region that hasn't been added
            if sig_region_start is not inf:
                if sig_region_end < sig_region_start:
                    raise ValueError("Region end cannot be less than region start")
                sig_region_pval = region_sig_calculator(
                    cur,chrom,sig_region_start,sig_region_end
                )
                values = ['"'+chrom+'"', sig_region_start, sig_region_end, sig_region_pval]
                sig_table_inserter(values)
