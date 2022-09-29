# stdlib
import re
from typing import Callable
import logging 
import logging.config 
import yaml
# local
from ...HopsDb import HopsDb
# outside
import pandas as pd
import scipy.stats as scistat

# with open('logging.config.yaml', 'r') as f: #pylint:disable=W1514
#     config = yaml.safe_load(f.read())
#     logging.config.dictConfig(config)

# logger = logging.getLogger(__name__)

class WithBackground(HopsDb):
    """A peak caller for yeast data. Requires background data"""

    @staticmethod
    def poisson_pval_factory(total_bg_hops:int, total_expr_hops:int, pseudocount:float) -> Callable[[int,int,dict], float]:
        """_summary_

        Args:
            total_bg_hops (int): _description_
            total_expr_hops (int): _description_
            pseudocount (float): _description_

        Returns:
            Callable[[int,int,dict], float]: _description_
        """

        # set constants
        total_bg_hops = float(total_bg_hops)
        total_expr_hops = float(total_expr_hops)
        hop_ratio = total_expr_hops / total_bg_hops
        def pval(bg_hops,expr_hops, **kwargs):
            """With the total_bg_hops and total_expr_hops set by the factory 
            function wrapper, this return function acts as a pvalue calculator.

            Args:
                bg_hops (int): True to scale the random variable by the hop_ratio 
                expr_hops (int): True to scale the random variable by the hop_ratio
                kwargs(dict): Recognizes the following key:value pairs:
                    - scale_expr = True to scale the random variable by the
                    hop_ratio. Default is False.
                    - scale_bg = True to scale lambda by the hop ratio. 
                    Default is false.
                    - ttaa = the number of ttaa sites in the region. 
                    Default is 1, and therefore has no effect if it is not 
                    passed in by the user

            Returns:
                float: A pvalue for the region
            """
            # extract values from kwargs and set defaults if dne
            scale_expr = kwargs.get("scale_background", False)
            scale_bg = kwargs.get("scale_experiment", False)
            num_ttaa = kwargs.get("ttaa", 1)

            if scale_expr and scale_bg:
                ValueError("Both scale_expr and scale_bg cannot be True")
            
            # if scale_bg is true, multiple bg_hops by the hop_ratio
            mu = (bg_hops * hop_ratio) / num_ttaa if scale_bg else total_bg_hops / num_ttaa
            
            # scale experiment
            # NOTE that the scaling in the original scripts FLIPPED the ratio 
            # from every other scaling (in the bg scaling and yeast)
            # I don't know if this is intended -- I don't understand why it would 
            # be this way -- so I used the hop_ratio defined above. BUT
            # this needs to be checked and quite possibly corrected back to the 
            # original scaling factor included below
            # 1-scistat.poisson.cdf(
            #    ((float(total_background_hops)/total_experiment_hops)*\
            #        num_exp_hops+pseudocounts),
            #        lambda_f*max(num_TTA                    value = 1-scistat.poisson.cdf(
            #           ((float(total_background_hops)/total_experiment_hops)*\
            #                num_exp_hops+pseudocounts),
            #                lambda_f*max(num_TTAAs,1)+pseudocounts)As,1)+pseudocounts)
            # if scale_expr is true, multiply the pseudocounted expr_hops by the hop_ratio
            x  = hop_ratio*(expr_hops + pseudocount) \
                if scale_expr else expr_hops + pseudocount

            # usage: scistat.poisson.cdf(x,mu)
            # where X are the experimental hops,
            # and mu is the background hops * (total_expr_hops/total_background_hops)
            return 1-scistat.poisson.cdf(x, mu)
        return pval

    @staticmethod
    def hypergeometric_pval_factory(total_bg_hops:int, total_expr_hops:int) -> Callable[[int,int], float]:
        """_summary_

        Args:
            total_bg_hops (int): _description_
            total_expr_hops (int): _description_

        Returns:
            Callable[[int,int], float]: _description_
        """
        # usage:
        # scistat.hypergeom.cdf(x,M,n,N)
        # where x is observed number of type I events
        # (white balls in draw) (experiment hops at locus)
        # M is total number of balls (total number of hops)
        # n is total number of white balls (total number of expeirment hops)
        # N is the number of balls drawn (total hops at a locus)
        M = total_bg_hops + total_expr_hops
        n = total_expr_hops
        
        def pval(bg_hops, expr_hops):
            x = expr_hops - 1
            N = bg_hops + expr_hops
            return 1-scistat.hypergeom.cdf(x,M,n,N)
        return pval

    def create_aggregate_view(self, regions_tbl:str) -> bool:
        """For each background and experiment table, create a view where the hops are 
        aggregated over the regions provided in the regions table

        Args:
            regions_tbl (str): name of the regions table to use to aggregate the hops

        Returns:
            bool: True if successful
        """

        # get table lists
        tbl_list = self.list_tables(self.con)
        background_and_expr_tbls = [x for x in tbl_list \
            if re.search(r"^background|^experiment",x)]
        available_regions = [x for x in tbl_list if x.startswith("regions_")]
        # check if region_tbl is in the regions table list
        if regions_tbl not in self.list_tables(self.con):
            AttributeError(f"No table named {regions_tbl} in database. "\
                           f"Current tables are: {','.join(available_regions)}")
        # create the views
        cur = self.con.cursor()
        for tbl in background_and_expr_tbls:

            viewname = regions_tbl + '_' +  tbl
            drop_view_sql = f'DROP VIEW IF EXISTS "main"."{viewname}";'
            self.db_execute(cur,drop_view_sql)
            
            # this sql is best read from the inside out starting at the second 
            # FROM statement.
            # first, the qbed (background or experiment) is left joined to the 
            # regions table. Note that this could be done either way, but right 
            # now, regions is on the left. Then the hops from the 
            # background/experiment qbed in the left/right bounds 
            # of the regions (inclusive [regions.start,regions.end] __on the start__ 
            # value of the qbed) are counted, which creates the 'hops_table'.
            # This is then left joined (regions table on the left, hops table  
            # on the right) back to the regions table on the chr, start, end 
            # fields, and if a region has hops, NULL is replaced by 0. The result 
            # is a table where all of the regions are present, and have either 0 
            # or the number of hops in that region
            # this needs to be done this way because there is no FULL OUTER JOIN
            # function in sql
            create_view_sql =  f"CREATE VIEW {viewname} AS "\
            f"SELECT regions_table.chr, regions_table.start, regions_table.end, "\
            f"hops_table.tmp_hops AS hops, hops_table.sample AS sample " \
            f"FROM {regions_tbl} as regions_table " \
            f"LEFT JOIN (SELECT chr, start, end, COUNT(*) as tmp_hops, sample " \
            f"FROM(SELECT r.chr AS chr, r.start AS start, r.end AS end, "\
                        f"r.strand AS strand, x.start AS insert_pos, "\
                        f"x.sample AS sample " \
                f"FROM {regions_tbl} as r LEFT JOIN {tbl} as x " \
	            f"WHERE (x.chr = r.chr AND x.start "\
                        f"BETWEEN r.start AND r.end)) as t " \
            f"GROUP BY sample, chr, start, end ORDER BY sample, chr, " \
            f"start ASC, end ASC) as hops_table " \
            f"ON regions_table.chr = hops_table.chr AND " \
            f"regions_table.start = hops_table.start AND " \
            f"regions_table.end = regions_table.end " \
            f"WHERE hops IS NOT NULL " \
            f"ORDER BY sample,regions_table.chr,regions_table.start,regions_table.end"
            # note the regions_table.end = regions_table.end does nothing. need to fix that

            self.db_execute(cur,create_view_sql)
            self.con.commit()

        return True

    def peak_caller(self, regions_tbl:str, background_tbl:str, 
                     expr_tbl:str, poisson_pseudocount:float = 0.2) -> None:
        """Call Peaks and add the result to the database.

        Args:
            regions_tbl (str): Name of the table from which extract regions over which to calculate significance
            background_tbl (str): Name of the background table
            expr_tbl (str): Name of the experiment (TF) table
            poisson_pseudocount (float, optional): peudocount to add to poisson pvalue calculation. Defaults to 0.2.

        Raises:
            AttributeError: _description_
        """
        # get the viewname
        hop_views = {'background': regions_tbl + '_' + background_tbl,
                     'experiment': regions_tbl + '_' + expr_tbl}
        # check that the views are in the database
        tbl_list = self.list_tables(self.con)
        for t in hop_views.values():
            if t not in tbl_list:
                raise AttributeError(f"No view called {t} exists. "\
                    f"Try calling region_aggregate() and then resubmit")

        # instantiate pvalue functions
        poisson_pval = self.poisson_pval_factory(
            self.get_total_hops(background_tbl),
            self.get_total_hops(expr_tbl),
            poisson_pseudocount)
        hypergeom_pval = self.hypergeometric_pval_factory(
            self.get_total_hops(background_tbl),
            self.get_total_hops(expr_tbl))
        # below I went through a number of options to do this. the first was to 
        # do a full outer join
        # sqlite does not have FULL OUTER JOIN. This emulates one.
        # see https://www.sqlitetutorial.net/sqlite-full-outer-join/
        # join_sql = " ".join(["SELECT b.chr AS chr, b.start AS start, b.end AS end,b.hops AS bg_hops, e.hops AS expr_hops",
        #                      f"FROM {hop_views['background']} as b",
        #                      f"LEFT JOIN {hop_views['experiment']} as e",
        #                      "UNION ALL",
        #                      "SELECT b.chr AS chr, b.start AS start, b.end AS end,b.hops AS bg_hops, e.hops AS expr_hops",
        #                      f"FROM {hop_views['experiment']} as e",
        #                      f"LEFT JOIN {hop_views['background']} as b",
        #                      "WHERE b.hops IS NULL;"])
        #
        # next I just reconfigured the view table so that a 0 would be entered 
        # if there was a null value. But, that didn't work with the sample names
        # so this won't work
        # join_sql = f"SELECT l.chr AS chr, l.start AS start, l.end AS end, "\
        #                   f"l.hops AS bg_hops, r.hops AS expr_hops " \
        #       f"FROM {hop_views['background']} as l "\
        #       f"LEFT JOIN {hop_views['experiment']} AS r " \
        #       f"ON l.chr = r.chr AND l.start = r.start AND l.end = r.end"

        # finally, I decided that the areas where the experiment is null actually 
        # don't matter -- no reason to save them, anyway.
        join_sql = " ".join(["SELECT e.chr AS chr, e.start AS start, e.end AS end,",
                             "IFNULL(b.hops,0) AS bg_hops, e.hops AS expr_hops,",
                             "e.sample as sample",
                             f"FROM {hop_views['experiment']} as e",
                             f"LEFT JOIN {hop_views['background']} as b USING(chr,start,end)",
                             "ORDER BY sample,chr,start,end"])
        
        # get the tables
        quant_df = pd.read_sql_query(join_sql, self.con)
        
        # calculate pvalues
        quant_df['poisson_pval'] = \
            [poisson_pval(quant_df.loc[index,'bg_hops'], 
                          quant_df.loc[index,'expr_hops']) for \
                            index in range(len(quant_df))]
        quant_df['hypergeom_pval'] = \
            [hypergeom_pval(quant_df.loc[index, 'bg_hops'], 
                            quant_df.loc[index, 'expr_hops']) for \
                                    index in range(len(quant_df))]

        # create the tablename
        sig_tablename = regions_tbl + '_' + background_tbl + \
                        "_" + expr_tbl + "_sig"

        # send the table to the database
        quant_df.to_sql(sig_tablename,
            con=self.con,
            if_exists='replace',
            index=False)
        
        # index the table
        index_col_string = self.index_col_string+',"sample"'
        self.index_table(sig_tablename, index_col_string)