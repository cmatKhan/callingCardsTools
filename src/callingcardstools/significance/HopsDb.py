import os
import sys
import re
import logging
from math import inf,floor,ceil

import pysqlite3 as sqlite3
import pandas as pd
import numpy as np
import scipy.stats as scistat

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

# TODO move all sql statements into db_execute or wrap in the logging error handling

# discovered the 'row_factory' attribute in writing this. The connection 
# row factory is set to sqlite3.Row. This is another interesting one:
# see https://docs.python.org/3/library/sqlite3.html#connection-objects
# self.con.row_factory = lambda cursor, row: row[0]

class HopsDb():
    
    # set attributes
    require_regions = True
    require_chr_map = True
    
    db_loc = ""
    con = ""

    chr_map_table = "chr_map"

    pk = "id"

    standard_chr_format = 'ucsc'

    required_fields = {
        'qbed': ['chr', 'start', 'end', 'depth', 'strand', 'annotation', 'sample'],
        'bed6': ["chr", "start", "end", "name", "score", "strand"],
        chr_map_table: [standard_chr_format, 'seqlength']
    }

    index_col_string = ','.join(['"chr"', '"start" ASC', '"end" ASC', '"strand"'])

    def __init__(self, db_path):
        self.open(db_path)

    def __del__(self):
        try:
            # save changes
            self.con.commit()
            self.con.close()
        except AttributeError:
            pass

    @staticmethod
    def db_execute(cur,sql):
        try:
            return cur.execute(sql)
        except sqlite3.OperationalError as exc:                                 #pylint: disable=E1101
            msg = f"Could not execute sql: {sql}. Error: {exc}"
            logging.critical(msg+" ", exc_info=(sys.exc_info()))                #pylint: disable=W1201
            raise

    @staticmethod
    def list_tables(con):

        sql = ' '.join(["SELECT name",
                        "FROM sqlite_master", 
                        "WHERE type='table' or type='view';"])
        
        cur = con.cursor()

        try:
            table_list = cur.execute(sql).fetchall()
        except sqlite3.OperationalError as exc:                                 #pylint: disable=E1101
            msg = f"Could not execute sql: {sql}. Error: {exc}"
            logging.critical(msg+" ", exc_info=(sys.exc_info()))                #pylint: disable=W1201
            raise

        return [x[0] for x in table_list if x[0]  != "sqlite_sequence"]
    
    @staticmethod
    def list_fields(con, tablename):

        cur = con.cursor()
        sql = F"SELECT * FROM {tablename}"
        try:
            res = cur.execute(sql)
        except sqlite3.OperationalError as exc:                                 #pylint: disable=E1101
            msg = f"Could not execute sql: {sql}. Error: {exc}"
            logging.critical(msg+" ", exc_info=(sys.exc_info()))                #pylint: disable=W1201
            raise
        
        return [x[0] for x in res.description]
    

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

    @staticmethod
    def poisson_pval_factory(total_bg_hops, total_expr_hops, pseudocount):
        """Create a pvalue calculator. 

        Args:
            total_bg_hops (_type_): _description_
            total_expr_hops (_type_): _description_
            pseudocount (_type_): _description_

        Returns:
            function: A pvalue calculator which takes 
            (bg_hops, expr_hops, **kwargs) and returns a p-value
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
    def hypergeometric_pval_factory(total_bg_hops, total_expr_hops):
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
     
    def count_hops_factory(self,tbl):

        aggregate_sql = \
            f"SELECT COUNT(*) as hops " \
            f"FROM {tbl} " \
            f"WHERE (chr = %s AND start BETWEEN %s AND %s)"
        
        cur = self.con.cursor()
        
        def hops(chrom,region_start,region_end):
            res = self.db_execute(
                cur,
                aggregate_sql %(chrom,tbl,region_start, region_end))
                
            return res.fetchall()
        return hops
    

    def get_total_hops(self, tbl):
        if not re.match("^background|^experiment", tbl):
            raise AttributeError(f"get_total_hops expects to act on tables " \
                f"with the prefix background_ or experiment_ in the name. "
                f"Therefore, {tbl} is invalid")
        
        return pd.read_sql_query(f"SELECT COUNT(*) as total FROM {tbl}", 
                                 self.con).total[-1]
    
    def is_open(self):
        is_open = False
        try:
            self.con.cursor()
            is_open = True
        except Exception: #pylint: disable=W0703
            pass
        return is_open

    def open(self, db_path = None):
        
        # args is going to be either empty, or the first value is a path 
        # to the database
        # kwargs are keyword arguments for validate, right now region and chr_map

        if not db_path:
            # this allows for for an empty instantiation, or to open the 
            # database connection after it has been manually closed
            try:
                self.con.open()
                # see https://docs.python.org/3/library/sqlite3.html#connection-objects
                #con.row_factory = lambda cursor, row: row[0]
                self.con.row_factory = sqlite3.Row
            except AttributeError:
                pass
        else:
            # check path
            if not os.path.exists(db_path) and db_path != ':memory:':
                print(f"file DNE: {db_path}. Creating new database")
            # close the current con, if one is open
            try:      
                self.con.close()
            except AttributeError:
                pass

            # if memory, create a database in memory
            if db_path == ":memory:":
                self.db_loc = ":memory:"
                self.con = sqlite3.connect(db_path)                             #pylint: disable=E1101
                # see https://docs.python.org/3/library/sqlite3.html#connection-objects
                #con.row_factory = lambda cursor, row: row[0]
                self.con.row_factory = sqlite3.Row                              #pylint: disable=E1101
            # else, check the path and validate
            else:
                self.db_loc = db_path
                self.con = sqlite3.connect(db_path)                             #pylint: disable=E1101
                # see https://docs.python.org/3/library/sqlite3.html#connection-objects
                #con.row_factory = lambda cursor, row: row[0]
                self.con.row_factory = sqlite3.Row                              #pylint: disable=E1101
                self.validate()

    def close(self):
        self.con.close()
    
    def validate(self, con = None):

        if not con:
            con = self.con
        
        table_list = self.list_tables(con)

        if len(table_list) == 0:
            print("Database is empty")
            return True
        else:
            print("Checking table column names...")
            for table in table_list:
                fields = self.list_fields(con, table)
                if re.match(r"^background|^experiment", table):
                    if len(np.setdiff1d(self.required_fields['qbed'], fields)) > 0:
                        raise AttributeError(f"table {table} must have at " \
                            f"least the following fields "\
                                f"{','.join(self.required_fields['qbed'])}")
                if re.match(r"^regions", table):
                    # ensure that at least the bed6 fields exist
                    if len(np.setdiff1d(self.required_fields['bed6'], fields)):
                        raise AttributeError(f"table {table} does not have "\
                            f"the expected fields: "\
                                f"{','.join(self.required_fields['bed6'])}")

            if self.chr_map_table in table_list: 
                print("Standardizing chr names...")
                self.standardize_chr_format(con)
            
            print("Current database tables are valid")
            return True
    
    # TODO return inserter, updater, etc as an 'overloaded' function (via 
    # internal factory function)
    def new_table(self,tablename,col_dict):
        # drop the table if it exists
        drop_sql = f"""DROP TABLE IF EXISTS {tablename}""" 
        # create the table if it exists
        parsed_col_dict = ",".join([" ".join(['"'+k.strip()+'"',v.strip()]) \
            for k,v in col_dict.items()])
        create_sql = f"CREATE TABLE {tablename}("\
            f'"{self.pk}" INTEGER NOT NULL PRIMARY KEY, {parsed_col_dict})'
        
        cur = self.con.cursor()
        for sql in [drop_sql, create_sql]:
            self.db_execute(cur,sql) 
    
        def inserter(values, columns = col_dict.keys()):
            insert_sql = f"INSERT INTO {tablename} ({','.join(columns)}) "\
                         f"VALUES ({','.join([str(x) for x in values])})"
            self.db_execute(cur,insert_sql)
            self.con.commit()
        return inserter

    def index_table(self,tablename):
        cur = self.con.cursor()
        # create index on table
        index_sql = f"CREATE INDEX {tablename + '_index'} "\
                    f"ON {tablename} ({self.index_col_string})"
        self.db_execute(cur, index_sql)

    def add_frame(self, df, table_type, tablename_suffix=None):

        table_types_dict = {
            'bed6': ['regions'],
            'qbed': ['background', 'experiment', 'ttaa'],
            self.chr_map_table: [self.chr_map_table]
        }
    
        table_format = ""
        for field_format,types in table_types_dict.items():
            if table_type in types:
                table_format = field_format
    
        if table_format not in table_types_dict:
            expected_tbl_list = [x for sublist in table_types_dict.values() \
                for x in sublist]
            raise KeyError(f"table_type {table_type} not recognized. "\
                f"Available options are {','.join(expected_tbl_list)}")
        
        if len(np.setdiff1d(self.required_fields[table_format], df.columns)) > 0:
            raise AttributeError(f"The format for this type_type is "\
                f"{table_format}. The columns of the dataframe must therefore "\
                f"be a subset of {self.required_fields[table_format]}") 

        if table_type == self.chr_map_table:
            tablename = table_type
        else:
            tablename = '_'.join([table_type, 
                                  tablename_suffix.removeprefix("_")]).strip()

        # assign a cursor to the database at cur
        cur = self.con.cursor()
        # drop the table if it exists
        drop_sql = f"""DROP TABLE IF EXISTS {tablename}"""
        try:
            cur.execute(drop_sql)
        except sqlite3.OperationalError as exc:                                 #pylint: disable=E1101
            msg = f"Could not drop table with sql: {drop_sql}. Error: {exc}"
            logging.critical(msg+" ", exc_info=(sys.exc_info()))                #pylint: disable=W1201
            raise

        # create the table
        #create_tbl_sql = \
        #    ' '.join([f"CREATE TABLE {tablename}",
        #              f"({','.join([f'{self.pk} INTEGER'] + list(df.columns))},",
        #              f"PRIMARY KEY({self.pk} AUTOINCREMENT));"])
        
        col_dict = {k:"" for k in list(df.columns)}
        
        # create table
        self.new_table(tablename,col_dict)
        # create index
        if table_format in [x for x in self.required_fields if x != 'chr_map']:
            self.index_table(tablename) 
        # add the data
        df.to_sql(tablename,
                con=self.con,
                if_exists='append',
                index=False)
        
        return self.validate()
  
    def standardize_chr_format(self, con):

        swap_chr_format_sql = \
        " ".join(["UPDATE %s",
                 "SET chr = r.%s",
                 f"FROM (SELECT %s,%s FROM {self.chr_map_table}) AS r",
                 "WHERE %s.chr = r.%s;"])
        
        unique_chrnames_sql = "SELECT DISTINCT(chr) FROM %s"

        chr_map = pd.read_sql_query(
            f'select * from {self.chr_map_table}', 
            con).to_dict('list')
        chr_map = {k:v for k,v in chr_map.items() if k != 'id'}

        cur = self.con.cursor()

        for table in [x for x in self.list_tables(con) if x != self.chr_map_table]:
            if 'chr' in self.list_fields(con, table):
                res = cur.execute(unique_chrnames_sql % table).fetchall()
                current_unique_chr_names = [x[0] for x in res]
                # instantiate sentinel
                curr_chrom_format = -1
                # loop over colnames in chr_map_df. HALT if the current naming 
                # convention is discovered
                for format, map_chr_names in chr_map.items():
                    # check if all levels in the current chr_map_df[naming_convention]
                    # contain the `df[chrom_colname]`.unique() levels
                    if len(np.setdiff1d(current_unique_chr_names, map_chr_names)) == 0:
                        curr_chrom_format = format
                        break
                # if the current chromosome naming convention is not intuited above,
                # raise an error
                if curr_chrom_format == -1:
                    raise AttributeError("Chromosome names are not "+\
                        "recognized in a recognized format. Unique " +\
                            "chr names which cause error are: %s." 
                            %",".join(current_unique_chr_names))

                # if the current names are already in the requested format, return
                if not curr_chrom_format == self.standard_chr_format:
                    sql = swap_chr_format_sql %(table, 
                    self.standard_chr_format, curr_chrom_format, 
                    self.standard_chr_format, table, curr_chrom_format)
                    # execute swap
                    self.db_execute(cur,sql)

    def region_aggregate(self, regions_tbl):
        """Create view where hops are aggregated by chrom, pos and strand

        :param con: sqlite3 connection to database
        :type con: sqlie3 connection object
        :param tablename: name of the table from which to create the view
        :type tablename: str
        :param view_name: name of the view
        :type view_name: str
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
            f"IFNULL(hops_table.tmp_hops,0) AS hops " \
            f"FROM {regions_tbl} as regions_table " \
            f"LEFT JOIN (SELECT chr, start, end, COUNT(*) as tmp_hops " \
            f"FROM(SELECT r.chr AS chr, r.start AS start, r.end AS end, "\
                        f"r.strand AS strand, x.start AS insert_pos, "\
                        f"x.sample AS sample " \
                f"FROM {regions_tbl} as r LEFT JOIN {tbl} as x " \
	            f"WHERE (x.chr = r.chr AND x.start "\
                        f"BETWEEN r.start AND r.end)) as t " \
            f"GROUP BY sample, chr, start, end ORDER BY sample, chr, " \
            f"start ASC, end ASC) as hops_table " \
            f"ON regions_table.chr = hops_table.chr AND "\
            f"regions_table.start = hops_table.start AND " \
            f"regions_table.end = regions_table.end"

            self.db_execute(cur,create_view_sql)

        return True
     
    def region_score(self, regions_tbl, background_tbl, 
                     expr_tbl, poisson_pseudocount = 0.2,
                     return_frame = False):

        hop_views = {'background': regions_tbl + '_' + background_tbl,
                     'experiment': regions_tbl + '_' + expr_tbl}

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
        # join background and experiment hops tables
        # note that since the hops views are created such that all regions 
        # have a value, either 0 or the number of hops, this amounts to a 
        # full join -- both the regions and background will have the same number 
        # records (with different hops of course)
        # see aggregate_hops function for more details
        join_sql = f"SELECT l.chr AS chr, l.start AS start, l.end AS end, "\
                          f"l.hops AS bg_hops, r.hops AS expr_hops " \
              f"FROM {hop_views['background']} as l "\
              f"LEFT JOIN {hop_views['experiment']} AS r " \
              f"ON l.chr = r.chr AND l.start = r.start AND l.end = r.end"
        
        # get the tables
        quant_df = pd.read_sql_query(join_sql, self.con)

        quant_df['poisson_pval'] = \
            [poisson_pval(quant_df.loc[index,'bg_hops'], 
                          quant_df.loc[index,'expr_hops']) for \
                            index in range(len(quant_df))]
        
        quant_df['hypergeom_pval'] = \
            [hypergeom_pval(quant_df.loc[index, 'bg_hops'], 
                            quant_df.loc[index, 'expr_hops']) for \
                                    index in range(len(quant_df))]
        
        if return_frame:
            return quant_df
        else:
            sig_tablename = regions_tbl + '_' + background_tbl + \
                            "_" + expr_tbl + "_sig"
            quant_df.to_sql(sig_tablename,
                con=self.con,
                if_exists='replace',
                index=False)
    
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

    def region_significance_calculator_bf(self,ttaa_tbl,experiment_tbl,
                                          background_region_size, 
                                          poisson_pseudocount):
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


    def range_score_macslike_bf(self, experiment_tbl, ttaa_tbl,
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
        region_sig_calculator = self.region_significance_calculator_bf(
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