import os
import sys
import sqlite3
import re
import logging

import pandas as pd
import numpy as np

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

class HopsDb():
    require_regions = True
    require_chr_map = True
    db_loc = ""
    con = ""

    chr_map_tablename = "chr_map"

    id_colname = "id"

    qbed_colnames = ['chr', 'start', 'end', 'depth', 'strand', 'annotation', 'sample']

    index_col_string = '("chrom", "start" ASC, "end" ASC, "strand");'

    table_type_dict = {
        
        'regions': {
        'index_col_string': ["chr", "start", "end", "name", "score", "strand"]},
        
        'background': {
        'index_col_string': qbed_colnames},

        'experiment': {
        'index_col_string': qbed_colnames}
    
    }

    standard_chr_format = 'ucsc'

    def __init__(self, *args, **kwargs):
        self.open(args, kwargs)

    def __del__(self):
        self.con.close()

    @staticmethod
    def list_tables(con):
        list_tables_sql = """SELECT name FROM sqlite_master WHERE type='table';"""

        con.cursor()
        cur = con.cursor()

        try:
            table_list = cur.execute(list_tables_sql).fetchall()
        except sqlite3.Error as exc:
            raise f"Error listing tables with sql: {list_tables_sql}" from exc
        
        return table_list
    
    @staticmethod
    def list_fields(con, tablename):
        # cite https://stackoverflow.com/a/7831685
        cur = con.cursor()
        cur.execute(f"select * from {tablename}")
        return [field[0] for field in cur.description]

    def open(self, *args, **kwargs):
        
        # args is going to be either empty, or the first value is a path 
        # to the database
        # kwargs are keyword arguments for validate, right now region and chr_map

        if not args[0]:
            # this allows for for an empty instantiation, or to open the 
            # database connection after it has been manually closed
            try:
                self.con.open()
            except AttributeError:
                pass
        else:
            db_path = args[0]
            # check path
            if not os.path.exists(db_path) and db_path != ':memory:':
                raise FileNotFoundError(f"file DNE: {db_path}")
            # close the current con, if one is open
            try:      
                self.con.close()
            except AttributeError:
                pass

            # open the connection and validate the database
            con = sqlite3.connect(db_path)
            self.validate(con, kwargs)

            self.db_loc = db_path
            self.con = con

    def close(self):
        self.con.close()

    def validate(self, con, **kwargs):

        chr_map_path = kwargs.get("chr_map", None)
        regions_path = kwargs.get('regions', None)

        self.require_regions = True if chr_map_path else False
        self.require_chr_map = True if regions_path else False
        
        table_list = self.list_tables(con)
        
        if self.required_tbls['chr_map'] not in table_list and self.require_chr_map:
            if chr_map_path:
                self.create_table(chr_map_path, 'chr_map')
            else:
                raise AttributeError("You must provide a chr_map frame or set chr_map to False")
        elif self.require_tbls['regions'] not in table_list and self.require_regions:
            if regions_path:
                self.create_table(regions_path, 'regions_path')
            else:
                raise AttributeError("You must provide a region frame or set regions to False")

        for table in table_list:
            fields = self.list_fields(con, table)
            if re.match(r"^background|^hops", table):
                if len(np.setdiff1d(fields, self.qbed_colnames)) > 0:
                    raise AttributeError(f"table {table} must have at " \
                        f"least the following fields {','.join(self.qbed_colnames)}")
            if re.match(r"^regions", table):
                # ensure that at least chr, start, stop, strand are included in the regions
                if len(np.setdiff1d(fields, self.qbed_colnames[0:3] +[self.qbed_colnames[4]]))
                
        self.standardize_chr_format(con)

        if self.require_regions == True:
            self.region_aggregate()        
  
    def standardize_chr_format(self, con):

        swap_chr_format_sql = \
        " ".join(["UPDATE %s",
                 "SET chr = %s",
                 "FROM (SELECT %s,%s FROM chr_map) AS r",
                 "WHERE %s.chr = r.%s;"])
        
        unique_chrnames_sql = "SELECT DISTINCT(chr) FROM %s"

        chr_map = pd.read_sql_query('select * from chr_map', con).to_dict('list')

        # see https://docs.python.org/3/library/sqlite3.html#connection-objects
        self.con.row_factory = lambda cursor, row: row[0]
        cur = self.con.cursor()

        for table in self.list_tables(con):
            if 'chr' in self.list_fields(con, table):
                current_unique_chr_names = \
                    cur.execute(unique_chrnames_sql % table).fetchall()
                # instantiate sentinel
                curr_chrom_format = -1
                # loop over colnames in chr_map_df. HALT if the current naming convention
                # is discovered
                for format, map_chr_names in chr_map:
                    # check if all levels in the current chr_map_df[naming_convention]
                    # contain the `df[chrom_colname]`.unique() levels
                    if len(np.setdiff1d(current_unique_chr_names, map_chr_names)) == 0:
                        curr_chrom_format = format
                        break
                # if the current chromosome naming convention is not intuited above,
                # raise an error
                if curr_chrom_format == -1:
                    raise AttributeError("Chromosome names are not "+\
                        "recognized in a recognized format. Unique chr names which cause "+\
                            " error are: %s." %",".join(current_unique_chr_names))

                # if the current names are already in the requested format, return
                if not curr_chrom_format == self.standard_chr_format:
                    cur.execute(swap_chr_format_sql % (table, self.standard_chr_format, 
                                                     'chr', self.standard_chr_format, 
                                                     table, 'chr'))
        # reset row_factory                  
        # see https://docs.python.org/3/library/sqlite3.html#connection-objects
        self.con.row_factory = None

    def create_table(self, df, table_type, tablename_suffix):

        if not table_type in self.table_type_dict:
            raise KeyError(f"table_type {table_type} not recognized. Available options are {self.table_type_dict.keys()}")
        
        tablename = table_type + "_" + tablename_suffix.removeprefix("_")

        # assign a cursor to the database at cur
        cur = self.con.cursor()
        # drop the table if it exists
        drop_sql = f"""DROP TABLE IF EXISTS {tablename}"""
        try:
            cur.execute(drop_sql)
        except sqlite3.OperationalError as exc:
            msg = f"Could not drop table with sql: {drop_sql}. Error: {exc}"
            logging.critical(msg+" ", exc_info=(sys.exc_info()))
            raise

        # create the table
        create_tbl_sql = f"""CREATE TABLE {tablename} 
                            ({','.join(['id INTEGER'] + list(df.columns))}, 
                            PRIMARY KEY({self.id_colname} AUTOINCREMENT));""" 
        try:
            cur.execute(create_tbl_sql)
        except sqlite3.OperationalError as exc:
            cur.execute(drop_sql)
            msg = f'Could not create table with sql:\n"{create_tbl_sql}".\nTable has been dropped.\n Error: {exc}.'
            logging.critical(msg+" ", exc_info=(sys.exc_info()))
            raise

        # create index on table
        index_sql = f"""CREATE INDEX {tablename + '_index'} ON {tablename} {self.index_col_string}"""
        try:
            cur.execute(index_sql)
        except sqlite3.OperationalError as exc:
            cur.execute(drop_sql)
            msg = f"""could not create index with sql:\n{index_sql}.\n Table has been dropped.\nError: {exc}:"""
            logging.critical(msg+" ", exc_info=(sys.exc_info()))
            raise
        
        # add the data
        df.to_sql(tablename,
                con=self.con,
                if_exists='append',
                index=False)

    def region_aggregate(self):
        """Create view where hops are aggregated by chrom, pos and strand

        :param con: sqlite3 connection to database
        :type con: sqlie3 connection object
        :param tablename: name of the table from which to create the view
        :type tablename: str
        :param view_name: name of the view
        :type view_name: str
        """
        cur = self.con.cursor()
        for tbl in self.list_tables(self.con):

            if re.match(r"^hops|^background", tbl):

                viewname = "region_" + tbl 
                drop_view_sql = f'DROP VIEW IF EXISTS "main"."{self.region_hop_view}";'

                cur.execute(drop_view_sql)

                create_view_sql =  (f"CREATE VIEW {viewname} AS "
                                    f"SELECT chr, start, end, COUNT(reads) as hops "
                                    f"FROM(SELECT r.chr AS chr, r.start AS start, "
                                            f"r.end AS end, r.strand AS strand, "
                                            f"x.start AS insert_pos, x.reads AS reads "
                                        f"FROM {tbl} as x "
                                        f"LEFT JOIN regions as r "
                                        f"WHERE (x.chr = p.chr AND "
                                                f"x.start BETWEEN "
                                                    f"r.start AND r.end)) as t "
                                    f"GROUP BY sample, chr, start, end "
                                    f"ORDER BY sample, chr, start ASC, end ASC;")

                try:
                    cur.execute(create_view_sql)
                except sqlite3.OperationalError as exc:
                    cur.execute(drop_view_sql)
                    msg = f"could not create view with sql: {create_view_sql}. view has been dropped. Error: {exc}:"
                    logging.critical(msg+" ", exc_info=(sys.exc_info()))
                    raise
    
    def append_enrichment_score_with_background(self, background_tbl, expr_tbl):

        # see https://docs.python.org/3/library/sqlite3.html#connection-objects
        # this means cur will be a row iterator that will be parse-able as a 
        # dict, eg row['chr']
        self.con.row_factory = sqlite3.Row
        cur = self.con.cursor()

        # get some hop numbers
        total_background_hops = pd.read_sql_query(f"SELECT COUNT(*) as total FROM {background_tbl}", self.con).total
        total_expr_hops       = pd.read_sql_query(f"SELECT COUNT(*) as total FROM {expr_tbl}", self.con).total
        hop_ratio             = float(total_expr_hops) / float(total_background_hops)
        bg_plus_expr_hops     = total_background_hops + total_expr_hops
        # get the tables
        promoter_bg_df   = pd.read_sql_query(f"SELECT * FROM {background_tbl}", self.con)
        promoter_expr_df = pd.read_sql_query(f"SELECT * FROM {expr_tbl}", self.con)
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
            
        return q