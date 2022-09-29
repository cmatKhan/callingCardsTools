import os
import sys
import re
import logging 
from typing import Callable, Iterator
from math import inf,floor,ceil

import pysqlite3 as sqlite3
from pysqlite3.dbapi2 import ProgrammingError
import pandas as pd
import numpy as np

logging.getLogger(__name__).addHandler(logging.NullHandler())

# TODO move all sql statements into db_execute or wrap in the logging error handling

# discovered the 'row_factory' attribute in writing this. The connection 
# row factory is set to sqlite3.Row. This is another interesting one:
# see https://docs.python.org/3/library/sqlite3.html#connection-objects
# self.con.row_factory = lambda cursor, row: row[0]

class HopsDb():
    """An object to aid in creating,modifying and using calling cards data from a database backend"""
    
    _db_loc = ""
    _con = ""

    _chr_map_table = "chr_map"

    _pk = "id"

    _standard_chr_format = 'ucsc'

    _required_fields = {
        'qbed': ['chr', 'start', 'end', 'depth', 'strand', 'annotation', 'sample'],
        'bed6': ["chr", "start", "end", "name", "score", "strand"],
        _chr_map_table: [_standard_chr_format, 'seqlength']
    }
    _required_fields['qbed_dtypes'] = {k:v for k,v in zip(_required_fields['qbed'], 
                                       ['str','int','int','int','str','int','int'])}
    _required_fields['bed6_dtypes'] = {k:v for k,v in zip(_required_fields['qbed'], 
                                       ['str','int','int','str','float','str'])}

    _index_col_string = ','.join(['"chr"', '"start" ASC', '"end" ASC', '"strand"'])

    def __init__(self, db_path:str) -> None:
        """_summary_

        Args:
            db_path (str): either a path to a database file -- if one DNE, 
            then a sqlite db will be created at that location -- or :memory: for 
            an in memory database
        """
        self.open(db_path)

    def __del__(self):
        try:
            # save changes
            self.con.commit()
            self.con.close()
        except AttributeError:
            pass
        except ProgrammingError:
            pass

    @property
    def db_loc(self):
        """filepath to the database (sqlite) or address of database"""
        return self._db_loc
    @db_loc.setter
    def db_loc(self, new_db_loc):
        self._db_loc = new_db_loc
    
    @property
    def con(self):
        """connection to the database"""
        return self._con
    @con.setter
    def con(self, new_con):
        self._con = new_con

    @property
    def chr_map_table(self):
        """name of the table which stores the mapping between chr naming conventions"""
        return self._chr_map_table
    
    @property
    def pk(self): #pylint: disable=C0103
        """name to use for the primary field key in tables in the database"""
        return self._pk
    
    @property
    def standard_chr_format(self):
        """A field in the chr_map_table to use as the standard chr naming format for all tables"""
        return self._standard_chr_format

    @property
    def required_fields(self):
        """A dict which describes the required fields for various table formats. Currently defined: qBed, bed6, chr_map"""
        return self._required_fields
 
    @property
    def index_col_string(self):
        """string to use to create the index on indexed tables in the database"""
        return self._index_col_string

    @staticmethod
    def zero_index(df: pd.DataFrame) -> pd.DataFrame:
        """_summary_

        Args:
            df (pandas.DataFrame): _description_

        Returns:
            pandas.DataFrame: _description_
        """
        df[['start', 'end']] = \
            df.apply(lambda row: [row['start']-1, row['end']-1], 
            axis = 1, result_type='expand')
        
        return df

    @staticmethod
    def db_execute(cur,sql:str) -> int:
        """_summary_

        Args:
            cur (_type_): _description_
            sql (str): _description_

        Returns:
            int: _description_
        """
        try:
            return cur.execute(sql)
        except sqlite3.OperationalError as exc:                                 #pylint: disable=E1101
            msg = f"Could not execute sql: {sql}. Error: {exc}"
            logging.critical(msg+" ", exc_info=(sys.exc_info()))                #pylint: disable=W1201
            raise

    @staticmethod
    def list_tables(con: sqlite3.Connection) -> list: #pylint: disable=E1101
        """_summary_

        Args:
            con (sqlite3.Connection): _description_

        Returns:
            list: _description_
        """

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
    def list_fields(con, tablename:str)-> list:
        """_summary_

        Args:
            con (_type_): _description_
            tablename (str): _description_

        Returns:
            list: _description_
        """

        cur = con.cursor()
        sql = F"SELECT * FROM {tablename}"
        try:
            res = cur.execute(sql)
        except sqlite3.OperationalError as exc:                                 #pylint: disable=E1101
            msg = f"Could not execute sql: {sql}. Error: {exc}"
            logging.critical(msg+" ", exc_info=(sys.exc_info()))                #pylint: disable=W1201
            raise
        
        return [x[0] for x in res.description]
    
     
    def count_hops_factory(self,tbl:str) -> Callable[[],int]:
        """count number of rows (hops) in a given region of a table

        Args:
            tbl (str): name of the table on which to calculate hops

        Raises:
            AttributeError: _description_

        Returns:
            Callable[[],int]: _description_
        """

        if not re.match("^background|^experiment", tbl):
            raise AttributeError(f"count hops expects to act on tables " \
                f"with the prefix background_ or experiment_ in the name. "
                f"Therefore, {tbl} is invalid")

        aggregate_sql = \
            f"SELECT COUNT(*) as hops " \
            f"FROM {tbl} " \
            f"WHERE (chr = %s AND start BETWEEN %s AND %s)"
        
        cur = self.con.cursor()
        
        def hops(chrom,region_start,region_end):
            # TODO make this return an int value
            res = self.db_execute(
                cur,
                aggregate_sql %(chrom,tbl,region_start, region_end))
                
            return res.fetchall()
        return hops
    

    def get_total_hops(self, tbl:str) -> int:
        """Get the total number of hops for a given table.

        Total hops is defined as the number of rows in given background/experiment 
        table 

        Args:
            tbl (str): name of the table

        Raises:
            AttributeError: raised if the tablename does not start with background or experiment

        Returns:
            int: the number of hops (rows) in the given table
        """
        if not re.match("^background|^experiment", tbl):
            raise AttributeError(f"get_total_hops expects to act on tables " \
                f"with the prefix background_ or experiment_ in the name. "
                f"Therefore, {tbl} is invalid")
        
        return int(pd.read_sql_query(f"SELECT COUNT(*) as total FROM {tbl}", 
                                 self.con).total)
    
    def is_open(self) -> bool:
        """check if the database connection is open

        Returns:
            bool: True if the database connection is open
        """
        is_open = False
        try:
            self.con.cursor()
            is_open = True
        except Exception: #pylint: disable=W0703
            pass
        return is_open

    def open(self, db_path:str = None) -> None:
        """Open a database connection

        Args:
            db_path (str, optional): Either a path or :memory: -- a string used 
            to open a database connection. Defaults to None.
        """
        
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
                self.con = sqlite3.connect(db_path)                            #pylint: disable=E1101
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
        """Close the database connection"""
        self.con.close()
    
    def validate(self, con: sqlite3.Connection = None) -> bool: #pylint:disable=E1101
        """A function to validate the database for expected structure/tables,etc

        Args:
            con (sqlite3.Connection, optional): connection to a database. Defaults to None, which will validate self.con.

        Raises:
            AttributeError: raised if an expectation on a table is not met

        Returns:
            bool: True if successful
        """

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
 
            print("Current database tables are valid")
            return True
    
    # TODO return inserter, updater, etc as an 'overloaded' function (via 
    # internal factory function)
    def new_table(self,tablename:str,col_dict:dict) -> Callable[[list,list],int]:
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
    
        def inserter(values:list, columns:list = col_dict.keys()) -> int:
            """_summary_

            Args:
                values (list): _description_
                columns (list, optional): _description_. Defaults to col_dict.keys().

            Returns:
                int: _description_
            """
            insert_sql = f"INSERT INTO {tablename} ({','.join(columns)}) "\
                         f"VALUES ({','.join([str(x) for x in values])})"
            self.db_execute(cur,insert_sql)
            self.con.commit()
        return inserter

    def index_table(self,tablename:str,index_col_string:list = None) -> int:
        """_summary_

        Args:
            tablename (str): _description_

        Returns:
            int: _description_
        """
        if not index_col_string:
            index_col_string = self.index_col_string
        cur = self.con.cursor()
        # create index on table
        index_sql = f"CREATE INDEX {tablename + '_index'} "\
                    f"ON {tablename} ({index_col_string})"
        self.db_execute(cur, index_sql)
        self.con.commit()

    def add_frame(self, df: pd.DataFrame, table_type:str, tablename_suffix:str=None, drop:bool=False) -> bool:
        """Add a pandas dataframe to the database

        Args:
            df (pandas.Dataframe): A pandas dataframe of the table
            table_type (str): One of regions,background,experiment,ttaa,chr_map. 
            The choice dictates what fields are expected. See the class attribute required_fields for more details.
            tablename_suffix (str): A suffix to add to the tablename, eg the TF in the case of 
            experimental hops. If you expect multiple tables of a similar type (eg multiple background 
            or more likely multiple experimental hops tables for mutiple TFs, use this argument)
            drop (bool): whether to drop an existing table. Default is False.

        Raises:
            OperationalError: raised if the upload to the database is unsuccessful
        
        Returns:
            bool: True if the table upload is successful
        """

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
        if drop:
            try:
                # drop the table if it exists
                drop_sql = f"""DROP TABLE IF EXISTS {tablename}"""
                cur.execute(drop_sql)
                self.con.commit()
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
        if tablename not in self.list_tables(self.con):
            self.new_table(tablename,col_dict)
            # create index
            if table_format in [x for x in self.required_fields if x != 'chr_map']:
                self.index_table(tablename) 
        # add the data
        df.to_sql(tablename,
                con=self.con,
                if_exists='append',
                index=False)
        # standardize the chromosomes if there is a column called chr
        if self.chr_map_table in self.list_tables(self.con): 
            print("Standardizing chr names...")
            self.standardize_chr_format(tablename)
        # do some more validating
        return self.validate()
  
    def standardize_chr_format(self,table:str) -> None: #pylint: disable=E1101
        """Use the chr_map table to standardize all 'chr' columns to the same naming format

        Args:
            table (str): name of the table to standardize

        Raises:
            AttributeError: raised if the chromosomes are not fully described by a chr_map field
        """

        swap_chr_format_sql = \
        " ".join(["UPDATE %s",
                 "SET chr = r.%s",
                 f"FROM (SELECT %s,%s FROM {self.chr_map_table}) AS r",
                 "WHERE %s.chr = r.%s;"])
        
        unique_chrnames_sql = "SELECT DISTINCT(chr) FROM %s"

        chr_map = pd.read_sql_query(
            f'select * from {self.chr_map_table}', 
            self.con).to_dict('list')
        chr_map = {k:v for k,v in chr_map.items() if k != 'id'}

        cur = self.con.cursor()

        if 'chr' in self.list_fields(self.con, table):
            res = cur.execute(unique_chrnames_sql % table).fetchall()
            current_unique_chr_names = [str(x[0]) for x in res]
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
                self.con.commit()