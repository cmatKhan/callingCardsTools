# pylint: disable=W1201,W1203,C0103
import os
import sys
import re
import logging
from typing import Callable, Literal

import pysqlite3 as sqlite
from pysqlite3.dbapi2 import ProgrammingError #pylint: disable=E0611
import pandas as pd
import numpy as np

__all__ = ['DatabaseApi']

logging.getLogger(__name__).addHandler(logging.NullHandler())

# https://github.com/pandas-dev/pandas/issues/14553#issuecomment-778599042
# class for SqlUpsert. Also consider SQLalchemy or that sql package i sent to 
# daniel

# discovered the 'row_factory' attribute in writing this. The connection
# row factory is set to sqlite.Row. This is another interesting one:
# see https://docs.python.org/3/library/sqlite.html#connection-objects
# self.con.row_factory = lambda cursor, row: row[0]

class DatabaseApi():
    """An object to aid in creating,modifying and using calling cards data from a database backend"""

    _db_loc = ""
    _db_timeout=1200
    _con = ""

    _chr_map_table = "chr_map"

    _pk = "id"

    _standard_chr_format = 'ucsc'

    _required_fields = {
        'qbed': ['chr', 'start', 'end', 'depth', 'strand', 'batch_id'],
        'bed6': ["chr", "start", "end", "name", "score", "strand"],
        'bed3': ["chr", "start", "end"],
        _chr_map_table: [_standard_chr_format, 'seqlength']
    }
    _required_fields['qbed_dtypes'] = {k: v for k, v in zip(_required_fields['qbed'],
                                       ['str', 'int', 'int', 'int', 'str', 'int', 'int'])}
    _required_fields['bed6_dtypes'] = {k: v for k, v in zip(_required_fields['qbed'],
                                       ['str', 'int', 'int', 'str', 'float', 'str'])}
    
    _table_types_dict = {
        'bed3': {'regions'},
        'qbed': {'background', 'experiment', 'ttaa'},
        _chr_map_table: [_chr_map_table]
    }

    _index_col_string_dict = {
        'bed3': ','.join(['"chr"', '"start" ASC', '"end" ASC']),
        'bed6': ','.join(['"chr"', '"start" ASC', '"end" ASC', '"strand"']),
        'qbed': ','.join(['"chr"', '"start" ASC', '"end" ASC', '"strand"', '"batch_id"'])
    }

    def __init__(self, db_path: str, validate: bool = True) -> None:
        """Constructor

        Args:
            db_path (str): either a path to a database file -- if one DNE, 
            then a sqlite db will be created at that location -- or :memory: for 
            an in memory database
            validate (bool, optional): Set to False to skip database validation 
            after opening from file. Defaults to True.
        """
        self.open(db_path, validate)

    def __del__(self):
        try:
            # attempt commit prior to deleting
            logging.info("Attempting to commit any changes before deleting hopsdb object")
            self.con.commit()
        except AttributeError as err:
            logging.info(f"connection not valid -- nothing commited: {err}")
        except ProgrammingError as err:
            logging.info(f"connection not valid -- nothing commited: {err}")
        # note that this is error handled in the close method
        self.con.close()

    @property
    def db_loc(self):
        """filepath to the database (sqlite) or address of database. Default to 20 minutes"""
        return self._db_loc

    @db_loc.setter
    def db_loc(self, new_db_loc):
        self._db_loc = new_db_loc
    
    @property
    def db_timeout(self):
        """Set the connection timeout limit before an error is thrown. see sqlite3.connection docs"""
        return self._db_timeout
    @db_timeout.setter
    def db_timeout(self, num_secs):
        self._db_timeout = num_secs

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
    def pk(self):
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
    def table_types_dict(self):
        """A dictionary with keys as table_format defined in keys of required_fields 
        and values a list of type_types. Note that the setter adds keys and values, 
        and does not remove any. Hence the default init set is always there"""
        return self._table_types_dict

    @table_types_dict.setter
    def table_types_dict(self, new_dict):
        # iterate over the keys and values in new_dict
        for k,v in new_dict.items():
            # if a given new key already exists, extract the values
            curr_tables = self.table_types_dict.get(k, None)
            if curr_tables:
                # if they key exists, then add the new tables and update the 
                # types_dict attribute
                self._table_types_dict[k] = curr_tables.add(v)
            else:
                # else, add an entirely key key:value pair
                self._table_types_dict[k] = set(v)

    @property
    def index_col_string_dict(self):
        """string to use to create the index on indexed tables in the database"""
        return self._index_col_string_dict

    @staticmethod
    def db_execute(cur: sqlite.Cursor, sql: str) -> int: #pylint: disable=E1101
        """A convenience wrapper to execute a given command. Has error handling. 

        Args:
            cur (sqlite.Cursor): A cursor to use to execute the statement
            sql (str): sql statement to execute

        Returns:
            int: Number of rows affected
        """
        try:
            return cur.execute(sql)
        except sqlite.OperationalError as exc:  # pylint: disable=E1101
            msg = f"Could not execute sql: {sql}. Error: {exc}"
            logging.critical(msg+" ", exc_info=(sys.exc_info()))
            raise

    @staticmethod
    def list_tables(con: sqlite.Connection) -> list:  # pylint: disable=E1101
        """Given a sqlite connection, list all tables

        Args:
            con (sqlite.Connection): Connection to a sqlite database

        Returns:
            list: A list of the tables in the database
        """

        sql = ' '.join(["SELECT name",
                        "FROM sqlite_master",
                        "WHERE type='table' or type='view';"])

        cur = con.cursor()

        try:
            table_list = cur.execute(sql).fetchall()
        except sqlite.OperationalError as exc:  # pylint: disable=E1101
            msg = f"Could not execute sql: {sql}. Error: {exc}"
            logging.critical(msg+" ", exc_info=(sys.exc_info()))
            raise

        return [x[0] for x in table_list if x[0] != "sqlite_sequence"]

    @staticmethod
    def list_fields(con: sqlite.Connection, tablename: str) -> list: #pylint: disable=E1101
        """Given a connection and a tablename, list all fields in the table

        Args:
            con (sqlite.Connection): connection to a sqlite database
            tablename (str): name of a table in the database

        Returns:
            list: A list of the fields for a given table
        """

        cur = con.cursor()
        sql = F"SELECT * FROM {tablename}"
        try:
            res = cur.execute(sql)
        except sqlite.OperationalError as exc:  # pylint: disable=E1101
            msg = f"Could not execute sql: {sql}. Error: {exc}"
            logging.critical(msg+" ", exc_info=(sys.exc_info()))
            raise

        return [x[0] for x in res.description]
    
    def get_batch_id(self,batch:str,tf:str,replicate:str) -> int:
        """_summary_

        Args:
            batch (str): _description_
            tf (str): _description_
            replicate (str): _description_

        Raises:
            ValueError: _description_

        Returns:
            int: _description_
        """

        sql = f"SELECT id FROM batch WHERE batch = '{batch}' "\
            f"AND tf = '{tf}' AND replicate = '{replicate}'"
        
        cur = self.con.cursor()
        
        res = self.db_execute(cur,sql).fetchall()

        if len(res) > 1:
            raise ValueError('More than 1 record returned')
        
        return res[0]['id']

    def count_hops_factory(self, tbl: str) -> Callable[[], int]:
        """count number of rows (hops) in a given region of a table

        Args:
            tbl (str): name of the table on which to calculate hops

        Raises:
            AttributeError: _description_

        Returns:
            Callable[[],int]: _description_
        """

        if not re.match("^background|^experiment", tbl):
            raise AttributeError(f"count hops expects to act on tables "
                                 f"with the prefix background_ or experiment_ in the name. "
                                 f"Therefore, {tbl} is invalid")

        aggregate_sql = \
            f"SELECT COUNT(*) as hops " \
            f"FROM {tbl} " \
            f"WHERE (chr = %s AND start BETWEEN %s AND %s)"

        cur = self.con.cursor()

        def hops(chrom, region_start, region_end):
            # TODO make this return an int value
            res = self.db_execute(
                cur,
                aggregate_sql % (chrom, tbl, region_start, region_end))

            return res.fetchall()
        return hops

    def get_total_hops(self, tbl: str, batch_id: str = None) -> int:
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
            raise AttributeError(f"get_total_hops expects to act on tables "
                                 f"with the prefix background_ or experiment_ in the name. "
                                 f"Therefore, {tbl} is invalid")

        sql = f"SELECT COUNT(*) as total FROM {tbl}"

        if batch_id:
            sql = " ".join([sql, f"WHERE batch_id = '{batch_id}'"])

        return int(pd.read_sql_query(sql, self.con).total)

    def is_open(self) -> bool:
        """check if the database connection is open

        Returns:
            bool: True if the database connection is open
        """
        is_open = False
        try:
            self.con.cursor()
            is_open = True
        except Exception:  # pylint: disable=W0703
            pass
        return is_open

    def open(self, db_path: str = None, validate: bool = False) -> None:
        """Open a database connection

        Args:
            db_path (str, optional): Either a path or :memory: -- a string used 
            to open a database connection. Defaults to None.
            validate (bool, optional): Whether or not to validate the database 
            after opening. Defaults to False
        """

        # args is going to be either empty, or the first value is a path
        # to the database
        # kwargs are keyword arguments for validate, right now region and chr_map

        if not db_path:
            # this allows for for an empty instantiation, or to open the
            # database connection after it has been manually closed
            try:
                self.con.open()
                # see https://docs.python.org/3/library/sqlite.html#connection-objects
                #con.row_factory = lambda cursor, row: row[0]
                self.con.row_factory = sqlite.Row
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
                self.con = sqlite.connect(db_path,timeout=self.db_timeout)  # pylint: disable=E1101
                # see https://docs.python.org/3/library/sqlite.html#connection-objects
                #con.row_factory = lambda cursor, row: row[0]
                self.con.row_factory = sqlite.Row  # pylint: disable=E1101
            # else, check the path and validate
            else:
                self.db_loc = db_path
                self.con = sqlite.connect(db_path,timeout=self.db_timeout)  # pylint: disable=E1101
                # see https://docs.python.org/3/library/sqlite.html#connection-objects
                #con.row_factory = lambda cursor, row: row[0]
                self.con.row_factory = sqlite.Row  # pylint: disable=E1101
                if validate:
                    self.validate()

    def close(self):
        """Close the database connection"""
        try:
            self.con.close()
        except AttributeError as err:
            logging.info(f"Failed attempt to close database. Likely not a valid DB. Error: {err}")
        except ProgrammingError as err:
            logging.info(f"Failed attempt to close database. Likely not a valid DB. Error: {err}")

    def validate(self, con: sqlite.Connection = None) -> bool:  # pylint:disable=E1101
        """A function to validate the database for expected structure/tables,etc

        Args:
            con (sqlite.Connection, optional): connection to a database. Defaults to None, which will validate self.con.

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
                    if set(self.required_fields['qbed']) - set(fields) != set():
                        raise AttributeError(f"table {table} must have at "
                                             f"least the following fields "
                                             f"{','.join(self.required_fields['qbed'])}")
                # if the table is strictly a regions table -- not accumulated over
                # background or experiment
                if re.match(r"^regions", table) and not len(re.findall(r"background|experiment", table)) > 0:
                    # ensure that at least the bed6 fields exist
                    if set(self.required_fields['bed3']) - set(fields) != set():
                        raise AttributeError(f"table {table} does not have "
                                             f"the expected fields: "
                                             f"{','.join(self.required_fields['bed3'])}")

            print("Current database tables are valid")
            return True

    # TODO return inserter, updater, etc as an 'overloaded' function (via
    # internal factory function)
    def new_table(self, tablename: str, col_dict: dict, 
                  fk_tablelist:list = None, on_delete_cascade:bool = True, 
                  clean:bool = False) -> Callable[[list, list], int]:
        """Create a new table in the database
        
        Args:
            tablename (str): name of the new table
            col_dict (dict): dictionary where the keys are field names and values are data types
            fk_tablelist (list): list of tables to which to foreign key, Note that the 
            field name table_id for each table in the fk_tablelist must exist, and the table must 
            already exist in the database 
            on_delete_cascade (bool): whether to delete the referencing record when the referenced record is deleted. 
            True by default.
            clean (bool): whether to drop an existing table of the same name and create a new (empty) table. 
            Defaults to False

        Returns:
            Callable[[list,list], int]: A function to insert values into the new table. 
            Note that this is in development -- there are some problems with appropriately 
            quoting columns and character values
        """
        if not tablename in self.list_tables(self.con) or clean:
            # drop the table if it exists
            drop_sql = f"""DROP TABLE IF EXISTS {tablename}"""
            # create the table if it exists
            parsed_col_dict = ",".join([" ".join(['"'+k.strip()+'"', v.strip()])
                                        for k, v in col_dict.items()])
            # create the foreign key constraint if a fk_tablelist is passed
            fk_constraints = []
            if fk_tablelist:
                for table in fk_tablelist:
                    if table not in self.list_tables(self.con):
                        raise AttributeError(f"Tables in the foreign key list must exist in the database. {table} DNE.")
                    elif table+"_id" not in col_dict:
                        raise AttributeError(f"{table+'_id'} does not exist as a field in col_dict. Foreign key columns must have format <fk_table>_id and exist in col_dict")
                    else:
                        fk_str = f"FOREIGN KEY ({table}_id) REFERENCES {table}(id)"
                        if on_delete_cascade:
                            fk_str = fk_str + ' ON DELETE CASCADE'
                        fk_constraints.append(fk_str)
            # turn the fk_constraints list into a sql statement, or an empty string
            # if there is on fk_tablelist
            fk_sql = ","+",".join(fk_constraints) if len(fk_constraints)>0 else ""
            
            # construct the create table sql
            create_sql = "".join([f"CREATE TABLE {tablename} (",
            ",".join([f'"{self.pk}" INTEGER PRIMARY KEY AUTOINCREMENT',parsed_col_dict]),
            fk_sql, ")"])
            # execute
            cur = self.con.cursor()
            for sql in [drop_sql, create_sql]:
                self.db_execute(cur, sql)
                self.con.commit()
        
        def inserter(values: list, columns: list = col_dict.keys()) -> int:
            """new_table() both creates a table, and returns a callable function 
            which takes a list of values and the columns into which to insert. By 
            default, the columns are all columns other than the primary key. Note 
            that this function will commit the updates itself.

            Args:
                values (list): list of values to insert into a set of columns
                columns (list, optional): list of columns into which to insert. 
                Defaults to col_dict.keys().

            Returns:
                int: Number of rows affected (I think)
            """
            # TODO this won't work as stands -- need some way of handling data type 
            # and appropriate quoting
            insert_sql = f"INSERT INTO {tablename} ({','.join(columns)}) "\
                        f"VALUES ({','.join([str(x) for x in values])})"
            # update_sql = " ".join([f"UPDATE {tablename}",])
            logging.debug(f"attempting to insert: {insert_sql}")
            self.db_execute(cur, insert_sql)
            self.con.commit()
        return inserter
    
    def index_table(self, tablename: str, index_col_string: list) -> int:
        """_summary_

        Args:
            tablename (str): _description_

        Returns:
            int: _description_
        """
        cur = self.con.cursor()
        # create index on table
        index_sql = f"CREATE INDEX IF NOT EXISTS {tablename + '_index'} "\
                    f"ON {tablename} ({index_col_string})"
        self.db_execute(cur, index_sql)
        self.con.commit()

    def add_frame(self, df: pd.DataFrame, table_format:str, 
                  tablename:str = None , table_type: str = None, 
                  tablename_suffix: str = None, drop: bool = False, 
                  fk_tablelist:list = None) -> bool:
        """Add a pandas dataframe to the database

        Args:
            df (pandas.Dataframe): A pandas dataframe of the table
            table_format (str): Format for the table. Must be a key in the dict attribute
            required_fields. This argument will determine what required fields are checked 
            to ensure that they exist in the input df.
            tablename (str, optional): The name of the new table in the DB. However, if you're 
            adding an expected type, eg background,experiment,regions,or ttaa, use the two 
            opens below to allow the function to build the appropriate tablename
            table_type (str): One of regions,background,experiment,ttaa,chr_map. 
            The choice dictates what fields are expected. See the class attribute required_fields for more details.
            tablename_suffix (str, optional): A suffix to add to the tablename, eg the TF in the case of 
            experimental hops. If you expect multiple tables of a similar type (eg multiple background 
            or more likely multiple experimental hops tables for mutiple TFs, use this argument)
            drop (bool, optional): whether to drop an existing table. Default is False.
            fk_table (str, optional): Name of a table onto which this table should 
            be keyed. Defaults to None.
            fk_keys (list, optional): List of columns

        Raises:
            OperationalError: raised if the upload to the database is unsuccessful
            Keyerror: if a passed table_type does not correspond to a field in 
            self.required_fields
            AttributeError: if the fields in the passed data do not conform to 
            the expectation for that table_format in self.required_fields

        Returns:
            bool: True if the table upload is successful
        """
        if not table_format:
            for field_format, types in self.table_types_dict.items():
                if table_type in types:
                    table_format = field_format

            if table_format not in self.table_types_dict:
                expected_tbl_list = [x for sublist in self.table_types_dict.values()
                                    for x in sublist]
                raise KeyError(f"table_type {table_type} not recognized. "
                            f"Available options are {','.join(expected_tbl_list)}")
        # verify that the expected columns exist
        if set(self.required_fields[table_format]) - set(df.columns) != set():
            raise AttributeError(f"The format for this type_type is "
                                 f"{table_format}. The columns of the dataframe must therefore "
                                 f"be a subset of {self.required_fields[table_format]}")

        if not tablename:
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
            except sqlite.OperationalError as exc:  # pylint: disable=E1101
                msg = f"Could not drop table with sql: {drop_sql}. Error: {exc}"
                logging.critical(msg+" ", exc_info=(sys.exc_info()))  # pylint: disable=W1201
                raise

        col_dict = {k: "" for k in list(df.columns)}

        # create table - do this to corretly set primary key, foreign keys and 
        # data types, etc.
        if tablename not in self.list_tables(self.con):
            self.new_table(tablename, col_dict, fk_tablelist)
            # create index
            if table_format != 'chr_map':
                self.index_table(
                    tablename,
                    self.index_col_string_dict[table_format])
        # add the data
        df.to_sql(tablename,
                  con=self.con,
                  if_exists='append',
                  index=False)
        # standardize the chromosomes if there is a column called chr
        if self.chr_map_table in self.list_tables(self.con) and 'chr' in df.columns:
            print("Standardizing chr names...")
            self.standardize_chr_format(tablename)
        
        return True

    def standardize_chr_format(self, table: str) -> None:  # pylint: disable=E1101
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
        chr_map = {k: v for k, v in chr_map.items() if k != 'id'}

        cur = self.con.cursor()

        if 'chr' in self.list_fields(self.con, table):
            res = cur.execute(unique_chrnames_sql % table).fetchall()
            current_unique_chr_names = {str(x[0]) for x in res}
            # instantiate sentinel
            curr_chrom_format = -1
            # loop over colnames in chr_map_df. HALT if the current naming
            # convention is discovered
            for format, map_chr_names in chr_map.items():
                # check if all levels in the current chr_map_df[naming_convention]
                # contain the `df[chrom_colname]`.unique() levels
                if current_unique_chr_names - set(map_chr_names) == set():
                    curr_chrom_format = format
                    break
            # if the current chromosome naming convention is not intuited above,
            # raise an error
            if curr_chrom_format == -1:
                raise AttributeError("Chromosome names are not " +
                                     "recognized in a recognized format. Unique " +
                                     "chr names which cause error are: %s."
                                     % ",".join(current_unique_chr_names))

            # if the current names are already in the requested format, return
            if not curr_chrom_format == self.standard_chr_format:
                sql = swap_chr_format_sql % (table,
                                             self.standard_chr_format, curr_chrom_format,
                                             self.standard_chr_format, table, curr_chrom_format)
                # execute swap
                self.db_execute(cur, sql)
                self.con.commit() 