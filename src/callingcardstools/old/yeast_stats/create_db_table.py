# standard lib
import sqlite3
import sys
import logging

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

def create_db_table(con, df, tbl_details_dict):
    """Create a table in the calling cards database

    :param con: sqlite3 connection object
    :type con: sqlite3 connection
    :param df: table to add to the database
    :type df: pandas DataFrame
    :param tbl_details_dict: a dictionary which has at least the attributes
    :type tbl_details_dict: dictionary
    """
    # assign a cursor to the database at cur
    cur = con.cursor()
    # drop the table if it exists
    drop_sql = "DROP TABLE IF EXISTS %s" %tbl_details_dict['tablename']
    try:
        cur.execute(drop_sql)
    except sqlite3.OperationalError as exc:
        msg = f"Could not drop table with sql: {drop_sql}. Error: {exc}"
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise
    # create the table
    create_tbl_sql = "CREATE TABLE %s (%s, PRIMARY KEY(%s AUTOINCREMENT));" \
        %(tbl_details_dict['tablename'], ",".join(["id INTEGER"] +
          list(df.columns)), "id")
    try:
        cur.execute(create_tbl_sql)
    except sqlite3.OperationalError as exc:
        cur.execute(drop_sql)
        msg = 'Could not create table with sql:\n"%s".\nTable has been dropped.\n \
            Error: %s.' %(create_tbl_sql, exc)
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise

    # create index on table
    index_sql = 'CREATE INDEX %s ON %s %s' %(tbl_details_dict['index_name'],
                                             tbl_details_dict['tablename'],
                                          tbl_details_dict['index_col_string'])
    try:
        cur.execute(index_sql)
    except sqlite3.OperationalError as e:
        cur.execute(drop_sql)
        msg = "could not create index with sql:\n%s.\n"%index_sql +\
        "Table has been dropped.\nError: %s:" %e
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise

    df.to_sql(tbl_details_dict['tablename'],
              con=con,
              if_exists='append',
              index=False)
