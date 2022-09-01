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

def create_hop_view(con, tablename, view_name):
    """Create view where hops are aggregated by chrom, pos and strand

    :param con: sqlite3 connection to database
    :type con: sqlie3 connection object
    :param tablename: name of the table from which to create the view
    :type tablename: str
    :param view_name: name of the view
    :type view_name: str
    """

    cur = con.cursor()

    drop_view_sql = 'DROP VIEW IF EXISTS "main"."%s";' %view_name

    cur.execute(drop_view_sql)

    create_view_sql = 'CREATE VIEW %s AS '%view_name +\
               'SELECT chrom, chromStart,chromEnd, COUNT(reads) as hops '+\
               'FROM(SELECT p.chrom AS chrom, p.chromStart AS chromStart, '+\
                            'p.chromEnd AS chromEnd, p.strand AS strand, '+\
                            'x.chromStart AS insert_pos, x.reads AS reads '+\
	                 'FROM %s as x '%tablename +\
	                 'LEFT JOIN promoters as p '+\
	                 'WHERE (x.chrom = p.chrom AND '+\
                        'x.chromStart BETWEEN p.chromStart AND p.chromEnd)) as t '+\
                'GROUP BY chrom, chromStart, chromEnd '+\
                'ORDER BY chrom, chromStart ASC, chromEnd ASC;'

    try:
        cur.execute(create_view_sql)
    except sqlite3.OperationalError as exc:
        cur.execute(drop_view_sql)
        msg = "could not create view with sql:\n%s.\n"%create_view_sql +\
        "View has been dropped.\nError: %s:" %exc
        logging.critical(msg+" ", exc_info=(sys.exc_info()))
        raise
