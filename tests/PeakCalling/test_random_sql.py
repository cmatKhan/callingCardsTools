
import os
import sqlite3
import pandas as pd
from pandas.testing import assert_frame_equal


def test_sig_hop_aggregate(tmp_path):

    create_hop_tbl_sql = {
        'drop': """DROP TABLE IF EXISTS qbed;""",
        'create': """CREATE TABLE "qbed" (
                            "chr"	TEXT,
                            "start"	INTEGER,
                            "end"	INTEGER,
                            "hops"	INTEGER,
                            "sig"	INTEGER
                        );""",
        'fill': """insert into qbed
                        values
                        ('chr1', 1,2, 10, 1),
                        ('chr1', 2,4, 10, 1),
                        ('chr1', 4,6, 5, 0),
                        ('chr2', 6,10, 10, 1);"""
    }

    db_path = os.path.join(tmp_path, "hopsdb.sqlite")

    con = sqlite3.connect(db_path)  # pylint: disable=E1101

    cur = con.cursor()

    for v in create_hop_tbl_sql.values():
        cur.execute(v)

    tablename = 'qbed'
    sig_hop_agg_sql = f"""SELECT DISTINCT(chr) as chr,
                         MIN(start) as start, 
                         MAX(end) as stop,
                         SUM(hops) as hops
                        FROM (
                            SELECT t1.*, SUM(group_flag) over (ORDER BY chr,start) as grp
                            FROM (
                                SELECT *,
                                CASE
                                    WHEN LAG(sig) OVER (ORDER BY start) = sig then NULL
                                    ELSE 1
                                END AS group_flag
                                FROM {tablename}
                            ) t1
                        ) t2
                        GROUP BY chr,grp
                        ORDER BY start"""

    query = cur.execute(sig_hop_agg_sql)
    cols = [column[0] for column in query.description]
    actual = pd.DataFrame.from_records(data=query.fetchall(), columns=cols)

    expected = pd.DataFrame(
        {'chr': ['chr1', 'chr1', 'chr2'], 'start': [1, 4, 6],
            'stop': [4, 6, 10], 'hops': [20, 5, 10]}
    )

    assert_frame_equal(actual, expected)


def test_rtree():
    con = sqlite3.connect(":memory:")  # pylint: disable=E1101
    # see https://docs.python.org/3/library/sqlite3.html#connection-objects
    # set row factory to sqlite3.Row to return iterable dictionary
    # objects per row where the keys are field names
    con.row_factory = lambda cursor, row: row[0]
    cur = con.cursor()

    rtree_create = """CREATE VIRTUAL TABLE demo_index USING rtree(
    id,              -- Integer primary key
    minX, maxX,      -- Minimum and maximum X coordinate
    minY, maxY       -- Minimum and maximum Y coordinate
    );"""

    rtree_fill = """INSERT INTO demo_index VALUES
  (28215, -80.781227, -80.604706, 35.208813, 35.297367),
  (28216, -80.957283, -80.840599, 35.235920, 35.367825),
  (28217, -80.960869, -80.869431, 35.133682, 35.208233),
  (28226, -80.878983, -80.778275, 35.060287, 35.154446),
  (28227, -80.745544, -80.555382, 35.130215, 35.236916),
  (28244, -80.844208, -80.841988, 35.223728, 35.225471),
  (28262, -80.809074, -80.682938, 35.276207, 35.377747),
  (28269, -80.851471, -80.735718, 35.272560, 35.407925),
  (28270, -80.794983, -80.728966, 35.059872, 35.161823),
  (28273, -80.994766, -80.875259, 35.074734, 35.172836),
  (28277, -80.876793, -80.767586, 35.001709, 35.101063),
  (28278, -81.058029, -80.956375, 35.044701, 35.223812),
  (28280, -80.844208, -80.841972, 35.225468, 35.227203),
  (28281, -80.844208, -80.841972, 35.0, 38.227203),
  (28282, -80.846382, -80.844193, 35.223972, 35.225655);"""

    rtree_query = ("SELECT id FROM demo_index "
                   "WHERE maxY>=35.0 AND minY<=35.0;")

    cur.execute(rtree_create)
    cur.execute(rtree_fill)
    x = cur.execute(rtree_query).fetchall()

    # not necessary here, just for demo and posterity purposes
    # reset row_factory
    con.row_factory = None

    assert x[0] == 28281
