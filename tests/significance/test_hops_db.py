import sqlite3
import os
from .conftests import *
import pandas as pd
from pandas.testing import assert_frame_equal

def test_sig_hop_aggregate(tmp_path, create_hop_tbl_sql):

    db_path = os.path.join(tmp_path, "hopsdb.sqlite")

    con = sqlite3.connect(db_path)

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
    actual = pd.DataFrame.from_records(data = query.fetchall(), columns = cols)

    expected = pd.DataFrame(
        {'chr': ['chr1', 'chr1', 'chr2'], 'start': [1,4,6], 'stop': [4,6,10], 'hops': [20,5,10]}
    )

    assert_frame_equal(actual, expected)
    
    