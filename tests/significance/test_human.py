
# pylint:disable=E1123,C0121,E1124
import os

import cProfile
import pysqlite3 as sqlite3
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal

#from callingcardstools.database_managers.yeast import HopsDb as yeast_hopsdb
from .conftests import *


def test_add_ttaa(human_hops_data):

    save_db = True

    hops_data = human_hops_data

    hops_db = DatabaseApi(":memory:")

    chr_map_df = pd.read_csv(hops_data['chr_map'])

    hops_db.add_frame(chr_map_df, "chr_map")

    qbed_colnames = ['chr', 'start', 'end',
                     'depth', 'strand', 'annotation', 'sample']

    ttaa_df = pd.read_csv(
        hops_data['ttaa'],
        sep="\t",
        names=qbed_colnames)

    ttaa_df['sample'] = ['hg38'] * len(ttaa_df)
    assert hops_db.add_frame(ttaa_df, 'ttaa', "hg38") == True

    expr_df = pd.read_csv(
        hops_data['experiment'],
        sep="\t",
        names=qbed_colnames)
    expr_df['sample'] = ['test1'] * len(expr_df)
    expr_df = expr_df[expr_df['chr'].isin(
        ['chr'+str(x) for x in range(1, 22)] + ["chrX", "chrY", "chrM"])]

    assert hops_db.add_frame(expr_df, 'experiment', 'test') == True

    if save_db:
        db_disk = "temp/human_hopsdb.sqlite"
        if os.path.exists(db_disk):
            os.remove(db_disk)
        hops_db.con.backup(sqlite3.Connection(db_disk))  # pylint:disable=E1101
    hops_db.close()


def test_region_score_human_bf():
    hops_db_path = 'temp/human_hopsdb.sqlite'
    hops_db = DatabaseApi(hops_db_path)
    if not os.path.exists("temp"):
        os.mkdir("temp")
    #memory_profile = open("temp/range_score_macslike_bf_memory_profile.txt")
    with cProfile.Profile() as pr:
        # @profile(stream = memory_profile)
        hops_db.range_score_macslike_bf('experiment_test', 'ttaa_hg38', 1e-4)
    pr.dump_stats("temp/range_score_macslike_bf_time_profile.txt")

    assert 2 == 2


def test_range_score_macslike(human_hopsdb):
    pass
