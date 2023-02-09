# pylint:disable=E1123,C0121,E1124
import os
from io import StringIO
from sys import exc_info

import cProfile
import pysqlite3 as sqlite3
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal

from callingcardstools.Database.yeast import HopsDb as _yeast_hopsdb
from callingcardstools.Resources import Resources
from .conftests import *


def test_hopsdb_constructor(tmp_path):
    tmp_file = os.path.join(tmp_path, "test.sqlite")
    # test instantiating at various locations
    for loc in [":memory:", tmp_file]:
        hops_db = _yeast_hopsdb(loc)
        assert hops_db.is_open() == True
        hops_db.close()
        assert hops_db.is_open() == False
        assert hops_db.db_loc == loc


def test_yeast_add_chr_map():
    hops_db = _yeast_hopsdb(":memory:")

    r = Resources()

    chr_map_df = pd.read_csv(StringIO(r.yeast_chr_map))

    assert hops_db.add_frame(chr_map_df, "chr_map",
                             tablename="chr_map") == True


def test_yeast_add_regions(yeast_hops_data):
    hops_db = _yeast_hopsdb(":memory:")  # pylint:disable=E1102

    r = Resources()

    chr_map_df = pd.read_csv(StringIO(r.yeast_chr_map))

    assert hops_db.add_frame(
        chr_map_df,
        "chr_map",
        tablename="chr_map") is True

    regions_df = pd.read_csv(
        yeast_hops_data["regions"],
        sep="\t",
        names=hops_db.required_fields['bed6'])

    assert hops_db.add_frame(regions_df, 'bed6',
                             table_type="regions",
                             tablename_suffix='test') == True


def test_yeast_add_qbed_and_standardize_seqnames(yeast_hops_data):
    hops_db = _yeast_hopsdb(":memory:")

    r = Resources()

    chr_map_df = pd.read_csv(StringIO(r.yeast_chr_map))

    hops_db.add_frame(chr_map_df, "chr_map", tablename="chr_map")

    background_df = pd.read_csv(
        yeast_hops_data["background"],
        sep="\t",
        names=hops_db.required_fields['qbed'])

    assert hops_db.add_frame(background_df, 'qbed',
                             table_type='background',
                             tablename_suffix='test') is True

    expr_df = pd.read_csv(
        yeast_hops_data['experiment'],
        sep="\t",
        names=hops_db.required_fields['qbed'])

    assert hops_db.add_frame(expr_df, 'qbed',
                             table_type='experiment',
                             tablename_suffix='test') is True


def test_region_score(yeast_hopsdb):

    with pytest.raises(KeyError):
        yeast_hopsdb.peak_caller()

    actual = yeast_hopsdb.peak_caller(
        regions='regions_test',
        background='background_dSir4',
        experiment='experiment_test')

    save_db = False
    if save_db:
        db_disk = "temp/hopsdb.sqlite"
        if os.path.exists(db_disk):
            os.remove(db_disk)
        yeast_hopsdb.con.backup(sqlite3.Connection(db_disk))

    assert actual.shape == (6, 8)