#pylint:disable=E1123,C0121,E1124
import os
from sys import exc_info

from memory_profiler import profile
import cProfile
import pysqlite3 as sqlite3
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal

from callingcardstools.database_managers.yeast import HopsDb as _yeast_hopsdb
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

def test_yeast_add_chr_map(yeast_hops_data):
    hops_db = _yeast_hopsdb(":memory:")

    chr_map_df = pd.read_csv(yeast_hops_data['chr_map'])

    assert hops_db.add_frame(chr_map_df, "chr_map", tablename="chr_map") == True

def test_yeast_add_regions(yeast_hops_data):
    hops_db = _yeast_hopsdb(":memory:") #pylint:disable=E1102

    chr_map_df = pd.read_csv(yeast_hops_data['chr_map'])

    assert hops_db.add_frame(chr_map_df, "chr_map", tablename="chr_map") == True

    regions_df = pd.read_csv(
        yeast_hops_data["regions"], 
        sep = "\t", 
        names = hops_db.required_fields['bed6'])

    assert hops_db.add_frame(regions_df, 'bed6', 
    table_type="regions", tablename_suffix='test') == True

def test_yeast_add_qbed_and_standardize_seqnames(yeast_hops_data):
    hops_db = _yeast_hopsdb(":memory:")

    chr_map_df = pd.read_csv(yeast_hops_data['chr_map'])

    hops_db.add_frame(chr_map_df, "chr_map", tablename="chr_map")
     
    background_df = pd.read_csv(
        yeast_hops_data["background"], 
        sep = "\t", 
        names = hops_db.required_fields['qbed'] )
    
    assert hops_db.add_frame(background_df, 'qbed', 
                             table_type='background', 
                             tablename_suffix = 'test') == True

    expr_df = pd.read_csv(
        yeast_hops_data['experiment'], 
        sep = "\t", 
        names = hops_db.required_fields['qbed'])

    assert hops_db.add_frame(expr_df, 'qbed', 
    table_type='experiment', tablename_suffix = 'test') == True


def test_aggregate_hops(yeast_hopsdb):

    save_db = False
    
    with pytest.raises(AttributeError):
        yeast_hopsdb.create_aggregate_view("regions_upstream_700")
    
    actual = yeast_hopsdb.create_aggregate_view('regions_test')

    if save_db:
        db_disk = "temp/hopsdb.sqlite"
        if os.path.exists(db_disk):
            os.remove(db_disk)
        yeast_hopsdb.con.backup(sqlite3.connect(db_disk))

    assert actual == True

def test_join_sql(yeast_hopsdb):
    actual_sql = yeast_hopsdb.regions_background_expr_sql(regions='regions_test', 
        background = 'background_dSir4', 
        experiment='experiment_test')

    expected_sql = \
        " ".join(["SELECT e.chr AS chr, e.start AS start, e.end AS end,",
                         "IFNULL(b.hops,0) AS bg_hops, e.hops AS expr_hops, e.sample as sample",
                  "FROM regions_test_experiment_test as e",
                  "LEFT JOIN regions_test_background_dSir4 as b",
                  "USING(chr,start,end,name)",
                  "ORDER BY sample,chr,start,end"])

    assert actual_sql == expected_sql


def test_region_score(yeast_hopsdb):
    yeast_hopsdb.create_aggregate_view("regions_test")

    with pytest.raises(KeyError):
        yeast_hopsdb.peak_caller()    

    actual = yeast_hopsdb.peak_caller(
        regions='regions_test', 
        background = 'background_dSir4', 
        experiment='experiment_test')

    save_db = False
    if save_db:
        db_disk = "temp/hopsdb.sqlite"
        if os.path.exists(db_disk):
            os.remove(db_disk)
        yeast_hopsdb.con.backup(sqlite3.Connection(db_disk))

    assert actual == True

def test_create_batch_and_qc_tables(yeast_hopsdb):

    qc_set = {'batch', 'qc_manual_review', 
    'qc_r1_to_r2_tf', 'qc_r2_to_r1_tf', 'qc_alignment', 'qc_hops'}

    assert len(qc_set-set(yeast_hopsdb.list_tables(yeast_hopsdb.con)))==0

    yeast_hopsdb.add_batch_qc(
        'run_6073', 
        '/mnt/scratch/calling_cards/sequence/run_6073/run_6073_barcode_details.json', 
        '/mnt/scratch/calling_cards/sequence/run_6073/cctools_split/id_bc_map.tsv')
    
    assert 2==2

def test_summarize_tf_to_r1_transposon(yeast_hopsdb):

    yeast_hopsdb.add_batch_qc(
        'run_5423_aro80', 
        'tests/test_data/yeast/run_5423_barcode_details.json', 
        'temp/yeast/run_5423_id_bc_map.tsv')
    
    assert 2==2