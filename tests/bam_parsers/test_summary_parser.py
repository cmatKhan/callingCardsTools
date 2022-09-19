import os
from .conftests import *
import pandas as pd
from pandas.testing import assert_frame_equal

def test_summary_parser_constructor(human_sp):
    assert isinstance(human_sp.summary, pd.DataFrame)
    assert human_sp.summary.shape[0] == 6
    assert human_sp.summary.shape[1] > 12

def test_to_qBed(human_sp):
    expected = {'chr': pd.Series(['chr1', 'chr1', 'chr1'], dtype = str), 
                'insert_start': pd.Series(['30191', '32303', '32303'], dtype = str),
                'insert_stop':  pd.Series([30195, 32307, 32307], dtype = str),
                'depth' : pd.Series([2,1,1], dtype = int),
                'strand': pd.Series(['+', '+', '-'], dtype = str),
                'annotation': pd.Series(['TTAA', 'TTAA', 'TTAA'], dtype = str)}

    expected_df = pd.DataFrame(expected)

    assert_frame_equal(human_sp.to_qBed('insert_seq'), expected_df, check_dtype = False)

def test_write_split(human_sp, tmpdir):
    human_sp.set_filter("mapq < 1000")
    path = os.path.join(tmpdir, "split")
    human_sp.write_split(path, 'insert_seq')
    assert os.listdir(tmpdir) == ['split_TTAA.qbed']