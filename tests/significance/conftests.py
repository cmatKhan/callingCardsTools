from io import StringIO
import pandas as pd

from ..conftests import *
from callingcardstools.Resources import Resources
from callingcardstools.database_managers.yeast import HopsDb as yeast_HopsDb


@pytest.fixture
def yeast_hopsdb(yeast_hops_data):
    """construct a hopsdb for yeast testing

    Args:
        yeast_hops_data (_type_): _description_

    Returns:
        _type_: _description_
    """

    hops_db = yeast_HopsDb(":memory:")

    r = Resources()

    chr_map_df = pd.read_csv(StringIO(r.yeast_chr_map))

    hops_db.add_frame(chr_map_df, "chr_map", tablename="chr_map")

    regions_df = pd.read_csv(yeast_hops_data["regions"],
                             sep="\t",
                             names=hops_db.required_fields['bed6'])

    hops_db.add_frame(
        regions_df,
        'bed6',
        table_type="regions",
        tablename_suffix='test')

    background_df = pd.read_csv(yeast_hops_data["background"],
                                sep="\t",
                                names=hops_db.required_fields['qbed'])

    hops_db.add_frame(background_df, 'qbed',
                      table_type='background',
                      tablename_suffix='dSir4')

    expr_df = pd.read_csv(yeast_hops_data['experiment'],
                          sep="\t",
                          names=hops_db.required_fields['qbed'])

    hops_db.add_frame(
        expr_df, 'qbed',
        table_type='experiment', 
        tablename_suffix='test')

    hops_db.create_aggregate_view('background_dSir4','regions_test')
    hops_db.create_aggregate_view('experiment_test','regions_test')

    return hops_db
