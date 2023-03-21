from callingcardstools.PeakCalling.yeast import *  #pylint:disable=W0164,W0401 # noqa

import pandas as pd
from scipy import stats as scistat

from ..conftests import *  #pylint:disable=W0401,W0614 # noqa


# def test_calculate_pvalues(yeast_qbed):

#     func_input = {
#         'qbed_df': pd.read_csv(
#             yeast_qbed,
#             sep='\t',
#             names=['chr', 'start', 'end', 'depth', 'strand', 'annotation']),
#         'regions_sample': 'yiming',
#         'background_sample': 'adh1',
#         'poisson_pseudocount': 0.1}

#     calculate_pvalues(**func_input)

def test_poisson_pval_factory():

    # set up parameters
    pseudocount = 0.2
    total_bg_hops = 103912
    total_expr_hops = 1181

    bg_hops = 58
    expr_hops = 1

    # create the poisson pval calc function
    pval_calculator = poisson_pval_factory( # noqa
        total_bg_hops,
        total_expr_hops,
        pseudocount)

    hop_ratio = total_expr_hops / total_bg_hops
    mu = (bg_hops * hop_ratio)+pseudocount
    x = expr_hops + pseudocount
    expected = 1-scistat.poisson.cdf(x, mu)

    actual = pval_calculator(bg_hops, expr_hops)

    assert expected == actual


def test_hypergeometric_pval_factory():

    # set up parameters
    total_bg_hops = 103912
    total_expr_hops = 1181

    bg_hops = 58
    expr_hops = 1

    # create the poisson pval calc function
    pval_calculator = hypergeometric_pval_factory( # noqa
        total_bg_hops,
        total_expr_hops)

    expected = 1-scistat.hypergeom.cdf(
        expr_hops-1,
        total_bg_hops + total_expr_hops,
        total_expr_hops,
        bg_hops + expr_hops)

    actual = pval_calculator(bg_hops, expr_hops)

    assert expected == actual


def test_peakcalling(yeast_qbed):

    func_input = {
        'qbed_df': pd.read_csv(
            yeast_qbed,
            sep='\t',
            names=['chr', 'start', 'end', 'depth', 'strand', 'annotation']),
        'regions_sample': 'yiming',
        'background_sample': 'adh1',
        'poisson_pseudocount': 0.1}

    call_peaks(**func_input) # noqa