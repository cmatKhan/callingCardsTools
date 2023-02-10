from callingcardstools.PeakCalling.yeast import call_peaks

import pandas as pd

from ..conftests import *


def test_peakcalling(yeast_qbed):

    func_input = {
        'qbed_df': pd.read_csv(
            yeast_qbed,
            sep='\t',
            names=['chr', 'start', 'end', 'depth', 'strand', 'annotation']),
        'regions_sample': 'yiming',
        'background_sample': 'adh1',
        'poisson_pseudocount': 0.1}

    call_peaks(**func_input)
