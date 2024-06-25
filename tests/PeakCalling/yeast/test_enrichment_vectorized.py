import numpy as np
import pandas as pd
from pandas import Series

from callingcardstools.PeakCalling.yeast.enrichment_vectorized import (
    enrichment_vectorized,
)


def test_enrichment_vectorized():
    # Create sample data
    total_background_hops = Series([103922] * 4)
    total_experiment_hops = Series([7000] * 4)
    background_hops = Series([1, 0, 5, 2])
    experiment_hops = Series([0, 10, 10, 1])
    # make sure this is the same as the function default
    pseudocount = 0.1

    # Expected result calculation
    numerator = experiment_hops / total_experiment_hops
    denominator = (background_hops + pseudocount) / total_background_hops
    expected_enrichment = numerator / denominator
    valid = np.logical_and(
        ~np.isinf(expected_enrichment), ~np.isnan(expected_enrichment)
    )
    expected_result = Series(-1, index=experiment_hops.index)
    expected_result[valid] = numerator[valid] / denominator[valid]

    # Call the function
    result = enrichment_vectorized(
        total_background_hops, total_experiment_hops, background_hops, experiment_hops
    )

    # Assert the result matches the expected result
    pd.testing.assert_series_equal(result, expected_result)
