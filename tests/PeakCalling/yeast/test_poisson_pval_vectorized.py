import pandas as pd
import pytest
from pandas import Series
from scipy.stats import poisson

from callingcardstools.PeakCalling.yeast.poisson_pval_vectorized import (
    poisson_pval_vectorized,
)


def test_poisson_pval_vectorized():
    # Create sample data
    total_background_hops = Series([103922] * 4)
    total_experiment_hops = Series([7000] * 4)
    background_hops = Series([1, 0, 5, 2])
    experiment_hops = Series([0, 10, 10, 1])
    # make sure this is the same as the function default
    pseudocount = 0.1

    # Expected result calculation
    hop_ratio = (total_experiment_hops / total_background_hops).astype("float")
    mu = ((background_hops + pseudocount) * hop_ratio).astype("float")
    x = experiment_hops.astype("float")

    expected_pval = pd.Series(1 - poisson.cdf(x, mu)) + poisson.pmf(x, mu)

    # Call the function
    result = poisson_pval_vectorized(
        total_background_hops,
        total_experiment_hops,
        background_hops,
        experiment_hops,
    )

    # Assert the result matches the expected result
    pd.testing.assert_series_equal(result, expected_pval)


def test_poisson_pval_vectorized_invalid_inputs():
    # Test for input length mismatch
    with pytest.raises(ValueError, match="All input Series must be the same length."):
        poisson_pval_vectorized(
            Series([100, 200]),
            Series([10, 20, 30]),
            Series([5, 10, 15]),
            Series([2, 4, 6]),
        )

    # Test for negative values
    with pytest.raises(
        ValueError, match="total_background_hops must be a non-negative integer."
    ):
        poisson_pval_vectorized(
            Series([-100, 200, 300]),
            Series([10, 20, 30]),
            Series([5, 10, 15]),
            Series([2, 4, 6]),
        )

    with pytest.raises(
        ValueError, match="total_experiment_hops must be a non-negative integer"
    ):
        poisson_pval_vectorized(
            Series([100, 200, 300]),
            Series([-10, 20, 30]),
            Series([5, 10, 15]),
            Series([2, 4, 6]),
        )


if __name__ == "__main__":
    pytest.main()
