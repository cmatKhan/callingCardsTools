import pandas as pd
import pytest
from pandas import Series
from scipy.stats import hypergeom

from callingcardstools.PeakCalling.yeast.hypergeom_pval_vectorized import (
    hypergeom_pval_vectorized,
)


def test_hypergeom_pval_vectorized():
    total_background_hops = Series([103922] * 4)
    total_experiment_hops = Series([7000] * 4)
    background_hops = Series([1, 0, 5, 2])
    experiment_hops = Series([0, 10, 10, 1])

    # Expected result calculation
    M = total_background_hops + total_experiment_hops
    n = total_experiment_hops
    N = background_hops + experiment_hops
    x = experiment_hops - 1

    valid = (M >= 1) & (N >= 1)
    expected_pval = Series(1, index=total_background_hops.index)
    expected_pval[valid] = 1 - hypergeom.cdf(x[valid], M[valid], n[valid], N[valid])

    # Call the function
    result = hypergeom_pval_vectorized(
        total_background_hops,
        total_experiment_hops,
        background_hops,
        experiment_hops,
    )

    # Assert the result matches the expected result
    pd.testing.assert_series_equal(result, expected_pval)


def test_hypergeom_pval_vectorized_invalid_inputs():
    # Test for input length mismatch
    with pytest.raises(ValueError, match="All input Series must be the same length."):
        hypergeom_pval_vectorized(
            Series([100, 200]),
            Series([10, 20, 30]),
            Series([5, 10, 15]),
            Series([2, 4, 6]),
        )

    # Test for negative values
    with pytest.raises(
        ValueError, match="total_background_hops must be a non-negative integer."
    ):
        hypergeom_pval_vectorized(
            Series([-100, 200, 300]),
            Series([10, 20, 30]),
            Series([5, 10, 15]),
            Series([2, 4, 6]),
        )

    with pytest.raises(
        ValueError, match="total_experiment_hops must be a non-negative integer"
    ):
        hypergeom_pval_vectorized(
            Series([100, 200, 300]),
            Series([-10, 20, 30]),
            Series([5, 10, 15]),
            Series([2, 4, 6]),
        )

    with pytest.raises(
        ValueError, match="background_hops must be a non-negative integer"
    ):
        hypergeom_pval_vectorized(
            Series([100, 200, 300]),
            Series([10, 20, 30]),
            Series([-5, 10, 15]),
            Series([2, 4, 6]),
        )

    with pytest.raises(
        ValueError, match="experiment_hops must be a non-negative integer"
    ):
        hypergeom_pval_vectorized(
            Series([100, 200, 300]),
            Series([10, 20, 30]),
            Series([5, 10, 15]),
            Series([-2, 4, 6]),
        )


if __name__ == "__main__":
    pytest.main()
