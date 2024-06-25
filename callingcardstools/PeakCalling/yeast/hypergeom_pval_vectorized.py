import logging

from numpy.typing import NDArray
from pandas import Series
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)


def hypergeom_pval_vectorized(
    total_background_hops: Series,
    total_experiment_hops: Series,
    background_hops: Series,
    experiment_hops: Series,
) -> NDArray:
    """
    Compute the hypergeometric p-value for the given hops counts.

    :param total_background_hops: a pandas Series (column of a dataframe)
        of total number of hops in the background.
    :type total_background_hops: Series[int64]
    :param total_experiment_hops: a pandas Series (column of a dataframe)
        of total number of hops in the experiment.
    :type total_experiment_hops: Series[int64]
    :param background_hops: a pandas Series (column of a dataframe)
        of number of hops in the background by promoter region.
    :type background_hops: Series[int64]
    :param experiment_hops: a pandas Series (column of a dataframe)
        of number of hops in the experiment by promoter region.
    :type experiment_hops: Series[int64]
    :return: A pandas Series of length equal to the input Series with the
        hypergeometric p-value for each row. If either of the `total hop`
        input Series is 0, the hypergeometric p-value is undefined and
        the output Series will have a value of 1 for that row.
    :rtype: NDArray[float]

    .. note:: This function is vectorized, so it can be applied to
        pandas Series (columns of dataframes) to compute the
        hypergeometric p-value for each row.

    :raises ValueError: If any of the input Series contain negative values,
        are not dtype int64 or the input Series are not all the same length.

    :Example:

    >>> import pandas as pd
    >>> total_background_hops = pd.Series([100, 200, 300])
    >>> total_experiment_hops = pd.Series([10, 20, 30])
    >>> background_hops = pd.Series([5, 10, 15])
    >>> experiment_hops = pd.Series([2, 4, 6])
    >>> vectorized_hypergeom_pval(
    ...     total_background_hops,
    ...     total_experiment_hops,
    ...     background_hops,
    ...     experiment_hops)
    0    0.122360
    1    0.027644
    2    0.006972
    dtype: float64
    """
    # check input
    if (
        not len(total_background_hops)
        == len(total_experiment_hops)
        == len(background_hops)
        == len(experiment_hops)
    ):
        raise ValueError("All input Series must be the same length.")
    if total_background_hops.min() < 0 or total_background_hops.dtype != "int64":
        raise ValueError(("total_background_hops must " "be a non-negative integer."))
    if total_experiment_hops.min() < 0 or total_background_hops.dtype != "int64":
        raise ValueError(("total_experiment_hops must " "be a non-negative integer"))
    if background_hops.min() < 0 or background_hops.dtype != "int64":
        raise ValueError("background_hops must be a non-negative integer")
    if experiment_hops.min() < 0 or experiment_hops.dtype != "int64":
        raise ValueError("experiment_hops must be a non-negative integer")

    # calculate hypergeometric p-values
    M = total_background_hops + total_experiment_hops
    n = total_experiment_hops
    N = background_hops + experiment_hops
    x = experiment_hops - 1

    # Handling edge cases
    valid = (M >= 1) & (N >= 1)
    pval = Series(1, index=total_background_hops.index)
    pval[valid] = 1 - hypergeom.cdf(x[valid], M[valid], n[valid], N[valid])

    return pval
