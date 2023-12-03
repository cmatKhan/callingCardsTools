import logging
from pandas import Series
from numpy import int64
from scipy.stats import poisson

logger = logging.getLogger(__name__)


def poisson_pval_vectorized(total_background_hops: Series,
                            total_experiment_hops: Series,
                            background_hops: Series,
                            experiment_hops: Series,
                            pseudocount: float = 1e-10) -> Series:
    """
    Compute the Poisson p-value for the given hops counts.

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
    :param pseudocount: , defaults to 1e-10
    :type pseudocount: float, optional
    :return: a pandas Series of length equal to the input Series with the
        Poisson p-value for each row.
    :rtype: Series[float]

    .. note:: This function is vectorized, so it can be applied to
        pandas Series (columns of dataframes) to compute the
        Poisson p-value for each row.

    :raises ValueError: If any of the input Series contain negative values or
        the input Series are not all the same length.

    :Example:

    >>> import pandas as pd
    >>> total_background_hops = pd.Series([100, 200, 300])
    >>> total_experiment_hops = pd.Series([10, 20, 30])
    >>> background_hops = pd.Series([5, 10, 15])
    >>> experiment_hops = pd.Series([2, 4, 6])
    >>> vectorized_poisson_pval(
    ...     total_background_hops,
    ...     total_experiment_hops,
    ...     background_hops,
    ...     experiment_hops)
    array([0.01438768, 0.00365985, 0.00092599])
    """
    # check input
    if not len(total_background_hops) == len(total_experiment_hops) == \
            len(background_hops) == len(experiment_hops):
        raise ValueError('All input Series must be the same length.')
    if total_background_hops.min() < 0 \
            or total_background_hops.dtype != 'int64':
        raise ValueError(('total_background_hops must '
                          'be a non-negative integer.'))
    if total_experiment_hops.min() < 0 \
            or total_background_hops.dtype != 'int64':
        raise ValueError(('total_experiment_hops must '
                          'be a non-negative integer'))
    
    # cast to `float` b/c of scipy
    hop_ratio = (total_experiment_hops 
                 / (total_background_hops + pseudocount)).astype('float')
    mu = ((background_hops * hop_ratio)
          + pseudocount).astype('float')
    x = (experiment_hops + pseudocount).astype('float')

    return 1 - poisson.cdf(x, mu)
