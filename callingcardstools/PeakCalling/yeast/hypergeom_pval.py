import logging
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)

def hypergeom_pval(total_background_hops: int,
                   total_experiment_hops: int,
                   background_hops: int,
                   experiment_hops: int) -> float:
    """
    Compute the hypergeometric p-value for the given hops counts.

    :param total_background_hops: Total number of hops in the background.
    :type total_background_hops: int
    :param total_experiment_hops: Total number of hops in the experiment.
    :type total_experiment_hops: int
    :param background_hops: Number of hops in the background promoter region.
    :type background_hops: int
    :param experiment_hops: Number of hops in the experiment promoter region.
    :type experiment_hops: int
    :return: The hypergeometric p-value.
    :rtype: float
    """
    # check input
    if total_background_hops < 0 or not isinstance(total_background_hops, int):
        raise ValueError(('total_background_hops must '
                          'be a non-negative integer'))
    if total_experiment_hops < 0 or not isinstance(total_experiment_hops, int):
        raise ValueError(('total_experiment_hops must '
                          'be a non-negative integer'))
    if background_hops < 0 or not isinstance(background_hops, int):
        raise ValueError('background_hops must be a non-negative integer')
    if experiment_hops < 0 or not isinstance(experiment_hops, int):
        raise ValueError('experiment_hops must be a non-negative integer')

    # total number of objects (hops) in the bag (promoter region)
    M = total_background_hops + total_experiment_hops
    # if M is 0, the hypergeometric distribution is undefined
    # (there is no bag), so return 1
    if M < 1:
        return 1
    # number of 'success' objects in the population (experiment hops)
    n = total_experiment_hops
    # sample size (total hops drawn drawn from the bag)
    N = background_hops + experiment_hops
    # if N is 0, the hypergeometric distribution is undefined
    # (there is no sample), so return 1
    if N < 1:
        return 1
    # number of 'success' objects (experiment hops).
    # since we're interested in the chance of drawing a number of
    # experiment hops equal to or greater than the observed number,
    # we subtract 1 from the observed number in the CDF calculation.
    # Subtracting the result from 1 yields the right tailed p-value
    x = max((experiment_hops - 1), 0)

    pval = 1 - hypergeom.cdf(x, M, n, N)

    return pval
