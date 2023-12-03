import logging
from scipy.stats import poisson

logger = logging.getLogger(__name__)

def poisson_pval(total_background_hops: int,
                 total_experiment_hops: int,
                 background_hops: int,
                 experiment_hops: int,
                 pseudocount: float = 1e-10) -> float:
    """
    Compute the Poisson p-value for the given hops counts.

    :param total_background_hops: Total number of hops in the background.
    :type total_background_hops: int
    :param total_experiment_hops: Total number of hops in the experiment.
    :type total_experiment_hops: int
    :param background_hops: Number of hops in the background promoter region.
    :type background_hops: int
    :param experiment_hops: Number of hops in the experiment promoter region.
    :type experiment_hops: int
    :param pseudocount: A small constant to avoid division by zero.
        Default is 0.2.
    :type pseudocount: float
    :return: The Poisson p-value.
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

    hop_ratio = total_experiment_hops / (total_background_hops+pseudocount)
    # expected number of hops in the promoter region
    mu = (background_hops * hop_ratio) + pseudocount
    # random variable -- observed hops in the promoter region
    x = experiment_hops + pseudocount

    pval = 1 - poisson.cdf(x, mu)

    return pval
