import logging

logger = logging.getLogger(__name__)


def enrichment(total_background_hops: int,
               total_experiment_hops: int,
               background_hops: int,
               experiment_hops: int,
               pseudocount: float = 1e-10):
    """
    Compute the Calling Cards effect (enrichment) for the given hops counts.

    :param total_background_hops: Total number of hops in the background.
    :type total_background_hops: int
    :param total_experiment_hops: Total number of hops in the experiment.
    :type total_experiment_hops: int
    :param background_hops: Number of hops in the background region.
    :type background_hops: int
    :param experiment_hops: Number of hops in the experiment region.
    :type experiment_hops: int
    :param pseudocount: A small constant to avoid division by zero.
        Default is 0.2.
    :type pseudocount: float
    :return: The Calling Cards effect (enrichment) value.
    :rtype: float
    """

    numerator = (experiment_hops / (total_experiment_hops + pseudocount))
    denominator = (background_hops / (total_background_hops + pseudocount))

    return (numerator / (denominator+pseudocount))
