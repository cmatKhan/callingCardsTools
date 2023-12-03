import logging
from pandas import Series

logger = logging.getLogger(__name__)


def enrichment_vectorized(total_background_hops: Series,
                          total_experiment_hops: Series,
                          background_hops: Series,
                          experiment_hops: Series,
                          pseudocount: float = 1e-10) -> Series:
    """
    Compute the Calling Cards effect (enrichment) for the given hops counts.

    :param total_background_hops: a pandas Series (column of a dataframe)
        of total number of hops in the background.
    :type total_background_hops: Series
    :param total_experiment_hops: a pandas Series (column of a dataframe)
        of total number of hops in the experiment.
    :type total_experiment_hops: Series
    :param background_hops: a pandas Series (column of a dataframe)
        of number of hops in the background by promoter region.
    :type background_hops: Series
    :param experiment_hops: a pandas Series (column of a dataframe)
        of number of hops in the experiment by promoter region.
    :type experiment_hops: Series
    :param pseudocount: , defaults to 1e-10
    :type pseudocount: float, optional
    :return: a pandas Series of length equal to the input Series with the
        Calling Cards effect (enrichment) value for each row.
    :rtype: Series
    """
    numerator = (experiment_hops / (total_experiment_hops + pseudocount))
    denominator = (background_hops / (total_background_hops + pseudocount))

    return numerator / (denominator + pseudocount)
