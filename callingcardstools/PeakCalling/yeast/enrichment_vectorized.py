import logging

import numpy as np
from pandas import Series

logger = logging.getLogger(__name__)


def enrichment_vectorized(
    total_background_hops: Series,
    total_experiment_hops: Series,
    background_hops: Series,
    experiment_hops: Series,
    pseudocount: float = 0.1,
    **kwargs
) -> Series:
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
    :param pseudocount: Added to the background hops to avoid division by zero,
    :type pseudocount: float, optional
    :param kwargs: Additional keyword arguments. None are currently used

    :return: a pandas Series of length equal to the input Series with the
        Calling Cards effect (enrichment) value for each row.
    :rtype: Series
    """
    # raise an error if any one of the 4 input Series is not a Series
    if not all(
        isinstance(x, Series)
        for x in [
            total_background_hops,
            total_experiment_hops,
            background_hops,
            experiment_hops,
        ]
    ):
        raise ValueError(
            "`total_background_hops`, `total_experiment_hops`, ",
            "`background_hops` and `experiment_hops` must all ",
            "be pandas Series. At least one is not.",
        )
    # validate that all input Series are the same length
    if (
        not len(total_background_hops)
        == len(total_experiment_hops)
        == len(background_hops)
        == len(experiment_hops)
    ):
        raise ValueError("All input Series must be the same length.")

    # validate that pseudocount is numeric (int or float). Cast to float if int
    if not isinstance(pseudocount, (int, float)):
        raise ValueError("pseudocount must be a number.")
    if isinstance(pseudocount, int):
        logger.warning("pseudocount is an integer. It will be cast to a float.")
        pseudocount = float(pseudocount)

    # NOTE: the total_experiment_hops and total_background_hops must be > 0 based on
    # input data verification. See read_in_experiment_data()
    # and read_in_background_data() in read_in_data.py
    numerator = experiment_hops / total_experiment_hops

    # Add a small pseudocount to background_hops to avoid division by zero in the
    # enrichment calculation below
    denominator = (background_hops + pseudocount) / total_background_hops

    enrichment = numerator / denominator

    # Check for invalid values
    if (enrichment < 0).any():
        raise ValueError("Enrichment values must be non-negative.")
    if enrichment.isnull().any():
        raise ValueError("Enrichment values must not be NaN.")
    if np.isinf(enrichment).any():
        raise ValueError("Enrichment values must not be infinite.")

    return enrichment
