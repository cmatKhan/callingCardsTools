import logging
import numpy as np

logger = logging.getLogger(__name__)


def create_partitions(vector_length, equal_parts=100):
    """
    Splits a vector of a specified length into nearly equal partitions.

    This function creates a partition vector where each partition is of equal
    size, except the last partition which may be smaller depending on the
    vector length and the number of equal parts specified. Each element in
    the partition vector represents the partition number.

    Args:
        vector_length (int): The total length of the vector to be partitioned.
        equal_parts (int, optional): The number of equal parts to divide the
            vector. Defaults to 100.

    Returns:
        numpy.ndarray: An array where each element represents the partition
            number for each element in the original vector.

    Examples:
        >>> create_partitions(10, 3)
        array([1, 1, 1, 2, 2, 2, 3, 3, 3, 3])
    """
    quotient, remainder = divmod(vector_length, equal_parts)
    return np.concatenate([np.repeat(np.arange(1, quotient + 1), equal_parts),
                           np.repeat(quotient + 1, remainder)])