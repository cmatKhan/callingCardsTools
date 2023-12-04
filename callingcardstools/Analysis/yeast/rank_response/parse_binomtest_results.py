import logging
from scipy.stats._result_classes import BinomTestResult

logger = logging.getLogger(__name__)


def parse_binomtest_results(binomtest_obj: BinomTestResult, **kwargs):
    """
    Parses the results of a binomtest into a tuple of floats.

    This function takes the results of a binomtest and returns a tuple of
    floats containing the response ratio, p-value, and confidence interval
    bounds.

    Args:
        binomtest_obj (scipy.stats.BinomTestResult): The results of a binomtest
            for a single rank bin.
        Additional keyword arguments: Additional keyword arguments are passed
            to the proportional_ci method of the binomtest object.

    Returns:
        tuple: A tuple of floats containing the response ratio, p-value, and
            confidence interval bounds.

    Example:
        >>> parse_binomtest_results(binomtest(1, 2, 0.5, alternative='greater')
        (0.5, 0.75, 0.2, 0.8)
    """
    return (binomtest_obj.statistic,
            binomtest_obj.pvalue,
            binomtest_obj.proportion_ci(
                confidence_level=kwargs.get('confidence_level', 0.95),
                method=kwargs.get('method', 'exact')).low,
            binomtest_obj.proportion_ci(
                confidence_level=kwargs.get('confidence_level', 0.95),
                method=kwargs.get('method', 'exact')).high)
