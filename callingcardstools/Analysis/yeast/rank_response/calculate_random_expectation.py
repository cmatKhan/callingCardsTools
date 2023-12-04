import logging
import pandas as pd

logger = logging.getLogger(__name__)


def calculate_random_expectation(df):
    """
    Calculates the random expectation of responsiveness in a DataFrame.

    This function takes a DataFrame that contains a 'responsive' column with
    boolean values.  It calculates the proportion of rows marked as responsive
    and unresponsive, and then computes the expected random proportion of
    responsiveness.

    Args:
        df (pd.DataFrame): A DataFrame containing at least a 'responsive'
        column with boolean values.

    Returns:
        pd.DataFrame: A DataFrame with columns 'unresponsive', 'responsive',
        and 'random', where 'unresponsive' and 'responsive' are counts of each
        category, and 'random' is the proportion of responsive rows over the
        total number of rows.

    Example:
        >>> df = pd.DataFrame({'responsive': [True, False, True, False]})
        >>> calculate_random_expectation(df)
        # Returns a DataFrame with the counts of responsive and unresponsive
        # rows and the proportion
        # of responsive rows.
    """
    # Calculate counts for responsive and unresponsive
    counts = df['responsive'].value_counts()

    # Create the DataFrame
    random_expectation_df = pd.DataFrame({
        'unresponsive': [counts.get(False, 0)],
        'responsive': [counts.get(True, 0)]
    })

    # Calculate the 'random' column
    total = random_expectation_df.sum(axis=1)
    random_expectation_df['random'] = \
        random_expectation_df['responsive'] / total

    return random_expectation_df
