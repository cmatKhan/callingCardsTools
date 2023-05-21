# pylint:disable=W1203
# standard library
import os
import logging
from typing import Union

# outside library
import pandas as pd

logger = logging.getLogger(__name__)

__all__ = ['SummaryParser']


class SummaryParser():
    """
    Class to parse summary data with provided grouping and ordering parameters.
    Able to convert this data into qBED format, a variant of the 
    BED format.
    """

    _query_string = "status == 0"

    _summary_columns = {'id': str, 'status': int, 'mapq': int, 'flag': int, 'chr': str,
                        'strand': str, 'five_prime': str, 'insert_start': str,
                        'insert_stop': str, 'insert_seq': str, 'depth': int}

    _grouping_fields = {'chr', 'insert_start', 'insert_stop', 'strand'}

    _qbed_col_order = \
        ['chr', 'start', 'end', 'depth', 'strand']

    _summary = None

    def __init__(self, summary: Union[str, pd.DataFrame]) -> None:
        """
        Initialize SummaryParser with given summary data.

        Args:
            summary (Union[str, pd.DataFrame]): Either a path to a CSV 
                file or an existing pandas DataFrame.
        """
        self.summary = summary

    @property
    def query_string(self):
        """
        Query string for filtering summary data. Default is "status == 0".
        """
        return self._query_string

    @query_string.setter
    def query_string(self, query_string: str):
        self._query_string = query_string

    @property
    def summary(self):
        """
        The summary data in DataFrame format.
        """
        return self._summary

    @summary.setter
    def summary(self, summary: Union[str, pd.DataFrame]):
        # check input
        if isinstance(summary, str):
            # check genome and index paths
            if not os.path.exists(summary):
                raise FileNotFoundError(f"Input file DNE: {summary}")
            summary = pd.read_csv(summary, dtype=self.summary_columns)
        elif isinstance(summary, pd.DataFrame):
            logger.info(f'passed a dataframe to SummaryParser')
        else:
            raise IOError(f'{summary} is not a data type recognized ' +
                          'as a summary by SummaryParser')

        if 'depth' not in summary.columns:
            summary['depth'] = 1

        self._verify(summary)

        self._summary = summary

    @property
    def summary_columns(self):
        """
        The expected structure (column names and data types) of 
        the summary data.
        """
        return self._summary_columns

    @summary_columns.setter
    def summary_columns(self, col_list: list):
        self._summary_columns = col_list

    @property
    def grouping_fields(self):
        """
        The set of fields to be used for grouping data in summary.
        """
        return self._grouping_fields

    @grouping_fields.setter
    def grouping_fields(self, new_grouping_fields: dict):
        self.grouping_fields = new_grouping_fields

    @property
    def qbed_col_order(self):
        """
        Order of columns to be used when generating a DataFrame in qBED format.
        """
        return self._qbed_col_order

    @qbed_col_order.setter
    def qbed_col_order(self, new_col_order: list):
        self._qbed_col_order = new_col_order

    def _verify(self, summary: pd.DataFrame) -> None:
        """
        Verifies that the provided summary DataFrame matches the 
        expected structure.

        Args:
            summary (pd.DataFrame): Summary data as a DataFrame.

        Raises:
            ValueError: Raised when the structure of the summary data does 
                not match the expected structure.
        """
        if not len(set(self.summary_columns.keys()) - set(summary.columns)) == 0:
            raise ValueError(
                f"The expected summary columns are "
                f"{','.join(self.summary_columns)} in that order")

    def to_qbed(self) -> pd.DataFrame:
        """
        Converts the summary data into a DataFrame in qBED format. It uses 
        the query string to filter data, groups by the defined grouping fields, 
        and orders columns as defined in qbed_col_order.

        Returns:
            pd.DataFrame: A DataFrame in qBED format.
        """

        local_grouping_fields = self.grouping_fields

        return self.summary\
            .query(self.query_string)[['chr', 'insert_start', 'insert_stop', 'depth', 'strand']]\
            .groupby(list(local_grouping_fields))['depth']\
            .agg(['sum'])\
            .reset_index()\
            .rename(columns={'sum': 'depth', 'insert_start': 'start', 'insert_stop': 'end'})[self.qbed_col_order]

    def write_qbed(self, output_path: str) -> None:
        """
        Writes the qBED-formatted DataFrame to a text file at the given path.

        Args:
            output_path (str): The path to the file where the output 
                should be written.
        """
        if not output_path[-4:] in ['.tsv', 'txt']:
            logger.warning(
                f"output path {output_path} does not end with tsv or txt")
        self.to_qbed().to_csv(output_path,
                              sep="\t",
                              header=None,
                              index=False)
