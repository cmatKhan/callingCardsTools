# standard library
import os
import logging 

# outside library
import pandas as pd

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ['SummaryParser']

class SummaryParser():

    _query_string = "status == 0"

    _summary_columns = {'id':str, 'status':int, 'mapq':int, 'flag':int, 'chr':str, 
                       'strand':str, 'five_prime':str, 'insert_start':str, 
                       'insert_stop':str, 'insert_seq':str, 'depth': int}

    _grouping_fields = {'chr', 'insert_start', 'insert_stop', 'strand'}

    _qbed_col_order = \
        ['chr', 'insert_start', 'insert_stop', 'depth', 'strand', 'annotation']
    
    _summary = None

    def __init__(self,summary_csv_path:str)->None:
        """_summary_

        Args:
            summary_csv_path (str): _description_
        """
        self.summary = summary_csv_path
    
        
        # switcher = {
        #     'status': self._status_filter(*args, **kwargs),
        #     'mapq': self._mapq_filter(*args, **kwargs),
        #     'region': self._region_filter(*args, **kwargs),
        #     'default': self.filterError(method)
        # }

        # switcher.get(method, "default")
    
    # def filterError(self, method):
    #     raise NotImplementedError(f"No filter method matches {method}")

    # def _status_filter(self, query_string):
    #     self.filter_string = query_string


    # def _mapq_filter(self):
    #     raise NotImplementedError
    
    # def _region_filter(self):
    #     raise NotImplementedError
	
    @property
    def query_string(self):
        """_summary_"""
        return self._query_string
    @query_string.setter
    def query_string(self, query_string:str):
        self._query_string = query_string

    @property
    def summary(self):
        """_summary_"""
        return self._summary
    @summary.setter
    def summary(self, summary_csv_path:str):
        # check genome and index paths
        for input_path in [summary_csv_path]:
            if not os.path.exists(input_path):
                raise FileNotFoundError(f"Input file DNE: {input_path}")
        summary = pd.read_csv(summary_csv_path, dtype = self.summary_columns)
        if 'depth' not in summary.columns:
            summary['depth'] = 1

        self._verify(summary)

        self._summary = summary
    
    @property
    def summary_columns(self):
        """_summary_"""
        return self._summary_columns
    @summary_columns.setter
    def summary_columns(self, col_list:list):
        self._summary_columns = col_list
    
    @property
    def grouping_fields(self):
        """_summary_"""
        return self._grouping_fields
    @grouping_fields.setter
    def grouping_fields(self, new_grouping_fields:dict):
        self.grouping_fields = new_grouping_fields
    
    @property
    def qbed_col_order(self):
        """_summary_"""
        return self._qbed_col_order
    @qbed_col_order.setter
    def qbed_col_order(self,new_col_order:list):
        self._qbed_col_order = new_col_order

    def _verify(self, summary:pd.DataFrame) -> None:
        """_summary_

        Args:
            summary (pd.DataFrame): _description_

        Raises:
            ValueError: _description_
        """
        if not len(set(self.summary_columns.keys()) - set(summary.columns)) == 0:
            raise ValueError(
                f"The expected summary columns are "\
                    f"{','.join(self.summary_columns)} in that order")

    def to_qbed(self, annotation:str) -> pd.DataFrame:
        """_summary_

        Args:
            annotation (str): _description_

        Raises:
            AttributeError: _description_

        Returns:
            pd.DataFrame: _description_
        """

        annote_list = ['tf', 'insert_seq', 'bc']
        if annotation not in annote_list:
            raise AttributeError(f"annotation must be one of {','.join(annote_list)}")

        local_grouping_fields = self.grouping_fields
        local_grouping_fields.add(annotation)

        return self.summary\
            .query(self.query_string)\
            [['chr', 'insert_start', 'insert_stop', 'depth', 'strand',annotation]]\
            .groupby(list(local_grouping_fields))['depth']\
        .agg(['sum'])\
        .reset_index()\
        .rename(columns={'sum':'depth', annotation: 'annotation'})[self.qbed_col_order]

    def write_split(self, prefix:str, annotation:str) -> None:
        """_summary_

        Args:
            prefix (str): _description_
            annotation (str): _description_
        """
        qbed = self.to_qbed(annotation).groupby('annotation')

        # output each grouped sheet with the name of the group in the filename
        for grouping,group_df in qbed:
            grouping = 'other' if grouping == '*' else grouping
            group_name = grouping if type(grouping) is str else "_".join(grouping)
            group_df.to_csv(prefix + "_" + group_name + ".qbed",
                            sep = "\t",
                            header = None,
                            index = False)
