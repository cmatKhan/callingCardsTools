# standard library
import sys
import os
import logging
# outside library
import pandas as pd

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

class SummaryParser:

    query_string = "status == 0"

    summary_columns = ['id', 'bc', 'status', 'mapq', 'flag', 'chr', 
                       'strand', 'five_prime', 'insert_start', 
                       'insert_stop', 'insert_seq', 'tf', 'depth']

    grouping_fields = {'chr', 'insert_start', 'insert_stop', 'strand'}

    qbed_col_order = ['chr', 'insert_start', 'insert_stop', 'depth', 'strand', 'annotation']

    def __init__(self,summary_csv_path):
        self.set_summary(summary_csv_path)
    
    def set_filter(self, query_string):
        self.query_string = query_string
        
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

    def set_summary(self, summary_csv_path):
        # check genome and index paths
        for input_path in [summary_csv_path]:
            if not os.path.exists(input_path):
                raise FileNotFoundError(f"Input file DNE: {input_path}")
        summary = pd.read_csv(summary_csv_path)
        if 'depth' not in summary.columns:
            summary['depth'] = 1

        self._check_summary(summary)

        self.summary = summary

    def _check_summary(self, summary):
        if summary.columns != self.summary_columns:
            raise ValueError(f"The expected summary columns are {','.join(self.summary_columns)} in that order")

    def set_summary_columns(self, col_list):
        self.summary_columns = col_list


    def set_grouping_fields(self, grouping_fields):
        self.grouping_fields = grouping_fields

    def to_qBed(self, annotation):
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

    def write_split(self, prefix, annotation):
        qbed = self.to_qBed(annotation).groupby('annotation')

        # output each grouped sheet with the name of the group in the filename
        for grouping,group_df in qbed:
            grouping = 'other' if grouping == '*' else grouping
            group_name = grouping if type(grouping) is str else "_".join(grouping)
            group_df.to_csv(prefix + "_" + group_name + ".qbed",
                            sep = "\t",
                            header = None,
                            index = False)
