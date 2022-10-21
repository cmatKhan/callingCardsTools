# stdlib
import re
from typing import Callable, Literal
import logging
# local
from ...DatabaseApi import DatabaseApi
from .peak_calling import call_peaks_with_background
# outside
import pandas as pd
import numpy as np

logging.getLogger(__name__).addHandler(logging.NullHandler())

class HopsDb(DatabaseApi):
    """A database manager for yeast calling cards experiments which includes background hops data"""
    
    def add_regions(self, regions_df: pd.DataFrame, *args, **kwargs) -> bool:
        """Add a regions table to the database. note that this function is a convenience function 
        which wraps DatabaseApi.add_table. It does some checking of the regions dataframe and then 
        calls DatabaseApi.add_table with the regions_df, "regions" set for the table_type, and 'bed3' 
        set for the table_format. You can pass additional arguments to DatabaseApi.add_table with either 
        positional or keyword arguments.

        Args:
            regions_df (pd.DataFrame): _description_
            *args (list, optional): additional arguments to DatabaseApi.add_table
            **kwargs (dict, optional): additional keyword arguments to DatabaseApi.add_table

        Raises:
            AttributeError: _description_
            ValueError: _description_

        Returns:
            bool: _description_
        """
        if len(np.setdiff1d(self.required_fields['bed3'], regions_df.columns)) > 0:
            raise AttributeError(f"regions df must have at least the headers {self.required_fields['bed3']}")
        elif len(regions_df) > len(regions_df.groupby(self.required_fields['bed3']).size()):
            raise ValueError(f"regions_df has duplicate entries in the combination of {self.required_fields['bed6']}. "\
                f"These regions should be unique as calling cards hop counting is unstranded. "\
                    f"Use an annotation table if you need to annotate the same region with multiple affiliations.")
        else:
            self.add_frame(regions_df,'regions','bed3', *args,**kwargs)
    
    def add_annotation_fk(self, df: pd.DataFrame, regions_fk_table:str) -> bool:
        """Add an annotation table. Optionally key to a regions table

        Args:
            df (pd.DataFrame): A dataframe with little restriction 
            to fields. However, the assumption is that you will want these 
            formatted in such a way that they can be used to annotate your 
            significant regions. Suggested minimum columns: chr,start,stop, which 
            can be used to join/overlap with the bed6 and qBed tables. Note that 
            adding a table with the same tablename_suffix will overwrite the current 
            table with the same suffix if one exists.
            tablename_suffix (str): suffix to be appended to annotation_ to form 
            the new tablename. Eg,
            regions_fk_table (str): Name of a regions table to which to key. The foreign key 
            will be <fk_tablename>_id

        Returns:
            bool: True if successful
        """
        # check that the fk table exists
        if not regions_fk_table in self.list_tables(self.con):
            raise AttributeError(f"{regions_fk_table} not in the database!")
        # create the fk_fieldname and ensure that it exists in the database
        fk_field = regions_fk_table+"_id"
        if not fk_field in df.columns:
            raise KeyError(f"{fk_field} is not a field in the input dataframe")
        # create the tablename
        tablename = 'annotation'+f'_{regions_fk_table}' 
        # add this table and the fk to the required_fields dict
        self.required_fields = {tablename: [fk_field]} #pylint:disable=W0201
        # add the frame
        self.add_frame(df,
                       table_type   = 'annotation',
                       table_format = 'annotation',
                       tablename    = tablename,
                       fk_table     = regions_fk_table,
                       fk_keys      = [fk_field])

    def create_aggregate_view(self, regions_tbl:str) -> bool:
        """For each background and experiment table, create a view where the hops are 
        aggregated over the regions provided in the regions table

        Args:
            regions_tbl (str): name of the regions table to use to aggregate the hops

        Returns:
            bool: True if successful
        """
        if not regions_tbl in self.list_tables(self.con):
            raise AttributeError(f'No table named {regions_tbl}')
        # get table lists
        tbl_list = self.list_tables(self.con)
        background_and_expr_tbls = [x for x in tbl_list \
            if re.search(r"^background|^experiment",x)]
        available_regions = [x for x in tbl_list if x.startswith("regions_")]
        # check if region_tbl is in the regions table list
        if regions_tbl not in self.list_tables(self.con):
            AttributeError(f"No table named {regions_tbl} in database. "\
                           f"Current tables are: {','.join(available_regions)}")
        # create the views
        cur = self.con.cursor()
        for tbl in background_and_expr_tbls:

            viewname = regions_tbl + '_' +  tbl
            drop_view_sql = f'DROP VIEW IF EXISTS "main"."{viewname}";'
            self.db_execute(cur,drop_view_sql)

            create_view_sql = " ".join([f"CREATE VIEW {viewname} AS",
                                         "SELECT r.chr as chr,r.start as start,r.end as end, count(*) as hops, sample, name",
                                        f"FROM {regions_tbl} as r",
                                        f"LEFT JOIN {tbl} as x",
                                        "WHERE (x.chr = r.chr AND x.start BETWEEN r.start AND r.end)",
                                        "GROUP BY sample, name, r.chr,r.start,r.end"])

            self.db_execute(cur,create_view_sql)
            self.con.commit()

        return True
    
    def regions_background_expr_sql(self, regions:str, background:str, experiment:str)->str:
        """_summary_

        Note that the keyword arguments must remain the same as the expected 
        keyword arguments in peak_caller

        Args:
            regions (str): _description_
            background (str): _description_
            experiment (str): _description_

        Raises:
            AttributeError: _description_

        Returns:
            str: _description_
        """
        # get the viewname
        hop_views = {'background': regions + '_' + background,
                     'experiment': regions + '_' + experiment}

        # check that the views are in the database
        tbl_list = self.list_tables(self.con)
        for t in hop_views.values():
            if t not in tbl_list:
                raise AttributeError(f"No view called {t} exists. "\
                    f"Try calling region_aggregate() and then resubmit")

        # finally, I decided that the areas where the experiment is null actually 
        # don't matter -- no reason to save them, anyway.
        sql = " ".join(["SELECT e.chr AS chr, e.start AS start, e.end AS end,",
                             "IFNULL(b.hops,0) AS bg_hops, e.hops AS expr_hops,",
                             "e.sample as sample",
                             f"FROM {hop_views['experiment']} as e",
                             f"LEFT JOIN {hop_views['background']} as b USING(chr,start,end,name)",
                             "ORDER BY sample,chr,start,end"])
        return sql
    

    def consolidate_tables(self,table_class:Literal['background','experiment']) -> bool:
        # take all tables which start with the prefix 'background_' and collapse them 
        # into a single table indexed by the sample column
        raise NotImplementedError

    def peak_caller(self,replicate_handling:Literal['separate','sum']='separate',poisson_pseudocount:float = 0.2, *args, **kwargs) -> None:
        """Call Peaks and add the result to the database.

        Args:
            replicate_handling (str, ['separate','sum'], optional): How to handle replicates. 'sum' 
            will add hops by region for all replicates. 'separate' will call peaks for each 
            replicate separately. Defaults to separate.
            *args (list): additional positional arguments -- currently unused
            poisson_pseudocount (float, optional): peudocount to add to poisson pvalue calculation. Defaults to 0.2.
            **kwargs (dict): Use the following combination of keyword arguments 
            to direct the peak_caller method:
                regions='regions_tbl',background='background_tbl',experiment='experiment_tbl' -- Call peaks with background. 

        Raises:
            AttributeError: _description_
        """
        if {'regions','background','experiment'} == set(kwargs):
            join_sql = self.regions_background_expr_sql(**kwargs)
            quant_df = pd.read_sql_query(join_sql, self.con)
            total_hops_dict = \
                {'background': self.get_total_hops(kwargs.get('background'))}
            if replicate_handling == 'separate':
                # group by sample (replicates)
                grouped_df = quant_df.groupby('sample')
                # add total hops for each replicate to the total_hop_dict
                for group,df in grouped_df:
                    total_hops_dict[group] = \
                        self.get_total_hops(kwargs.get('experiment'),group)
                # create the tablename
                sig_tablename = kwargs.get('regions') + '_' + kwargs.get('background')+ \
                                "_" + kwargs.get('experiment') + "_sig"
            else:
                NotImplementedError(f'Replicate handling method: '\
                    f'{replicate_handling} is not implemented')
            # call peaks
            output_df = call_peaks_with_background(grouped_df, total_hops_dict,poisson_pseudocount)
        else:
            raise KeyError(f'The combination of tables {kwargs} '\
                'does not match any expected combinations -- cannot find an '
                'appropriate peak calling method')
            
        # send the table to the database
        output_df.to_sql(sig_tablename,
            con=self.con,
            if_exists='replace',
            index=False)
        
        # index the table
        index_col_string = self.index_col_string_dict['qbed']+',"sample"'
        self.index_table(sig_tablename, index_col_string)
    
        return True