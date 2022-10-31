from typing import Callable

import pandas as pd
import scipy.stats as scistat

def poisson_pval_factory(total_bg_hops:int, total_expr_hops:int, pseudocount:float) -> Callable[[int,int,dict], float]:
    """_summary_

    Args:
        total_bg_hops (int): _description_
        total_expr_hops (int): _description_
        pseudocount (float): _description_

    Returns:
        Callable[[int,int,dict], float]: _description_
    """

    # set constants
    total_bg_hops = float(total_bg_hops)
    total_expr_hops = float(total_expr_hops)
    hop_ratio = total_expr_hops / total_bg_hops
    def pval(bg_hops,expr_hops):
        """With the total_bg_hops and total_expr_hops set by the factory 
        function wrapper, this return function acts as a pvalue calculator.

        Args:
            bg_hops (int): True to scale the random variable by the hop_ratio 
            expr_hops (int): True to scale the random variable by the hop_ratio

        Returns:
            float: A pvalue for the region
        """
        # usage: scistat.poisson.cdf(x,mu)
        # where X are the experimental hops,
        # and mu is the background hops * (total_expr_hops/total_background_hops)
        mu = (bg_hops * hop_ratio)+pseudocount
        x  = expr_hops + pseudocount
        return 1-scistat.poisson.cdf(x, mu)
    return pval


def add_pvalues(poisson_pval:Callable[[int,int],float], 
                hypergeom_pval:Callable[[int,int],float]) -> Callable[[pd.Series],float]:
    """Calculate any number of pvalues for a row of a given dataframe

    Args:
        poisson_pval (Callable[[int,int],float]): a poisson pvalue calculator function
        hypergeom_pval (Callable[[int,int],float]): a hypergeometric pvalue calculator function

    Returns:
        Callable[[pd.Series],float]: A function which returns any number of pvalues
    """
    def calculate_pvalues(row:pd.Series):
        poisson = poisson_pval(row['bg_hops'],row['expr_hops'])
        hypergeom = hypergeom_pval(row['bg_hops'], row['expr_hops'])
        
        return poisson,hypergeom
    return calculate_pvalues

def hypergeometric_pval_factory(total_bg_hops:int, total_expr_hops:int) -> Callable[[int,int], float]:
    """_summary_

    Args:
        total_bg_hops (int): _description_
        total_expr_hops (int): _description_

    Returns:
        Callable[[int,int], float]: _description_
    """
    # usage:
    # scistat.hypergeom.cdf(x,M,n,N)
    # where x is observed number of type I events
    # (white balls in draw) (experiment hops at locus)
    # M is total number of balls (total number of hops)
    # n is total number of white balls (total number of expeirment hops)
    # N is the number of balls drawn (total hops at a locus)
    M = total_bg_hops + total_expr_hops
    n = total_expr_hops
    
    def pval(bg_hops, expr_hops):
        x = expr_hops - 1
        N = bg_hops + expr_hops
        return 1-scistat.hypergeom.cdf(x,M,n,N)
    return pval

def call_peaks_with_background(grouped_df: pd.DataFrame, 
                               total_hops_dict:dict,
                               poisson_pseudocount:float,
                               group_field:str="batch_id") -> pd.DataFrame:
    """_summary_

    Args:
        grouped_df (pd.DataFrame): _description_
        total_hops_dict (dict): _description_
        poisson_pseudocount (float): _description_
        group_field (str, optional): _description_. Defaults to "batch_id".

    Raises:
        KeyError: _description_
        AttributeError: _description_

    Returns:
        pd.DataFrame: _description_
    """
    
    
    try:
        if not set(grouped_df.groups) - set(total_hops_dict) == set():
            raise KeyError("All group names must be keys in the total_hops_dict")
    except AttributeError as err:
        raise AttributeError("The dataframe must be grouped, even if only a group of 1") from err
    
    output_df_list = []
    for group,df in grouped_df:
        # currently, don't collapse replicates
        # instantiate pvalue functions
        poisson_pval = poisson_pval_factory(
            total_hops_dict['background'],
            total_hops_dict[group],
            poisson_pseudocount)

        hypergeom_pval = hypergeometric_pval_factory(
            total_hops_dict['background'],
            total_hops_dict[group])

        pvalue_appender = add_pvalues(poisson_pval,hypergeom_pval)
        
        df[['poisson_pval', 'hypergeom_pval']] = \
            df.apply(lambda row: pvalue_appender(row), 
                        axis = 1, result_type = "expand")
        
        df[group_field] = group
        
        output_df_list.append(df)
    
    output_df = pd.concat(output_df_list,axis=0)

    return output_df