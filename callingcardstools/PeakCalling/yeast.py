# pylint:disable=W0640,W0108,C0114,W1203,W0127
from typing import Callable
import logging

import numpy as np
import sqlalchemy
import pandas as pd
import scipy.stats as scistat

from callingcardstools.Resources.PackageResources import Resources

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ['poisson_pval_factory',
           'hypergeometric_pval_factory',
           'evaluate_hops_constructor',
           'call_peaks']


def poisson_pval_factory(
        total_bg_hops: int,
        total_expr_hops: int,
        pseudocount: float) -> Callable[[int, int, dict], float]:
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
    # bg_hops / total_bg_hops -> fraction of total hops within region
    # in background (control) when multiplying by total_expr_hops,
    # gives E[expr_hops] in 'background model' goal: compare actual expr
    # hops to poisson with mean = E[expr_hops]
    # (expr_hops/total_expr_hops) / (bg_hops/total_bg_hops) -- this gives ratio
    # of fraction of reads (fold change) in expr to background

    def pval(bg_hops, expr_hops):
        """With the total_bg_hops and total_expr_hops set by the factory
         function wrapper, this return function acts as a pvalue calculator.

        Args:
            bg_hops (int): True to scale the random variable by the hop_ratio
             expr_hops (int): True to scale the random variable by
             the hop_ratio

        Returns:
            dict: keys pval and fold_change
        """
        # usage: scistat.poisson.cdf(x,mu)
        # where X are the experimental hops,
        # and mu is:
        # background hops * (total_expr_hops/total_background_hops)
        mu = (bg_hops * hop_ratio)+pseudocount
        x = expr_hops + pseudocount

        return 1-scistat.poisson.cdf(x, mu)

    return pval


def hypergeometric_pval_factory(
        total_bg_hops: int,
        total_expr_hops: int) -> Callable[[int, int], float]:
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
        return 1-scistat.hypergeom.cdf(x, M, n, N)
    return pval


def evaluate_hops_constructor(
        qbed_df: pd.DataFrame,
        total_bg_hops: int,
        total_expr_hops: int,
        pseudocount: float) -> Callable[[pd.Series], float]:
    """A factory function that returns a function that calculates 
    expression hops, log2 fold change over background, and poisson and 
    hypergeometric pvalues for a given range (chr, start, end) 

    Args:
        qbed_df (pd.DataFrame): _description_
        total_bg_hops (int): _description_
        total_expr_hops (int): _description_
        pseudocount (float): _description_

        Returns:
            Callable[[pd.Series], float]: _description_
    """
    qbed_df = qbed_df
    total_bg_hops = total_bg_hops
    total_expr_hops = total_expr_hops
    pseudocount = pseudocount

    def calculate_pvalues(chr: str, start: int,
                          end: int, bg_hops: int) -> np.array:
        """with the constants set by the factory function wrapper, this 
        calculates expression hops, log2 fold change over background, and 
        poisson and hypergeometric pvalues for a given range (chr, start, 
        end) presumably derived from a regions of interest table

        Args:
            chr (str): _description_
            start (int): _description_
            end (int): _description_
            bg_hops (int): _description_

        Returns:
            np.array: order of items: 
            [expr_hops, log2fc, poisson_pval,hypergeom_pval]
        """

        expr_hops = len(qbed_df[(qbed_df['chr'] == chr) &
                                (qbed_df['start'] >= start) &
                                (qbed_df['end'] <= end)].index)

        poisson_pval = poisson_pval_factory(
            total_bg_hops,
            total_expr_hops,
            pseudocount)

        hypergeom_pval = hypergeometric_pval_factory(
            total_bg_hops,
            total_expr_hops)

        # add psuedocount to avoid divide by zero errors
        log2fc = np.log2(((expr_hops/total_expr_hops)+pseudocount) /
                         ((bg_hops/total_bg_hops)+pseudocount))

        poisson = poisson_pval(bg_hops, expr_hops)

        hypergeom = hypergeom_pval(bg_hops, expr_hops)

        results_arr = np.array([expr_hops, log2fc, poisson, hypergeom])

        return results_arr

    return calculate_pvalues


def call_peaks(qbed_df: pd.DataFrame,
               regions_sample: str,
               background_sample: str,
               poisson_pseudocount: float) -> pd.DataFrame:
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
    if regions_sample not in {'yiming', 'not_orf'}:
        raise KeyError('regions_sample must be one of "yiming" or "not_orf"')
    if background_sample not in {'adh1', 'dsir4'}:
        raise KeyError('background_sample must be one of "adh1" or "dsir4"')

    total_hops_dict = {'total_expr_hops': len(qbed_df.index)}
    with Resources().yeast_resources.get('yeast_db') as yeast_db:
        engine = sqlalchemy.create_engine('sqlite:///{}'.format(yeast_db))

        # create total hop dictionary
        with engine.connect() as con:
            res = con.execute(f"select * from total_bg_hops "
                              f"where sample == '{background_sample}'")
            total_hops_dict.setdefault('total_bg_hops', res.all()[0][1])

            regions_df = pd.read_sql(f"select * from regions "
                                     f"where sample = '{regions_sample}'", con)
            background_agg = pd.read_sql(f"select * from region_background_agg "  # noqa
                                         f"where regions_sample = "
                                         f"'{regions_sample}' AND "
                                         f"background_sample = "
                                         f"'{background_sample}'", con)

        regions_full = regions_df.merge(
            background_agg, on=['chr', 'start', 'end'], how='left')
        regions_full.background_hops = regions_full.background_hops.fillna(0)
        regions_full.regions_sample = [regions_sample]*len(regions_full.index)
        regions_full.background_sample = \
            [background_sample]*len(regions_full.index)

        evaluate_hops = np.vectorize(evaluate_hops_constructor(
            qbed_df=qbed_df,
            **total_hops_dict,
            pseudocount=poisson_pseudocount), otypes=[np.double])

        effect_results = np.empty(len(regions_full.index), dtype=np.ndarray)
        for index, row in regions_full.iterrows():
            effect_results[index] = evaluate_hops(
                row.chr,
                row.start,
                row.end,
                row.background_hops)

        logging.debug(f"len effect_results: {len(effect_results)}, "
                      f"len regions_full: {len(regions_full.index)}")

        regions_full[['expr_hops', 'log2fc',
                      'poisson_pval', 'hypergeometric_pval']] = \
            pd.DataFrame.from_records(effect_results)

        return regions_full.loc[:, ['chr', 'start', 'end', 'background_hops',
                                    'expr_hops', 'log2fc', 'poisson_pval',
                                    'hypergeometric_pval']]\
            .rename({'background_hops': 'bg_hops'}, axis=1)
