import logging
import pandas as pd
from .read_in_data import read_in_data
from .validate_config import validate_config
from .rank_response_ratio_summarize import rank_response_ratio_summarize

logger = logging.getLogger(__name__)


def create_rank_response_table(config_dict: dict) -> (pd.DataFrame,
                                                      pd.DataFrame,
                                                      pd.DataFrame):
    """
    Create a rank repsonse table from a dictionary which contains the
    configuration parameters. See docs at
    https://cmatkhan.github.io/callingCardsTools/file_format_specs/yeast_rank_response/ # noqa
    for details

    Args:
        config_dict (dict): A dictionary containing the configuration
            parameters

    Returns:
        tuple: A tuple containing three DataFrames
            (see rank_response_summarize):
               1. The input with `repsonsive` genes labeled and the random
                    expectation column added
               2. The random expectation DataFrame
               3. The Rank Response DataFrame

    Raises:
        KeyError: if the configuration dictionary is missing any of the
            required keys
        FileExistsError: if the data files do not exist
        AttributeError: if there are NA values in the effect or pvalue columns
        ValueError: if there are incomplete cases in the data
    """
    # validate the configuration key/value pairs
    args = validate_config(config_dict)
    # read i the binding data
    try:
        binding_data = read_in_data(
            args['binding_data_path'],
            args['binding_identifier_col'],
            args['binding_effect_col'],
            args['binding_pvalue_col'],
            args['binding_source'],
            'binding')
    except (KeyError, FileExistsError, AttributeError) as exc:
        logger.error("Error reading in binding data: %s", exc)
        raise

    # read in the expression data
    try:
        expression_data = read_in_data(
            args['expression_data_path'],
            args['expression_identifier_col'],
            args['expression_effect_col'],
            args['expression_pvalue_col'],
            args['expression_source'],
            'expression')
    except (KeyError, FileExistsError, AttributeError) as exc:
        logger.error("Error reading in expression data: %s", exc)
        raise

    df = expression_data.merge(binding_data[['binding_effect',
                                             'binding_pvalue',
                                             'binding_source',
                                             'feature']],
                               how='inner',
                               on='feature')
    # test that there no incomplete cases. raise an error if there are
    if df.isnull().values.any():
        raise ValueError("There are incomplete cases in the data")

    logger.info('There are %s genes in the data after merging '
                'the %s binding data and '
                ' %s expression data',
                str(df.shape[0]),
                args['binding_source'],
                args['expression_source'])

    try:
        # the first two items in the return tuple aren't passed out of
        # this function, hence _, _
        _, _, rank_response_df = rank_response_ratio_summarize(
            df,
            effect_expression_thres=args['expression_effect_thres'],
            p_expression_thres=args['expression_pvalue_thres'],
            normalize=args['normalize'],
            bin_size=args['rank_bin_size'])
    except (KeyError) as exc:
        logger.error("Error summarizing data: %s", exc)
        raise

    return rank_response_df
