import logging
import warnings
import os
from .set_none_str_to_none import set_none_str_to_none

logger = logging.getLogger(__name__)


def validate_config(config: dict) -> dict:
    """
    Validate the yeast rank_response input configuration file.

    Args:
        config (dict): the configuration dictionary.

    Returns:
        dict: the validated configuration dictionary.

    Raises:
        KeyError: if the configuration is invalid due to either a missing
            key or an invalid value.
        TypeError: if the configuration is invalid due to an invalid type.
        FileNotFoundError: if the configuration is invalid due to a missing
    """
    # set default values if they are not in the config file
    config.setdefault('rank_by_effect', False)
    config.setdefault('rank_bin_size', 5)
    config.setdefault('normalize', False)
    config.setdefault('output_file', 'rank_response.csv')
    config.setdefault('compress', False)

    # this is used to check if a column is set to 'none' and replace it
    # with None. currently set for the expression_{effect/pvalue}_{col/thres}

    try:
        if not os.path.exists(config['binding_data_path']):
            raise FileNotFoundError(f"Binding data file "
                                    f"{config['binding_data_path']} "
                                    "does not exist")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_data_path' in config") from exc

    try:
        config['binding_source'] = str(config['binding_source'])
    except KeyError as exc:
        raise KeyError("Missing key 'binding_source' in config") from exc

    try:
        if not isinstance(config['binding_identifier_col'], str):
            raise TypeError("binding_identifier_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_identifier_col' in config") \
            from exc

    try:
        if not isinstance(config['binding_effect_col'], str):
            raise TypeError("binding_effect_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_effect_col' in config") \
            from exc

    try:
        if not isinstance(config['binding_pvalue_col'], str):
            raise TypeError("binding_pvalue_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'binding_pvalue_col' in config") \
            from exc

    try:
        if not isinstance(config['rank_by_effect'], bool):
            raise TypeError("rank_by_effect must be a boolean")
    except KeyError as exc:
        raise KeyError("Missing key 'rank_by_effect' in config") \
            from exc

    try:
        if not os.path.exists(config['expression_data_path']):
            raise FileNotFoundError(f"Expression data file "
                                    f"{config['expression_data_path']} "
                                    "does not exist")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_data_path' in config") \
            from exc
    
    try:
        config['expression_source'] = str(config['expression_source'])
    except KeyError as exc:
        raise KeyError("Missing key 'expression_source' in config") from exc

    try:
        if not isinstance(config['expression_identifier_col'], str):
            raise TypeError("expression_identifier_col must be a string")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_identifier_col' in config") \
            from exc

    try:
        if not isinstance(config['expression_effect_col'],
                          (str, type(None))):
            raise TypeError("expression_effect_col must be a string or None")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_effect_col' in config") \
            from exc

    try:
        if not isinstance(config['expression_effect_thres'],
                          (int, float, type(None))):
            raise TypeError("expression_effect_thres must be numeric or None")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_effect_thres' in config") \
            from exc

    for key in ['expression_effect_col', 'expression_effect_thres']:
        config[key] = set_none_str_to_none(config[key])

    if (config['expression_effect_col'] is None
            and config['expression_effect_thres'] is not None)\
        or (config['expression_effect_col'] is not None
            and config['expression_effect_thres'] is None):
        raise KeyError("expression_effect_thres must be None if "
                       "expression_effect_col is None")

    try:
        if not isinstance(config['expression_pvalue_col'], (str, type(None))):
            raise TypeError("expression_pvalue_col must be a string or None")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_pvalue_col' in config") \
            from exc

    try:
        if not isinstance(config['expression_pvalue_thres'],
                          (int, float, type(None))):
            raise TypeError("expression_pvalue_thres must be numeric or None")
    except KeyError as exc:
        raise KeyError("Missing key 'expression_pvalue_thres' in config") \
            from exc

    for key in ['expression_pvalue_col', 'expression_pvalue_thres']:
        config[key] = set_none_str_to_none(config[key])

    if (config['expression_pvalue_col'] is None
            and config['expression_pvalue_thres'] is not None)\
        or (config['expression_pvalue_col'] is not None
            and config['expression_pvalue_thres'] is None):
        raise KeyError("expression_pvalue_thres must be None if "
                       "expression_pvalue_col is None")

    if config['expression_pvalue_col'] is None \
            and config['expression_effect_col'] is None:
        raise KeyError("expression_pvalue_col and expression_effect_col "
                       "cannot both be None")

    try:
        if not isinstance(config['rank_bin_size'], int):
            raise TypeError("rank_bin_size must be an integer")
    except KeyError as exc:
        raise KeyError("Missing key 'rank_bin_size' in config") \
            from exc

    try:
        if not isinstance(config['normalize'], bool):
            raise TypeError("normalize must be a boolean")
    except KeyError as exc:
        raise KeyError("Missing key 'normalize' in config") \
            from exc

    try:
        if not isinstance(config['output_file'], str):
            raise TypeError("output_file must be a string")
        if os.path.exists(config['output_file']):
            warnings.warn("The output file already exists. "
                          "It will be overwritten.")
        if not config['output_file'].endswith('.csv'):
            warnings.warn("We suggest that the output file end with .csv")
    except KeyError as exc:
        raise KeyError("Missing key 'output_file' in config") \
            from exc

    try:
        if not isinstance(config['compress'], bool):
            raise TypeError("compress must be a boolean")
    except KeyError as exc:
        raise KeyError("Missing key 'compress' in config") \
            from exc

    return config