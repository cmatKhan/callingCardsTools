import logging
from typing import Literal
import os
from math import inf
import pandas as pd

logger = logging.getLogger(__name__)


def read_in_data(
        data_path: str,
        identifier_col: str,
        effect_col: str,
        pval_col: str,
        source: str,
        data_type: Literal['binding', 'expression']) -> pd.DataFrame:
    """
    Read in data from a file and return a dataframe with
    the following columns: gene_id, {binding/expression}_effect,
    {binding/expression}_pvalue, source

    Args:
        data_path (str): path to the data file
        identifier_col (str): name of the feature identifier column in the data
        effect_col (str): name of the effect column in the data
        pval_col (str): name of the pvalue column in the data
        source (str): source of the data
        data_type (str): type of data, either 'binding' or 'expression'

    Returns:
        pd.DataFrame: dataframe with the following columns:
            feature, {binding/expression}_effect, {binding/expression}_pvalue,
            source

    Raises:
        FileExistsError: if data_path does not exist
        KeyError: if identifier_col, effect_col, or pval_col is not in the
            data, or if the `identifier_col` is something other than `feature`
            and the column `feature` also exists in the data
        AttributeError: if there are NA values in the effect or pvalue columns
    """
    if not os.path.exists(data_path):
        raise FileExistsError(f"{data_path} does not exist")

    compressed = data_path.endswith('.gz')
    logger.debug("data compressed: %s", compressed)

    sep = '\t' if any(data_path.endswith(ext)
                      for ext in ['.tsv', '.txt', '.tsv.gz', 'txt.gz']) \
        else ','
    logger.debug("data separator: %s", sep)

    df = pd.read_csv(
        data_path,
        sep=sep,
        compression='gzip' if compressed else None)

    if identifier_col not in df.columns:
        raise KeyError(f"Column `{identifier_col}` is not in {data_path}")
    if 'feature' in df.columns and identifier_col != 'feature':
        raise KeyError(f"Column `feature` exists in the data, but is not the "
                       f"`identifier_col` {identifier_col}. Please rename the "
                       f"current `feature` column to avoid confusion.")

    try:
        effect_colname = data_type + '_effect'
        # Assuming df is your DataFrame and effect_col is a variable
        # indicating column name
        df[effect_colname] = df[effect_col] if effect_col else float('inf')

        # Check for NA values in the effect_colname
        if pd.isna(df[effect_colname]).any():
            raise AttributeError(f"NA values found in column {effect_colname}."
                                 " This must not be.")
    except KeyError as exc:
        raise KeyError(f"Column {effect_col} is not `none` and "
                       "does not exist in {data_path}") from exc

    try:
        pval_colname = data_type + '_pvalue'
        df[pval_colname] = df[pval_col] if pval_col else 0.0

        # Check for NA values in the pval_colname
        if pd.isna(df[pval_colname]).any():
            raise AttributeError(f"NA values found in column {pval_colname}. "
                                 "This must not be.")
    except KeyError as exc:
        raise KeyError(f"Column {pval_col} is not `none` and "
                       "does not exist in {data_path}") from exc

    source_colname = data_type + '_source'
    df[source_colname] = source

    df.rename(columns={identifier_col: 'feature'}, inplace=True)

    return df[['feature', effect_colname, pval_colname, source_colname]]
