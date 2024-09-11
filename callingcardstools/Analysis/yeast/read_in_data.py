import logging
import os
from typing import Callable, List, Literal

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

__all__ = ["read_in_data", "combine_data"]


def read_in_data(
    data_path: str,
    identifier_col: str,
    effect_col: str,
    pval_col: str,
    source: str,
    data_type: Literal["binding", "expression"],
) -> pd.DataFrame:
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

    compressed = data_path.endswith(".gz")
    logger.debug("data compressed: %s", compressed)

    sep = (
        "\t"
        if any(data_path.endswith(ext) for ext in [".tsv", ".txt", ".tsv.gz", "txt.gz"])
        else ","
    )
    logger.debug("data separator: %s", sep)

    df = pd.read_csv(data_path, sep=sep, compression="gzip" if compressed else None)

    if identifier_col not in df.columns:
        raise KeyError(f"Column `{identifier_col}` is not in {data_path}")
    if "feature" in df.columns and identifier_col != "feature":
        raise KeyError(
            f"Column `feature` exists in the data, but is not the "
            f"`identifier_col` {identifier_col}. Please rename the "
            f"current `feature` column to avoid confusion."
        )

    try:
        effect_colname = data_type + "_effect"
        # Assuming df is your DataFrame and effect_col is a variable
        # indicating column name
        df[effect_colname] = df[effect_col] if effect_col else float("inf")

        # Check for NA values in the effect_colname
        if pd.isna(df[effect_colname]).any():
            raise AttributeError(
                f"NA values found in column {effect_colname}." " This must not be."
            )
    except KeyError as exc:
        raise KeyError(
            f"Column {effect_col} is not `none` and " "does not exist in {data_path}"
        ) from exc

    try:
        pval_colname = data_type + "_pvalue"
        df[pval_colname] = df[pval_col] if pval_col else 0.0

        # Check for NA values in the pval_colname
        if pd.isna(df[pval_colname]).any():
            raise AttributeError(
                f"NA values found in column {pval_colname}. This must not be."
            )

        # Check for 'inf' values in the pval column
        if np.isinf(df[pval_colname]).any():
            # Remove inf values and check the largest remaining value
            non_inf_pvals = df[pval_colname].replace([np.inf, -np.inf], np.nan).dropna()

            # if the largest remaining p-value is 0, and there exists at least 1
            # negative value, then the pvalue column is logged.
            # convert positive infinity to negative infinity (this is for the raw
            # chipexo data)
            if (
                not non_inf_pvals.empty
                and non_inf_pvals.max() <= 0
                and non_inf_pvals.min() < 0
            ):
                # Convert positive infinity (resulting from log(0)) to negative infinity
                logger.warning(
                    "The pvalue column is logged, but there exist "
                    "positive infinity values. This is assumed to be an error. "
                    "Converting positive infinity to negative infinity."
                )
                df[pval_colname].replace(np.inf, -np.inf, inplace=True)

    except KeyError as exc:
        raise KeyError(
            f"Column {pval_col} is not `none` and does not exist in {data_path}"
        ) from exc

    source_colname = data_type + "_source"
    df[source_colname] = source

    df.rename(columns={identifier_col: "feature"}, inplace=True)

    return df[["feature", effect_colname, pval_colname, source_colname]]


def combine_pvals_detect_logged(pvals: pd.Series) -> float:
    """
    Detects if p-values are logged. If they are, converts them to linear scale,
    combines them using geometric mean, and then converts back to log scale.

    Args:
        pvals (pd.Series): Series of p-values, potentially logged.

    Returns:
        float: Combined p-value, logged if the original data was logged.
    """
    # check that pvals is a pandas Series
    if not isinstance(pvals, pd.Series):
        raise TypeError("pvals must be a pandas Series")
    # if pvals are ints, convert to floats
    if pvals.dtype == int:
        pvals = pvals.astype(float)
    # Detect if p-values are logged (log10 scale typically results in negative values)
    if pvals.max() <= 0:
        # If the p-values are likely log-scaled (log10 or log_e)
        logger.warning(
            "Detected log-scaled p-values, converting to linear scale to combine. "
            "They will be converted back to log scale after combination."
        )
        linear_pvals = 10**pvals  # Assuming log10, convert to linear scale
        combined_linear_pval = np.exp(
            np.mean(
                np.log(np.where(linear_pvals == 0, np.finfo(float).eps, linear_pvals))
            )
        )
        # Convert back to log scale after combination
        return np.log10(combined_linear_pval)
    else:
        # If p-values are already in linear scale, use the standard geometric mean
        combined_pval = np.exp(
            np.mean(np.log(np.where(pvals == 0, np.finfo(float).eps, pvals)))
        )
        return combined_pval


def combine_data(
    data_paths: List[str],
    identifier_col: str,
    effect_col: str,
    pval_col: str,
    source: str,
    data_type: Literal["binding", "expression"],
    combine_effect_func: Callable[[pd.Series], float] = np.mean,
    combine_pval_func: Callable[[pd.Series], float] = combine_pvals_detect_logged,
) -> pd.DataFrame:
    """
    Read in multiple data files and combine the effect and pvalue columns
    using specified functions (defaults are additive mean for effect and
    log mean for pvalue).

    Args:
        data_paths (List[str]): List of data file paths
        identifier_col (str): Name of the feature identifier column
        effect_col (str): Name of the effect column
        pval_col (str): Name of the pvalue column
        source (str): Source of the data
        data_type (str): Type of data, either 'binding' or 'expression'
        combine_effect_func (Callable): Function to combine effect columns
        combine_pval_func (Callable): Function to combine pvalue columns

    Returns:
        pd.DataFrame: Combined dataframe with averaged effect and pvalue
    """
    logger.info("combining data for data type {data_type} from {data_paths}")
    all_dfs = []

    for data_path in data_paths:
        df = read_in_data(
            data_path=data_path,
            identifier_col=identifier_col,
            effect_col=effect_col,
            pval_col=pval_col,
            source=source,
            data_type=data_type,
        )
        all_dfs.append(df)

    combined_df = (
        pd.concat(all_dfs)
        .groupby(["feature", f"{data_type}_source"])
        .agg(
            {
                f"{data_type}_effect": combine_effect_func,
                f"{data_type}_pvalue": combine_pval_func,
            }
        )
        .reset_index()
    )

    return combined_df
