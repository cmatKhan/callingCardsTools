import numpy as np
import pytest

from callingcardstools.Analysis.yeast import combine_data, read_in_data

# Sample test data for read_in_data
sample_data_1 = """gene_id,effect,pvalue
gene1,0.5,0.01
gene2,0.7,0.05
gene3,0.3,0.10
"""

sample_data_2 = """gene_id,effect,pvalue
gene1,0.4,0.02
gene2,0.8,0.04
gene3,0.2,0.09
"""


@pytest.fixture
def create_temp_csv(tmpdir):
    """
    Fixture to create temporary CSV files.
    """

    def _create_temp_csv(content: str, file_name: str):
        temp_file = tmpdir.join(file_name)
        temp_file.write(content)
        return str(temp_file)

    return _create_temp_csv


# Test read_in_data
def test_read_in_data(create_temp_csv):
    # Create a temp file with sample data
    data_path = create_temp_csv(sample_data_1, "test_data.csv")

    # Run the read_in_data function
    df = read_in_data(
        data_path=data_path,
        identifier_col="gene_id",
        effect_col="effect",
        pval_col="pvalue",
        source="test_source",
        data_type="binding",
    )

    # Check that the dataframe has the expected columns
    assert list(df.columns) == [
        "feature",
        "binding_effect",
        "binding_pvalue",
        "binding_source",
    ]

    # Check some data integrity
    assert df.loc[df["feature"] == "gene1", "binding_effect"].values[0] == 0.5
    assert df.loc[df["feature"] == "gene2", "binding_pvalue"].values[0] == 0.05


def test_read_in_data_file_not_exist():
    with pytest.raises(FileExistsError):
        read_in_data(
            data_path="non_existent_file.csv",
            identifier_col="gene_id",
            effect_col="effect",
            pval_col="pvalue",
            source="test_source",
            data_type="binding",
        )


def test_read_in_data_missing_column(create_temp_csv):
    # Create a temp file without the identifier column
    data_path = create_temp_csv(sample_data_1, "test_data.csv")

    with pytest.raises(KeyError):
        read_in_data(
            data_path=data_path,
            identifier_col="nonexistent_column",
            effect_col="effect",
            pval_col="pvalue",
            source="test_source",
            data_type="binding",
        )


def test_read_in_data_na_values(create_temp_csv):
    # Create a temp file with NA values in the effect column
    data_with_na = """gene_id,effect,pvalue
    gene1,,0.01
    gene2,0.7,0.05
    """
    data_path = create_temp_csv(data_with_na, "test_data.csv")

    with pytest.raises(AttributeError):
        read_in_data(
            data_path=data_path,
            identifier_col="gene_id",
            effect_col="effect",
            pval_col="pvalue",
            source="test_source",
            data_type="binding",
        )


# Test combine_data
def test_combine_data(create_temp_csv):
    # Create two temp files with sample data
    data_path_1 = create_temp_csv(sample_data_1, "test_data_1.csv")
    data_path_2 = create_temp_csv(sample_data_2, "test_data_2.csv")

    combined_df = combine_data(
        data_paths=[data_path_1, data_path_2],
        identifier_col="gene_id",
        effect_col="effect",
        pval_col="pvalue",
        source="test_source",
        data_type="binding",
    )

    # Check that the combined dataframe has the expected shape and values
    assert len(combined_df) == 3  # 3 unique features
    assert (
        combined_df.loc[combined_df["feature"] == "gene1", "binding_effect"].values[0]
        == 0.45
    )
    assert combined_df.loc[combined_df["feature"] == "gene1", "binding_pvalue"].values[
        0
    ] == pytest.approx(0.014142, rel=1e-5)


def test_combine_data_with_zero_pvalues(create_temp_csv):
    # Create a temp file with zero p-values
    data_with_zero_pvals = """gene_id,effect,pvalue
    gene1,0.5,0.0
    gene2,0.7,0.0
    """
    data_path_1 = create_temp_csv(data_with_zero_pvals, "test_data_1.csv")

    combined_df = combine_data(
        data_paths=[data_path_1],
        identifier_col="gene_id",
        effect_col="effect",
        pval_col="pvalue",
        source="test_source",
        data_type="binding",
    )

    # make sure that no erroneous whitespace was added as a result of the
    # string conversion
    combined_df["feature"] = combined_df["feature"].str.strip()

    # Check that zeros are handled properly
    assert (
        combined_df.loc[combined_df["feature"] == "gene1", "binding_pvalue"].values[0]
        > 0
    )  # Should be replaced by a small value


def test_combine_data_custom_functions(create_temp_csv):
    # Create a temp file with sample data
    data_path_1 = create_temp_csv(sample_data_1, "test_data_1.csv")
    data_path_2 = create_temp_csv(sample_data_2, "test_data_2.csv")

    # Custom combine function for effect and p-values
    custom_effect_func = lambda x: np.median(x)
    custom_pval_func = lambda pvals: np.max(pvals)

    combined_df = combine_data(
        data_paths=[data_path_1, data_path_2],
        identifier_col="gene_id",
        effect_col="effect",
        pval_col="pvalue",
        source="test_source",
        data_type="binding",
        combine_effect_func=custom_effect_func,
        combine_pval_func=custom_pval_func,
    )

    # Check that the custom function was applied
    assert (
        combined_df.loc[combined_df["feature"] == "gene1", "binding_effect"].values[0]
        == 0.45
    )  # Median of 0.5 and 0.4
    assert (
        combined_df.loc[combined_df["feature"] == "gene1", "binding_pvalue"].values[0]
        == 0.02
    )  # Max of 0.01 and 0.02
