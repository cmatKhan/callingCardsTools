import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from callingcardstools.Analysis.yeast import (
    bin_by_binding_rank,
    compute_rank_response,
    create_partitions,
    create_rank_response_table,
    label_responsive_genes,
    parse_binomtest_results,
    rank_response_main,
    read_in_data,
    validate_config,
)


def test_read_in_data_columns():
    test_data_directory = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "test_data/yeast/Analysis",
    )

    assert os.path.isdir(test_data_directory) is True

    test_data_path = os.path.join(
        test_data_directory, "hap5_exprid_17_yiming_adh1_promoter_sig.csv.gz"
    )

    test_identifier_col = "target_gene_id"
    test_effect_col = "callingcards_enrichment"  # Replace with actual column name
    test_pval_col = "poisson_pval"  # Replace with actual column name
    test_source = "callingcards"
    data_type = "binding"

    # Run the function under test
    df = read_in_data(
        test_data_path,
        test_identifier_col,
        test_effect_col,
        test_pval_col,
        test_source,
        data_type,
    )

    # Expected columns
    expected_columns = [
        "feature",
        f"{data_type}_effect",
        f"{data_type}_pvalue",
        f"{data_type}_source",
    ]

    # Check if all expected columns are present
    assert all(
        [col in df.columns for col in expected_columns]
    ), f"Dataframe columns do not match expected for type '{data_type}'"

    # Optionally, you can also check the types of columns or other properties


def test_create_partitions():
    # Test case with perfect division
    vector_length = 100
    equal_parts = 10
    expected = np.repeat(np.arange(1, 11), 10)
    result = create_partitions(vector_length, equal_parts)
    assert np.array_equal(result, expected), "Failed for perfect division"

    # Test case with a remainder
    vector_length = 103
    equal_parts = 10
    expected = np.concatenate([np.repeat(np.arange(1, 11), 10), np.repeat(11, 3)])
    result = create_partitions(vector_length, equal_parts)
    assert np.array_equal(result, expected), "Failed for division with remainder"

    # Test case to check the behavior with default equal_parts
    vector_length = 150
    # Expected behavior with default equal_parts
    expected = np.concatenate([np.repeat(1, 100), np.repeat(2, 50)])
    result = create_partitions(vector_length)
    assert np.array_equal(result, expected), "Failed for default equal_parts"


def test_label_responsive_genes():
    df = pd.DataFrame(
        {
            "expression_effect": [0.5, 0.7, 1.2, 0.2, 1.5, -1.5],
            "expression_pvalue": [0.01, 0.05, 0.2, 0.03, 0.04, 0.001],
        }
    )

    # Test without normalization
    expected_responsive = [True, True, False, False, False, False]
    result_df = label_responsive_genes(df, 0.6, 0.05)
    assert all(
        result_df["responsive"] == expected_responsive
    ), "Failed test without normalization"

    # Test with normalization (assuming normalization functionality is added)
    # Note: This test is placeholder and should be adjusted according to
    # the normalization implementation
    expected_responsive_normalized = [True, False, False, False, False, False]
    result_df_normalized = label_responsive_genes(df, 0.6, 0.05, 1)

    assert all(
        result_df_normalized["responsive"] == expected_responsive_normalized
    ), "Failed test with normalization"


def test_calculate_rank_response():
    df = pd.DataFrame(
        {
            "binding_effect": [1.2, 0.5, 0.8, 1.5, 0.4],
            "binding_pvalue": [0.05, 0.3, 0.004, 0.0002, 0.00006],
        }
    )
    # Manually create the expected DataFrame
    expected_df = df.copy()

    # Call the function
    bin_size = 2
    result_df = bin_by_binding_rank(df, bin_size, rank_by_binding_effect=True)

    expected_df = expected_df.sort_values(
        by=["binding_effect", "binding_pvalue"], ascending=[False, True]
    )

    # Assuming create_partitions works as expected
    expected_df["rank_bin"] = [2, 2, 4, 4, 6]

    # Assert the result is as expected
    pd.testing.assert_frame_equal(result_df, expected_df.reset_index(drop=True))

    result_df = bin_by_binding_rank(df, bin_size)

    expected_df = expected_df.sort_values(
        by=["binding_pvalue", "binding_effect"], ascending=[True, False]
    )
    # Assuming create_partitions works as expected
    expected_df["rank_bin"] = [2, 2, 4, 4, 6]

    # Assert the result is as expected
    pd.testing.assert_frame_equal(result_df, expected_df.reset_index(drop=True))


def test_compute_rank_response():
    df = pd.DataFrame(
        {
            "rank_bin": [2, 2, 4, 4, 6],
            "responsive": [True, False, False, False, True],
            "random": [0.5, 0.5, 0.5, 0.5, 0.5],
        }
    )

    # Call the function
    result_df = compute_rank_response(df)

    # Define the expected output
    # Note: The expected values here are placeholders and should be replaced
    # with the actual expected values based on your specific implementation of binom_test
    expected_df = pd.DataFrame(
        {
            "rank_bin": [2, 4, 6],
            "n_responsive_in_rank": [1, 0, 1],
            "random": [0.5, 0.5, 0.5],
            "n_successes": [1, 1, 2],
            # Replace with actual expected values
            "response_ratio": [0.5, 0.25, 0.333333],
            # Replace with actual expected values
            "pvalue": [1.0, 0.6250, 0.6875],
            # Replace with actual expected values
            "ci_lower": [0.012579, 0.006309, 0.043272],
            # Replace with actual expected values
            "ci_upper": [0.987421, 0.805880, 0.777222],
        }
    )

    # Specify the tolerance levels
    atol = 1e-5  # Absolute tolerance
    rtol = 1e-5  # Relative tolerance

    # Assert that the DataFrames are equal within the specified tolerance
    pd.testing.assert_frame_equal(result_df, expected_df, atol=atol, rtol=rtol)


def test_parse_binomtest_results():
    binom_test_res = binomtest(1, 2, 0.5, alternative="two-sided")
    actual = parse_binomtest_results(binom_test_res)
    assert actual == (0.5, 1.0, 0.01257911709367899, 0.987420882906321)


def test_validate_config(tmpdir):
    config_path = tmpdir.join("config.json")

    test_data_directory = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "test_data/yeast/Analysis",
    )

    assert os.path.isdir(test_data_directory) is True

    config = {
        "binding_data_path": [
            os.path.join(
                test_data_directory, "hap5_exprid_17_yiming_adh1_promoter_sig.csv.gz"
            )
        ],
        "binding_source": "cc_17",
        "binding_identifier_col": "target_gene_id",
        "binding_effect_col": "callingcards_enrichment",
        "binding_pvalue_col": "poisson_pval",
        "rank_by_binding_effect": False,
        "expression_data_path": [
            os.path.join(test_data_directory, "hap5_15min_mcisaac.csv.gz")
        ],
        "expression_source": "mcisaac_hap5_15",
        "expression_identifier_col": "gene_id",
        "expression_effect_col": "log2_shrunken_timecourses",
        "expression_effect_thres": 0.00,
        "expression_pvalue_col": None,
        "expression_pvalue_thres": None,
        "rank_bin_size": 5,
        "normalization_cutoff": -1,
    }

    # Write to the config file
    config_path.write(json.dumps(config), "w")

    # Assert that the config file exists
    assert os.path.exists(str(config_path))

    with open(config_path, "r", encoding="utf-8") as config_file:
        parsed_config = json.load(config_file)

    # Call the function
    actual = validate_config(parsed_config)

    assert isinstance(actual, dict)

    rank_response_df, _, _ = create_rank_response_table(actual)

    assert isinstance(rank_response_df, pd.DataFrame)

    args = argparse.Namespace(
        config=config_path,
        output_file=str(tmpdir.join("rank_response.csv")),
        compress=False,
    )

    rank_response_main(args)

    assert os.path.exists(str(tmpdir.join("rank_response.csv"))) is True
