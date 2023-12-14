import argparse
import os

import pandas as pd

from callingcardstools.PeakCalling.yeast.call_peaks import (call_peaks,
                                                            count_hops)
from callingcardstools.PeakCalling.yeast.call_peaks import \
    main as call_peaks_main


def count_hops_unstranded():
    # Setup example data
    promoter_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 200, 300, 400, 500],
        'end': [200, 300, 400, 500, 600],
        'name': ['promoter1', 'promoter2',
                 'promoter3', 'promoter4', 'promoter5'],
        'score': [100, 100, 100, 100, 100],
        'strand': ['+', '+', '+', '+', '+']
    })

    qbed_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [150, 250, 350, 450, 550],
        'end': [200, 300, 400, 500, 600],
        'strand': ['+', '-', '+', '-', '+'],
        'depth': [10, 100, 20, 200, 3]
    })

    # Expected output
    expected_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 200, 300, 400, 500],
        'end': [200, 300, 400, 500, 600],
        'name': ['promoter1', 'promoter2',
                 'promoter3', 'promoter4', 'promoter5'],
        'strand': ['+', '+', '+', '+', '+'],
        'experiment_hops': [1, 1, 1, 1, 1]
    })

    return promoter_df, qbed_df, expected_df


def count_hops1():
    # Setup example data
    promoter_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 200, 300, 400, 500],
        'end': [200, 300, 400, 500, 600],
        'name': ['promoter1', 'promoter2',
                 'promoter3', 'promoter4', 'promoter5'],
        'score': [100, 100, 100, 100, 100],
        'strand': ['+', '+', '+', '+', '+']
    })

    qbed_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [150, 250, 350, 450, 550],
        'end': [200, 300, 400, 500, 600],
        'strand': ['+', '-', '+', '-', '+'],
        'depth': [10, 100, 20, 200, 3]
    })

    # Expected output
    expected_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 200, 300, 400, 500],
        'end': [200, 300, 400, 500, 600],
        'name': ['promoter1', 'promoter2',
                 'promoter3', 'promoter4', 'promoter5'],
        'strand': ['+', '+', '+', '+', '+'],
        'experiment_hops': [1, 1, 1, 1, 1]
    })

    return promoter_df, qbed_df, expected_df


def count_hops_stranded():
    # Setup example data
    promoter_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 200, 300, 400, 500],
        'end': [200, 300, 400, 500, 600],
        'name': ['promoter1', 'promoter2',
                 'promoter3', 'promoter4', 'promoter5'],
        'score': [100, 100, 100, 100, 100],
        'strand': ['+', '+', '+', '+', '+']
    })

    qbed_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [150, 250, 350, 450, 550],
        'end': [200, 300, 400, 500, 600],
        'strand': ['+', '-', '+', '-', '+'],
        'depth': [10, 100, 20, 200, 3]
    })

    # Expected output
    expected_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1'],
        'start': [100, 300, 500],
        'end': [200, 400, 600],
        'name': ['promoter1', 'promoter3', 'promoter5'],
        'strand': ['+', '+', '+'],
        'experiment_hops': [1, 1, 1]
    })

    return promoter_df, qbed_df, expected_df


def count_hops_strand_testing():
    # Setup example data
    promoter_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 100, 300, 400, 500],
        'end': [200, 200, 400, 500, 600],
        'name': ['promoter1', 'promoter2',
                 'promoter3', 'promoter4', 'promoter5'],
        'score': [100, 100, 100, 100, 100],
        'strand': ['+', '-', '+', '+', '+']
    })

    qbed_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [150, 150, 350, 450, 550],
        'end': [200, 200, 400, 500, 600],
        'strand': ['+', '-', '+', '-', '+'],
        'depth': [10, 100, 20, 200, 3]
    })

    # Expected output
    expected_df_stranded = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 100, 300, 500],
        'end': [200, 200, 400, 600],
        'name': ['promoter1', 'promoter2', 'promoter3', 'promoter5'],
        'strand': ['+', '-', '+', '+'],
        'experiment_hops': [1, 1, 1, 1]
    })

    expected_df_unstranded = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 100, 300, 400, 500],
        'end': [200, 200, 400, 500, 600],
        'name': ['promoter1', 'promoter2', 'promoter3',
                 'promoter4', 'promoter5'],
        'strand': ['+', '-', '+', '+', '+'],
        'experiment_hops': [1, 1, 1, 1, 1]
    })

    return promoter_df, qbed_df, expected_df_stranded, expected_df_unstranded


def test_count_hops():

    promoter_df1, qbed_df1, expected_df1 = count_hops_unstranded()

    # Call the function
    result_df1 = count_hops(promoter_df1, qbed_df1, 'experiment_hops', False)

    # Assertions to check if the result is as expected
    pd.testing.assert_frame_equal(result_df1, expected_df1)

    promoter_df2, qbed_df2, expected_df2 = count_hops_stranded()

    # Call the function
    result_df2 = count_hops(promoter_df2, qbed_df2, 'experiment_hops', True)

    # Assertions to check if the result is as expected
    pd.testing.assert_frame_equal(result_df2, expected_df2)

    (promoter_df3,
     qbed_df3,
     expected_df3_stranded,
     expected_df3_unstranded) = count_hops_strand_testing()

    # Call the function
    result_df3 = count_hops(promoter_df3, qbed_df3, 'experiment_hops', True)

    # Assertions to check if the result is as expected
    pd.testing.assert_frame_equal(result_df3, expected_df3_stranded)

    # Call the function
    result_df4 = count_hops(promoter_df3, qbed_df3, 'experiment_hops', False)

    # Assertions to check if the result is as expected
    pd.testing.assert_frame_equal(result_df4, expected_df3_unstranded)


def test_call_peaks(tmp_path):

    chrmap_data = pd.DataFrame({
        'curr_chr_name_convention': ['chr1', 'chr2', 'chr3'],
        'new_chr_name_convention': ['chrI', 'chrII', 'chrIII'],
        'type': ['genomic', 'genomic', 'genomic']})
    promoter_data = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 100, 300, 400, 500],
        'end': [200, 200, 400, 500, 600],
        'name': ['promoter1', 'promoter2',
                 'promoter3', 'promoter4', 'promoter5'],
        'score': [100, 100, 100, 100, 100],
        'strand': ['+', '-', '+', '+', '+']
    })

    experiment_data = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [150, 150, 350, 450, 550],
        'end': [200, 200, 400, 500, 600],
        'depth': [10, 100, 20, 200, 3],
        'strand': ['+', '-', '+', '-', '+']
    })

    background_data = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [150, 150, 350, 450, 550],
        'end': [200, 200, 400, 500, 600],
        'depth': [10, 100, 20, 200, 3],
        'strand': ['+', '-', '+', '-', '+']
    })

    # Expected output
    expected_df_stranded = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 100, 300, 500],
        'end': [200, 200, 400, 600],
        'name': ['promoter1', 'promoter2', 'promoter3', 'promoter5'],
        'strand': ['+', '-', '+', '+'],
        'experiment_hops': [1, 1, 1, 1]
    })

    # Write mock data to temporary files
    experiment_data_path = tmp_path / "experiment.qbed"
    promoter_data_path = tmp_path / "promoter.bed"
    background_data_path = tmp_path / "background.qbed"
    chrmap_data_path = tmp_path / "chrmap.csv"

    experiment_data.to_csv(experiment_data_path, sep='\t', index=False)
    promoter_data.to_csv(promoter_data_path, sep='\t', index=False)
    background_data.to_csv(background_data_path, sep='\t', index=False)
    chrmap_data.to_csv(chrmap_data_path, index=False)

    # Call the function
    result_df = call_peaks(
        str(experiment_data_path),
        "curr_chr_name_convention",  # replace with actual convention
        str(promoter_data_path),
        "curr_chr_name_convention",    # replace with actual convention
        str(background_data_path),
        "curr_chr_name_convention",  # replace with actual convention
        str(chrmap_data_path),
        False,  # consider_strand
        "new_chr_name_convention"  # unified_chr_name_convention
    )

    # Assertions
    # Here you should write assertions that verify the DataFrame
    # returned by call_peaks is as expected. This will depend on the
    # structure and content of your mock data.
    assert isinstance(result_df, pd.DataFrame)
    # More specific assertions...


def test_with_data():
    test_data_directory = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        'test_data/yeast/Analysis')

    assert os.path.isdir(test_data_directory) is True

    args = argparse.Namespace(
        experiment_data_path=os.path.join(
            test_data_directory, 'hap5_expr17.qbed.gz'),
        experiment_orig_chr_convention='id',
        promoter_data_path=os.path.join(
            test_data_directory, 'yiming_promoters.bed.gz'),
        promoter_orig_chr_convention='id',
        background_data_path=os.path.join(
            test_data_directory, 'adh1_background.qbed.gz'),
        background_orig_chr_convention='ucsc',
        chrmap_data_path=os.path.join(test_data_directory, 'chrmap.csv.gz'),
        unified_chr_convention='ucsc',
        output_path='test_sig_output.csv',
        consider_strand=False,
        pseudocount=0.2,
        compress_output=False
    )

    call_peaks_main(args)

    assert True
