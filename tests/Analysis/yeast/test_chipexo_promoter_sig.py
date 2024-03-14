import argparse
import os

import pandas as pd

from callingcardstools.Analysis.yeast.chipexo_promoter_sig import \
    chipexo_promoter_sig
from callingcardstools.Analysis.yeast.chipexo_promoter_sig import \
    main as chipexo_main


def test_chipexo_promoter_sig(tmp_path):
    # Create DataFrames
    chipexo_df = pd.DataFrame({
        'chr': ['chr1'],
        'start': [150],
        'end': [151],
        'YPD_log2Fold': [2.0],
        'YPD_log2P': [0.05]
    })

    promoter_df = pd.DataFrame({
        'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'start': [100, 100, 300, 400, 500],
        'end': [200, 200, 400, 500, 600],
        'name': ['promoter1', 'promoter2',
                 'promoter3', 'promoter4', 'promoter5'],
        'score': [100, 100, 100, 100, 100],
        'strand': ['+', '-', '+', '+', '+']
    })

    chrmap_df = pd.DataFrame({
        'chr': ['chr1'],
        'id': ['1'],
        'type': ['genomic']
    })

    # Write DataFrames to temporary files
    chipexo_path = tmp_path / "chipexo.tsv"
    chipexo_df.to_csv(chipexo_path, index=False)

    promoter_path = tmp_path / "promoter.tsv"
    promoter_df.to_csv(promoter_path, sep='\t', index=False)

    chrmap_path = tmp_path / "chrmap.tsv"
    chrmap_df.to_csv(chrmap_path, index=False)

    # Call the function
    result = chipexo_promoter_sig(chipexo_path,
                                  'chr',
                                  promoter_path,
                                  'chr',
                                  chrmap_path,
                                  'id')

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame)


def test_chipexo_data(tmpdir):
    test_data_directory = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        'test_data/yeast/Analysis')

    assert os.path.isdir(test_data_directory) is True

    chipexo_data_path = os.path.join(test_data_directory,
                                     'chipexo_10352.csv.gz')
    chrmap_data_path = os.path.join(test_data_directory,
                                    'chrmap.csv.gz')
    promoter_data_path = os.path.join(test_data_directory,
                                      'yiming_promoters.bed.gz')
    
    output_path = os.path.join(tmpdir, 'chipexo_promoter_sig.csv')

    args = argparse.Namespace(
        chipexo_data_path=chipexo_data_path,
        chipexo_orig_chr_convention='ucsc',
        unified_chr_convention='ucsc',
        chrmap_data_path=chrmap_data_path,
        promoter_data_path=promoter_data_path,
        promoter_orig_chr_convention='ucsc',
        output_file=output_path,
        compress=False)

    chipexo_main(args)

    assert os.path.isfile(output_path)
