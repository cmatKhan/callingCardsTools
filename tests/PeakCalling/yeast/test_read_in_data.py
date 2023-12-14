import pandas as pd

from callingcardstools.PeakCalling.yeast import (read_in_background_data,
                                                 read_in_experiment_data,
                                                 read_in_promoter_data,
                                                 relabel_chr_column)

# Test for relabel_chr_column


def test_relabel_chr_column():
    # Setup
    data_df = pd.DataFrame({
        'chr': ['chr1', 'chr2', 'chr3'],
        'start': [1, 2, 3]})
    chrmap_df = pd.DataFrame({
        'curr_chr_name_convention': ['chr1', 'chr2', 'chr3'],
        'new_chr_name_convention': ['chrI', 'chrII', 'chrIII'],
        'type': ['genomic', 'genomic', 'genomic']})

    # Call function
    relabeled_df = relabel_chr_column(data_df, chrmap_df,
                                      'curr_chr_name_convention',
                                      'new_chr_name_convention')

    # Assertions
    assert list(relabeled_df.columns) == ['chr', 'start']
    assert list(relabeled_df['chr']) == ['chrI', 'chrII', 'chrIII']

# Test for read_in_experiment_data


def test_read_in_experiment_data(tmp_path):
    # Setup
    qbed_data = pd.DataFrame({'chr': ['chr1'],
                              'start': [1],
                              'end': [2],
                              'depth': [1],
                              'strand': ['+']})
    qbed_file = tmp_path / "test.qbed"
    qbed_data.to_csv(qbed_file, sep='\t', index=False, header=True)

    chrmap_df = pd.DataFrame({'curr_chr_name_convention': ['chr1'],
                              'new_chr_name_convention': ['chrI'],
                              'type': ['genomic']})

    # Call function
    experiment_df, experiment_total_hops = \
        read_in_experiment_data(str(qbed_file),
                                'curr_chr_name_convention',
                                'new_chr_name_convention',
                                chrmap_df)

    # Assertions
    assert list(experiment_df.columns) == ['chr', 'start',
                                           'end', 'depth', 'strand']
    assert experiment_total_hops == 1

# Test for read_in_promoter_data


def test_read_in_promoter_data(tmp_path):
    # Setup
    bed_data = pd.DataFrame({'chr': ['chr1'],
                             'start': [1],
                             'end': [2],
                             'name': ['test'],
                             'score': [1.0],
                             'strand': ['+']})

    bed_file = tmp_path / "test.bed"
    bed_data.to_csv(bed_file, sep='\t', index=False, header=True)

    chrmap_df = pd.DataFrame({'curr_chr_name_convention': ['chr1'],
                              'new_chr_name_convention': ['chrI'],
                              'type': ['genomic']})

    # Call function
    promoter_df = read_in_promoter_data(str(bed_file),
                                        'curr_chr_name_convention',
                                        'new_chr_name_convention',
                                        chrmap_df)

    # Assertions
    assert list(promoter_df.columns) == ['chr', 'start', 'end',
                                         'name', 'score', 'strand']
    assert len(promoter_df) == 1


def test_read_in_background_data(tmp_path):
    # Setup
    background_data = pd.DataFrame({
        'chr': ['chr1'],
        'start': [1],
        'end': [2],
        'depth': [1],
        'strand': ['+']})
    background_file = tmp_path / "background.qbed"
    background_data.to_csv(background_file, sep='\t', index=False, header=True)

    chrmap_df = pd.DataFrame({'curr_chr_name_convention': ['chr1'],
                              'new_chr_name_convention': ['chrI'],
                              'type': ['genomic']})

    # Call function
    background_df, background_total_hops = \
        read_in_background_data(str(background_file),
                                'curr_chr_name_convention',
                                'new_chr_name_convention',
                                chrmap_df)

    # Assertions
    assert list(background_df.columns) == ['chr', 'start', 'end',
                                           'depth', 'strand']
    assert background_total_hops == 1
