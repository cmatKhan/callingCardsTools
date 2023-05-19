
from argparse import Namespace
import os
import logging

import pandas as pd

from callingcardstools.Alignment.mammals import process_alignments
from callingcardstools.Alignment.mammals import combine_qc
from callingcardstools.Alignment.mammals.Qbed import Qbed
from callingcardstools.BarcodeParser.mammals.BarcodeQcCounter import BarcodeQcCounter


def test_parse_bam_plain(
        barcode_details_file,
        genome_file,
        bam_file,
        tmpdir,
        caplog):

    cur_dir = os.getcwd()

    os.chdir(tmpdir)

    caplog.set_level(logging.INFO)
    args = Namespace(input=bam_file,
                     barcode_details=barcode_details_file,
                     genome=genome_file,
                     suffix="",
                     pickle=False,
                     mapq_threshold=10)
    process_alignments.main(args)
    # todo -- test that the number of hops recorded in the qbed
    # is the same as the number of reads in the passing bam
    # note this is 3 b/c the summarize QC is not implemented at the moment
    assert len(os.listdir(tmpdir)) == 6

    os.chdir(cur_dir)


def test_parse_bam_pickle(
        barcode_details_file,
        genome_file,
        bam_file,
        tmpdir,
        caplog):

    cur_dir = os.getcwd()

    os.chdir(tmpdir)

    caplog.set_level(logging.INFO)
    args = Namespace(input=bam_file,
                     barcode_details=barcode_details_file,
                     genome=genome_file,
                     suffix="",
                     pickle=True,
                     mapq_threshold=10)
    process_alignments.main(args)
    # todo -- test that the number of hops recorded in the qbed
    # is the same as the number of reads in the passing bam
    # note this is 3 b/c the summarize QC is not implemented at the moment
    assert len(os.listdir(tmpdir)) == 4

    os.chdir(cur_dir)


def test_combine(parsed_bam_pickled, tmpdir, caplog):

    cur_dir = os.getcwd()

    os.chdir(tmpdir)

    # caplog.set_level(logging.INFO)

    qbed = Qbed(parsed_bam_pickled.get('qbed_pickle'))
    barcode_qc = BarcodeQcCounter(parsed_bam_pickled.get('barcode_qc_pickle'))

    qbed.write('single', '', False)
    barcode_qc.write('single', '', False)

    args = Namespace(qbed_list=[parsed_bam_pickled.get('qbed_pickle')]*3,
                     barcodeQcCounter_list=[
                         parsed_bam_pickled.get('barcode_qc_pickle')]*3,
                     filename="combined",
                     suffix="",
                     pickle=False)

    combine_qc.main(args)
    # todo -- test that the number of hops recorded in the qbed
    # is the same as the number of reads in the passing bam
    # note this is 3 b/c the summarize QC is not implemented at the moment
    assert len(os.listdir(tmpdir)) == 8

    suffixes = ['_aln_summary.tsv',
                '_barcode_qc.tsv',
                '.qbed',
                '_srt_count.tsv']

    single_filenames = ['single'+s for s in suffixes]
    single_paths = [os.path.join(tmpdir, f) for f in single_filenames]

    combined_filenames = ['combined'+s for s in suffixes]
    combined_paths = [os.path.join(tmpdir, f) for f in combined_filenames]

    single_dict = {k: pd.read_csv(v, sep='\t', header=None) for
                   k, v in zip(single_filenames, single_paths)}
    combined_dict = {k: pd.read_csv(v, sep='\t', header=None) for
                     k, v in zip(combined_filenames, combined_paths)}

    assert all(single_dict['single_aln_summary.tsv'].iloc[:, 1]*3 ==
               combined_dict['combined_aln_summary.tsv'].iloc[:, 1])

    assert [int(x)*3 for x
            in single_dict['single_barcode_qc.tsv'].iloc[1:, 4]] == \
        [int(x) for x in combined_dict['combined_barcode_qc.tsv'].iloc[1:, 4]]

    assert [int(x)*3 for x
            in single_dict['single.qbed'].iloc[1:, 4]] == \
        [int(x) for x in combined_dict['combined.qbed'].iloc[1:, 4]]

    assert all(single_dict['single_srt_count.tsv'].iloc[:, 1] ==
               combined_dict['combined_srt_count.tsv'].iloc[:, 1])

    os.chdir(cur_dir)
