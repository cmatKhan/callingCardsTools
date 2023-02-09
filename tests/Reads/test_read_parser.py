
from _pytest._code.code import ExceptionRepr
from argparse import Namespace
import os
from callingcardstools.Reads.ReadParser import ReadParser
from .conftests import *
from Bio import SeqIO

from callingcardstools.Reads.split_fastq import split_fastq


def test_barcode_extractor_open(yeast_barcode_details, yeast_reads):
    rp = ReadParser(yeast_barcode_details,
                    yeast_reads['r1'], yeast_reads['r2'])

    assert rp.r1_path == yeast_reads['r1']
    assert rp.r2_path == yeast_reads['r2']
    assert rp.barcode_details_json == yeast_barcode_details

    rp.r1_path = yeast_reads['r1_gz']
    rp.r2_path = yeast_reads['r2_gz']

    assert rp.r1_path == yeast_reads['r1_gz']
    assert rp.r2_path == yeast_reads['r2_gz']

    rp.open()

    assert isinstance(rp.r1_handle, SeqIO.QualityIO.FastqPhredIterator)
    assert isinstance(rp.r2_handle, SeqIO.QualityIO.FastqPhredIterator)

    with pytest.raises(StopIteration):
        i = 0
        while True:
            rp.next()
            i += 1

    assert i == 100000


def test_barcode_extractor_process(yeast_barcode_details, yeast_reads):
    rp = ReadParser(yeast_barcode_details,
                    yeast_reads['r1'], yeast_reads['r2'])

    rp.open()

    rp.next()

    actual = rp.parse()

    assert actual['status']['passing'] == False
    assert actual['status']['details']['tf']['name'] == 'MTH1'


def test_split_fastq(tmpdir):

    ns = Namespace(
        read1='tests/test_data/yeast/run_6177/r1_100_subsample.fq.gz',
        read2='tests/test_data/yeast/run_6177/r2_100_subsample.fq.gz',
        barcode_details='tests/test_data/yeast/run_6177/barcode_details.json',
        split_key='tf',
        split_suffix='split',
        output_prefix=tmpdir)

    split_fastq(ns)

    assert len(os.listdir(tmpdir)) == 17
