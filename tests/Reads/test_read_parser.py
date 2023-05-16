from argparse import Namespace
import os
from callingcardstools.Reads.ReadParser import ReadParser
from callingcardstools.BarcodeParser import BarcodeParser
from callingcardstools.BarcodeParser.yeast.BarcodeQcCounter import BarcodeQcCounter
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


def test_split_fastq(tmpdir, yeast_barcodeqccounter_data):

    pickle_qc = False
    verbose_qc = True

    r1, r2, barcode_details = yeast_barcodeqccounter_data

    ns = Namespace(
        read1=r1,  # 'tests/test_data/yeast/run_6177/r1_100_subsample.fq.gz',
        read2=r2,  # 'tests/test_data/yeast/run_6177/r2_100_subsample.fq.gz',
        # 'tests/test_data/yeast/run_6177/barcode_details.json',
        barcode_details=barcode_details,
        split_key='tf',
        split_suffix='split',
        output_dirpath=tmpdir,
        verbose_qc=verbose_qc,
        pickle_qc=pickle_qc)

    split_fastq(ns)

    if pickle_qc and verbose_qc:
        assert len(os.listdir(tmpdir)) == 19
    elif pickle_qc:
        assert len(os.listdir(tmpdir)) == 17
    elif verbose_qc:
        assert len(os.listdir(tmpdir)) == 19
    else:
        assert len(os.listdir(tmpdir)) == 18


def test_combine_split_qc(tmpdir, yeast_barcodeqccounter_data, yeast_split_qc_files):
    bp = BarcodeParser(yeast_barcodeqccounter_data[2])
    # construct the input to the BarcodeQcCounter summarize method
    component_dict = {k: [] for k in ['tf', 'r1_primer', 'r2_transposon']}
    r1_primer_start = bp.barcode_dict['r1']['primer']['index'][0]
    r1_primer_end = bp.barcode_dict['r1']['primer']['index'][1]
    r2_transposon_start = r1_primer_end + \
        bp.barcode_dict['r2']['transposon']['index'][0]
    r2_transposon_end = (r2_transposon_start +
                         bp.barcode_dict['r2']['transposon']['index'][1] -
                         bp.barcode_dict['r2']['transposon']['index'][0])
    for k, v in bp.barcode_dict['components']['tf']['map'].items():
        r1_primer_seq = k[r1_primer_start:r1_primer_end]
        r2_transposon_seq = k[r2_transposon_start:r2_transposon_end]
        component_dict['tf'].append(v)
        component_dict['r1_primer'].append(r1_primer_seq)
        component_dict['r2_transposon'].append(r2_transposon_seq)

    r1_bc_obj_list = [BarcodeQcCounter(x) for x in yeast_split_qc_files
                      .get('split')]

    r1_bc_combined = BarcodeQcCounter.combine(r1_bc_obj_list)

    r1_bc_combined.write(component_dict=component_dict,
                         output_dirpath=tmpdir)

    assert 42 == 42
