from .conftests import *
from Bio import SeqIO

def test_read_in_reads(yeast_reads):
    r1 = next(SeqIO.parse(yeast_reads['r1'], format="fastq"))
    r2 = next(SeqIO.parse(yeast_reads['r2'], format="fastq"))

    assert len(r1) == len(r2)


def test_barcode_extractor(tmp_path, yeast_be, yeast_reads):
    
    be = yeast_be

    r1 = next(SeqIO.parse(yeast_reads['r1'], format="fastq"))
    r2 = next(SeqIO.parse(yeast_reads['r2'], format="fastq"))

    assert len(r1) == len(r2)
    
    actual = be.process_read(r1,r2)
    expected_r1_primer = "TCAGT"
    expected_r1_transposon = "AATTCACTACGTCAACA"
    expected_r2_primer = 'CCCGTTGG'
    expected_r2_restriction = 'TCGAGCGCCCGG'
    expected_bc = expected_r1_primer+expected_r1_transposon+\
        expected_r2_primer+expected_r2_restriction
    
    assert actual['bc'] == expected_bc

    assert len(actual['r1']) == \
        (len(r1) - len(expected_r1_primer) - len(expected_r1_transposon))

    assert actual['r1'].seq == \
        r1.seq[len(expected_r1_primer)+len(expected_r1_transposon):]
    
    assert len(actual['r2']) == \
        (len(r2) - len(expected_r2_primer) - len(expected_r2_restriction))
    
    assert actual['r2'].seq == \
        r2.seq[len(expected_r2_primer)+len(expected_r2_restriction):]