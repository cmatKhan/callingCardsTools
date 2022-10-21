
from _pytest._code.code import ExceptionRepr
from callingcardstools.ReadParser import ReadParser
from .conftests import *
from Bio import SeqIO

def test_barcode_extractor_open(yeast_barcode_details, yeast_reads):
    rp = ReadParser(yeast_barcode_details, yeast_reads['r1'], yeast_reads['r2'])
    
    assert rp.r1_path == yeast_reads['r1']
    assert rp.r2_path == yeast_reads['r2']
    assert rp.barcode_details_json == yeast_barcode_details
    
    rp.r1_path = yeast_reads['r1_gz']
    rp.r2_path = yeast_reads['r2_gz']
    
    assert rp.r1_path == yeast_reads['r1_gz']
    assert rp.r2_path == yeast_reads['r2_gz']

    rp.open()

    assert isinstance(rp.r1_handle,SeqIO.QualityIO.FastqPhredIterator)
    assert isinstance(rp.r2_handle,SeqIO.QualityIO.FastqPhredIterator)

    with pytest.raises(StopIteration):
        i=0
        while True:
            rp.next()
            i+=1
    
    assert i == 100000

def test_barcode_extractor_process(yeast_barcode_details, yeast_reads):
    rp = ReadParser(yeast_barcode_details, yeast_reads['r1'], yeast_reads['r2'])

    rp.open()

    rp.next()

    actual = rp.parse()

    assert actual['bc_check']['pass'] == False
    assert actual['bc_check']['tf'] == 'MIG2_5'

    # actual = be.process_read(r1,r2)
    # expected_r1_primer = "TCAGT"
    # expected_r1_transposon = "AATTCACTACGTCAACA"
    # expected_r2_primer = 'CCCGTTGG'
    # expected_r2_restriction = 'TCGAGCGCCCGG'
    # expected_bc = expected_r1_primer+expected_r1_transposon+\
    #     expected_r2_primer+expected_r2_restriction
    
    # assert actual['bc'] == expected_bc

    # assert len(actual['r1']) == \
    #     (len(r1) - len(expected_r1_primer) - len(expected_r1_transposon))

    # assert actual['r1'].seq == \
    #     r1.seq[len(expected_r1_primer)+len(expected_r1_transposon):]
    
    # assert len(actual['r2']) == \
    #     (len(r2) - len(expected_r2_primer) - len(expected_r2_restriction))
    
    # assert actual['r2'].seq == \
    #     r2.seq[len(expected_r2_primer)+len(expected_r2_restriction):]