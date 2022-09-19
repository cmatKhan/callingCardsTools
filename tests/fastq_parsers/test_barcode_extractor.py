from .conftests import *
from Bio import SeqIO

def test_barcode_extractor(yeast_be_full, yeast_be_concise):
    assert yeast_be_full.barcode_dict == yeast_be_concise.barcode_dict

def test_barcode_extractor(tmp_path, yeast_be_full, yeast_be):
    
    be = yeast_be

    r1 = SeqIO.parse(yeast_reads['r1'])
    r2 = SeqIO.parse(yeast_reads['r2'])

    processed_r1, processed_r1_bc = be.process_read(r1, 'r1')
    processed_r2, processed_r2_bc = be.process_read(r2, 'r2')