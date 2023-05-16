import pytest
import pathlib

from ..conftests import *


@pytest.fixture
def yeast_alignment_data():
    """paths to chr1 fasta (mitra lab yeast genome) and dal80 bam

    Returns:P
        tuple: chr1_fasta, dal80_bam
    """
    chr1_fasta = pathlib\
        .Path("tests/Alignment/data/yeast/chr2.fasta")\
        .resolve()
    dal80_bam = pathlib\
        .Path("tests/Alignment/data/yeast/run_6177_T1_DAL80_chr2.bam")\
        .resolve()
    barcode_details = pathlib\
        .Path("tests/Alignment/data/yeast/run_6177_barcode_details.json")\
        
    return chr1_fasta, dal80_bam, barcode_details
