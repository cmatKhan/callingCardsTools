"""Set up bam_parser fixtures"""
import pytest
import pysam
from callingcardstools.bam_parsers.BarcodeParser import BarcodeParser
from callingcardstools.bam_parsers.ReadTagger import ReadTagger

# yeast fixtures ---------------------------------------------------------------

@pytest.fixture
def yeast_barcode_details():
    """Path to the test yeast barcode details json file

    Returns:
        Str: Path to the yeast barcode details json
    """
    barcode_details = "tests/test_data/yeast/barcode_details.json"
    return barcode_details

@pytest.fixture
def yeast_barcode_dict():
    """expected result of parsing the yeast barcode details json

    Returns:
        dict: the expected barcode dict from the yeast barcode details json file
    """
    expected = {
        'indicies':{
            'r1_primer_bc': [0, 5], 
            'transposon_seq': [5, 22],
            'r2_primer_bc': [22, 30], 
            'restriction_site': [30, 34]},
        'components': {
            'r1_primer_bc': ['TCAGT', 'GCCTG', 'ATTTG', 'TTGGT', 'CTCGG'],
            'transposon_seq': ['AATTCACTACGTCAACA'],
            'r2_primer_bc': ['CCCGTTGG', 'GGCGGCAG', 'GGGGGGGT', 'GGGGGTAG', 'TCGTCAGT'],
            'restriction_site': ['TCGA', 'GCGC', 'CCGG']},
        'tf_map': {
            'bc_components': ['r1_primer_bc', 'r2_primer_bc'],
            'tf': ['MIG2', 'CAT8', 'GLN3', 'ARO80', 'CBF1']},
        'insert_seqs': ['*'],
        'match_allowance': {
            'transposon_seq': 1,
            'r1_primer_bc': 0,
            'r2_primer_bc': 0,
            'restriction_site': 1,
            'max': 1},
        'tf_dict':{
            'TCAGTCCCGTTGG': 'MIG2',
            'GCCTGGGCGGCAG': 'CAT8',
            'ATTTGGGGGGGGT': 'GLN3',
            'TTGGTGGGGGTAG': 'ARO80',
            'CTCGGTCGTCAGT': 'CBF1'
        }
    }

    return expected

@pytest.fixture
def valid_mig2_yeast_barcode():
    """Valid yeast barcode for a TF in the tf_map

    Returns:
        str: a full yeast barcode which corresponds to current yeast barcode 
        details json
    """
    return "TCAGTAATTCACTACGTCAACACCCGTTGGTCGA"

@pytest.fixture
def mig2_tseq_dist1_yeast_barcode():
    """cat8 barcode with edit distance 1 in transposon_seq

    Returns:
        str: a full yeast barcode which does not meet expectations
        details json
    """
    return "TCAGTNATTCACTACGTCAACACCCGTTGGTCGA"

@pytest.fixture
def invalid_yeast_barcode():
    """Invalid yeast barcode, but created from MIG2. Edit distance from MIG2 is 
    5 with an edit distance of 3 in the barcode

    Returns:
        str: a full yeast barcode which does not meet expectations
        details json
    """
    return "NCANTAATTCACNACGTCANCACCCGNTGGTCGA"

@pytest.fixture
def yeast_bp():
    """create a barcode parser for yeast

    Returns:
        BarcodeParser: a barcode parser obj with the test barcode details for yeast
    """
    bp = BarcodeParser("tests/test_data/yeast/barcode_details.json")
    return bp

@pytest.fixture
def yeast_readtagger():
    """A read tagger object for yeast data

    Returns:
        ReadTagger: A default read tagger with data from test_data
    """

    setup = {
        "fasta_path": "tests/test_data/yeast/NC_001133.9_1_60000.fasta",
        "fasta_index_path": "tests//test_data/yeast/NC_001133.9_1_60000.fasta.fai",
        "barcode_length":  34,
        "insertion_length": 1
    }
    
    rt = ReadTagger(**setup)

    return rt

@pytest.fixture
def yeast_bamfile():
    """An open AlignmentFile object of yeast untagged yeast alignments"""
    bampath = "tests/test_data/yeast/untagged.bam"
    return pysam.AlignmentFile(bampath, "rb")

# Human fixtures ---------------------------------------------------------------

@pytest.fixture
def human_barcode_dict():
    """expected result of parsing the human barcode details json

    Returns:
        dict: the expected barcode dict from the human barcode details json file
    """
     
    expected = {
        "indicies": {
            "om_pb": [0, 3],
            "pb_lrt1": [3, 28],
            "srt": [28, 32],
            "pb_lrt2": [32, 38]
        },
        "components": {
            "om_pb":   ["TAG"],
            "pb_lrt1": ["CGTCAATTTTACGCAGACTATCTTT"],
            "srt":     ["CTAG", "CAAC", "CTGA", "GCAT", "GTAC", "CACA", "TGAC", "GTCA",
                        "CGAT", "CTCT", "GAAG", "TCGA", "CATG", "GTTG", "CTTC", "GCTA",
                        "GAGA", "GTGT", "CGTA", "TGGT", "GGAA", "ACAC", "TCAG", "TTGG",
                        "CAGT", "TTTT"],
            "pb_lrt2": ["GGTTAA"]
        },
        "insert_seqs": ["TTAA"],
        "match_allowance" : {
            "om_pb": 0,
            "pb_lrt1": 0,
            "srt": 0,
            "pb_lrt2": 0,
            "max": 0
        }
    }

    return expected

@pytest.fixture
def human_barcode_details():
    """Path to the test yeast barcode details json file

    Returns:
        Str: Path to the yeast barcode details json
    """
    barcode_details =  "tests/test_data/human/barcode_details.json"
    return barcode_details

@pytest.fixture
def valid_human_barcode():
    """Valid yeast barcode for a TF in the tf_map

    Returns:
        str: a full yeast barcode which corresponds to current yeast barcode 
        details json
    """
    return "TAGCGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA"

@pytest.fixture
def invalid_human_barcode():
    """Invalid human barcode -- first char of each component set to N, except 
    pb_lrt2 where the first and last char are set to N

    Returns:
        str: a full human barcode which is off by 1 in all but pb_lrt2 which 
        is off by 2
        details json
    """
    return "NAGNGTCAATTTTACGCAGACTATCTTTNTAGNGTTAA"

@pytest.fixture
def human_bp():
    """create a barcode parser for human

    Returns:
        BarcodeParser: a barcode parser obj with the test barcode details for human
    """
    bp = BarcodeParser("tests/test_data/human/barcode_details.json")
    return bp

@pytest.fixture
def human_readtagger():
    """A read tagger object for human data

    Returns:
        ReadTagger: A default read tagger with data from test_data
    """

    setup = {
        "fasta_path": "tests/test_data/human/chr1_grch38_1_8793000.fasta",
        "fasta_index_path": "tests/test_data/human/chr1_grch38_1_8793000.fasta.fai",
        "barcode_length":  38,
        "insertion_length": 4
    }

    rt = ReadTagger(**setup)

    return rt

@pytest.fixture
def human_bamfile():
    """an open AlignmentFile object for human untagged alignments"""
    bampath = "tests/test_data/human/untagged.bam"
    return pysam.AlignmentFile(bampath, "rb")