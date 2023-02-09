"""Set up bam_parser fixtures"""
import random
import pytest
import pathlib
import pandas as pd
from callingcardstools.BarcodeParser.BarcodeParser import BarcodeParser
from callingcardstools.Alignment.AlignmentTagger import AlignmentTagger

# yeast fixtures --------------------------------------------------------------


@pytest.fixture
def yeast_barcode_details():
    """Path to the test yeast barcode details json file

    Returns:
        Str: Path to the yeast barcode details json
    """
    barcode_details = pathlib\
        .Path("tests/test_data/yeast/barcode_details.json")\
        .resolve()
    return barcode_details


@pytest.fixture
def yeast_bp():
    """create a barcode parser for yeast

    Returns:
        BarcodeParser: a barcode parser obj with the test barcode
         details for yeast
    """
    bp = BarcodeParser("tests/test_data/yeast/barcode_details.json")
    return bp


@pytest.fixture
def yeast_readtagger():
    """A read tagger object for yeast data

    Returns:
        AlignmentTagger: A default read tagger with data from test_data
    """

    setup = {
        "barcode_details_json": "tests/test_data/yeast/barcode_details.json",
        "fasta_path": "tests/test_data/yeast/NC_001133.9_1_60000.fasta",
    }

    rt = AlignmentTagger(**setup)

    return rt


@pytest.fixture
def yeast_bamfile():
    """An open AlignmentFile object of yeast untagged yeast alignments"""
    bampath = pathlib.Path("tests/test_data/yeast/untagged.bam").resolve()
    return bampath


@pytest.fixture
def yeast_fasta():
    """An open AlignmentFile object of yeast untagged yeast alignments"""
    fasta = pathlib\
        .Path("tests/test_data/yeast/NC_001133.9_1_60000.fasta")\
        .resolve()
    return fasta


@pytest.fixture
def yeast_hops_data():
    """A set of region, background and hop files to test hopdb

    Returns:
        dict: a dict with keys regions, background and experiment storing 
        paths to corresponding bed/qbed files 
    """
    file_dict = {
        "chr_map": "src/callingcardstools/resources/yeast/chr_map.csv",
        "regions": "tests/test_data/yeast/regions.bed",
        "background": "tests/test_data/yeast/background.qbed",
        "experiment": "tests/test_data/yeast/experiment.qbed"
    }

    return file_dict


# Mouse fixtures --------------------------------------------------------------

@pytest.fixture
def mouse_barcode_details():
    """Path to the test yeast barcode details json file

    Returns:
        Str: Path to the yeast barcode details json
    """
    barcode_details = "tests/test_data/mouse/barcode_details.json"
    return barcode_details


@pytest.fixture(autouse=True)
def mouse_barcodes():
    d = {
        "dist0": "TAGCGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA",
        "pb_dist1": "TGGCGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA",
        'lrt_dist1': "TAGCGTCAANTTTACGCAGACTATCTTTCTAGGGTTAA",
        "error": "TAGCGTCAATTTTACGCAGACTATCTTTCTCTGGTAAT"
    }
    return d


@pytest.fixture
def background_hops():
    background_dict = {
        'chr': pd.Series(['chr1']*5 + ['chr2']*4, dtype=str),
        'start': pd.Series([1, 3, 5, 7, 10]+[1, 3, 5, 10], dtype=int),
        'end': pd.Series([x+1 for x in [1, 3, 5, 7, 10]+[1, 3, 5, 10]], dtype=int),  # noqa
        'depth': pd.Series([random.randrange(1, 300, 1) for x in range(9)], dtype=int),  # noqa
        'strand': pd.Series(['+', '+', '-', '+', '-']+['+', '-', '+', '-'], dtype=str)  # noqa
    }
