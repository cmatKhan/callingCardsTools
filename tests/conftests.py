"""Set up bam_parser fixtures"""
import random
import pytest
import pysam
import pandas as pd
from callingcardstools.BarcodeParser import BarcodeParser
from callingcardstools.AlignmentTagger import AlignmentTagger
from callingcardstools.SummaryParser import SummaryParser
from callingcardstools.ReadParser import ReadParser

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
    bampath = "tests/test_data/yeast/untagged.bam"
    return pysam.AlignmentFile(bampath, "rb")


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


# Human fixtures --------------------------------------------------------------


@pytest.fixture
def human_barcode_details():
    """Path to the test yeast barcode details json file

    Returns:
        Str: Path to the yeast barcode details json
    """
    barcode_details = "tests/test_data/human/barcode_details.json"
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
        AlignmentTagger: A default read tagger with data from test_data
    """

    setup = {
        "fasta_path": "tests/test_data/human/chr1_grch38_1_8793000.fasta",
        "fasta_index_path": "tests/test_data/human/chr1_grch38_1_8793000.fasta.fai",
        "barcode_length":  38,
        "insertion_length": 4
    }

    rt = AlignmentTagger(**setup)

    return rt


@pytest.fixture
def human_bamfile():
    """an open AlignmentFile object for human untagged alignments"""
    bampath = "tests/test_data/human/untagged.bam"
    return pysam.AlignmentFile(bampath, "rb")


@pytest.fixture
def human_sp():
    summary_path = "tests/test_data/human/summary.csv"
    sp = SummaryParser(summary_path)
    return (sp)


@pytest.fixture
def background_hops():
    background_dict = {
        'chr': pd.Series(['chr1']*5 + ['chr2']*4, dtype=str),
        'start': pd.Series([1, 3, 5, 7, 10]+[1, 3, 5, 10], dtype=int),
        'end': pd.Series([x+1 for x in [1, 3, 5, 7, 10]+[1, 3, 5, 10]], dtype=int),
        'depth': pd.Series([random.randrange(1, 300, 1) for x in range(9)], dtype=int),
        'strand': pd.Series(['+', '+', '-', '+', '-']+['+', '-', '+', '-'], dtype=str)
    }


@pytest.fixture
def human_hops_data():
    """A set of region, background and hop files to test hopdb

    Returns:
        dict: a dict with keys regions, background and experiment storing 
        paths to corresponding bed/qbed files 
    """
    file_dict = {
        "chr_map": "src/callingcardstools/resources/human/chr_map.csv",
        "ttaa": "temp/human/TTAA_hg38.qbed",
        "experiment": "temp/human/TAG_AY53-1_50k_downsampled_human_map_sort_final.qbed"
    }

    return file_dict


@pytest.fixture
def human_hopsdb(human_hops_data):

    hops_data = human_hops_data

    hops_db = DatabaseApi(":memory:")

    chr_map_df = pd.read_csv(hops_data['chr_map'])

    hops_db.add_frame(chr_map_df, "chr_map")

    qbed_colnames = ['chr', 'start', 'end',
                     'depth', 'strand', 'annotation', 'sample']

    ttaa_df = pd.read_csv(
        hops_data['ttaa'],
        sep="\t",
        names=qbed_colnames
    )
    hops_db.add_frame(ttaa_df, 'ttaa')

    expr_df = pd.read_csv(
        hops_data['experiment'],
        sep="\t",
        names=qbed_colnames)

    hops_db.add_frame(expr_df, 'experiment', 'ay53-1_50k')

    return hops_db
