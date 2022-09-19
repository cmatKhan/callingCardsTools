"""Set up bam_parser fixtures"""
import random
import pytest
import pysam
import pandas as pd
from callingcardstools.bam_parsers.BarcodeParser import BarcodeParser
from callingcardstools.bam_parsers.ReadTagger import ReadTagger
from callingcardstools.bam_parsers.SummaryParser import SummaryParser
from callingcardstools.fastq_parsers.BarcodeExtractor import BarcodeExtractor
from callingcardstools.significance.HopsDb import HopsDb

# yeast fixtures ---------------------------------------------------------------

@pytest.fixture
def yeast_reads():
    """paths relative to project root to yeast r1 and r2 fastq files

    Returns:
        dict: key, value where keys are r1, r2 and values are paths to fastq files
    """
    return {'r1':"tests/test_data/yeast/r1.fastq", 
            'r2':"tests/test_data/yeast/r2.fastq"}

@pytest.fixture
def yeast_be_full():
    be = BarcodeExtractor("tests/test_data/yeast/read_barcode_details_full.json")
    return be

@pytest.fixture
def yeast_be_concise():
    be = BarcodeExtractor("tests/test_data/read_barcode_details_concise.json")
    return be

@pytest.fixture
def yeast_be_trim():
    be = BarcodeExtractor("tests/test_data/read_barcode_details_trim.json")
    return be

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
            'r2_trans_bc': [22, 30], 
            'restriction_site': [30, 42]},
        'components': {
            'r1_primer_bc': ['TCAGT', 'GCCTG', 'ATTTG', 'TTGGT', 'CTCGG'],
            'transposon_seq': ['AATTCACTACGTCAACA'],
            'r2_trans_bc': ['CCCGTTGG', 'GGCGGCAG', 'GGGGGGGT', 'GGGGGTAG', 'TCGTCAGT'],
            'restriction_site': {"seq": ['TCGAGCGCCCGG', 'TCGAGCGC', 'TCGA'], 
                                 "name": ["Hpall", "HinP1I", "TaqAI"]}
            },
        'tf_map': {
            'bc_components': ['r1_primer_bc', 'r2_trans_bc'],
            'tf': ['MIG2', 'CAT8', 'GLN3', 'ARO80', 'CBF1']},
        'insert_seqs': ['*'],
        'match_allowance': {
            'transposon_seq': 1,
            'r1_primer_bc': 0,
            'r2_trans_bc': 0,
            'restriction_site': 1,
            'max': 1},
        'tf_dict':{
            'TCAGTCCCGTTGG': 'MIG2',
            'GCCTGGGCGGCAG': 'CAT8',
            'ATTTGGGGGGGGT': 'GLN3',
            'TTGGTGGGGGTAG': 'ARO80',
            'CTCGGTCGTCAGT': 'CBF1'
        },
        'length': 42,
        'insert_length': 1
    }

    return expected

@pytest.fixture
def valid_mig2_yeast_barcode():
    """Valid yeast barcode for a TF in the tf_map

    Returns:
        str: a full yeast barcode which corresponds to current yeast barcode 
        details json
    """
    return "TCAGTAATTCACTACGTCAACACCCGTTGGTCGANNNNNNNN"

@pytest.fixture
def valid_mig2_yeast_barcode_Hall():
    """Valid yeast barcode for a TF in the tf_map with the Hall restriction site

    Returns:
        str: a full yeast barcode which corresponds to current yeast barcode 
        details json
    """
    return "TCAGTAATTCACTACGTCAACACCCGTTGGTCGAGCGCCCGG"

@pytest.fixture
def valid_mig2_yeast_barcode_HinP1I():
    """Valid yeast barcode for a TF in the tf_map with the HpinP1I restriction site

    Returns:
        str: a full yeast barcode which corresponds to current yeast barcode 
        details json
    """
    return "TCAGTAATTCACTACGTCAACACCCGTTGGNNTCGAGCGCNN"
                                            
@pytest.fixture
def mig2_tseq_dist1_yeast_barcode():
    """cat8 barcode with edit distance 1 in transposon_seq

    Returns:
        str: a full yeast barcode which does not meet expectations
        details json
    """
    return "TCAGTNATTCACTACGTCAACACCCGTTGGTCGANNNNNNNN"

@pytest.fixture
def invalid_yeast_barcode():
    """Invalid yeast barcode, but created from MIG2. Edit distance from MIG2 is 
    5 with an edit distance of 3 in the barcode

    Returns:
        str: a full yeast barcode which does not meet expectations
        details json
    """
    return "NCANTAATTCACNACGTCANCACCCGNTGGTNGANNNNNNNN"

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

@pytest.fixture
def yeast_hopsdb(yeast_hops_data):

    hops_db = HopsDb(":memory:")

    chr_map_df = pd.read_csv(yeast_hops_data['chr_map'])

    hops_db.add_table(chr_map_df, "chr_map")

    regions_df = pd.read_csv(yeast_hops_data["regions"], sep = "\t", names = ['chr', 'start', 'end','name', 'value', 'strand'])

    hops_db.add_table(regions_df, 'regions', 'upstream_700')

    qbed_colnames = ['chr', 'start', 'end', 'depth', 'strand', 'annotation', 'sample']
    background_df = pd.read_csv(yeast_hops_data["background"], sep = "\t", names = qbed_colnames )

    hops_db.add_table(background_df, 'background', 'sir4')

    expr_df = pd.read_csv(yeast_hops_data['experiment'], sep = "\t", names = qbed_colnames)

    hops_db.add_table(expr_df, 'experiment', 'test')

    return hops_db

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
            "om_pb": 1,
            "pb_lrt1": 0,
            "srt": 1,
            "pb_lrt2": 1,
            "max": 2
        },
        "length": 38,
        "insert_length": 4
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

@pytest.fixture
def human_sp():
    summary_path = "tests/test_data/human/summary.csv"
    sp = SummaryParser(summary_path)
    return(sp)

@pytest.fixture
def background_hops():
    background_dict = {
        'chr': pd.Series(['chr1']*5 + ['chr2']*4, dtype = str),
        'start': pd.Series([1,3,5,7,10]+[1,3,5,10],dtype = int ),
        'end': pd.Series([x+1 for x in [1,3,5,7,10]+[1,3,5,10]], dtype = int),
        'depth': pd.Series([random.randrange(1,300,1) for x in range(9)], dtype = int),
        'strand': pd.Series(['+','+','-','+','-']+['+','-','+','-'], dtype = str)
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

    hops_db = HopsDb(":memory:")

    chr_map_df = pd.read_csv(hops_data['chr_map'])

    hops_db.add_table(chr_map_df, "chr_map")

    qbed_colnames = ['chr', 'start', 'end', 'depth', 'strand', 'annotation', 'sample']
    
    ttaa_df = pd.read_csv(
        hops_data['ttaa'],
        sep = "\t",
        names = qbed_colnames
    )
    hops_db.add_table(ttaa_df, 'ttaa')

    expr_df = pd.read_csv(
        hops_data['experiment'], 
        sep = "\t", 
        names = qbed_colnames)

    hops_db.add_table(expr_df, 'experiment', 'ay53-1_50k')

    return hops_db