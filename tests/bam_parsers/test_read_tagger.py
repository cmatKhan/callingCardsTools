from .conftests import *
from callingcardstools.bam_parsers import ReadTagger
from callingcardstools.bam_parsers import BarcodeParser
from callingcardstools.bam_parsers import StatusFlags

def test_human_constructor(human_readtagger):
    assert human_readtagger.is_open() == True

def test_yeast_constructor(yeast_readtagger):
    assert yeast_readtagger.is_open() == True

YEAST_RT_SETUP = {
    "fasta_path": "tests/test_data/yeast/NC_001133.9_1_60000.fasta",
    "fasta_index_path": "tests//test_data/yeast/NC_001133.9_1_60000.fasta.fai",
    "barcode_length":  34,
    "insertion_length": 1
}

YEAST_ALN = "tests/test_data/yeast/untagged.bam"
YEAST_BARCODE_DETAILS = "tests/test_data/yeast/barcode_details.json"


HUMAN_RT_SETUP = {
    "genome_path": "tests/test_data/human/chr1_grch38_1_8793000.fasta",
    "genome_index_path": "tests/test_data/human/chr1_grch38_1_8793000.fasta.fai",
    "barcode_length":  38,
    "insertion_length": 4
}

HUMAN_ALN = ""
HUMAN_BARCODE_DETAILS = ""

def test_ReadTagger_yeast():


    input_bamfile = pysam.AlignmentFile(YEAST_ALN, "rb")

    rt = ReadTagger(**YEAST_RT_SETUP)

    bp = BarcodeParser(YEAST_BARCODE_DETAILS)

    read_summary = []
    for read in input_bamfile.fetch():
        tagged_read = rt.tag_read(read)
        barcode_details = bp.barcode_check(tagged_read.get_tag("RG"))
        tagged_read.set_tag("XF", barcode_details['tf'])

        status_code = 0
        if not barcode_details['pass']:
            status_code += StatusFlags['BARCODE']
        if tagged_read.mapping_quality < mapq_threshold:
            status_code += StatusFlags['MAPQ']
        try:
            if tagged_read.get_tag("XZ") == bp.get_insert_seq():
                status_code += StatusFlags['INSERT_SEQ']
        except AttributeError as exc:
            print(f"insert sequence not found in Barcode Parser. {exc}")

        read_summary += {"id": tagged_read.query_name(),
                        "bc": tagged_read.get_tag("RG"),
                        "status": status_code,
                        "mapq": tagged_read.mapping_quality,
                        "three_prime": tagged_read.get("XT"),
                        "flag": tagged_read.flag,
                        "insert_start": tagged_read.get_tag("XI"),
                        "insert_stop": tagged_read.get_tag("XE"),
                        "insert_seq": tagged_read.get_tag("XZ")}
        
        assert 2==2
    
