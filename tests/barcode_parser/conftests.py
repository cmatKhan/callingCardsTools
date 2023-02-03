# pylint:disable=C0114
import pytest
from ..conftests import *  # pylint:disable=W0401,W0614 # noqa


@pytest.fixture
def yeast_barcode_dict():
    """expected result of parsing the yeast barcode details json

    Returns:
        dict: the expected barcode dict from the yeast
         barcode details json file
    """
    expected = {
        "r1": {
            "primer": {"trim": True,
                       "index": [0, 5]},

            "transposon": {"trim": True,
                           "index": [5, 22]}
        },
        "r2": {
            "transposon": {"trim": True,
                           "index": [0, 8]},

            "restriction": {"trim": True,
                            "index": [8, 20]}
        },
        "components": {
            "r1_transposon":  {"map": ["AATTCACTACGTCAACA"],
                               "bam_tag": "RT"},
            "r2_restriction": {"map": {"TCGAGCGCCCGG": "Hpall",
                                       "TCGAGCGC": "HinP1I",
                                       "TCGA": "TaqAI"},
                               "match_type": "greedy",
                               "require": False,
                               "bam_tag": "RS"},
            "tf":             {"components": ["r1_primer", "r2_transposon"],
                               "map": {"TGATAACCTGTTT": "RIM101",
                                       "TACTCCAGAGGGG": "RTG1",
                                       "CCTGCACGCACGC": "ERT1",
                                       "TCGTCTACGGCGT": "MTH1",
                                       "AACGCTTTTTCTA": "MET32"},
                               "bam_tag": "TF",
                               "match_allowance": 1}
        },
        "batch": "run_1234"
    }
    return expected
