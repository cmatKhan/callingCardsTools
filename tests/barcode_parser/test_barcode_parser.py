import pytest
from math import inf as infinity

from callingcardstools.BarcodeParser import BarcodeParser
from .conftests import *

def test_constructor_yeast(yeast_barcode_details, yeast_barcode_dict):

    yeast_bp = BarcodeParser(yeast_barcode_details)

    assert yeast_bp.barcode_dict == yeast_barcode_dict

def test_yeast_barcode_check1(yeast_barcode_details):

    yeast_bp = BarcodeParser(yeast_barcode_details)

    component_dict = {'r1_primer':'TGATA',
                      'r1_transposon': 'AATTCACTACGTCAACA',
                      'r2_transposon': 'ACCTGTTT',
                      'r2_restriction': 'TCGANCNCGTAA'}
    expected = {
        'passing': True,
        'details': {
            'r1_transposon': {'name':'AATTCACTACGTCAACA','dist':0, 'bam_tag':'RT'},
            'r2_restriction': {'name':'TaqAI','dist':0,'bam_tag':'RS'},
            'tf': {'name':'RIM101','dist':0, 'bam_tag': 'TF'}
        }
    }

    actual = yeast_bp.component_check(component_dict)

    assert actual == expected


def test_yeast_barcode_check2(yeast_barcode_details):

    yeast_bp = BarcodeParser(yeast_barcode_details)

    component_dict = {'r1_primer':'TGATA',
                      'r1_transposon': 'AATTCACTACGTCAACA',
                      'r2_transposon': 'ACCTGTTT',
                      'r2_restriction': 'TCGGGCGCCCGG'}
    expected = {
        'passing': True,
        'details': {
            'r1_transposon': {'name':'AATTCACTACGTCAACA','dist':0, 'bam_tag':'RT'},
            'r2_restriction': {'name':'*','dist':infinity,'bam_tag':'RS'},
            'tf': {'name':'RIM101','dist':0, 'bam_tag': 'TF'}
        }
    }

    actual = yeast_bp.component_check(component_dict)

    assert actual == expected


def test_yeast_barcode_check3(yeast_barcode_details):

    yeast_bp = BarcodeParser(yeast_barcode_details)

    component_dict = {'r1_primer':'TGATA',
                      'r1_transposon': 'AATTCACTACGTCAACA',
                      'r2_transposon': 'ACCTGCTT',
                      'r2_restriction': 'TCGAGCGCCCGG'}
    expected = {
        'passing': True,
        'details': {
            'r1_transposon': {'name':'AATTCACTACGTCAACA','dist':0, 'bam_tag':'RT'},
            'r2_restriction': {'name':'Hpall','dist':0,'bam_tag':'RS'},
            'tf': {'name':'RIM101','dist':1, 'bam_tag': 'TF'}
        }
    }

    actual = yeast_bp.component_check(component_dict)

    assert actual == expected

def test_get_tagged_components(yeast_barcode_details):

    yeast_bp = BarcodeParser(yeast_barcode_details)
    
    actual = yeast_bp.tagged_components
    expected = {'r1_transposon': 'RT', 'r2_restriction': 'RS', 'tf': 'TF'}

    assert actual == expected

def test_yeast_barcode_check4(yeast_barcode_details):

    yeast_bp = BarcodeParser(yeast_barcode_details)

    component_dict = {'r1_primer':'TGATA',
                      'r1_transposon': 'AATTCACTACGTCAACA',
                      'r2_transposon': 'ACCTGCTT',
                      'r2_restriction': 'TCGAGCGCCCGG'}
    expected = {
        'passing': True,
        'details': {
            'r1_transposon': {'name':'AATTCACTACGTCAACA','dist':0, 'bam_tag':'RT'},
            'r2_restriction': {'name':'Hpall','dist':0,'bam_tag':'RS'},
            'tf': {'name':'RIM101','dist':1, 'bam_tag': 'TF'}
        }
    }

    actual = yeast_bp.component_check(component_dict)

    assert actual == expected

def test_constructor_human(human_bp, human_barcode_dict):
    """test"""

    assert human_bp.barcode_dict == human_barcode_dict

def test_barcode_length_getter(yeast_bp, yeast_barcode_dict):
    """test barcode length getter"""
    assert yeast_bp.barcode_length == yeast_barcode_dict['length']

def test_insert_length_getter(human_bp, human_barcode_dict):
    """test barcode length getter"""
    assert human_bp.insert_length == human_barcode_dict['insert_length']

def test_inexact_component_edit_distance_with_diffs(yeast_bp, invalid_yeast_barcode):
    """test with barcode which is not exact"""
    
    expected = {
        'transposon_seq': 2,
        'r1_primer_bc': 2,
        'r2_trans_bc': 1,
        'restriction_site': 1}
        
    
    with pytest.raises(TypeError):
        yeast_bp.component_edit_distance()
    
    actual = yeast_bp.component_edit_distance(invalid_yeast_barcode)

    assert actual == expected
    
def test_exact_component_edit_distance(yeast_bp, valid_mig2_yeast_barcode):
    """test with barcode which is exact"""
    
    expected = {
        'transposon_seq': 0,
        'r1_primer_bc': 0,
        'r2_trans_bc': 0,
        'restriction_site': 0}
    
    with pytest.raises(TypeError):
        yeast_bp.component_edit_distance()
    
    barcode_summary = yeast_bp.component_edit_distance(valid_mig2_yeast_barcode)

    assert barcode_summary == expected

def test_edit_distance(yeast_bp, valid_mig2_yeast_barcode):

    # check that string matches itself with 0 edit dist
    assert yeast_bp.edit_dist(valid_mig2_yeast_barcode, valid_mig2_yeast_barcode) == 0
    
    # check that the known number of mismatches produces the correct result
    num_gs = sum([n=='G' for n in valid_mig2_yeast_barcode])
    barcode_with_mismatches = valid_mig2_yeast_barcode.replace('G', 'N', num_gs)
    assert yeast_bp.edit_dist(valid_mig2_yeast_barcode, barcode_with_mismatches) == num_gs

def test_get_insert_seq(human_bp, human_barcode_dict):
    assert human_bp.insert_seqs == human_barcode_dict["insert_seqs"]

def test_get_restriction_enzyme(yeast_bp, valid_mig2_yeast_barcode_Hall, valid_mig2_yeast_barcode_HinP1I):
    assert 'Hpall' == yeast_bp.get_restriction_enzyme(valid_mig2_yeast_barcode_Hall)
    assert "HinP1I" == yeast_bp.get_restriction_enzyme(valid_mig2_yeast_barcode_HinP1I)

def test_get_tf_barcode(yeast_bp, yeast_barcode_dict, valid_mig2_yeast_barcode):
    """test with barcode which is exact"""
    
    with pytest.raises(TypeError):
        yeast_bp.get_tf()

    expected = 'TCAGTCCCGTTGG'

    actual = yeast_bp.get_tf_barcode(valid_mig2_yeast_barcode) 

    assert actual == expected

def test_get_tf_named(yeast_bp, valid_mig2_yeast_barcode):
    """test with barcode which is exact"""
    
    with pytest.raises(TypeError):
        yeast_bp.get_tf()
    
    expected = 'MIG2'

    actual = yeast_bp.get_tf(valid_mig2_yeast_barcode) 

    assert actual == expected

def test_get_tf_unnamed(human_bp,valid_human_barcode):
    """test with barcode which is exact"""
    
    with pytest.raises(TypeError):
        human_bp.get_tf()

    expected = '*'

    actual = human_bp.get_tf(valid_human_barcode)

    assert actual == expected

def test_get_tf_best_match(yeast_bp, invalid_yeast_barcode):
    """test with barcode which is exact"""
    
    with pytest.raises(TypeError):
        yeast_bp.get_tf()
    
    expected = 'MIG2_3'

    actual = yeast_bp.get_tf(invalid_yeast_barcode) 

    assert actual == expected

def test_min_edit_dist(yeast_bp):
    """test with barcode which is exact"""
    
    transposon_seq_exact = 'AATTCACTACGTCAACA'
    transposon_seq_1_diff = 'NATTCACTACGTCAACA'
    transposon_seq_4_diff = 'AANTCANTACGTNAACN'

    assert yeast_bp.min_edit_dist('transposon_seq', transposon_seq_exact) == 0
    assert yeast_bp.min_edit_dist('transposon_seq', transposon_seq_1_diff) == 1
    assert yeast_bp.min_edit_dist('transposon_seq', transposon_seq_4_diff) == 4

def test_tf_dict(yeast_bp, yeast_barcode_dict):
    """test with barcode which is not exact"""

    assert yeast_bp.barcode_dict['tf_dict'] == yeast_barcode_dict['tf_dict']

def test_barcode_check(yeast_bp, valid_mig2_yeast_barcode_HinP1I):
    expected = {"pass": True, "tf": "MIG2", "restriction_enzyme": "HinP1I"}
    actual = yeast_bp.barcode_check(valid_mig2_yeast_barcode_HinP1I)

    assert actual == expected
