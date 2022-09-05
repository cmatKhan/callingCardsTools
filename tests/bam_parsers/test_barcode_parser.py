import pytest
from .conftests import *

def test_BarcodeParser_constructor_yeast(yeast_bp, yeast_barcode_dict):
    """test""" 

    assert yeast_bp.barcode_dict == yeast_barcode_dict

def test_BarcodeParser_constructor_human(human_bp, human_barcode_dict):
    """test"""

    assert human_bp.barcode_dict == human_barcode_dict

def test_inexact_component_edit_distance_with_diffs(yeast_bp, invalid_yeast_barcode):
    """test with barcode which is not exact"""
    
    expected = {
        'transposon_seq': 2,
        'r1_primer_bc': 2,
        'r2_primer_bc': 1,
        'restriction_site': 0}
        
    
    with pytest.raises(AttributeError):
        yeast_bp.component_edit_distance()
    
    yeast_bp.set_barcode(invalid_yeast_barcode)
    actual = yeast_bp.component_edit_distance()

    assert actual == expected
    
def test_exact_component_edit_distance(yeast_bp, valid_mig2_yeast_barcode):
    """test with barcode which is exact"""
    
    expected = {
        'transposon_seq': 0,
        'r1_primer_bc': 0,
        'r2_primer_bc': 0,
        'restriction_site': 0}
    
    with pytest.raises(AttributeError):
        yeast_bp.component_edit_distance()
    
    yeast_bp.set_barcode(valid_mig2_yeast_barcode)
    barcode_summary = yeast_bp.component_edit_distance()

    assert barcode_summary == expected

def test_edit_distance(yeast_bp, valid_mig2_yeast_barcode):

    # check that string matches itself with 0 edit dist
    assert yeast_bp.edit_dist(valid_mig2_yeast_barcode, valid_mig2_yeast_barcode) == 0
    
    # check that the known number of mismatches produces the correct result
    num_gs = sum([n=='G' for n in valid_mig2_yeast_barcode])
    barcode_with_mismatches = valid_mig2_yeast_barcode.replace('G', 'N', num_gs)
    assert yeast_bp.edit_dist(valid_mig2_yeast_barcode, barcode_with_mismatches) == num_gs

def test_get_insert_seq(human_bp, human_barcode_dict):
    insert_seqs = human_bp.get_insert_seqs() 
    assert insert_seqs == human_barcode_dict["insert_seqs"]

def test_get_tf_barcode(yeast_bp, yeast_barcode_dict, valid_mig2_yeast_barcode):
    """test with barcode which is exact"""
    
    with pytest.raises(AttributeError):
        yeast_bp.get_tf()
    
    yeast_bp.set_barcode(valid_mig2_yeast_barcode)

    expected = 'TCAGTCCCGTTGG'

    actual = yeast_bp.get_tf_barcode() 

    assert actual == expected

def test_get_tf_named(yeast_bp, valid_mig2_yeast_barcode):
    """test with barcode which is exact"""
    
    with pytest.raises(AttributeError):
        yeast_bp.get_tf()
    
    yeast_bp.set_barcode(valid_mig2_yeast_barcode)

    expected = 'MIG2'

    actual = yeast_bp.get_tf() 

    assert actual == expected

def test_get_tf_unnamed(human_bp,valid_human_barcode):
    """test with barcode which is exact"""
    
    with pytest.raises(AttributeError):
        human_bp.get_tf()
    
    human_bp.set_barcode(valid_human_barcode)

    expected = '*'

    actual = human_bp.get_tf()

    assert actual == expected

def test_get_tf_best_match(yeast_bp, invalid_yeast_barcode):
    """test with barcode which is exact"""
    
    with pytest.raises(AttributeError):
        yeast_bp.get_tf()
    
    yeast_bp.set_barcode(invalid_yeast_barcode)

    expected = 'MIG2_3'

    actual = yeast_bp.get_tf() 

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
