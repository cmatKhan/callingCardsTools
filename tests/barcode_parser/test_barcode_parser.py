from math import inf as infinity

from callingcardstools.BarcodeParser import BarcodeParser
from .conftests import *


def test_constructor_yeast(yeast_barcode_details, yeast_barcode_dict):

    yeast_bp = BarcodeParser(yeast_barcode_details)

    assert yeast_bp.barcode_dict == yeast_barcode_dict


def test_yeast_barcode_check1(yeast_barcode_details):

    yeast_bp = BarcodeParser(yeast_barcode_details)

    component_dict = {'r1_primer': 'TGATA',
                      'r1_transposon': 'AATTCACTACGTCAACA',
                      'r2_transposon': 'ACCTGTTT',
                      'r2_restriction': 'TCGANCNCGTAA'}

    expected = {
        'passing': True,
        'details': {
            'r1_transposon': {'query': 'AATTCACTACGTCAACA',
                              'name': 'AATTCACTACGTCAACA',
                              'dist': 0, 'bam_tag': 'RT'},
            'r2_restriction': {'name': 'TaqAI',
                               'dist': 0,
                               'bam_tag': 'RS'},
            'tf': {'query': 'TGATAACCTGTTT',
                   'name': 'RIM101',
                   'dist': 0,
                   'bam_tag': 'TF'}
        }
    }

    actual = yeast_bp.component_check(component_dict)

    assert actual == expected


def test_yeast_barcode_check2(yeast_barcode_details):

    yeast_bp = BarcodeParser(yeast_barcode_details)

    component_dict = {'r1_primer': 'TGATA',
                      'r1_transposon': 'AATTCACTACGTCAACA',
                      'r2_transposon': 'ACCTGTTT',
                      'r2_restriction': 'TCGGGCGCCCGG'}
    expected = {
        'passing': True,
        'details': {
            'r1_transposon': {'query': 'AATTCACTACGTCAACA', 'name': 'AATTCACTACGTCAACA', 'dist': 0, 'bam_tag': 'RT'},
            'r2_restriction': {'name': '*', 'dist': infinity, 'bam_tag': 'RS'},
            'tf': {'query': 'TGATAACCTGTTT', 'name': 'RIM101', 'dist': 0, 'bam_tag': 'TF'}
        }
    }

    actual = yeast_bp.component_check(component_dict)

    assert actual == expected


def test_yeast_barcode_check3(yeast_barcode_details):

    yeast_bp = BarcodeParser(yeast_barcode_details)

    component_dict = {'r1_primer': 'TGATA',
                      'r1_transposon': 'AATTCACTACGTCAACA',
                      'r2_transposon': 'ACCTGCTT',
                      'r2_restriction': 'TCGAGCGCCCGG'}
    expected = {
        'passing': True,
        'details': {
            'r1_transposon': {'query': 'AATTCACTACGTCAACA',
                              'name': 'AATTCACTACGTCAACA',
                              'dist': 0,
                              'bam_tag': 'RT'},
            'r2_restriction': {'name': 'Hpall', 'dist': 0, 'bam_tag': 'RS'},
            'tf': {'query': 'TGATAACCTGCTT',
                   'name': 'RIM101',
                   'dist': 1,  # match allowance set to 1 for tf in test barcode details # noqa
                   'bam_tag': 'TF'}
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

    component_dict = {'r1_primer': 'TGATA',
                      'r1_transposon': 'AATTCACTACGTCAACA',
                      'r2_transposon': 'ACCTGCTT',
                      'r2_restriction': 'TCGAGCGCCCGG'}
    expected = {
        'passing': True,
        'details': {
            'r1_transposon': {'query': 'AATTCACTACGTCAACA',
                              'name': 'AATTCACTACGTCAACA',
                              'dist': 0,
                              'bam_tag': 'RT'},
            'r2_restriction': {'name': 'Hpall',
                               'dist': 0,
                               'bam_tag': 'RS'},
            'tf': {'query': 'TGATAACCTGCTT',
                   'name': 'RIM101',
                   'dist': 1,
                   'bam_tag': 'TF'}
        }
    }

    actual = yeast_bp.component_check(component_dict)

    assert actual == expected

def test_constructor(mouse_barcode_details):
	bp = BarcodeParser(mouse_barcode_details)
	assert bp.barcode_dict['tf'] == ''

def test_barcode_breakdown(mouse_barcode_details, mouse_barcodes):

	bp = BarcodeParser(mouse_barcode_details)
	assert bp.decompose_barcode(mouse_barcodes.get('dist0')).get('passing') == True
	
	pb_dist1 = bp.decompose_barcode(mouse_barcodes.get('pb_dist1'))
	assert pb_dist1.get('passing') == False
	assert pb_dist1.get('details').get('r1_pb').get('dist') == 1
    
	# note: turns out not actually an error -- this was part of debugging
	# which needs to be removed. It ended up showing that current behavior 
	# is correct
	error_bc = bp.decompose_barcode(mouse_barcodes.get('error'))
	assert error_bc.get('details').get('r1_lrt2').get('dist') == 2

def test_annotation_tag_list(mouse_barcode_details):
	bp = BarcodeParser(mouse_barcode_details)
	assert bp.annotation_tags == ['ST']
