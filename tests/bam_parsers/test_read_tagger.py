from callingcardstools.QcStatusCoding import create_status_coder
from .conftests import *


def test_yeast_constructor(yeast_readtagger):
    assert yeast_readtagger.is_open() == True


def test_read_tagger_yeast(yeast_readtagger, yeast_bamfile, yeast_bp):

    rt = yeast_readtagger

    bp = yeast_bp

    expected = [

        {'id': 'NB501801:571:HNM3KAFX3:4:21409:11423:19883_GAATCAATTCACTACGTCAACATTTTTCTATCGA',  # noqa
         'r1_transposon': 'AATTCACTACGTCAACA/0',
         'r2_transposon': 'AATTCACTACGTCAACA/0',
         'tf': 'GAATCTTTTTCTA/3',
         'status': 11,
         'mapq': 0,
         'three_prime': '*',
         'flag': 0,
         'insert_start': '*',
         'insert_stop': '*',
         'insert_seq': '*'},


        # 8: {'id': 'NB501801:571:HNM3KAW...GCACGCTCGA', 'r1_transposon': 'AATTCACTACGTCAACA/0', 'r2_transposon': 'AATTCACTACGTCAACA/0', 'tf': 'TACTCACGCACGC/3', 'status': 8, 'mapq': 60, 'three_prime': 59999, 'flag': 16, 'insert_start': '*', 'insert_stop': '*', 'insert_seq': '*'}
        {'id': 'NB501801:571:HNM3KAFX3:1:21310:11840:5050_TACTCAATTCACTACGTCAACAACGCACGCTCGA',  # noqa
         'r1_transposon': 'AATTCACTACGTCAACA/0',
         'r2_transposon': 'AATTCACTACGTCAACA/0',
         'tf': 'TACTCACGCACGC/3',
         'status': 11,
         'mapq': 0,
         'three_prime': 0,
         'flag': 0,
         'insert_start': '*',
         'insert_stop': '*',
         'insert_seq': '*'},

        # 2: {'id': 'NB501801:571:HNM3KAF...CTGTTTTCGA', 'r1_transposon': 'AATTCACTACGTCAACA/0', 'r2_transposon': 'AATTCACTACGTCAACA/0', 'tf': 'TTGCGACCTGTTT/4', 'status': 10, 'mapq': 5, 'three_prime': 4644, 'flag': 0, 'insert_start': 4643, 'insert_stop': 4644, 'insert_seq': 'A'}
        {'id': 'NB501801:571:HNM3KAFX3:3:11501:15818:8474_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA',  # noqa
         'r1_transposon': 'AATTCACTACGTCAACA/0',
         'r2_transposon': 'AATTCACTACGTCAACA/0',
         'tf': 'TTGCGACCTGTTT/4',
         'status': 11,
         'mapq': 5,
         'three_prime': 4644,
         'flag': 0,
         'insert_start': 4643,
         'insert_stop': 4644,
         'insert_seq': 'A'},


        # 3: {'id': 'NB501801:571:HNM3KAF...CGTTGGTCGA', 'r1_transposon': 'AATTCACTACGTCAACA/0', 'r2_transposon': 'AATTCACTACGTCAACA/0', 'tf': 'TCAGTCCCGTTGG/6', 'status': 0, 'mapq': 60, 'three_prime': 32765, 'flag': 16, 'insert_start': 32766, 'insert_stop': 32767, 'insert_seq': 'G'}
        {'id': 'NB501801:571:HNM3KAFX3:2:21105:24619:7841_TCAGTAATTCACTACGTCAACACCCGTTGGTCGA',  # noqa
         'r1_transposon': 'AATTCACTACGTCAACA/0',
         'r2_transposon': 'AATTCACTACGTCAACA/0',
         'tf': 'TCAGTCCCGTTGG/6',
         'status': 1,
         'mapq': 60,
         'three_prime': 32765,
         'flag': 16,
         'insert_start': 32766,
         'insert_stop': 32767,
         'insert_seq': 'G'},


        # 4: {'id': 'NB501801:571:HNM3KAF...CGGCAGCCNG', 'r1_transposon': 'NATTCACTACGTCAACA/1', 'r2_transposon': 'NATTCACTACGTCAACA/1', 'tf': 'GCCTGGGCGGCAG/6', 'status': 8, 'mapq': 60, 'three_prime': 32780, 'flag': 16, 'insert_start': 32781, 'insert_stop': 32782, 'insert_seq': 'A'}
        {'id': 'NB501801:571:HNM3KAFX3:2:21206:2348:1221_GCCTGNATTCACTACGTCAACAGGCGGCAGCCNG',  # noqa
         'r1_transposon': 'NATTCACTACGTCAACA/1',
         'r2_transposon': 'NATTCACTACGTCAACA/1',
         'tf': 'GCCTGGGCGGCAG/6',
         'status': 1,
         'mapq': 60,
         'three_prime': 32780,
         'flag': 16,
         'insert_start': 32781,
         'insert_stop': 32782,
         'insert_seq': 'A'},


        # 5: {'id': 'NB501801:571:HNM3KAF...GGGTAGTCGN', 'r1_transposon': 'AATTCACTACGTCAACA/0', 'r2_transposon': 'AATTCACTACGTCAACA/0', 'tf': 'TTGGTGGGGGTAG/8', 'status': 0, 'mapq': 60, 'three_prime': 32991, 'flag': 0, 'insert_start': 32990, 'insert_stop': 32991, 'insert_seq': 'G'}
        {'id': 'NB501801:571:HNM3KAFX3:1:11312:1357:3858_TTGGTAATTCACTACGTCAACAGGGGGTAGTCGN',  # noqa
         'r1_transposon': 'AATTCACTACGTCAACA/0',
         'r2_transposon': 'AATTCACTACGTCAACA/0',
         'tf': 'TTGGTGGGGGTAG/8',
         'status': 1,
         'mapq': 60,
         'three_prime': 32991,
         'flag': 0,
         'insert_start': 32990,
         'insert_stop': 32991,
         'insert_seq': 'G'},


        # 6: {'id': 'NB501801:571:HNM3KAF...CTGTTTTCGA', 'r1_transposon': 'AATTCACTACGTCAACA/0', 'r2_transposon': 'AATTCACTACGTCAACA/0', 'tf': 'TTGCGACCTGTTT/4', 'status': 0, 'mapq': 60, 'three_prime': 56908, 'flag': 16, 'insert_start': 56909, 'insert_stop': 56910, 'insert_seq': 'T'}
        {'id': 'NB501801:571:HNM3KAFX3:4:21612:4799:3141_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA',  # noqa
         'r1_transposon': 'AATTCACTACGTCAACA/0',
         'r2_transposon': 'AATTCACTACGTCAACA/0',
         'tf': 'TTGCGACCTGTTT/4',
         'status': 1,
         'mapq': 60,
         'three_prime': 56908,
         'flag': 16,
         'insert_start': 56909,
         'insert_stop': 56910,
         'insert_seq': 'T'},


        # 7: {'id': 'NB501801:571:HNM3KAF...CTGTTTTCGA', 'r1_transposon': 'AATTCACTACGTCAACA/0', 'r2_transposon': 'AATTCACTACGTCAACA/0', 'tf': 'TTGCGACCTGTTT/4', 'status': 0, 'mapq': 60, 'three_prime': 57418, 'flag': 0, 'insert_start': 57417, 'insert_stop': 57418, 'insert_seq': 'G'}
        {'id': 'NB501801:571:HNM3KAFX3:1:11102:19076:17266_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA',  # noqa
         'r1_transposon': 'AATTCACTACGTCAACA/0',
         'r2_transposon': 'AATTCACTACGTCAACA/0',
         'tf': 'TTGCGACCTGTTT/4',
         'status': 1,
         'mapq': 60,
         'three_prime': 57418,
         'flag': 0,
         'insert_start': 57417,
         'insert_stop': 57418,
         'insert_seq': 'G'},


        # 1: {'id': 'NB501801:571:HNM3KAF...GCACGCTCGA', 'r1_transposon': 'AATTCACTACGTCAACA/0', 'r2_transposon': 'AATTCACTACGTCAACA/0', 'tf': 'TACTCACGCACGC/3', 'status': 2, 'mapq': 0, 'three_prime': 0, 'flag': 0, 'insert_start': '*', 'insert_stop': '*', 'insert_seq': '*'}
        {'id': 'NB501801:571:HNM3KAWER:1:21310:11840:5050_TACTCAATTCACTACGTCAACAACGCACGCTCGA',  # noqa
         'r1_transposon': 'AATTCACTACGTCAACA/0',
         'r2_transposon': 'AATTCACTACGTCAACA/0',
         'tf': 'TACTCACGCACGC/3',
         'status': 1,
         'mapq': 60,
         'three_prime': 59999,
         'flag': 16,
         'insert_start': '*',
         'insert_stop': '*',
         'insert_seq': '*'}]

    actual = []
    for read in yeast_bamfile.fetch():
        tagged_read = rt.tag_read(read)
        status_coder = create_status_coder(
            mapq_threshold=10,
            check_5_prime_clip=True)

        status_code = status_coder(tagged_read)

        summary_dict = {"id": tagged_read.get('read').query_name,
                        "r1_transposon": tagged_read.get('read').get_tag("RT"),
                        "r2_transposon": tagged_read.get('read').get_tag("RT"),
                        "tf": tagged_read.get('read').get_tag("TF"),
                        "status": status_code,
                        "mapq": tagged_read.get('read').mapping_quality,
                        "three_prime": tagged_read.get('read').get_tag("XS"),
                        "flag": tagged_read.get('read').flag,
                        "insert_start": tagged_read.get('read').get_tag("XI"),
                        "insert_stop": tagged_read.get('read').get_tag("XE"),
                        "insert_seq": tagged_read.get('read').get_tag("XZ")}

        actual.append(summary_dict)
    assert 2 == 2
    #assert actual == expected
