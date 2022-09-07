from .conftests import *
from callingcardstools.bam_parsers import StatusFlags

def check_tags(read, expected, tag_list):
    """check tags against an expected dictionary

    Args:
        read (AlignedSegment): A pysam AlignedSegment object
        expected (dict): A dictionary where read.query_name are keys pointing 
        to dicts of expected tags, eg {RG: asdf, XI 1234, ...}
        tag_list (list): A list of tags to check 
    """
    for tag in tag_list:
        assert read.get_tag(tag) == expected[read.query_name][tag]

def test_human_constructor(human_readtagger):
    assert human_readtagger.is_open() == True

def test_yeast_constructor(yeast_readtagger):
    assert yeast_readtagger.is_open() == True


def test_tag_read(yeast_readtagger, 
                  yeast_bamfile):
    
    expected = {
        "NB501801:571:HNM3KAFX3:2:21105:24619:7841_TCAGTAATTCACTACGTCAACACCCGTTGGTCGA": 
        {"RG": "TCAGTAATTCACTACGTCAACACCCGTTGGTCGA", "XS": 32765, "XI":32766, "XE":32767, "XZ": "G"},

        "NB501801:571:HNM3KAFX3:1:11312:1357:3858_TTGGTAATTCACTACGTCAACAGGGGGTAGTCGN": 
        {"RG": "TTGGTAATTCACTACGTCAACAGGGGGTAGTCGN", "XS":32991, "XI": 32990, "XE": 32991, "XZ": "G"},

        "NB501801:571:HNM3KAFX3:4:21612:4799:3141_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA": 
        {"RG": "TTGCGAATTCACTACGTCAACAACCTGTTTTCGA", "XS": 56908, "XI": 56909, "XE": 56910, "XZ": "T"},

        "NB501801:571:HNM3KAFX3:1:11102:19076:17266_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA": 
        {"RG": "TTGCGAATTCACTACGTCAACAACCTGTTTTCGA", "XS": 57418, "XI":57417, "XE": 57418, "XZ": "G"},

        "NB501801:571:HNM3KAFX3:2:21206:2348:1221_GCCTGNATTCACTACGTCAACAGGCGGCAGCCNG": 
        {"RG": "GCCTGNATTCACTACGTCAACAGGCGGCAGCCNG", "XS":32780, "XI": 32781, "XE":32782, "XZ": "A"},

        "NB501801:571:HNM3KAFX3:3:11501:15818:8474_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA": 
        {"RG": "TTGCGAATTCACTACGTCAACAACCTGTTTTCGA", "XS":4644, "XI": 4643, "XE":4644, "XZ": "A"},

        "NB501801:571:HNM3KAFX3:1:21310:11840:5050_TACTCAATTCACTACGTCAACAACGCACGCTCGA": 
        {"RG": "TACTCAATTCACTACGTCAACAACGCACGCTCGA", "XS": 0, "XI": '*', "XE": '*', "XZ": '*'},

        "NB501801:571:HNM3KAWER:1:21310:11840:5050_TACTCAATTCACTACGTCAACAACGCACGCTCGA": 
        {"RG": "TACTCAATTCACTACGTCAACAACGCACGCTCGA", "XS": 59999, "XI": "*", "XE": "*", "XZ": "*"},

        "NB501801:571:HNM3KAFX3:4:21409:11423:19883_GAATCAATTCACTACGTCAACATTTTTCTATCGA":
        {"RG": "GAATCAATTCACTACGTCAACATTTTTCTATCGA", "XS": '*', "XI": '*', "XE": '*', "XZ": '*'}}
    
    tag_list = ['RG', 'XS', 'XI', 'XE', 'XZ']

    rt = yeast_readtagger
    for read in yeast_bamfile.fetch():
        tagged_read = rt.tag_read(read)
        check_tags(tagged_read, expected, tag_list)
        

def test_read_tagger_yeast(yeast_readtagger, yeast_bamfile, yeast_bp):

    rt = yeast_readtagger

    bp = yeast_bp

    mapq_threshold = 10

    expected_tf_dict = {
        'NB501801:571:HNM3KAFX3:4:21409:11423:19883_GAATCAATTCACTACGTCAACATTTTTCTATCGA': 'MIG2_9',
        'NB501801:571:HNM3KAFX3:1:21310:11840:5050_TACTCAATTCACTACGTCAACAACGCACGCTCGA': 'MIG2_7',
        'NB501801:571:HNM3KAFX3:3:11501:15818:8474_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA': 'MIG2_7',
        'NB501801:571:HNM3KAFX3:2:21105:24619:7841_TCAGTAATTCACTACGTCAACACCCGTTGGTCGA': 'MIG2',
        'NB501801:571:HNM3KAFX3:2:21206:2348:1221_GCCTGNATTCACTACGTCAACAGGCGGCAGCCNG': "CAT8",
        'NB501801:571:HNM3KAFX3:1:11312:1357:3858_TTGGTAATTCACTACGTCAACAGGGGGTAGTCGN': "ARO80",
        'NB501801:571:HNM3KAFX3:4:21612:4799:3141_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA': 'MIG2_7',
        'NB501801:571:HNM3KAFX3:1:11102:19076:17266_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA': 'MIG2_7',
        'NB501801:571:HNM3KAWER:1:21310:11840:5050_TACTCAATTCACTACGTCAACAACGCACGCTCGA': 'MIG2_7'
    }

    expected = [
        {'id': 'NB501801:571:HNM3KAFX3:4:21409:11423:19883_GAATCAATTCACTACGTCAACATTTTTCTATCGA', 
        'bc': 'GAATCAATTCACTACGTCAACATTTTTCTATCGA', 
        'status': 3, 
        'mapq': 0, 
        'three_prime': '*', 
        'flag': 0, 
        'insert_start': '*', 
        'insert_stop': '*', 
        'insert_seq': '*'}, 
        {'id': 'NB501801:571:HNM3KAFX3:1:21310:11840:5050_TACTCAATTCACTACGTCAACAACGCACGCTCGA', 
        'bc': 'TACTCAATTCACTACGTCAACAACGCACGCTCGA', 
        'status': 3, 
        'mapq': 0, 
        'three_prime': 0, 
        'flag': 0, 
        'insert_start': '*', 
        'insert_stop': '*', 
        'insert_seq': '*'}, 
        {'id': 'NB501801:571:HNM3KAFX3:3:11501:15818:8474_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA', 
        'bc': 'TTGCGAATTCACTACGTCAACAACCTGTTTTCGA', 
        'status': 3, 
        'mapq': 5, 
        'three_prime': 4644, 
        'flag': 0, 
        'insert_start': 4643, 
        'insert_stop': 4644, 
        'insert_seq': 'A'},
        {'id': 'NB501801:571:HNM3KAFX3:2:21105:24619:7841_TCAGTAATTCACTACGTCAACACCCGTTGGTCGA', 
        'bc': 'TCAGTAATTCACTACGTCAACACCCGTTGGTCGA', 
        'status': 0, 
        'mapq': 60, 
        'three_prime': 32765, 
        'flag': 16, 
        'insert_start': 32766, 
        'insert_stop': 32767, 
        'insert_seq': 'G'}, 
        {'id': 'NB501801:571:HNM3KAFX3:2:21206:2348:1221_GCCTGNATTCACTACGTCAACAGGCGGCAGCCNG', 
        'bc': 'GCCTGNATTCACTACGTCAACAGGCGGCAGCCNG', 
        'status': 1,
        'mapq': 60, 
        'three_prime': 32780, 
        'flag': 16, 
        'insert_start': 32781, 
        'insert_stop': 32782, 
        'insert_seq': 'A'}, 
        {'id': 'NB501801:571:HNM3KAFX3:1:11312:1357:3858_TTGGTAATTCACTACGTCAACAGGGGGTAGTCGN', 
        'bc': 'TTGGTAATTCACTACGTCAACAGGGGGTAGTCGN', 
        'status': 0, 
        'mapq': 60, 
        'three_prime': 32991, 
        'flag': 0, 
        'insert_start': 32990, 
        'insert_stop': 32991, 
        'insert_seq': 'G'},
         {'id': 'NB501801:571:HNM3KAFX3:4:21612:4799:3141_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA', 
         'bc': 'TTGCGAATTCACTACGTCAACAACCTGTTTTCGA', 
         'status': 1, 
         'mapq': 60, 
         'three_prime': 56908, 
         'flag': 16, 
         'insert_start': 56909, 
         'insert_stop': 56910, 
         'insert_seq': 'T'}, 
         {'id': 'NB501801:571:HNM3KAFX3:1:11102:19076:17266_TTGCGAATTCACTACGTCAACAACCTGTTTTCGA', 
         'bc': 'TTGCGAATTCACTACGTCAACAACCTGTTTTCGA', 
         'status': 1, 
         'mapq': 60, 
         'three_prime': 57418, 
         'flag': 0, 
         'insert_start': 57417, 
         'insert_stop': 57418, 
         'insert_seq': 'G'}, 
         {'id': 'NB501801:571:HNM3KAWER:1:21310:11840:5050_TACTCAATTCACTACGTCAACAACGCACGCTCGA', 
         'bc': 'TACTCAATTCACTACGTCAACAACGCACGCTCGA', 
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
        bp.set_barcode(tagged_read.get_tag("RG"))
        barcode_details = bp.barcode_check()
        tagged_read.set_tag("XF", barcode_details['tf'])

        # check the TF tagging
        assert tagged_read.get_tag("XF") == expected_tf_dict[tagged_read.query_name]


        status_code = 0
        if not barcode_details['pass']:
            status_code += StatusFlags.BARCODE.flag()
        if tagged_read.mapping_quality < mapq_threshold:
            status_code +=  StatusFlags.MAPQ.flag()
        try:
            if not bp.get_insert_seqs() == ["*"]:
                if tagged_read.get_tag("XZ") not in bp.get_insert_seqs():
                    status_code += StatusFlags.INSERT_SEQ.flag()
        except AttributeError as exc:
            print(f"insert sequence not found in Barcode Parser. {exc}")

        actual.append({"id": tagged_read.query_name,
                        "bc": tagged_read.get_tag("RG"),
                        "status": status_code,
                        "mapq": tagged_read.mapping_quality,
                        "three_prime": tagged_read.get_tag("XS"),
                        "flag": tagged_read.flag,
                        "insert_start": tagged_read.get_tag("XI"),
                        "insert_stop": tagged_read.get_tag("XE"),
                        "insert_seq": tagged_read.get_tag("XZ")})
    
    assert actual == expected

def test_read_tagger_human(human_readtagger, human_bamfile, human_bp):
    rt = human_readtagger

    bp = human_bp

    mapq_threshold = 10

    expected = [
        {'id': 'MN00200:583:000H3LYJT:1:13109:24273:8517_TAGCGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA', 
        'bc': 'TAGCGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA', 
        'status': 0, 
        'mapq': 60, 
        'three_prime': 30195, 
        'flag': 0, 
        'insert_start': 30191, 
        'insert_stop': 30195, 
        'insert_seq': 'TTAA'}, 
        {'id': 'MN00200:583:000H3LYJT:1:12110:14874:16649_TAGCGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA', 
        'bc': 'TAGCGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA', 
        'status': 4, 
        'mapq': 50, 
        'three_prime': 32249, 
        'flag': 0, 
        'insert_start': 32245, 
        'insert_stop': 32249, 
        'insert_seq': 'GGGT'}, 
        {'id': 'MN00200:583:000H3LYJT:1:21108:15538:6862_AAGCGTCAATTTTACGCAGACTATCTTTTTTTGGTTAA', 
        'bc': 'AAGCGTCAATTTTACGCAGACTATCTTTTTTTGGTTAA', 
        'status': 0, 
        'mapq': 30, 
        'three_prime': 32307, 
        'flag': 0, 
        'insert_start': 32303, 
        'insert_stop': 32307, 
        'insert_seq': 'TTAA'}, 
        {'id': 'MN00200:583:000H3LYJT:1:21105:10495:3666_ACGCGTCAATTTTACGCAGACTATCTTTGGTTGGTTAA', 
        'bc': 'ACGCGTCAATTTTACGCAGACTATCTTTGGTTGGTTAA', 
        'status': 1, 
        'mapq': 10, 
        'three_prime': 34621, 
        'flag': 0, 
        'insert_start': 34617, 
        'insert_stop': 34621, 
        'insert_seq': 'TTAA'}, 
        {'id': 'MN00200:583:000H3LYJT:1:11109:8966:8780_TAGCGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA', 
        'bc': 'TAGCGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA', 
        'status': 0, 
        'mapq': 40, 
        'three_prime': 35234, 
        'flag': 16, 
        'insert_start': 35235, 
        'insert_stop': 35239, 
        'insert_seq': 'TTAA'}, 
        {'id': 'MN00200:583:000H3LYJT:1:12103:17924:13875_TANCGTCAATTTTACGCAGACTATCTTTNTTTGGTTAN', 
        'bc': 'TANCGTCAATTTTACGCAGACTATCTTTNTTTGGTTAN', 
        'status': 1, 
        'mapq': 20, 
        'three_prime': 64714, 
        'flag': 16, 
        'insert_start': 64715, 
        'insert_stop': 64719, 
        'insert_seq': 'TTAA'}]
    
    actual = []
    for read in human_bamfile.fetch():
        tagged_read = rt.tag_read(read)
        bp.set_barcode(tagged_read.get_tag("RG"))
        barcode_details = bp.barcode_check()
        tagged_read.set_tag("XF", barcode_details['tf'])

        # check the TF tagging
        assert tagged_read.get_tag("XF") == "*"


        status_code = 0
        if not barcode_details['pass']:
            status_code += StatusFlags.BARCODE.flag()
        if tagged_read.mapping_quality < mapq_threshold:
            status_code += StatusFlags.MAPQ.flag()
        try:
            if not bp.get_insert_seqs() == ["*"]:
                if tagged_read.get_tag("XZ") not in bp.get_insert_seqs():
                    status_code += StatusFlags.INSERT_SEQ.flag()
        except AttributeError as exc:
            print(f"insert sequence not found in Barcode Parser. {exc}")

        actual.append({"id": tagged_read.query_name,
                        "bc": tagged_read.get_tag("RG"),
                        "status": status_code,
                        "mapq": tagged_read.mapping_quality,
                        "three_prime": tagged_read.get_tag("XS"),
                        "flag": tagged_read.flag,
                        "insert_start": tagged_read.get_tag("XI"),
                        "insert_stop": tagged_read.get_tag("XE"),
                        "insert_seq": tagged_read.get_tag("XZ")})
    
    assert actual == expected


    