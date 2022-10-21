from ..conftests import *

@pytest.fixture
def yeast_reads():
    """paths relative to project root to yeast r1 and r2 fastq files

    Returns:
        dict: key, value where keys are r1, r2 and values are paths to fastq files
    """
    return {
        'r1':"tests/test_data/yeast/r1.fastq",
        'r2':"tests/test_data/yeast/r1.fastq",
        'r1_gz':"tests/test_data/yeast/R1.fq.gz", 
        'r2_gz':"tests/test_data/yeast/R2.fq.gz"}