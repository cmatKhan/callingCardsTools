import pathlib
import os
import re
import pytest
from ..conftests import *


@pytest.fixture
def yeast_barcodeqccounter_data():
    """paths to r1 r2 and barcode details which identify it as a correct read

    Returns:
        tuple: r1, r2, barcode_details
    """
    r1 = pathlib\
        .Path("tests/Reads/data/yeast/r1.fq")\
        .resolve()
    r2 = pathlib\
        .Path("tests/Reads/data/yeast/r2.fq")\
        .resolve()
    barcode_details = pathlib\
        .Path("tests/Reads/data/yeast/run_6177_barcode_details.json")\
        .resolve()
    return r1, r2, barcode_details


@pytest.fixture
def yeast_reads():
    """paths relative to project root to yeast r1 and r2 fastq files

    Returns:
        dict: key, value where keys are r1, r2 and values are paths to fastq files
    """
    return {
        'r1': "tests/test_data/yeast/r1.fastq",
        'r2': "tests/test_data/yeast/r1.fastq",
        'r1_gz': "tests/test_data/yeast/R1.fq.gz",
        'r2_gz': "tests/test_data/yeast/R2.fq.gz"}


@pytest.fixture
def yeast_split_qc_files():
    """
    collect the paths of the qc files in data/yeast/split_fastq
    and return them as a dict
    """
    x = [pathlib.Path("tests/Reads/data/yeast/split_qc/" + x).resolve()
         for x in os.listdir("tests/Reads/data/yeast/split_qc")]

    return {'split': [path for path in x if re.search(r".pickle$", str(path))],
            'check': [path for path in x if re.search(r".csv$", str(path))]}
