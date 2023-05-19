import pathlib
import pytest


@pytest.fixture
def yeast_summary_data():
    path = pathlib.Path("tests/test_data/yeast/"
                        "ERT1_novoalign_mapped_sorted_tagged_summary.csv")\
        .resolve()
    return path
