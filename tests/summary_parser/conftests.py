from ..conftests import *
from callingcardstools.Database.yeast import HopsDb as yeast_HopsDb

@pytest.fixture
def yeast_summary():
	return "tests/test_data/yeast/ERT1_novoalign_mapped_sorted_tagged_summary.csv"