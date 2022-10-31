from .conftests import *
import pytest

from callingcardstools.SummaryParser import SummaryParser

def test_summary_parser_constructor(yeast_summary):

	sp = SummaryParser(yeast_summary)

	x = sp.summary

	df = sp.to_qbed('tf')