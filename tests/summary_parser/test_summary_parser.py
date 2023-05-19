from .conftests import *
import pytest

from callingcardstools.Alignment.SummaryParser import SummaryParser


def test_summary_parser_constructor(yeast_summary_data):

    sp = SummaryParser(yeast_summary_data.as_posix())

    x = sp.summary

    assert 2 == 2

    #df = sp.to_qbed('tf')
