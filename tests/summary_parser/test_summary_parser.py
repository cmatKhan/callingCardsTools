from .conftests import *
import pytest

from callingcardstools.SummaryParser import SummaryParser

def test_summary_parser_constructor(yeast_summary):

	sp = SummaryParser(yeast_summary)

	x = sp.summary

	df = sp.to_qbed('tf')


# def test_hops_db_tmp():

# 	from io import StringIO
# 	import glob
# 	import os
# 	import re
# 	# the peak callers all inherit from DatabaseApi, a class which offers an interface 
# 	# to a sqlite database to store Calling Cards data
# 	from callingcardstools.database_managers.yeast import HopsDb
# 	from callingcardstools.PackageResources import Resources
# 	import pandas as pd

# 	# This object allows retrieval of package resources
# 	cc_resources = Resources()

# 	# create a database either in memory or at a specified location
# 	#yeast_db = hopsdb("/home/oguzkhan/Desktop/cc_metadata/hops_db.sqlite")
# 	yeast_db = HopsDb("/home/oguzkhan/projects/rank_response_shiny/data/qc_db_v2.sqlite")

# 	regions_list = ['regions_yiming']
# 	background_list = ['background_adh1']
# 	experiment_list = ['experiment_GZF3']

# 	for region_tbl in regions_list:
# 		print(f"region: {region_tbl}")
# 		for background_tbl in background_list:
# 			print(f"background: {background_tbl}")
# 			for experiment_tbl in experiment_list:
# 				print(f"experiment: {experiment_tbl}")
# 				yeast_db.peak_caller(regions = region_tbl,
# 										background = background_tbl, 
# 										experiment = experiment_tbl,
# 										if_exists='replace')