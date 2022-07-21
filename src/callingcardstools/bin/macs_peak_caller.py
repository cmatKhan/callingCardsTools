#!/usr/bin/env python

"""This module roughly follows the algorithm used by MACS to call
ChIP-Seq peaks, but applies it to calling card data.  The main peak calling
function is passed a experiment frame, a background frame, and an TTAA_frame, 
all in ccf format.  It then builds interval trees containing
all of the background and experiment hops and all of the TTAAs.  Next, 
it scans through the genome with a window of window_size and step size of step_size 
and looks for regions that have significantly more experiment hops than background 
hops (poisson w/ pvalue_cutoff).  It merges consecutively enriched windows and 
computes the center of the peak.  Next it computes lambda, the number 
of insertions per TTAA expected from the background distribution by taking the max of 
lambda_bg, lamda_1, lamda_5, lambda_10.  It then computes a p-value based 
on the expected number of hops = lamda * number of TTAAs in peak * number of hops
in peak.  Finally, it returns a frame that has Chr,Start,End,Center,experiment_hops,
fraction_experiment,Background Hops,Fraction Background,poisson_pvalue

There is also a background free version of this algorithm that computes significance 
based on the number of hops in the neighboring region.

"""

#stdlib
import argparse
import os
import sys
# local import
from callingcardstools.macs_style import with_background
from callingcardstools.macs_style import background_free
#third party
import pandas as pd

# TODO figure out how to set this in config; make part of user input
CCF_COLNAMES = ["Chr","Start","End","Reads","Strand","Barcode"]
TTAA_COLNAMES = ["Chr","Start","End"]

def parse_args(args=None):
	"""parse command line arguments

	NOTE that any 'file' argument should have the word 'file' in the long 
	  name of the argument. Anything with 'file' in the long name will be 
	  checked for existence in main.

	Args:
		args (str, optional): Command line arguments. Defaults to None.

	Returns:
		dict: dictionary of command line arguments + default arguments
	"""
	Description = "A MACS mimic for calling regions significantly enriched for "+\
		"calling cards transpositions. If a background file is provided, "+\
			"it will be used to calculate the statistics. If a background "+\
				"file is not provided, then a background free algorithm "+\
					"will be used instead."
	Epilog = "write me!" # TODO write me

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument('-e', '--experiment_ccf_file', 
						required=True,
						help='experiment CCF file')
	parser.add_argument('-t', '--ttaa_file', 
						required=True,
						help='TTAA file')
	parser.add_argument('-a', '--annotation_file',
						required=True,
						help='annotation filename (full path)')
	parser.add_argument('-b', '--background_ccf_file', 
						required=False, 
						default=False,
						help='background filename (full path)')
	parser.add_argument('-pc', '--peak_pval_cutoff', 
						required=False, 
						default=1e-4,
						help='pvalue cutoff for significant peaks')
	parser.add_argument('--window_size',
						required=False,
						default=1000,
						help='window size')
	parser.add_argument('--lam_win_size',
						required=False,
						default=100000,
						help='window size for lambda computation in "+\
							"background free peak-caller')
	parser.add_argument('--step_size',
						required=False,
						default=500,
						help='step size',)
	# TODO check where pseudocount is added and update help
	parser.add_argument('--pseudocount',
						required=False,
						default=0.2,
						help='pseudocount to add to poisson(?) pvalue calcuations')

	return parser.parse_args(args)

def main(args=None):
	args = parse_args(args)

	# Any input with 'file' in the long name will be checked 
	# for existence raise FileNotFoundError if DNE.
	input_path_list = [v for k,v in args.__dict__.items() if "file" in k]
	for input_path in input_path_list:
		if not os.path.exists(input_path):
			raise FileNotFoundError("Input file DNE: %s" %input_path)

	experiment_ccf_df = pd.read_csv(args.experiment_ccf_file,
	                                sep="\t",
						            names = CCF_COLNAMES)
	# TODO check field expectations
	background_ccf_df = pd.read_csv(args.background_ccf_file,
	                                sep="\t",
						            names = CCF_COLNAMES)
	# TODO check field expectations

	ttaa_df = pd.read_csv(args.ttaa_file,
	                      sep="\t",
						  names = TTAA_COLNAMES)
	# TODO check field expectations

	output_filename = \
		os.path.splitext(args.experiment_ccf_file)[0] + "_peaks.tsv"

	output_dict = {}
	
	if args.background_ccf_file:
		try:
			with_background.call_peaks_and_annotate(
				experiment_ccf_df,
				background_ccf_df,
				ttaa_df,
				output_filename,
				args.annotation_file,
				float(args.peak_pval_cutoff),
				int(args.window_size),
				int(args.step_size),
				float(args.pseudocount))

		except:
			raise
	else:
		try:
			background_free.call_peaks_and_annotate(
				args.exp_file,
				args.output_file,
				args.TTAA_file,
				args.annotation_file,
				float(args.pvaluecutoff),
				float(args.peak_finder_pvalue),
				int(args.window),
				int(args.lam_win_size),
				int(args.step),
				float(args.pseudocounts))
		except:
			raise
	
	return 0

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))

