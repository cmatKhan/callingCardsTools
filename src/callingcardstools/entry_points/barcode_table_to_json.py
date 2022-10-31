import json
import os
import sys
import argparse

import pandas as pd

def parse_args(args=None):
    """parse command line arguments for tagBam

    Args:
        args (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    Description = ""
    Epilog = ""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument('-b', 
	                    '--barcode_json', 
						help = 'path to a valid barcode dict json', 
						required=True)
    parser.add_argument('-t',
	                    '--barcode_table', 
						help = 'old pipeline barcode table', 
						required=True)
    parser.add_argument('-r',
	                    '--batch',
						help = 'batch name, eg the run number like run_1234',
						required=True)

    return parser.parse_args(args)

def main(args=None):
	args = parse_args(args)
	
	if not os.path.exists(args.barcode_table):
		raise FileNotFoundError(f'{args.barcode_table} DNE')
	
	if not os.path.exists(args.barcode_json):
		raise FileNotFoundError(f'{args.barcode_json} DNE')
	
	# open json, read in as dict
	with open(args.barcode_json, 'r') as f1: #pylint:disable=W1514
		barcode_dict = json.load(f1)

		df = pd.read_csv(args.barcode_table, sep = "\t", names = ['tf', 'r1', 'r2'])

		tf_map = {x[1]:x[0] for x in \
			df.assign(bc = lambda x: x['r1']+x['r2'])[['tf','bc']]\
				.to_dict('tight')['data']} 

		barcode_dict['components']['tf']['map'] = tf_map
		barcode_dict['batch'] = args.batch
		with open(f"{args.run_number}_barcode_details.json", 'w') as f2:
			json_object = json.dumps(barcode_dict, indent=4)
			f2.write(json_object)


if __name__ == "__main__":
    sys.exit(main())