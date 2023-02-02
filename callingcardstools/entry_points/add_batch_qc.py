#pylint:disable=w1201,w1203
import os
import sys
import argparse
import logging

from callingcardstools.utils import database_switcher

def parse_args(args=None)->dict:
    """parse command line arguments for tagbam

    args:
        args (_type_, optional): _description_. defaults to none.

    returns:
        _type_: _description_
    """
    description = "add barcode information from reads to the database"
    epilog = ""

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('-o', 
	                    '--organism', 
						help = 'either yeast or mammal', 
						required=True)
    parser.add_argument('-d',
	                    '--database', 
						help = 'path to the cc database', 
						required=True)
    parser.add_argument('-i',
	                    '--id_to_barcode_map',
						type=str, 
						nargs='+',
						help = 'path to list of id_to_barcode_map',
						required=True)
    parser.add_argument('-j',
	                    '--barcode_details',
						help = "path to the barcode_details json for this run",
						required=True)
    parser.add_argument('-r',
	                    '--batch',
						help = "name of the batch, eg the run 1234, for the database entry",
						required=True)
    parser.add_argument("-l",
	                    "--log_level",
                        help="the desired log level (default warning).",
                        choices=("critical", "error", "warning", "info", "debug"),
                        default="warning")

    return parser.parse_args(args)

def main(args=None)->None:
	"""_summary_

	args:
		args (_type_, optional): _description_. defaults to none.

	raises:
		filenotfounderror: _description_
	"""
	log_config = {
		"version":1,
		"root":{
			"handlers" : ["console"],
			"level": f"{args.log_level}"
		},
		"handlers":{
			"console":{
				"formatter": "std_out",
				"class": "logging.streamhandler",
				"level": "info"
			}
		},
		"formatters":{
			"std_out": {
				"format": "%(asctime)s : %(levelname)s : %(module)s : %(funcname)s : %(lineno)d : (process details : (%(process)d, %(processname)s), thread details : (%(thread)d, %(threadname)s))\nlog : %(message)s",
				"datefmt":"%d-%m-%y %i:%m:%s"
			}
		},
	}
	logging.config.dictConfig(log_config)
    
	# check inputs
	logging.info('checking input...')
	input_path_list = [args.database, args.barcode_details,args.id_to_barcode_map]
	for input_path in input_path_list:
		if not os.path.exists(input_path):
			raise FileNotFoundError(f"input file dne: {input_path}")
    
	logging.info('connecting to database...')
	db = database_switcher(args.organism,args.database)
	
	logging.info(f'adding batch {args.batch} data...') 
	db.add_batch_qc(args.batch,args.id_to_barcode_map,args.barcode_details)	

if __name__ == "__main__":
    sys.exit(main())
