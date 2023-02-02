#pylint:disable=w1201,w1203
from typing import literal
import os
import sys
import argparse
import logging

from pysqlite3.dbapi2 import connection

from callingcardstools.database_managers.yeast import hopsdb as yeast_db

def database_manager(organism:literal['yeast','mammal'], db_path:str) -> connection:
	"""_summary_

	args:
		organism (literal[&#39;yeast&#39;,&#39;mammal&#39;]): _description_
		db_path (str): _description_

	raises:
		filenotfounderror: _description_

	returns:
		connection: _description_
	"""
    
	switcher = {
		'yeast': yeast_db(db_path)
	}

	out = switcher.get(organism)

	return out

def parse_args(args=none)->dict:
    """parse command line arguments for tagbam

    args:
        args (_type_, optional): _description_. defaults to none.

    returns:
        _type_: _description_
    """
    description = "add barcode information from reads to the database"
    epilog = ""

    parser = argparse.argumentparser(description=description, epilog=epilog)

    parser.add_argument('-o', 
	                    '--organism', 
						help = 'either yeast or mammal', 
						required=true)
    parser.add_argument('-d',
	                    '--database', 
						help = 'path to the cc database', 
						required=true)
    parser.add_argument('-i',
	                    '--id_to_barcode_map',
						type=str, 
						nargs='+',
						help = 'path to list of id_to_barcode_map',
						required=true)
    parser.add_argument('-j',
	                    '--barcode_details',
						help = "path to the barcode_details json for this run",
						required=true)
    parser.add_argument('-r',
	                    '--batch',
						help = "name of the batch, eg the run 1234, for the database entry",
						required=true)
    parser.add_argument("-l",
	                    "--log_level",
                        help="the desired log level (default warning).",
                        choices=("critical", "error", "warning", "info", "debug"),
                        default="warning")

    return parser.parse_args(args)

def main(args=none)->none:
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
	input_path_list = [args.database, 
	                   args.hops_path]
	for input_path in input_path_list:
		if not os.path.exists(input_path):
			raise filenotfounderror("input file dne: %s" %input_path)
    
	logging.info('connecting to database...')
	db = database_manager(args.organism,args.database)
	
	logging.info(f'adding hops {args.hops} data...') 
	yeast_db.add_frame(
		df,
		'qbed',
		table_type='experiment',
		tablename_suffix = args.tablename_suffix, 
		fk_tablelist=['batch'])

if __name__ == "__main__":
    sys.exit(main())
