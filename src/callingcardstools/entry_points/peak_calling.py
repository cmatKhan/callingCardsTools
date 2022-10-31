
#pylint:disable=W1203
import logging
import sys
import os
import argparse

from callingcardstools.utils import database_switcher

def parse_args(args=None):
    """parse command line arguments for yeast peak calling

    Args:
        args (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    Description = "parse alignment summary into qbed file"
    Epilog = ""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument('-o', 
	                    '--organism', 
						help = 'either yeast or mammal', 
						required=True)
    parser.add_argument('-d',
	                    '--database', 
						help = 'path to the cc database', 
						required=True)
    parser.add_argument('-r',
	                    '--regions_tbl', 
						help = 'Name of the regions table in the db', 
						required=True)
    parser.add_argument('-b',
	                    '--background_tbl', 
						help = 'Name of the background table in the db', 
						required=True)
    parser.add_argument('-e',
	                    '--experiment_tbl', 
						help = 'Name of the experiment table', 
						required=True)
    parser.add_argument("-l",
	                    "--log_level",
                        help="the desired log level (default warning).",
                        choices=("critical", "error", "warning", "info", "debug"),
                        default="warning")

    return parser.parse_args(args)


def main(args=None):

    args = parse_args(args)

    log_config = {
        "version": 1,
        "root": {
            "handlers": ["console"],
            "level": f"{args.log_level}"
        },
        "handlers": {
            "console": {
                "formatter": "std_out",
                "class": "logging.streamhandler",
                "level": "info"
            }
        },
        "formatters": {
            "std_out": {
                "format": "%(asctime)s : %(levelname)s : %(module)s : "+\
                    "%(funcname)s : %(lineno)d : (process details : "\
                        "(%(process)d, %(processname)s), thread details : "\
                            "(%(thread)d, %(threadname)s))\nlog : %(message)s",
                "datefmt": "%d-%m-%y %i:%m:%s"
            }
        },
    }
    logging.config.dictConfig(log_config)

    # Check inputs
    logging.info('checking input...')
    input_path_list = [args.database]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" % input_path)
 
    logging.info('connecting to database...')
    db = database_switcher(args.organism,args.database)
    
    logging.info(f'adding batch {args.batch} data...')
    
	# wonder if I could just pass **args to generalize?
    db.peak_caller(
		regions    = args.region_tbl,
		background = args.background_tbl, 
		experiment = args.experiment_tbl)	
    
if __name__ == '__main__':
    sys.exit(main())