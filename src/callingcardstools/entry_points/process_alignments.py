#pylint:disable=W1203
import logging
from logging.config import dictConfig
from re import U
import sys
import os
import argparse

from callingcardstools.utils import database_switcher
from callingcardstools import BarcodeParser
from callingcardstools import SummaryParser
from callingcardstools import tag_bam

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
    parser.add_argument("-a",
	                    "--alignment",
                        help="path to the input bam file",
						required=True)
    parser.add_argument("-g",
	                    "--genome",
                        help=" ".join(["Path to a genome .fasta file.",
						"Note that an index .fai file must exist in the same path"]),
                        required=True)
    parser.add_argument("-j",
		                "--barcode_details",
                        help = "Path to the barcode details json file",
						required=True)
    parser.add_argument("-i",
		                "--id_to_barcode_map",
                        help = "Path to the id_to_bc_map from split_fastq for a given run",
						required=True)
    parser.add_argument('-d',
	                    '--database', 
						help = 'path to the cc database', 
						required=True)
    parser.add_argument('-o', 
	                    '--organism', 
						help = 'One of either "yeast" or "mammal"',
						choices = ['yeast', 'mammal'],
						required=True)
    parser.add_argument('-t',
	                    '--tf', 
						help = 'Name of the tf for this alignment', 
						required=True)
    parser.add_argument('-n',
	                    '--replicate', 
						help = 'if the tf is replicated within the batch, then you must specify the replicate identifier', 
						required=True,
                        default='none')
    parser.add_argument('-r',
	                    '--regions_tblname', 
						help = 'Name of the regions table in the db', 
						required=True)
    parser.add_argument('-b',
	                    '--background_tblname', 
						help = 'Name of the background table in the db', 
						required=True)
    parser.add_argument("-q",
		                "--mapq_threshold",
                        help = "", 
						default = 10,
						type=int)
    parser.add_argument("-l",
	                    "--log_level",
                        help="the desired log level (default warning).",
                        choices=("critical", "error", "warning", "info", "debug"),
                        default="warning")

    return parser.parse_args(args)


def main(args=None):
    """_summary_

    Args:
        args (_type_, optional): _description_. Defaults to None.

    Raises:
        FileNotFoundError: _description_
    """
    args = parse_args(args)

    log_config = {
        "version": 1,
        "root": {
            "handlers": ["console"],
            "level": f"{args.log_level.upper()}"
        },
        "handlers": {
            "console": {
                "formatter": "std_out",
                "class": "logging.StreamHandler"
            }
        },
        "formatters": {
            "std_out": {
                "format": "%(asctime)s : %(module)s : %(funcName)s : line: %(lineno)d\n"+\
                    "process details : %(process)d, %(processName)s\n"+\
                        "thread details : %(thread)d, %(threadName)s\n"+\
                                "%(levelname)s : %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S"
            }
        }
    }
    dictConfig(log_config)

    # Check inputs
    input_path_list = [args.alignment,
                       args.genome,
                       args.genome + '.fai',
                       args.barcode_details,
					   args.id_to_barcode_map]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    barcode_details = BarcodeParser(args.barcode_details)

    # loop over the reads in the bam file and add the read group (header and tag)
    # and the XI and XZ tags
    aln_summary_df = tag_bam(
		args.alignment,
		args.genome,
		args.barcode_details,
		int(args.mapq_threshold))

    logging.info('connecting to database...')
    db = database_switcher(args.organism,args.database)
    
    logging.info(f'adding batch {barcode_details.barcode_dict["batch"]} data...') 
    db.add_batch_qc(
        barcode_details,
        args.id_to_barcode_map)
 
    sp = SummaryParser(aln_summary_df)
    qbed_df = sp.to_qbed()
    
    try:
        batch_id = db.get_batch_id(
            barcode_details.barcode_dict['batch'],
            args.tf,
            args.replicate)
    except ValueError as exc:
        logging.critical(f"get_batch_id returned more than 1 record: {exc}")
        raise

    qbed_df = qbed_df.assign(batch_id = [batch_id]*len(qbed_df))
    
    db.add_frame(
        qbed_df,
        'qbed',
        table_type='experiment',
        tablename_suffix = args.tf,
        drop = True,
        fk_tablelist=['batch'])
	
    experiment_tblname = 'experiment_'+args.tf

    if args.organism == 'yeast':
        db.create_aggregate_view('experiment_'+args.tf,args.regions_tblname)
	
    db.peak_caller(
		regions    = args.regions_tblname,
		background = args.background_tblname,
		experiment = experiment_tblname)

if __name__ == "__main__":
    sys.exit(main())