#pylint:disable=W1203
import logging
import sys
import os
import argparse

from callingcardstools.utils import database_switcher
from callingcardstools import SummaryParser

def parse_args(args=None):
    """parse command line arguments for tagBam

    Args:
        args (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    Description = "parse alignment summary into qbed file"
    Epilog = ""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument('-a',
                        '--aln_summary',
                        help='Path to the alignment summary tsv',
                        required=True)
    parser.add_argument('-b',
                        '--barcode_details',
                        help='Path to the barcode_details json',
                        required=True)
    parser.add_argument('regions_tbl')

    return parser.parse_args(args)


def main(args=None):

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
    logging.config.dictConfig(log_config)

    # Check inputs
    logging.info('checking input...')
    input_path_list = [args.aln_summary,
                       args.barcode_details]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" % input_path)
 
    logging.info('connecting to database...')
    db = database_switcher(args.organism,args.database)
    
    logging.info(f'adding batch {args.batch} data...') 
    db.add_batch_qc(args.batch,args.id_to_barcode_map,args.barcode_details)	
    
    sp = SummaryParser(args.aln_summary)
    df = sp.to_qbed()
    
    try:
        batch_id = db.get_batch_id(args.batch,args.tf,args.replicate)
    except ValueError as exc:
        logging.critical(f"get_batch_id returned more than 1 record: {exc}")
        raise

    df[['batch_id']] = [batch_id]*len(df)
    
    db.add_frame(
        df,
        'qbed',
        table_type='experiment',
        tablename_suffix = args.tf, 
        fk_tablelist=['batch'])

    if args.organism == 'yeast':
        db.create_aggregate_view('experiment_'+args.tf,args.regions_tbl)

if __name__ == '__main__':
    sys.exit(main())
