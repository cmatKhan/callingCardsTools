import json
import os
import argparse

import pandas as pd

from callingcardstools.Resources import Resources


def parse_args(
        subparser: argparse.ArgumentParser,
        script_desc: str,
        common_args: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """This is intended to be used as a subparser for a parent parser passed 
    from __main__.py. This function adds the arguments for this script to the
    subparser and returns the subparser to be used in __main__.py

    Args:
        subparser (argparse.ArgumentParser): See __main__.py -- this is the 
        subparser for the parent parser in __main__.py
        script_desc (str): Description of this script, which is set in 
        __main__.py. The description is set in __main__.py so that all of 
        the script descriptions are together in one spot and it is easier to 
        write a unified cmd line interface
        common_args (argparse.ArgumentParser): These are the common arguments 
        for all scripts in callingCardsTools, for instance logging level

    Returns:
        argparse.ArgumentParser: The subparser with the this additional 
        cmd line tool added to it -- intended to be gathered in __main__.py 
        to create a unified cmd line interface for the package
    """

    parser = subparser.add_parser(
        'barcode_table_to_json',
        help=script_desc,
        prog='barcode_table_to_json',
        parents=[common_args]
    )

    parser.set_defaults(func=main)

    parser.add_argument('-t',
                        '--barcode_table',
                        help='old pipeline barcode table',
                        required=True)
    parser.add_argument('-r',
                        '--batch',
                        help='batch name, eg the run number like run_1234',
                        required=True)

    return subparser


def main(args: argparse.Namespace) -> None:
    """This function takes a barcode table from the old pipeline and converts
    it to a barcode_details.json file for use in the new pipeline

    Args:
        args (argparse.Namespace): The cmd line arguments passed from
        __main__.py

    Raises:
        FileNotFoundError: If the barcode table does not exist
        TypeError: If the barcode_details.json file does not exist in the
        yeast PackageResources
    """

    if not os.path.exists(args.barcode_table):
        raise FileNotFoundError(f'{args.barcode_table} Does Not Exist!')

    try:
        barcode_dict = json.loads(Resources().yeast.get('barcode_details'))
    except TypeError as exc:
        raise TypeError('barcode_details.json not found in '
                       'yeast PackageResources. Please post an issue at '
                       'https://github.com/cmatKhan/callingCardsTools/issues') from exc  # noqa

    df = pd.read_csv(args.barcode_table, sep="\t",
                     names=['tf', 'r1', 'r2'])

    tf_map = {x[1]: x[0] for x in
              df.assign(bc=lambda x: x['r1']+x['r2'])[['tf', 'bc']]
              .to_dict('tight')['data']}

    barcode_dict['components']['tf']['map'] = tf_map
    barcode_dict['batch'] = args.batch
    with open(f"{args.batch}_barcode_details.json", 'w', encoding='utf-8') as f2: # noqa
        json_object = json.dumps(barcode_dict, indent=4)
        f2.write(json_object)
