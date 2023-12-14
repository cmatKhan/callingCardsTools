import argparse
import json
import logging

from .create_rank_response_table import create_rank_response_table
from .validate_config import validate_config

logger = logging.getLogger(__name__)

__all__ = ['parse_args', 'main']


def parse_args(
        subparser: argparse.ArgumentParser,
        script_desc: str,
        common_args: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Parse the command line arguments.

    :param subparser: the subparser object.
    :type subparser: argparse.ArgumentParser
    :param script_desc: the description of the script.
    :type script_desc: str
    :param common_args: the common arguments.
    :type common_args: argparse.ArgumentParser
    :return: the parser.
    :rtype: argparse.ArgumentParser
    """

    parser = subparser.add_parser(
        'yeast_rank_response',
        help=script_desc,
        prog='yeast_rank_response',
        parents=[common_args]
    )

    parser.set_defaults(func=main)

    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help="Path to the configuration json file. "
        "For details, see "
        "https://cmatkhan.github.io/callingCardsTools/file_format_specs/yeast_rank_response/"  # noqa
    )

    return subparser


def main(args: argparse.Namespace):
    # Load the JSON configuration file
    with open(args.config, 'r', encoding='utf-8') as config_file:
        config_dict = json.load(config_file)

    try:
        config_dict = validate_config(config_dict)
    except (KeyError, TypeError, FileNotFoundError) as exc:
        logger.error("Error in configuration file: %s", exc)
        raise

    rank_response_df = create_rank_response_table(config_dict)

    compression = 'gzip' if config_dict.get('compress', False) else None
    if compression and not config_dict['output_file'].endswith('.gz'):
        config_dict['output_file'] += '.gz'
    rank_response_df.to_csv(config_dict['output_file'],
                            compression=compression,
                            index=False)
