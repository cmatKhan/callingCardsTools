from typing import Literal
import logging
from difflib import SequenceMatcher
from pysqlite3.dbapi2 import Connection
from callingcardstools.database_managers.yeast import HopsDb as yeast_db

__all__ = ['remove_suffix', 'convert_logger_level', 'database_switcher']

def database_switcher(organism:Literal['yeast','mammal'], db_path:str) -> Connection:
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

def convert_logger_level(level:str)-> int:
    """Convert string logger level, eg 'info' to corresponding int
	Cite: Christian Tremblay https://github.com/ChristianTremblay/BAC0

    Args:
        level (str): one of info, debug, warning, error, critical

    Raises:
        ValueError: if the level is not recognized

    Returns:
        int: the integer corresponding to the input level, eg 'info' returns 20
    """
    if not level:
        return None
    _valid_levels = [
        logging.DEBUG,
        logging.INFO,
        logging.WARNING,
        logging.ERROR,
        logging.CRITICAL,
    ]
    if level in _valid_levels:
        return level
    if level.lower() == "info":
        return logging.INFO
    elif level.lower() == "debug":
        return logging.DEBUG
    elif level.lower() == "warning":
        return logging.WARNING
    elif level.lower() == "error":
        return logging.ERROR
    elif level.lower() == "critical":
        return logging.CRITICAL
    raise ValueError(f"Wrong log level use one of the following : {_valid_levels}")

def remove_suffix(input_string: str, suffix: str) -> str:
	"""a mimic of python 3.9 removesuffix

	Args:
		input_string (str): _description_
		suffix (str): _description_

	Returns:
		str: _description_
	"""
	if suffix and input_string.endswith(suffix):
		return input_string[:-len(suffix)]
	return input_string

def get_best_match(query, corpus, step=4, flex=3, case_sensitive=False, verbose=False):
    """Return best matching substring of corpus.

    cite: https://stackoverflow.com/a/36132391

    Parameters
    ----------
    query : str
    corpus : str
    step : int
        Step size of first match-value scan through corpus. Can be thought of
        as a sort of "scan resolution". Should not exceed length of query.
    flex : int
        Max. left/right substring position adjustment value. Should not
        exceed length of query / 2.

    Outputs
    -------
    output0 : str
        Best matching substring.
    output1 : float
        Match ratio of best matching substring. 1 is perfect match.
    """

    def _match(a, b):
        """Compact alias for SequenceMatcher."""
        return SequenceMatcher(None, a, b).ratio()

    def scan_corpus(step):
        """Return list of match values from corpus-wide scan."""
        match_values = []

        m = 0
        while m + qlen - step <= len(corpus):
            match_values.append(_match(query, corpus[m : m-1+qlen]))
            if verbose:
                print(query, "-", corpus[m: m + qlen], _match(query, corpus[m: m + qlen]))
            m += step

        return match_values

    def index_max(v):
        """Return index of max value."""
        return max(range(len(v)), key=v.__getitem__)

    def adjust_left_right_positions():
        """Return left/right positions for best string match."""
        # bp_* is synonym for 'Best Position Left/Right' and are adjusted 
        # to optimize bmv_*
        p_l, bp_l = [pos] * 2
        p_r, bp_r = [pos + qlen] * 2

        # bmv_* are declared here in case they are untouched in optimization
        bmv_l = match_values[p_l // step]
        bmv_r = match_values[p_l // step]

        for f in range(flex):
            ll = _match(query, corpus[p_l - f: p_r])
            if ll > bmv_l:
                bmv_l = ll
                bp_l = p_l - f

            lr = _match(query, corpus[p_l + f: p_r])
            if lr > bmv_l:
                bmv_l = lr
                bp_l = p_l + f

            rl = _match(query, corpus[p_l: p_r - f])
            if rl > bmv_r:
                bmv_r = rl
                bp_r = p_r - f

            rr = _match(query, corpus[p_l: p_r + f])
            if rr > bmv_r:
                bmv_r = rr
                bp_r = p_r + f

            if verbose:
                print("\n" + str(f))
                print("ll: -- value: %f -- snippet: %s" % (ll, corpus[p_l - f: p_r]))
                print("lr: -- value: %f -- snippet: %s" % (lr, corpus[p_l + f: p_r]))
                print("rl: -- value: %f -- snippet: %s" % (rl, corpus[p_l: p_r - f]))
                print("rr: -- value: %f -- snippet: %s" % (rl, corpus[p_l: p_r + f]))

        return bp_l, bp_r, _match(query, corpus[bp_l : bp_r])

    if not case_sensitive:
        query = query.lower()
        corpus = corpus.lower()

    qlen = len(query)

    if flex >= qlen/2:
        print("Warning: flex exceeds length of query / 2. Setting to default.")
        flex = 3

    match_values = scan_corpus(step)
    pos = index_max(match_values) * step

    pos_left, pos_right, match_value = adjust_left_right_positions()

    return corpus[pos_left: pos_right].strip(), match_value