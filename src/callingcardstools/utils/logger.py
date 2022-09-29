import sys
import logging

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