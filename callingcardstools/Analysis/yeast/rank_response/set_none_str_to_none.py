import logging
import re

logger = logging.getLogger(__name__)


def set_none_str_to_none(
        value: (str, type(None), int, float)) -> (str, type(None), int, float):
    """
    Test whether a string matches 'none' in any case. Return None if it does.
    Otherwise, return the original value

    Args:
        value (str, type(None), int, float): the string to test, or a value
            already set to None.

    Returns:
        str, type(None), int, float: the original value if it does not
            match 'none' in any case, or None if it does.

    Raises:
        TypeError: if the value is not a string, None or base python numeric.
    """
    if not isinstance(value, (str, type(None), int, float)):
        raise TypeError("value must be a string, None or base python numeric")

    none_pattern = r'(?i)^none$'

    # if the value is a string
    if isinstance(value, str):
        # test whether it matches 'none' in any case and return None if it does
        if bool(re.match(none_pattern, value)):
            return None
        # else, return the original value
        else:
            return value
    # if the value is not a string, it must be None, so return it
    return value