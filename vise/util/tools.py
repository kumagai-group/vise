# -*- coding: utf-8 -*-

from distutils.util import strtobool
from typing import Callable, Any
from xml.etree.ElementTree import ParseError

from vise.util.logger import get_logger

logger = get_logger(__name__)


def parse_file(class_method_name: Callable, parsed_filename: str) -> Any:
    """Check filename and parse and return cls via __init__ or class method.

    Args:
         class_method_name (Callable):
            Method to parse the given file. E.g., CLASS.from_file
        parsed_filename (str):
            Parsed file name.

    Return:
         Return of the class method.
    """
    try:
        logger.info(f"Parsing {parsed_filename}...")
        return class_method_name(parsed_filename)
    except ParseError:
        logger.warning(f"Parsing {parsed_filename} failed.")
        raise ParseError
    except FileNotFoundError:
        logger.warning(f"File {parsed_filename} not found.")
        raise


def str2bool(string: str) -> bool:
    return bool(strtobool(string))


def is_str_digit(n: str) -> bool:
    """Check whether the given string is a digit or not.

    Args:
        n (str): The checked string.

    Returns:
        Bool.
    """
    try:
        float(n)
        return True
    except ValueError:
        return False


def is_str_int(n: str) -> bool:
    """Check whether the given string is an integer or not.

    Args:
        n (str): The checked string.

    Returns:
        Bool.
    """
    try:
        if int(n) - float(n) < 1e-5:
            return True
        else:
            return False
    except ValueError:
        return False