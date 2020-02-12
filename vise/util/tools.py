# -*- coding: utf-8 -*-
from typing import Callable, Any
from xml.etree.ElementTree import ParseError

from pydefect.util.tools import logger


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
