# -*- coding: utf-8 -*-

from distutils.util import strtobool
from typing import Callable, Any
from xml.etree.ElementTree import ParseError

from vise.util.logger import get_logger

logger = get_logger(__name__)


def str2bool(val: str) -> bool:
    """Convert a string representation of truth to true (1) or false (0).
    """
    val = val.lower()
    if val in ('t', 'true'):
        return True
    elif val in ('false'):
        return False
    else:
        raise ValueError("invalid truth value %r" % (val,))


def is_str_digit(n: str) -> bool:
    """Check whether the given string is a digit or not. """
    try:
        float(n)
        return True
    except ValueError:
        return False


def is_str_int(n: str, rounding_error=1e-7) -> bool:
    """Check whether the given string is an integer or not. """
    try:
        if int(n) - float(n) < rounding_error:
            return True
        else:
            return False
    except ValueError:
        return False
