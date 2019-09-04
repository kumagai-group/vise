# -*- coding: utf-8 -*-

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InvalidFileError(Exception):
    """Raised when the given file is invalid."""
    pass


class VaspNotConvergedError(Exception):
    pass