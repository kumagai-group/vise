# -*- coding: utf-8 -*-

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InvalidFileError(Exception):
    """Raised when the given file is invalid."""
    pass


class InvalidStructureError(Exception):
    """Raised when the structure is invalid."""
    pass


class VaspNotConvergedError(Exception):
    """Raised when the vasp calculation is not converged."""
    pass
