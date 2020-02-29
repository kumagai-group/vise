# -*- coding: utf-8 -*-

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InvalidFileError(Exception):
    """Raised when the given file is invalid."""
    pass


class InvalidStructureError(Exception):
    """Raised when the structure is invalid."""
    pass


class NoVaspCommandError(Exception):
    """Raised when the vasp command is not set."""
    pass


class VaspNotConvergedError(Exception):
    """Raised when the vasp calculation is not converged."""
    pass


class KptNotConvergedError(Exception):
    """Raised when the k-point set is not converged."""
    pass
