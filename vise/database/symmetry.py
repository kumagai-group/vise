# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import inv


def num_symmetry_operation(point_group: str) -> int:
    """ Return number of symmetry operations from Hermann–Mauguin notation.

    Args:
         point_group (str):
            Point group in Hermann–Mauguin notation.
            "." will be removed, e.g., ..6 -> 6.

    Returns:
         Number of symmetry operation.
    """
    d = {"1":     1,
         "-1":    2,
         "2":     2,
         "m":     2,
         "2/m":   4,
         "222":   4,
         "2mm":   4,
         "m2m":   4,
         "mm2":   4,
         "mmm":   8,
         "2/mmm": 8,
         "4":     4,
         "-4":    4,
         "4/m":   8,
         "422":   8,
         "4mm":   8,
         "-4m2":  8,
         "-42m":  8,
         "4/mmm": 16,
         "3":     3,
         "-3":    6,
         "32":    6,
         "3m":    6,
         "-3m":   12,
         "6":     6,
         "-6":    6,
         "6/m":   12,
         "622":   12,
         "6mm":   12,
         "-6m2":  12,
         "6/mmm": 24,
         "23":    12,
         "m3":    24,
         "432":   24,
         "-43m":  24,
         "m-3m":  48}

    # remove "." from the given point_group
    point_group = "".join([s for s in point_group if s != "."])

    return d[point_group]


def transmat_standard2primitive(centering: str) -> np.ndarray:
    """Transformation matrix from standardized cell to primitive cell

    Args:
        centering (str):
            Centering in one character.
    Return:
        Transformation matrix in numpy.ndarray.
     """

    if centering == "P":
        matrix = np.eye(3)
    elif centering == "A":
        matrix = np.array([[   1,    0,    0],
                           [   0,  1/2, -1/2],
                           [   0,  1/2,  1/2]])
    elif centering == "C":
        matrix = np.array([[ 1/2,  1/2,    0],
                           [-1/2,  1/2,    0],
                           [   0,    0,    1]])
    elif centering == "R":
        matrix = np.array([[ 2/3,  1/3,  1/3],
                           [-1/3,  1/3,  1/3],
                           [-1/3, -2/3,  1/3]])
#        matrix = np.array([[ 2/3, -1/3, -1/3],
#                           [ 1/3,  1/3, -2/3],
#                           [ 1/3,  1/3,  1/3]])
    elif centering == "I":
        matrix = np.array([[-1/2,  1/2,  1/2],
                           [ 1/2, -1/2,  1/2],
                           [ 1/2,  1/2, -1/2]])
    elif centering == "F":
        matrix = np.array([[   0,  1/2,  1/2],
                           [ 1/2,    0,  1/2],
                           [ 1/2,  1/2,    0]])
    else:
        raise ValueError(f"Centering {centering} is invalid")

    return matrix


def transmat_primitive2standard(centering: str) -> np.ndarray:
    """Transformation matrix from primitive cell to standardized cell

    Args:
        centering (str):
            Centering in one character.
    Return:
        Transformation matrix in numpy.ndarray.
    """
    matrix = inv(transmat_standard2primitive(centering))
    return matrix.astype(int)


