# -*- coding: utf-8 -*-

from itertools import groupby
import operator
from typing import List

import numpy as np
from pymatgen.core import Structure

from vise.util.logger import get_logger


logger = get_logger(__name__)


def create_symbol_list(structure: Structure):
    species = [str(s) for s in structure.species]
    # unique_justseen https://docs.python.org/ja/3/library/itertools.html
    # ["H", "H", "O", "O", "H"] -> ['H', 'O', 'H']
    return list(map(next, map(operator.itemgetter(1), groupby(species, None))))


def sanitize_matrix(matrix: List[int]) -> List[List[int]]:

    if len(matrix) == 9:
        return [matrix[:3], matrix[3:6], matrix[6:]]
    elif len(matrix) == 3:
        result = np.eye(3, dtype=int)
        for i in range(3):
            result[i, i] = matrix[i]
        return result.tolist()
    elif len(matrix) == 1:
        result = np.eye(3, dtype=int)
        for i in range(3):
            result[i, i] = matrix[0]
        return result.tolist()
    else:
        raise ValueError(f"Matrix element length {len(matrix)} is improper.")