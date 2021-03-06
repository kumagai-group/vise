# -*- coding: utf-8 -*-

from itertools import groupby
import operator

from pymatgen.core import Structure

from vise.util.logger import get_logger


logger = get_logger(__name__)


def create_symbol_list(structure: Structure):
    species = [str(s) for s in structure.species]
    # unique_justseen https://docs.python.org/ja/3/library/itertools.html
    # ["H", "H", "O", "O", "H"] -> ['H', 'O', 'H']
    return list(map(next, map(operator.itemgetter(1), groupby(species, None))))


