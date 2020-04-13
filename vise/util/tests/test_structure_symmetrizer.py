# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from vise.util.structure_symmetrizer import cell_to_structure
from pymatgen.core.structure import Structure


def test_cell_to_structure():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0]]
    cell = (lattice, coords, [1])

    expected = Structure(lattice=lattice, species="H", coords=coords)
    assert cell_to_structure(cell) == expected


