# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
import numpy as np

from vise.util.structure_symmetrizer import (
    cell_to_structure, StructureSymmetrizer)
from pymatgen.core.structure import Structure


def test_cell_to_structure():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0]]
    cell = (lattice, coords, [1])

    expected = Structure(lattice=lattice, species="H", coords=coords)
    assert cell_to_structure(cell) == expected


@pytest.fixture()
def symmetrizer_mc():
    lattice =[[ 6.0,-3.0, 0.0],
              [ 6.0, 3.0, 0.0],
              [-4.0, 0.0, 7.0]]
    coords = [[0.0, 0.0, 0.0]]
    structure = Structure(lattice=lattice, species=["H"], coords=coords)
    symmetrizer = StructureSymmetrizer(structure)
    return symmetrizer


def test_structure_symmetrizer_mc_primitive(symmetrizer_mc):
    matrix = symmetrizer_mc.primitive.lattice.matrix
    expected = np.array([[ 6.0,-3.0, 0.0],
                         [ 6.0, 3.0, 0.0],
                         [-4.0, 0.0, 7.0]])
    np.testing.assert_array_almost_equal(matrix, expected)


def test_structure_symmetrizer_mc_band_primitive(symmetrizer_mc):
    band_matrix = symmetrizer_mc.band_primitive.lattice.matrix
    expected = np.array([[ 6.0, 3.0, 0.0],
                         [-6.0, 3.0, 0.0],
                         [-4.0, 0.0, 7.0]])
    np.testing.assert_array_almost_equal(band_matrix, expected)


def test_structure_symmetrizer_mc_lattice_change(symmetrizer_mc):
    assert symmetrizer_mc.is_band_primitive_lattice_changed is True
    assert symmetrizer_mc.is_primitive_lattice_changed is False
    assert symmetrizer_mc.band_primitive_differ_primitive is True


@pytest.fixture()
def symmetrizer_bcc():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    structure = Structure(lattice=lattice, species=["H", "H"],
                          coords=coords)
    symmetrizer = StructureSymmetrizer(structure)
    return symmetrizer


def test_structure_symmetrizer_bcc_primitive(symmetrizer_bcc):
    matrix = symmetrizer_bcc.primitive.lattice.matrix
    expected = np.array([[-0.5,  0.5,  0.5],
                         [ 0.5, -0.5,  0.5],
                         [ 0.5,  0.5, -0.5]])
    np.testing.assert_array_almost_equal(matrix, expected)


def test_structure_symmetrizer_bcc_lattice_change(symmetrizer_bcc):
    assert symmetrizer_bcc.is_band_primitive_lattice_changed is True
    assert symmetrizer_bcc.is_primitive_lattice_changed is True
    assert symmetrizer_bcc.band_primitive_differ_primitive is False


