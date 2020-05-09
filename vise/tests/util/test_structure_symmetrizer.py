# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pymatgen.core.structure import Structure

from vise.util.structure_symmetrizer import (
    cell_to_structure, StructureSymmetrizer)


def test_cell_to_structure():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0]]
    cell = (lattice, coords, [1])

    expected = Structure(lattice=lattice, species="H", coords=coords)
    assert cell_to_structure(cell) == expected


@pytest.fixture()
def symmetrizer_mc():
    lattice =[[10.0,  0.0,  0.0],
              [ 0.0, 10.0,  0.0],
              [-2.0,  0.0, 10.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]
    structure = Structure(lattice=lattice, species=["H"] * 2, coords=coords)
    symmetrizer = StructureSymmetrizer(structure)
    return symmetrizer


def test_structure_symmetrizer_mc_primitive(symmetrizer_mc):
    matrix = symmetrizer_mc.primitive.lattice.matrix
    expected = np.array([[5., -5., 0.],
                         [5., 5., 0.],
                         [-2., 0., 10.]])
    np.testing.assert_array_almost_equal(matrix, expected)


def test_structure_symmetrizer_mc_conventional(symmetrizer_mc):
    matrix = symmetrizer_mc.conventional.lattice.matrix
    expected = np.array([[10.,  0.,  0.],
                         [ 0., 10.,  0.],
                         [-2.,  0., 10.]])
    np.testing.assert_array_almost_equal(matrix, expected)


def test_structure_symmetrizer_mc_band_primitive(symmetrizer_mc):
    band_matrix = symmetrizer_mc.band_primitive.lattice.matrix
    expected = np.array([[5., 5., 0.],
                         [-5., 5., 0.],
                         [-2., 0., 10.]])
    np.testing.assert_array_almost_equal(band_matrix, expected)


def test_structure_symmetrizer_mc_irreducible_kpoints(symmetrizer_mc):
    num_kpt_list = [2, 2, 2]
    kpt_shift = [0.5, 0.5, 0.5]
    actual = symmetrizer_mc.irreducible_kpoints(num_kpt_list, kpt_shift)
    expected = [([0.25, 0.25, 0.25], 4), ([0.75, 0.25, 0.25], 4)]
    assert actual == expected


def test_structure_symmetrizer_mc_lattice_change(symmetrizer_mc):
    assert symmetrizer_mc.is_band_primitive_lattice_changed is True
    assert symmetrizer_mc.is_primitive_lattice_changed is True
    assert symmetrizer_mc.band_primitive_differ_primitive is True


def test_sg_number(symmetrizer_mc):
    assert symmetrizer_mc.sg_number == 12


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


def test_structure_symmetrizer_bcc_irreducible_kpoints(symmetrizer_bcc):
    num_kpt_list = [2, 2, 2]
    kpt_shift = [0.0, 0.0, 0.0]
    actual = symmetrizer_bcc.irreducible_kpoints(num_kpt_list, kpt_shift)
    expected = [([0.0, 0.0, 0.0], 1), ([0.5, 0.0, 0.0], 3),
                ([0.5, 0.5, 0.0], 3), ([0.5, 0.5, 0.5], 1)]
    assert actual == expected


def test_structure_symmetrizer_sc_irreducible_kpoints():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0]]
    structure = Structure(lattice=lattice, species=["H"], coords=coords)
    symmetrizer = StructureSymmetrizer(structure)

    num_kpt_list = [2, 2, 2]
    kpt_shift = [0.5, 0.5, 0.5]
    actual = symmetrizer.irreducible_kpoints(num_kpt_list, kpt_shift)
    expected = [([0.25, 0.25, 0.25], 8)]
    assert actual == expected


def test_grouped_atom_indices(complex_ortho_structure):
    symmetrizer = StructureSymmetrizer(complex_ortho_structure)
    actual = symmetrizer.grouped_atom_indices()
    assert actual == {'H_a1': [0], 'He_m1': [1, 2], 'He_m2': [3, 4]}

"""
TODO: 
- Add grouped_atom_indices to structure symmetrizer.

"""



