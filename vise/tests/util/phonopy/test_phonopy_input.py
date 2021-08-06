# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from phonopy.interface.calculator import read_crystal_structure
from phonopy.structure.atoms import PhonopyAtoms
from vise.util.phonopy.phonopy_input import structure_to_phonopy_atoms
import numpy as np


def assert_same_phonopy_atoms(actual: PhonopyAtoms,
                              expected: PhonopyAtoms):
    assert (actual.get_cell() == expected.get_cell()).all()
    assert (actual.get_scaled_positions()
            == expected.get_scaled_positions()).all()
    assert actual.symbols == expected.symbols


def test_phonopy_atoms_behavior(sc_structure, tmpdir):
    print(tmpdir)
    tmpdir.chdir()
#    actual = structure_to_phonopy_atoms(sc_structure)
    sc_structure.to(fmt="poscar", filename="POSCAR")
    a, _ = read_crystal_structure("POSCAR")
    b = PhonopyAtoms(atoms=a)
    print(type(a.get_cell()))
    print(a.get_atomic_numbers())
    assert_same_phonopy_atoms(a, b)


def test_structure_to_phonopy_atoms(sc_structure):
    actual = structure_to_phonopy_atoms(sc_structure)
    expected = PhonopyAtoms(symbols=["H"],
                            cell=np.array([[1.0, 0.0, 0.0],
                                           [0.0, 1.0, 0.0],
                                           [0.0, 0.0, 1.0]]),
                            scaled_positions=np.array([[0.0, 0.0, 0.0]]))
    assert_same_phonopy_atoms(actual, expected)

#
# def test_make_phonopy_input(mc_structure, mc_structure_conv):
#     actual = make_phonopy_input(unitcell=mc_structure,
#                                 supercell_matrix=np.eye(3).tolist(),
#                                 conventional_base=True)
#     supercell_matrix = [[ 1., 1., 0.],
#                         [-1., 1., 0.],
#                         [ 0., 0., 1.]]
#     supercell = mc_structure * supercell_matrix
#     expected = PhonopyInput(unitcell=mc_structure,
#                             supercell=supercell,
#                             supercell_matrix=supercell_matrix)
#     assert actual == expected
#
#
# def test_make_phonopy_input_default(mc_structure, mc_structure_conv):
#     actual = make_phonopy_input(unitcell=mc_structure)
#     supercell_matrix = [[ 2., 2., 0.],
#                         [-2., 2., 0.],
#                         [ 0., 0., 2.]]
#     supercell = mc_structure * supercell_matrix
#     expected = PhonopyInput(unitcell=mc_structure,
#                             supercell=supercell,
#                             supercell_matrix=supercell_matrix)
#     assert actual == expected
#
#
# def test_make_phonopy_input_default_hexa():
#     structure = Structure(Lattice.hexagonal(1.0, 2.0), species=["H"],
#                           coords=[[0.0]*3])
#     actual = make_phonopy_input(unitcell=structure)
#     supercell_matrix = [[2, -1, 0], [2, 1, 0], [0, 0, 2]]
#     supercell = structure * supercell_matrix
#     expected = PhonopyInput(unitcell=structure,
#                             supercell=supercell,
#                             supercell_matrix=supercell_matrix)
#     assert actual == expected
