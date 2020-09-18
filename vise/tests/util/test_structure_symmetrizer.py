# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen import Element
from vise.tests.conftest import assert_msonable

from vise.util.structure_symmetrizer import (
    cell_to_structure, StructureSymmetrizer, Site, num_symmetry_operation,
    sum_frac_coords)
from vise.util.bravais_lattice import BravaisLattice


def test_cell_to_structure():
    lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    coords = [[0.0, 0.0, 0.0]]
    cell = (lattice, coords, [1])

    expected = Structure(lattice=lattice, species="H", coords=coords)
    assert cell_to_structure(cell) == expected


@pytest.fixture()
def symmetrizer_mc():
    lattice = [[10.0,  0.0,  0.0],
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


def test_grouped_atom_indices(complex_monoclinic_structure):
    symmetrizer = StructureSymmetrizer(complex_monoclinic_structure)
    actual = symmetrizer.grouped_atom_indices()
    assert actual == {'H1_a': [0], 'He1_m': [1, 2], 'He2_m': [3, 4]}


def test_sites(complex_monoclinic_structure):
    symmetrizer = StructureSymmetrizer(complex_monoclinic_structure)
    actual = symmetrizer.sites
    expected = {"H1": Site("H", "a", "2/m", [0]),
                "He1": Site("He", "m", "m", [1, 2]),
                "He2": Site("He", "m", "m", [3, 4])}
    assert actual == expected


def test_sites_for_bio2():
    bio2_structure = Structure.from_str("""Bi4 O8
    1.0
    4.880207 0.000000 -2.019673
             -0.667412 5.402283 -1.612690
                       -0.040184 -0.016599 6.808846
    Bi O
    4 8
    direct
    0.729106 0.250000 1.000000 Bi
    0.270894 0.750000 0.000000 Bi
    0.000000 0.000000 0.500000 Bi
    0.500000 0.500000 0.500000 Bi
    0.321325 0.021684 0.814642 O
    0.114894 0.650123 0.351314 O
    0.885106 0.349877 0.648686 O
    0.763580 0.849877 0.648686 O
    0.236420 0.150123 0.351314 O
    0.493318 0.521684 0.814642 O
    0.506682 0.478316 0.185358 O
    0.678675 0.978316 0.185358 O""", fmt="POSCAR")
    symmetrizer = StructureSymmetrizer(bio2_structure)
    actual = symmetrizer.sites
    expected = \
        {'Bi1': Site(element='Bi', wyckoff_letter='e', site_symmetry='2', equivalent_atoms=[0, 1]),
         'Bi2': Site(element='Bi', wyckoff_letter='c', site_symmetry='-1', equivalent_atoms=[2, 3]),
         'O1': Site(element='O', wyckoff_letter='f', site_symmetry='1', equivalent_atoms=[4, 9, 10, 11]),
         'O2': Site(element='O', wyckoff_letter='f', site_symmetry='1', equivalent_atoms=[5, 6, 7, 8])}
    assert actual == expected


def test_():
    structure = Structure.from_str(""" Na2 Cu2 O4
1.0
1.8088488050 -5.4059198794 0.0000000000
1.8088488050 5.4059198794 0.0000000000
0.0000000000 0.0000000000 5.3190514901
Na Cu O
2 2 4
direct
0.7123738820 0.2876261180 0.2500000000 Na
0.2876261180 0.7123738820 0.7500000000 Na
0.9970663194 0.0029336806 0.2500000000 Cu
0.0029336806 0.9970663194 0.7500000000 Cu
0.1138815169 0.8861184831 0.4939795157 O
0.8861184831 0.1138815169 0.5060204843 O
0.1138815169 0.8861184831 0.0060204843 O
0.8861184831 0.1138815169 0.9939795157 O 
""", fmt="POSCAR")
    symmetrizer = StructureSymmetrizer(structure)
    print(symmetrizer.primitive)
    sstructure = Structure.from_str(""" Na24 Cu24 O48
1.0
10.853093 0.000000 0.000000
0.000000 10.811840 0.000000
0.000000 0.000000 10.638103
Na Cu O
24 24 48
direct
0.166667 -0.212374 0.125000 Na
0.166667 -0.212374 0.625000 Na
0.500000 -0.212374 0.125000 Na
0.500000 -0.212374 0.625000 Na
0.833333 -0.212374 0.125000 Na
0.833333 -0.212374 0.625000 Na
0.333333 0.287626 0.125000 Na
0.333333 0.287626 0.625000 Na
0.666667 0.287626 0.125000 Na
0.666667 0.287626 0.625000 Na
1.000000 0.287626 0.125000 Na
1.000000 0.287626 0.625000 Na
0.166667 0.212374 0.375000 Na
0.166667 0.212374 0.875000 Na
0.500000 0.212374 0.375000 Na
0.500000 0.212374 0.875000 Na
0.833333 0.212374 0.375000 Na
0.833333 0.212374 0.875000 Na
0.333333 0.712374 0.375000 Na
0.333333 0.712374 0.875000 Na
0.666667 0.712374 0.375000 Na
0.666667 0.712374 0.875000 Na
1.000000 0.712374 0.375000 Na
1.000000 0.712374 0.875000 Na
0.166667 -0.497066 0.125000 Cu
0.166667 -0.497066 0.625000 Cu
0.500000 -0.497066 0.125000 Cu
0.500000 -0.497066 0.625000 Cu
0.833333 -0.497066 0.125000 Cu
0.833333 -0.497066 0.625000 Cu
0.333333 0.002934 0.125000 Cu
0.333333 0.002934 0.625000 Cu
0.666667 0.002934 0.125000 Cu
0.666667 0.002934 0.625000 Cu
1.000000 0.002934 0.125000 Cu
1.000000 0.002934 0.625000 Cu
0.166667 0.497066 0.375000 Cu
0.166667 0.497066 0.875000 Cu
0.500000 0.497066 0.375000 Cu
0.500000 0.497066 0.875000 Cu
0.833333 0.497066 0.375000 Cu
0.833333 0.497066 0.875000 Cu
0.333333 0.997066 0.375000 Cu
0.333333 0.997066 0.875000 Cu
0.666667 0.997066 0.375000 Cu
0.666667 0.997066 0.875000 Cu
1.000000 0.997066 0.375000 Cu
1.000000 0.997066 0.875000 Cu
0.166667 0.386118 0.246990 O
0.166667 0.386118 0.746990 O
0.500000 0.386118 0.246990 O
0.500000 0.386118 0.746990 O
0.833333 0.386118 0.246990 O
0.833333 0.386118 0.746990 O
0.333333 0.886118 0.246990 O
0.333333 0.886118 0.746990 O
0.666667 0.886118 0.246990 O
0.666667 0.886118 0.746990 O
1.000000 0.886118 0.246990 O
1.000000 0.886118 0.746990 O
0.166667 -0.386118 0.253010 O
0.166667 -0.386118 0.753010 O
0.500000 -0.386118 0.253010 O
0.500000 -0.386118 0.753010 O
0.833333 -0.386118 0.253010 O
0.833333 -0.386118 0.753010 O
0.333333 0.113882 0.253010 O
0.333333 0.113882 0.753010 O
0.666667 0.113882 0.253010 O
0.666667 0.113882 0.753010 O
1.000000 0.113882 0.253010 O
1.000000 0.113882 0.753010 O
0.166667 0.386118 0.003010 O
0.166667 0.386118 0.503010 O
0.500000 0.386118 0.003010 O
0.500000 0.386118 0.503010 O
0.833333 0.386118 0.003010 O
0.833333 0.386118 0.503010 O
0.333333 0.886118 0.003010 O
0.333333 0.886118 0.503010 O
0.666667 0.886118 0.003010 O
0.666667 0.886118 0.503010 O
1.000000 0.886118 0.003010 O
1.000000 0.886118 0.503010 O
0.166667 -0.386118 0.496990 O
0.166667 -0.386118 0.996990 O
0.500000 -0.386118 0.496990 O
0.500000 -0.386118 0.996990 O
0.833333 -0.386118 0.496990 O
0.833333 -0.386118 0.996990 O
0.333333 0.113882 0.496990 O
0.333333 0.113882 0.996990 O
0.666667 0.113882 0.496990 O
0.666667 0.113882 0.996990 O
1.000000 0.113882 0.496990 O
1.000000 0.113882 0.996990 O""", fmt="POSCAR")
    symmetrizer = StructureSymmetrizer(sstructure)
    print(symmetrizer.primitive)



def test_bravais_lattice():
    lattice =[[10.0,  0.0,  0.0],
              [ 0.0, 10.0,  0.0],
              [-2.0,  0.0, 10.0]]
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]
    structure = Structure(lattice=lattice, species=["H"] * 2, coords=coords)
    symmetrizer = StructureSymmetrizer(structure)
    assert symmetrizer.bravais == BravaisLattice.mC
    assert symmetrizer.centering == "C"


def test_species_order():
    lattice = Lattice.orthorhombic(5, 6, 7)
    coords = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],

        [0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.5],
    ]
    structure = Structure(lattice=lattice,
                          species=["H"] * 4 + ["He"] * 4,
                          coords=coords)
    supercell = structure * [1, 1, 2]
    actual = [e.specie for e in StructureSymmetrizer(supercell).conventional]
    expected = [Element.H] * 4 + [Element.He] * 4
    assert actual == expected


# def test_spglib_cyclic_behavior():
#     input_structure = Structure.from_str("""Ca4 Sc2 Sb2 O12
# 1.0
# 5.469740 0.000000 0.000000
# 0.000000 5.632621 0.000000
# -5.467914 0.000000 7.830251
# O
# 12
# direct
# 0.344251 0.197893 0.550856 O
# 0.344251 0.302107 0.050856 O
# 0.655749 0.802107 0.449144 O
# 0.655749 0.697893 0.949144 O
# 0.752694 0.208674 0.946720 O
# 0.752694 0.291326 0.446720 O
# 0.247306 0.791326 0.053280 O
# 0.247306 0.708674 0.553280 O
# 0.144882 0.034671 0.743977 O
# 0.144882 0.465329 0.243977 O
# 0.855118 0.965329 0.256023 O
# 0.855118 0.534671 0.756023 O
# """, fmt="POSCAR")

    # ss = StructureSymmetrizer(input_structure)
    # assert ss.primitive == input_structure


def test_sum_frac_coords():
    structure = Structure.from_str(""" 
1.0
1 0 0
0 1 0
0 0 1
H
2
direct
0.1 0.2 0.3 H
0.4 0.5 0.6 H""", fmt="POSCAR")
    assert sum_frac_coords(structure) == 2.1


@pytest.fixture
def site():
    return Site(element="H",
                wyckoff_letter="a",
                site_symmetry="m3m",
                equivalent_atoms=[0, 1, 2, 3, 5, 6])


def test_site_msonable(site):
    assert_msonable(site)


def test_site_pprint_equiv_atoms(site):
    assert site.pprint_equiv_atoms == "0..3 5 6"


def test_num_symmetry_operation():
    assert num_symmetry_operation("m-3m") == 48
    assert num_symmetry_operation(".m.") == 2



"""
TODO: 
- Add grouped_atom_indices to structure symmetrizer.

"""



