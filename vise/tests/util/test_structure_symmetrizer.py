# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pymatgen.core import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core import Element
from vise.tests.helpers.assertion import assert_msonable

from vise.util.structure_symmetrizer import (
    cell_to_structure, StructureSymmetrizer, Site, num_symmetry_operation,
    first_structure_is_primitive)
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


def test_repr(symmetrizer_mc):
    expected = """Symprec: 0.01
Angle tolerance: 5.0
Space group: C2/m
Is primitive: False
site    wyckoff    site sym    equiv sites
------  ---------  ----------  -------------
H1      a          2/m         0 1"""
    assert symmetrizer_mc.__repr__() == expected


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


def test_spglib_cyclic_behavior():
    input_structure = Structure.from_str("""Ca4 Sc2 Sb2 O12
1.0
5.469740 0.000000 0.000000
0.000000 5.632621 0.000000
-5.467914 0.000000 7.830251
O
12
direct
0.344251 0.197893 0.550856 O
0.344251 0.302107 0.050856 O
0.655749 0.802107 0.449144 O
0.655749 0.697893 0.949144 O
0.752694 0.208674 0.946720 O
0.752694 0.291326 0.446720 O
0.247306 0.791326 0.053280 O
0.247306 0.708674 0.553280 O
0.144882 0.034671 0.743977 O
0.144882 0.465329 0.243977 O
0.855118 0.965329 0.256023 O
0.855118 0.534671 0.756023 O
""", fmt="POSCAR")
    ss = StructureSymmetrizer(input_structure)
    sss = StructureSymmetrizer(ss.primitive)
    assert ss.primitive == sss.primitive


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


def test_is_first_primitive():
    structure = Structure.from_str("""Hg6
1.0
4.4679886286 2.5795944374 4.4083109087
-4.4679886286 2.5795944374 4.4083109087
0.0000000000 -5.1591888749 4.4083109087
H Hg
1 6 
direct
0.0 0.0 0.0 H
0.2500000000 0.6035253347 0.8964746653 Hg
0.8964746653 0.2500000000 0.6035253347 Hg
0.6035253347 0.8964746653 0.2500000000 Hg
0.7500000000 0.3964746653 0.1035253347 Hg
0.1035253347 0.7500000000 0.3964746653 Hg
0.3964746653 0.1035253347 0.7500000000 Hg""", fmt="poscar")
    structure2 = Structure.from_str("""Hg6
1.0
4.4679886286 2.5795944374 4.4083109087
-4.4679886286 2.5795944374 4.4083109087
0.0000000000 -5.1591888749 4.4083109087
H Hg
1 6 
direct
0.0 0.0 0.0 H
0.750000 0.103525 0.396475 Hg
0.396475 0.750000 0.103525 Hg
0.103525 0.396475 0.750000 Hg
0.250000 0.896475 0.603525 Hg
0.603525 0.250000 0.896475 Hg
0.896475 0.603525 0.250000 Hg
""", fmt="poscar")
    assert first_structure_is_primitive(structure, structure2) is True
    assert first_structure_is_primitive(structure2, structure) is False
