# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest
from pymatgen.core import Lattice, IStructure
from vise.util.centering import Centering
from vise.util.structure_symmetrizer import StructureSymmetrizer


def test_names():
    assert Centering.names_string() == "P, A, C, R, I, F"


@pytest.fixture(scope="session")
def rhombohedral():
    lattice = Lattice.rhombohedral(a=1, alpha=45)
    coords = [[0.0, 0.0, 0.0]]
    return IStructure(lattice=lattice, species=["H"], coords=coords)


@pytest.fixture(scope="session")
def a_centered_orthorhombic():
    lattice = Lattice([[1,  0, 0],
                       [0,  2, 3],
                       [0, -2, 3]])
    coords = [[0.5, 0.8, 0.8],
              [0.0, 0.3, 0.0],
              [0.0, 0.0, 0.3]]

    return IStructure(lattice=lattice, species=["H"] * 3, coords=coords)


@pytest.fixture(scope="session")
def c_centered_monoclinic():
    lattice = Lattice.monoclinic(3, 4, 5, 100)
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]
    return IStructure(lattice=lattice, species=["H", "H"], coords=coords)


@pytest.fixture(scope="session")
def bcc():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    return IStructure(lattice=lattice, species=["H"] * 2, coords=coords)


def test_a_centering(a_centered_orthorhombic):
    s = StructureSymmetrizer(a_centered_orthorhombic)
    primitive = s.primitive
    to_conventional = primitive * Centering.A.primitive_to_conv
    assert to_conventional.lattice == s.conventional.lattice
    assert Centering.A.conv_multiplicity == 2


def test_c_centering(c_centered_monoclinic):
    s = StructureSymmetrizer(c_centered_monoclinic)
    primitive = s.primitive
    to_conventional = primitive * Centering.C.primitive_to_conv
    assert to_conventional == s.conventional
    assert Centering.C.conv_multiplicity == 2


def test_r_centering(rhombohedral):
    s = StructureSymmetrizer(rhombohedral)
    rhombohedral = s.primitive
    to_conventional = rhombohedral * Centering.R.primitive_to_conv
    assert to_conventional == s.conventional
#    assert Centering.R.conv_multiplicity == 3


def test_f_centering(bcc):
    s = StructureSymmetrizer(bcc)
    primitive = s.primitive
    to_conventional = primitive * Centering.I.primitive_to_conv
    assert to_conventional == s.conventional
    assert Centering.I.conv_multiplicity == 2


"""
TODO
- 

DONE
"""