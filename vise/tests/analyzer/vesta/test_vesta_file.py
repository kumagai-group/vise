# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest

from pathlib import Path

from pymatgen.core import Lattice, DummySpecies
from pymatgen.core.structure import Structure

from vise.analyzer.vesta.vesta_file import VestaFile, Title, Cellp, \
    Struc, Bound, SBond, Vect, Style, ImportDensity, add_density, calc_isurfs

parent_dir = Path(__file__).parent


@pytest.fixture
def structure():
    return Structure(Lattice.cubic(3.9),
                     species=["Ba", DummySpecies()],
                     coords=[[0.5]*3, [0.0]*3])


def test_title(structure):
    actual = repr(Title(structure))
    expected = '''#VESTA_FORMAT_VERSION 3.5.0
TITLE
X1 Ba1'''
    assert actual == expected


def test_cellp(structure):
    actual = repr(Cellp(structure))
    expected = '''CELLP
3.900000 3.900000 3.900000 90.000000 90.000000 90.000000'''
    assert actual == expected


def test_struc(structure):
    actual = repr(Struc(structure))
    expected = '''STRUC
1 Ba Ba1 1.0 0.500000 0.500000 0.500000
 0.0 0.0 0.0 
2 XX XX2 1.0 0.000000 0.000000 0.000000
 0.0 0.0 0.0 
 0 0 0 0 0 '''
    assert actual == expected


def test_import_density():
    actual = repr(ImportDensity(volumetric_filename="CHG"))
    expected = """IMPORT_DENSITY 0.1
+1.000000 CHG"""
    assert actual == expected


def test_bound():
    actual = repr(Bound((0, 2, 0, 2, 0, 2)))
    expected = '''BOUND
0.000000 2.000000 0.000000 2.000000 0.000000 2.000000
 0 0 0 0 0 '''
    assert actual == expected


def test_sbond():
    s = Structure(Lattice.cubic(3.9), species=["Ba", "Ti"], coords=[[0.5]*3, [0.0]*3])
    actual = repr(SBond(s, bond_factor=1.2))
    expected = '''SBOND
1 Ba Ti 0.0  2.81  0  0  1  0  1
2 Ti Ba 0.0  2.81  0  0  1  0  1
 0 0 0 0 '''
    assert actual == expected


def test_vect():
    actual = repr(Vect({2: [0., 0., 0.], 3: [0., 0., -0.1]}))
    expected = '''VECTR
1 0.000000 0.000000 0.000000
2  0 0 0 0
 0 0 0 0 0 
2 0.000000 0.000000 -0.100000
3  0 0 0 0
 0 0 0 0 0 
 0 0 0 0 0 

VECTT
1 0.5 0 0 0 2
2 0.5 0 0 0 2
 0 0 0 0 0 '''
    assert actual == expected


def test_style(structure):
    actual = repr(Style(0.1, is_ionic=True))
    expected = '''STYLE
VECTS  1.0
SECTS  64  1
SECTP 
   1 0 0  0.00000E+00  0.00000E+00  -1.00000E+01  1.00000E+01
UCOLP 
 0  2  1.000   0   0   0
ATOMS  1  0  1
BONDP 
   1  16  0.1  1.000 127 127 127'''
    assert actual == expected


def test_vesta_file_to(structure, tmpdir):
    tmpdir.chdir()
    VestaFile(structure).write_file("test.vesta")
    actual = Path("test.vesta").read_text()
    expected = """#VESTA_FORMAT_VERSION 3.5.0
TITLE
X1 Ba1

CELLP
3.900000 3.900000 3.900000 90.000000 90.000000 90.000000

STRUC
1 Ba Ba1 1.0 0.500000 0.500000 0.500000
 0.0 0.0 0.0 
2 XX XX2 1.0 0.000000 0.000000 0.000000
 0.0 0.0 0.0 
 0 0 0 0 0 



SBOND

 0 0 0 0 

SITET
 1  Ba1  0.5 30 239 44 30 239 44 204 1
 2  X0+2  0.5 30 30 30 30 30 30 204 1
 0 0 0 0 



STYLE
VECTS  1.0
SECTS  64  1
SECTP 
   1 0 0  0.00000E+00  0.00000E+00  -1.00000E+01  1.00000E+01
UCOLP 
 0  2  1.000   0   0   0
ATOMS  1  0  1
BONDP 
   1  16  0.12  1.000 127 127 127"""
    assert actual == expected


def test_add_density(tmpdir):
    tmpdir.chdir()
    vesta_text = """#VESTA_FORMAT_VERSION 3.5.0
TITLE
X1 Ba1

CELLP
3.900000 3.900000 3.900000 90.000000 90.000000 90.000000

STRUC
1 Ba Ba1 1.0 0.500000 0.500000 0.500000
0.0 0.0 0.0
2 XX XX2 1.0 0.000000 0.000000 0.000000
0.0 0.0 0.0
0 0 0 0 0"""

    actual = add_density(original_vesta_text=vesta_text,
                         isurfs=[0.1, 0.2], volumetric_filename="PARCHG")

    expected = """#VESTA_FORMAT_VERSION 3.5.0
TITLE
X1 Ba1

IMPORT_DENSITY 0.1
+1.000000 PARCHG

CELLP
3.900000 3.900000 3.900000 90.000000 90.000000 90.000000

STRUC
1 Ba Ba1 1.0 0.500000 0.500000 0.500000
0.0 0.0 0.0
2 XX XX2 1.0 0.000000 0.000000 0.000000
0.0 0.0 0.0
0 0 0 0 0

ISURF
  1   1  0.1  0  0  255  50  50
  1   1  0.2  0  0  255  50  50
  0   0   0   0"""
    assert actual == expected


def test_calc_isurfs():
    au_volume_in_ang = 0.529177210903**3
    actual = calc_isurfs([1 / au_volume_in_ang], is_chg=True, volume=10.0)
    expected = [0.1]
    assert actual == expected

    actual = calc_isurfs([0.1], is_chg=False, volume=10.0)
    expected = [0.1]
    assert actual == expected
