# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

from pymatgen.core import Element
from vise.atom_energies.make_atom_vasp_set import is_target_element, \
    make_atom_poscar_dirs


def test_is_target_element():
    assert is_target_element(Element.H) is True
    assert is_target_element(Element.Pr) is False


def test_make_atom_poscar_dirs(tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    make_atom_poscar_dirs(Path.cwd(), ["H"])
    prior_info = Path("H") / "prior_info.yaml"
    expected = """incar:
  ISPIN: 2
  NELM: 300
  NUPDOWN: 1.0
is_cluster: true
"""

    assert prior_info.read_text() == expected

    poscar = Path("H") / "POSCAR"
    expected = """H1
1.0
  10.0000000000000000    0.0000000000000000    0.0000000000000000
   0.0000000000000000   10.0000000000000000    0.0000000000000000
   0.0000000000000000    0.0000000000000000   10.0000000000000000
H
1
direct
   0.5000000000000000    0.5000000000000000    0.5000000000000000 H
"""

    assert poscar.read_text() == expected


