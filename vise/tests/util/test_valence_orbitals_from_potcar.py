# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from pymatgen.io.vasp import Potcar
from vise.util.valence_orbitals_from_potcar import valence_orbitals_from_potcar


def test_valence_orbitals_from_potcar():
    potcar_single = Potcar(["Ba_sv"], "PBE")[0]
    assert valence_orbitals_from_potcar(potcar_single) == "(6s)2 (5p)6 (5s)2"
