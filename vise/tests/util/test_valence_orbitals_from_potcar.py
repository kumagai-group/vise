# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from pymatgen.io.vasp import Potcar
from vise.input_set.potcar_generator import generate_potcar
from vise.input_set.xc import Xc
from vise.util.valence_orbitals_from_potcar import valence_orbitals_from_potcar


def test_valence_orbitals_from_potcar():
    potcar_single = generate_potcar(["Mg"], xc=Xc.pbe)[0]
    assert valence_orbitals_from_potcar(potcar_single) == "(3s)2"
