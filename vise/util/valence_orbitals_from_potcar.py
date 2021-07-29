# -*- coding: utf-8 -*-
#  Copyright (c) 2021. Distributed under the terms of the MIT License.
from pymatgen.io.vasp import PotcarSingle


def valence_orbitals_from_potcar(potcar_single: PotcarSingle):

    result = []
    for config in potcar_single.electron_configuration:
        principal_num, orb_type, num_ele = config
        result.append(f"({principal_num}{orb_type}){num_ele}")

    return " ".join(result)


