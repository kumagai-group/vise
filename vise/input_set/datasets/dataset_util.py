# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from math import ceil
from pathlib import Path
from typing import Tuple

from monty.serialization import loadfn
from pymatgen import Composition
from pymatgen.io.vasp import Potcar

unoccupied_bands = loadfn(Path(__file__).parent / "unoccupied_bands.yaml")
# magmom = loadfn(Path(__file__).parent / "magmom.yaml")

# This incar_flags should be OrderedDict, but from python 3.6, dict uses
# order-preserving semantics. Besides, it does not affect vasp result.
incar_tags = loadfn(Path(__file__).parent / "incar_flags.yaml")
kpar_set = loadfn(Path(__file__).parent / "kpar_set.yaml")


def nbands(composition: Composition, potcar: Potcar) -> int:
    """Number of bands required for optical absorption, band structure, & DOS"""
    results = 0
    for c, p in zip(composition, potcar):
        num_atoms_per_element = composition[c]
        num_bands_per_atom = p.nelectrons / 2 + unoccupied_bands[str(c)]
        results += num_atoms_per_element * num_bands_per_atom

    return ceil(results)


def npar_kpar(num_kpoints: int,
              num_nodes: int) -> Tuple[int, int]:

    num_kpt_key = num_kpoints if num_kpoints in kpar_set else "None"
    if num_nodes == 2:
        kpar = kpar_set[num_kpt_key][1]
    elif num_nodes == 4:
        kpar = kpar_set[num_kpt_key][2]
    else:
        kpar = kpar_set[num_kpt_key][0]

    npar = 2 if kpar == 1 else 1

    return kpar, npar

