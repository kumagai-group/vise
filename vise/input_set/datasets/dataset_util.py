# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from math import ceil
from pathlib import Path
from typing import Tuple, List, Dict, Any, Optional
from copy import deepcopy

from monty.serialization import loadfn
from pymatgen import Composition, Element
from pymatgen.io.vasp import Potcar
from vise.util.enum import ExtendedEnum

unoccupied_bands = loadfn(Path(__file__).parent / "unoccupied_bands.yaml")

# This incar_flags should be OrderedDict, but from python 3.6, dict uses
# order-preserving semantics. Besides, it does not affect vasp result.
all_incar_flags: Dict[str, Any] = \
    loadfn(Path(__file__).parent / "incar_flags.yaml")


def has_f_elements(symbol_list: List[str]):
    return any([Element(el).Z > 56 for el in symbol_list])


class LDAU:
    def __init__(self,
                 symbol_list: List[str],
                 override_ldauu: Optional[Dict[str, float]] = None,
                 override_ldaul: Optional[Dict[str, float]] = None):

        self.symbol_list = symbol_list
        ldau_set = loadfn(Path(__file__).parent / "u_parameter_set.yaml")

        ldauu_set = ldau_set["LDAUU"]
        ldauu_set.update(override_ldauu or {})
        self.ldauu = [ldauu_set.get(el, 0) for el in symbol_list]

        ldaul_set = ldau_set["LDAUL"]
        ldaul_set.update(override_ldaul or {})
        self.ldaul = [ldaul_set.get(el, -1) for el in symbol_list]

    @property
    def lmaxmix(self):
        return 6 if has_f_elements(self.symbol_list) else 4

    @property
    def is_ldau_needed(self) -> bool:
        return set(self.ldauu) != {0}


class PotcarSet(ExtendedEnum):
    normal = "normal"
    mp_relax_set = "mp_relax_set"
    gw = "gw"

    def overridden_potcar_dict(self,
                               override_potcar_set: Optional[dict] = None
                               ) -> Dict[str, str]:
        result = deepcopy(self.potcar_dict())
        if override_potcar_set:
            result.update(override_potcar_set)

        return result  # e.g. {"Zr": "Zr_pv", ...}

    def potcar_dict(self) -> Dict[str, str]:
        _potcar_list = loadfn(Path(__file__).parent / "potcar_set.yaml")
        set_names = _potcar_list.pop("set_names")
        result = {}
        for name in set_names:
            result[name] = {}  # normal, mp_relax_set, gw

        for element, potcar_string in _potcar_list.items():
            potcars = potcar_string.split()

            def sanitize(val: str):
                if val == "---":
                    return potcars[0]
                elif val == "None":
                    return None
                else:
                    return val

            for index, potcar_single in enumerate(potcars):
                set_name = set_names[index]
                result[set_name][element] = sanitize(potcar_single)

        return result[self.value]


def num_bands(composition: Composition, potcar: Potcar) -> int:
    """Required for optical absorption, band structure, and DOS"""
    results = 0
    for element, potcar_single in zip(composition, potcar):
        num_atoms_per_element = composition[element]
        occupied_bands = potcar_single.nelectrons / 2
        num_bands_per_atom = occupied_bands + unoccupied_bands[str(element)]
        results += num_atoms_per_element * num_bands_per_atom

    return ceil(results)


def npar_kpar(num_kpoints: int, num_nodes: int) -> Tuple[int, int]:
    kpar_set = loadfn(Path(__file__).parent / "kpar_set.yaml")
    num_kpt_key = num_kpoints if num_kpoints in kpar_set else "None"
    if num_nodes == 2:
        kpar = kpar_set[num_kpt_key][1]
    elif num_nodes == 4:
        kpar = kpar_set[num_kpt_key][2]
    else:
        kpar = kpar_set[num_kpt_key][0]

    npar = 2 if kpar == 1 else 1

    return kpar, npar

