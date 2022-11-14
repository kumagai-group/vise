# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from math import ceil
from pathlib import Path
from typing import Tuple, List, Dict, Any

from monty.serialization import loadfn
from pymatgen.core import Composition, Element
from pymatgen.io.vasp import Potcar
from vise.defaults import defaults

unoccupied_bands = loadfn(Path(__file__).parent / "unoccupied_bands.yaml")

# This incar_flags should be OrderedDict, but from python 3.6, dict uses
# order-preserving semantics. Besides, it does not affect vasp result.
incar_categories: Dict[str, Any] = \
    loadfn(Path(__file__).parent / "incar_flags.yaml")
all_incar_flags: List[str] = sum(incar_categories.values(), [])


def has_f_elements(symbol_list: List[str]):
    return any([Element(el).Z > 56 for el in symbol_list])


class LDAU:
    def __init__(self,
                 symbol_list: List[str]):

        self.symbol_list = symbol_list
        ldau_set = loadfn(defaults.u_parameter_set_yaml_file)
        self.ldauu = [ldau_set["LDAUU"].get(el, 0) for el in symbol_list]
        self.ldaul = [ldau_set["LDAUL"].get(el, -1) for el in symbol_list]

    @property
    def lmaxmix(self):
        return 6 if has_f_elements(self.symbol_list) else 4

    @property
    def is_ldau_needed(self) -> bool:
        return set(self.ldauu) != {0}


def num_bands(composition: Composition, potcar: Potcar) -> int:
    """Required for optical absorption, band structure, and DOS"""
    results = 0
    for element, potcar_single in zip(composition, potcar):
        num_atoms_per_element = composition[element]
        occupied_bands = potcar_single.nelectrons / 2
        num_bands_per_atom = occupied_bands + unoccupied_bands[str(element)]
        results += num_atoms_per_element * num_bands_per_atom

    return ceil(results)


def calc_kpar(num_kpoints: int, num_cores: int,
              unused_core_ratio_threshold: float) -> int:
    divisors = [i for i in range(num_cores, 0, -1) if num_cores % i == 0]
    for d in divisors:
        kpt_by_core = round(num_kpoints / d, 5)
        frac_kpt_by_core = kpt_by_core % 1
        frac_kpt_by_core = 0.0 if not frac_kpt_by_core else 1 - frac_kpt_by_core
        unused_core_ratio = frac_kpt_by_core / ceil(kpt_by_core)
        if unused_core_ratio < unused_core_ratio_threshold:
            return d

    raise ValueError(f"The threshold for unused core ratio "
                     f"{unused_core_ratio_threshold} is not adequate.")


