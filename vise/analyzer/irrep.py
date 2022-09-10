# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import Dict, List

from monty.json import MSONable


@dataclass
class Irrep(MSONable):
    frac_coords: List[float]
    symbols: List[str]
    energies: List[float]
    degeneracies: List[int]

    @property
    def irrep_info_set(self):
        return zip(self.symbols, self.energies, self.degeneracies)


@dataclass
class Irreps(MSONable):
    sg_num: int
    # key is special point name. Gamma = GM
    irreps: Dict[str, Irrep]

    def __call__(self):
        return self.irreps



