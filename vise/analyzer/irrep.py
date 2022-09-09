# -*- coding: utf-8 -*-
#  Copyright (c) 2022. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import Dict, List, Tuple


COORDS_TYPE = Tuple[float, float, float]


@dataclass
class Irrep:
    frac_coords: COORDS_TYPE
    symbols: List[str]
    energies: List[float]
    degeneracies: List[int]

    @property
    def irrep_info_set(self):
        return zip(self.symbols, self.energies, self.degeneracies)


@dataclass
class Irreps:
    sg_num: int
    # key is special point name. Gamma = GM
    irreps: Dict[str, Irrep]



