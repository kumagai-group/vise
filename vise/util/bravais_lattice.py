# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path
from typing import List, Dict

from monty.serialization import loadfn
from vise.util.enum import ExtendedEnum
from vise.util.structure_symmetrizer import StructureSymmetrizer
hermann_mauguin_list: Dict[int, str] = \
    loadfn(Path(__file__).parent / "Hermann–Mauguin.yaml")


class BravaisLattice(ExtendedEnum):
    aP = "aP"
    mP = "mP"
    mC = "mC"
    oP = "oP"
    oF = "oF"
    oI = "oI"
    oA = "oA"
    oC = "oC"
    tP = "tP"
    tI = "tI"
    hR = "hR"
    hP = "hP"
    cP = "cP"
    cF = "cF"
    cI = "cI"

    @classmethod
    def from_sg_num(cls, sg_num: int) -> "BravaisLattice":
        cs = {"a": (1, 2),
              "m": (3, 15),
              "o": (16, 74),
              "t": (75, 142),
              "h": (143, 194),
              "c": (195, 230)}
        centering_capital_letter = hermann_mauguin_list[sg_num][0]
        for crystal_system, (initial_sg, final_sg) in cs.items():
            if initial_sg <= sg_num <= final_sg:
                bravais_letter = crystal_system + centering_capital_letter
                break
        else:
            raise ValueError

        return cls.from_string(bravais_letter)

    @classmethod
    def from_structure_symmetrizer(cls,
                                   structure_symmetrizer: StructureSymmetrizer):
        bravais_string = structure_symmetrizer.seekpath_data["bravais_lattice"]
        return cls.from_string(bravais_string)

    @property
    def kpt_centering(self) -> List[float]:
        if self in [self.oF, self.tI, self.cF, self.cI]:
            return [0.0, 0.0, 0.0]
        elif self is self.hP:
            return [0.0, 0.0, 0.5]
        else:
            return [0.5, 0.5, 0.5]

    @property
    def need_same_num_kpt(self) -> bool:
        return True if self in (self.oI, self.tI) else False
