# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import List
from vise.util.enum import ExtendedEnum

from pathlib import Path
from monty.serialization import loadfn

hermann_mauguin_list: List[str] = \
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
    def from_sg_number(cls, sg_num):
        cs = {"a": (1, 2),
              "m": (3, 15),
              "o": (16, 74),
              "t": (75, 142),
              "h": (143, 194),
              "c": (195, 230)}

        for crystal_system, (initial_sg, final_sg) in cs.items():
            if initial_sg <= sg_num <= final_sg:
                sg = hermann_mauguin_list[sg_num]
                centering = sg[0]
                break
        else:
            raise ValueError

        return cls(crystal_system + centering)

    @property
    def kpt_centering(self):
        if self in [self.oF, self.tI, self.cF, self.cI]:
            return [0.0, 0.0, 0.0]
        elif self is self.hP:
            return [0.0, 0.0, 0.5]
        else:
            return [0.5, 0.5, 0.5]

    @property
    def need_same_num_kpt(self):
        return True if self in (self.oI, self.tI) else False
