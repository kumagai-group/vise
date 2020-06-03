# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from monty.json import MSONable
from vise.util.enum import ExtendedEnum


class Xc(MSONable, ExtendedEnum):
    pbe = "pbe"
    pbesol = "pbesol"
    lda = "lda"
    scan = "scan"
    pbe0 = "pbe0"
    hse = "hse"

    @classmethod
    def from_string(cls, s: str):
        if s == "perdew-zunger81":
            return cls.lda
        return super().from_string(s)

    @property
    def potcar_functional(self) -> str:
        return "LDA" if self == Xc.lda else "PBE_54"

    @property
    def is_lda_or_gga(self):
        return self in (self.pbe, self.pbesol, self.lda)

    @property
    def is_hybrid_functional(self):
        return self in (self.pbe0, self.hse)

    @property
    def is_local_or_semilocal(self):
        return self.is_lda_or_gga or self is self.scan

    @property
    def is_nonlocal(self):
        return self.is_hybrid_functional
