# -*- coding: utf-8 -*-

from enum import unique, Enum


@unique
class Xc(Enum):
    """ Supported exchange-correlation treatment. """
    pbe = "pbe"
    pbesol = "pbesol"
    lda = "lda"
    scan = "scan"
    pbe0 = "pbe0"
    hse = "hse"

    def __str__(self):
        return self.name

    @classmethod
    def from_string(cls, s):

        for m in Xc:
            if m.value == s:
                return m
        if s == "perdew-zunger81":
            return Xc.lda
        raise AttributeError(f"Xc:{s} is not proper.\n "
                             f"Supported Xc:\n {cls.name_list()}")

    @classmethod
    def name_list(cls):
        return ', '.join([e.value for e in cls])

    @property
    def require_wavefunctions(self):
        return True if self in BEYOND_MGGA else False


LDA_OR_GGA = (Xc.pbe, Xc.pbesol, Xc.lda)
SEMILOCAL = (Xc.pbe, Xc.pbesol, Xc.lda, Xc.scan)
HYBRID_FUNCTIONAL = (Xc.pbe0, Xc.hse)
MGGA_OR_HYBRID = (Xc.pbe0, Xc.hse, Xc.scan)
BEYOND_MGGA = HYBRID_FUNCTIONAL
