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
    gw0 = "gw0"

    def __str__(self):
        return self.name

    @classmethod
    def from_string(cls, s):

        for m in Xc:
            if m.value == s:
                return m
        if s == "perdew-zunger81":
            return Xc.lda
        raise AttributeError("Xc: " + str(s) + " is not proper.\n" +
                             "Supported Xc:\n" + cls.name_list())

    @classmethod
    def name_list(cls):
        return ', '.join([e.value for e in cls])

    @property
    def require_wavefunctions(self):
        return True if self in BEYOND_DFT else False


LDA_OR_GGA = (Xc.pbe, Xc.pbesol, Xc.lda)
DFT_FUNCTIONAL = (Xc.pbe, Xc.pbesol, Xc.lda, Xc.scan)
HYBRID_FUNCTIONAL = (Xc.pbe0, Xc.hse)
BEYOND_GGA = (Xc.pbe0, Xc.hse, Xc.scan)
GW = (Xc.gw0, )
BEYOND_DFT = HYBRID_FUNCTIONAL + GW