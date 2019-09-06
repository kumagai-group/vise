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


class XcInputSet:

    incar_required_flags = ["ALGO", "NELM", "LWAVE"]
    incar_optional_flags = ["LDAU", "NKRED", "LHFCALC", "TIME"]
    incar_flags = incar_required_flags + incar_optional_flags

    def __init__(self, incar_settings: dict):
        self.incar_settings = incar_settings

    @classmethod
    def from_options(cls, xc: Xc, factor: int, hubbard_u: bool):
        incar_settings = dict()
        incar_settings["NELM"] = 100

        if xc in SEMILOCAL:
            incar_settings["ALGO"] = "N"
            incar_settings["LWAVE"] = False

        elif xc in HYBRID_FUNCTIONAL:
            incar_settings["ALGO"] = "D"
            incar_settings["LWAVE"] = True
            incar_settings["TIME"] = 0.5

            if factor > 1:
                incar_settings["NKRED"] = 1
        else:
            raise ValueError





