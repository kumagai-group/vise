from enum import unique, Enum
from typing import Optional

from pymatgen.io.vasp import Potcar, Kpoints
from pymatgen.core.structure import Structure

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


@unique
class Task(Enum):
    """ Supported tasks """
    structure_opt = "structure_opt"
    structure_opt_rough = "structure_opt_rough"
    structure_opt_tight = "structure_opt_tight"
    cluster_opt = "cluster_opt"
    phonon_force = "phonon_force"
    defect = "defect"
    band = "band"
    dos = "dos"
    dielectric_dfpt = "dielectric_dfpt"
    dielectric_finite_field = "dielectric_finite_field"
    dielectric_function = "dielectric_function"

    def __str__(self):
        return self.name

    @classmethod
    def from_string(cls, s):
        for m in Task:
            if m.value == s:
                return m
        raise AttributeError(f"Task: {s} is not proper.\n"
                             f"Supported Task:\n {cls.name_list()}")

    @classmethod
    def name_list(cls):
        return ', '.join([e.value for e in cls])


def check_keys(d: dict, required: set, optional: set):

    check_required = set(required)
    for key in d:
        check_required.remove(key)
        if key not in required | optional:
            raise KeyError(f"{key} does not belong.")

    if check_required:
        raise KeyError(f"{check_required} must be set.")

def nelect(potcar: Potcar, structure: Structure):
    """ Gets the default number of electrons for a given structure. """
    # if structure is not sorted this can cause problems, so must take
    # care to remove redundant symbols when counting electrons
    site_symbols = list(set(structure.site_symbols))
    nelect = 0
    for ps in potcar:
        if ps.element in site_symbols:
            site_symbols.remove(ps.element)
            nelect += structure.composition.element_composition[
                          ps.element] * ps.ZVAL

    return nelect - structure.charge


TASK_REQUIRED_FLAGS = {"LREAL", "ENCUT", "ISPIN", "ISIF", "EDIFF", "IBRION",
                       "ISMEAR", "NSW"}
TASK_OPTIONAL_FLAGS = {"EDIFFG", "POTIM", "ADDGRID", "KPAR"}

XC_REQUIRED_FLAGS = {"ALGO", "NELM", "LWAVE"}
XC_OPTIONAL_FLAGS = {"LDAU", "NKRED", "LHFCALC", "TIME"}

XC_TASK_REQUIRED_FLAGS = {"PREC", "NBANDS"}
XC_TASK_OPTIONAL_FLAGS = set()

COMMON_REQUIRED_FLAGS = {"LASPH", "LORBIT", "LCHARG"}
COMMON_OPTIONAL_FLAGS = {"SIGMA", "NELECT"}


class TaskInputSet:

    def __init__(self,
                 incar_settings: dict,
                 structure: Structure,
                 kpoints: Kpoints):

        check_keys(incar_settings, TASK_REQUIRED_FLAGS, TASK_OPTIONAL_FLAGS)
        self.incar_settings = incar_settings
        self.structure = structure

    @classmethod
    def from_options(cls,
                     task: Task,
                     structure: Structure,
                     kpt_mode: str,
                     kpt_density: float,
                     band_ref_dist: float,
                     band_gap: float,
                     standardize_structure: bool,
                     factor: int,
                     encut: float,
                     structure_opt_encut_factor: float,
                     is_magnetization: bool,
                     hubbard_u: Optional[bool]):

        # Be careful LREAL = F with large unitcell
        incar_settings = dict()
        incar_settings["NELM"] = 100


class XcInputSet:

    def __init__(self, incar_settings: dict):
        check_keys(incar_settings, XC_REQUIRED_FLAGS, XC_OPTIONAL_FLAGS)
        self.incar_settings = incar_settings

    @classmethod
    def from_options(cls,
                     xc: Xc,
                     structure: Structure,
                     factor: int,
                     hubbard_u: Optional[bool],
                     ldauu: Optional[dict],
                     ldaul: Optional[dict]):

        incar_settings = dict()
        incar_settings["NELM"] = 100

        if xc in SEMILOCAL:
            incar_settings["ALGO"] = "N"
            incar_settings["LWAVE"] = False

        elif xc in HYBRID_FUNCTIONAL:
            incar_settings["ALGO"] = "D"
            incar_settings["LWAVE"] = True
            incar_settings["TIME"] = 0.5
            incar_settings["LHFCALC"] = True
            incar_settings["PRECFOCK"] = "Fast"
            incar_settings["AEXX"] = 0.25
            if xc is Xc.hse:
                incar_settings["HFSCREEN"] = 0.208

            if factor > 1:
                incar_settings["NKRED"] = 1
        else:
            raise ValueError

        if hubbard_u is None:
            hubbard_u = True if xc in LDA_OR_GGA else False

        if hubbard_u:
            incar_settings["LDAU"] = True
            incar_settings["LDAUTYPE"] = 2
            incar_settings["LDAUPRINT"] = True
            # Bring'LDAUU', 'LDAUJ', 'LDAUL' from incar property
            # incar_settings["LMAXMIX"] = 4
            # if 'LMAXMIX' not in settings.keys():
            #     # contains f-electrons
            #     if any([el.Z > 56 for el in structure.composition]):
            #         incar['LMAXMIX'] = 6
            #     # contains d-electrons
            #     elif any([el.Z > 20 for el in

        return cls(incar_settings=incar_settings)


class XcTaskInputSet:

    # Add allowed xc and task combinations here.

    def __init__(self,
                 incar_settings: dict, potcar: Potcar):

        check_keys(incar_settings, XC_REQUIRED_FLAGS, XC_OPTIONAL_FLAGS)
        self.incar_settings = incar_settings
        self.potcar = potcar

    @classmethod
    def from_options(cls,
                     xc: Xc,
                     task: Task,
                     structure: Structure):

    # potcar here


class CommonInputSet:
    def __init__(self, incar_settings: dict):
        check_keys(incar_settings, XC_REQUIRED_FLAGS, XC_OPTIONAL_FLAGS)
        self.incar_settings = incar_settings

    @classmethod
    def from_options(cls,
                     potcar: Potcar,
                     structure: Structure):
        """ """
        incar_settings = {"LASPH": True,
                          "LORBIT": 12,
                          "LCHARG": False,
                          "SIGMA": 0.1,
                          "NELCT": nelect(potcar, structure)}

        return cls(incar_settings=incar_settings)



