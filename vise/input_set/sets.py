from copy import deepcopy
from enum import unique, Enum
from math import ceil
from pathlib import Path
from typing import Optional

from monty.serialization import loadfn
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Potcar, Kpoints
from vise.core.config import BAND_GAP_CRITERION
from vise.input_set.datasets.element_parameters import unoccupied_bands
from vise.input_set.kpoints import make_kpoints
from vise.util.logger import get_logger
from vise.util.structure_handler import find_spglib_primitive

logger = get_logger(__name__)

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


LATTICE_RELAX_TASK = (Task.structure_opt, Task.structure_opt_rough,
                      Task.structure_opt_tight)
SPECTRA_TASK = (Task.dos, Task.dielectric_function)

SET_DIR = Path(__file__).parent / "datasets"


def load_potcar_yaml(set_name: Optional[str] = "normal",
                     override_potcar_set: Optional[dict] = None) -> dict:
    """Load the yaml setting files for config and POTCAR list.

    Args:
        set_name (str):
            Potcar set name.
        override_potcar_set (dict):
            User specifying POTCAR set

    Return:
          Dictionary of potcar_set, like {"Zr": "Zr_pv", ...}
    """
    potcar_set = loadfn(SET_DIR / "potcar_set.yaml")
    try:
        potcar = potcar_set[set_name]
    except KeyError:
        logger.warning(f"The accepted potcar set name is {potcar_set.keys()}")
        raise

    if override_potcar_set:
        potcar.update(override_potcar_set)

    return potcar


def load_default_incar_settings(yaml_filename: str,
                                required_flags: set,
                                optional_flags: set,
                                key_name: str) -> dict:
    """Load the yaml setting files for config and POTCAR list.

    Args:
        yaml_filename (str):
            Yaml filename with xxx_incar_set.yaml.
        required_flags (set):
            Required INCAR flags.
        optional_flags (set):
            Optional INCAR flags.
        key_name (str):
            Key name such as "structure_opt" or "hse".

    Return:
          settings (dict):
    """
    incar_set = loadfn(SET_DIR / yaml_filename)
    settings = {}
    for f in required_flags:
        settings[f] = incar_set[f][key_name]

    for f in optional_flags:
        try:
            settings[f] = incar_set[f][key_name]
        except KeyError:
            continue

    return settings


def check_keys(d: dict, required: set, optional: set) -> bool:
    """Check if the required keys exist and other keys are included in optional

    Args:
        d (dict):
            Dictionary to be checked.
        required (set):
            Required key set.
        optional (set):
            Optional key set.

    Return:
        Simply raise KeyError if the condition is not satisfied.
        If succeed, True is returned.
    """
    check_required = set(required)
    for key in d:
        check_required.discard(key)
        if key not in required | optional:
            raise KeyError(f"{key} does not belong.")

    if check_required:
        raise KeyError(f"{check_required} must exist in {d}.")

    return True


def nelect(composition: Composition, potcar: Potcar, charge: int = 0) -> int:
    """Gets the default number of electrons for a given structure.


     """
    # if structure is not sorted this can cause problems, so must take
    # care to remove redundant symbols when counting electrons
    num_nelect = 0
    for pt in potcar:
        num_nelect += composition.element_composition[pt.element] * pt.ZVAL

    return int(num_nelect) - charge


def nbands(composition: Composition, potcar: Potcar) -> int:
    """
    Calculate the total number of bands required for the unoccupied related
    properties such as optical absorption, band structure, and DOS.

    Args:
        composition:
        potcar (Potcar):
    """
    num_bands = sum([composition[c] * (p.nelectrons / 2 + unoccupied_bands[str(c)])
                    for c, p in zip(composition, potcar)])

    return ceil(num_bands)


class TaskStructureKpoints:

    def __init__(self,
                 structure: Structure,
                 kpoints: Kpoints,
                 is_structure_changed: bool,
                 sg: int,
                 num_kpts: int):

        self.structure = structure
        self.kpoints = kpoints
        self.is_structure_changed = is_structure_changed
        self.sg = sg
        self.num_kpts = num_kpts

    @classmethod
    def from_options(cls,
                     task: Task,
                     original_structure: Structure,
                     standardize_structure: bool,
                     sort_structure: bool,
                     is_magnetization: bool,
                     kpt_mode: str,
                     kpt_density: float,
                     kpt_shift: list,
                     only_even: bool,
                     band_ref_dist: float,
                     factor: Optional[int],
                     symprec: float,
                     angle_tolerance: float):
        """

        See docstrings in the make_input method in InputSet class
        :return:
        """
        if sort_structure:
            structure = original_structure.get_sorted_structure()
        else:
            structure = original_structure.copy()

        is_structure_changed = False
        if task == Task.defect:
            kpt_mode = "manual_set"
        elif task == Task.cluster_opt:
            kpt_mode = "manual_set"
            kpt_density = 1e-5
            only_even = False
            kpt_shift = [0, 0, 0]
        elif task == Task.band:
            kpt_mode = "band"
        else:
            primitive_structure, is_structure_changed = \
                find_spglib_primitive(structure, symprec, angle_tolerance)

            if standardize_structure:
                structure = primitive_structure
            else:
                if is_structure_changed and kpt_mode != "manual_set":
                    logger.warning(
                        "Standardizaion is set to False and the given "
                        "structure is not a primitive cell. Thus, the "
                        "kpoint set is switched to manual_set.")
                    kpt_mode = "manual_set"
                else:
                    logger.info("The structure is a standardized primitive.")

        # - KPOINTS construction
        kpoints, structure, sg, num_kpts = \
            make_kpoints(mode=kpt_mode,
                         structure=structure,
                         kpts_density=kpt_density,
                         only_even=only_even,
                         ref_distance=band_ref_dist,
                         kpt_shift=kpt_shift,
                         factor=factor,
                         symprec=symprec,
                         angle_tolerance=angle_tolerance,
                         is_magnetization=is_magnetization)

        return cls(structure, kpoints, is_structure_changed, sg, num_kpts)


class XcTaskPotcar:

    def __init__(self,
                 potcar: Potcar,
                 max_enmax: int):

        self.potcar = potcar
        self.max_enmax = max_enmax

    @classmethod
    def from_options(cls,
                     xc: Xc,
                     task: Task,
                     symbol_set: tuple,
                     potcar_set_name: str = None,
                     override_potcar_set: dict = None):

        potcar_functional = "LDA" if xc == Xc.lda else "PBE"
        potcar_set_name = potcar_set_name or "normal"

        potcar_list = load_potcar_yaml(potcar_set_name, override_potcar_set)
        potcar_symbols = [potcar_list.get(el, el) for el in symbol_set]
        potcar = Potcar(potcar_symbols, functional=potcar_functional)

        max_enmax = max([p.enmax for p in potcar])

        return cls(potcar, max_enmax)


TASK_REQUIRED_FLAGS = {"LREAL", "ISPIN", "ISIF", "EDIFF", "IBRION", "ISMEAR",
                       "PREC"}
TASK_OPTIONAL_FLAGS = {"NSW", "EDIFFG", "POTIM", "ADDGRID", "KPAR", "ENCUT",
                       "NBANDS"}
TASK_FLAGS = TASK_REQUIRED_FLAGS | TASK_OPTIONAL_FLAGS

XC_REQUIRED_FLAGS = {"ALGO", "LWAVE"}
XC_OPTIONAL_FLAGS = {"LDAU", "LDAUTYPE", "LDAUPRINT", "LDAUU", "LDAUL",
                     "LMAXMIX", "NKRED", "LHFCALC", "PRECFOCK", "TIME",
                     "LHFSCREEN", "AEXX", "NKRED"}
XC_FLAGS = XC_REQUIRED_FLAGS | XC_OPTIONAL_FLAGS

XC_TASK_REQUIRED_FLAGS = set()
XC_TASK_OPTIONAL_FLAGS = set()
XC_TASK_FLAGS = XC_TASK_REQUIRED_FLAGS | XC_TASK_OPTIONAL_FLAGS

COMMON_REQUIRED_FLAGS = {"NELM", "LASPH", "LORBIT", "LCHARG", "SIGMA"}
COMMON_OPTIONAL_FLAGS = {"NELECT"}
COMMON_FLAGS = COMMON_REQUIRED_FLAGS | COMMON_OPTIONAL_FLAGS

ALL_FLAGS = set(sum(loadfn(SET_DIR / "incar_flags.yaml").values(), []))
OTHER_FLAGS = ALL_FLAGS - (TASK_FLAGS | XC_FLAGS | XC_TASK_FLAGS | COMMON_FLAGS)


class TaskIncarSettings:

    def __init__(self,
                 settings: dict):

        check_keys(settings, TASK_REQUIRED_FLAGS, TASK_OPTIONAL_FLAGS)
        self.settings = settings

    @classmethod
    def from_options(cls,
                     task: Task,
                     composition: Optional[Composition],
                     potcar: Potcar,
                     num_kpoints: int,
                     max_enmax: float,
                     is_magnetization: bool,
                     band_gap: float,
                     vbm_cbm: list,
                     npar_kpar: bool,
                     num_cores: list,
                     encut: Optional[float],
                     structure_opt_encut_factor: float):
        """
        See docstrings in the make_input method in InputSet class
        """
        band_gap = band_gap or vbm_cbm[1] - vbm_cbm[0] if vbm_cbm else None
        if band_gap:
            logger.info(f"Band gap: {band_gap} eV.")
            is_band_gap = band_gap > BAND_GAP_CRITERION
        else:
            is_band_gap = None

        required = deepcopy(TASK_REQUIRED_FLAGS)
        settings = \
            load_default_incar_settings(yaml_filename="task_incar_set.yaml",
                                        required_flags=required,
                                        optional_flags=TASK_OPTIONAL_FLAGS,
                                        key_name=str(task))

        settings["LREAL"] = "A" if task == Task.defect else "F"

        if is_band_gap is False or num_kpoints < 4:
            settings["ISMEAR"] = 0

        if task == Task.defect or is_magnetization:
            settings["ISPIN"] = 2

        if npar_kpar:
            settings["KPAR"] = 2

        if encut:
            settings["ENCUT"] = encut
        else:
            if task in LATTICE_RELAX_TASK:
                settings["ENCUT"] = max_enmax * structure_opt_encut_factor

        if task in SPECTRA_TASK:
            settings["NBANDS"] = nbands(composition, potcar)

        return cls(settings=settings)


class XcIncarSettings:

    def __init__(self, settings: dict):
        print(XC_OPTIONAL_FLAGS)
        check_keys(settings, XC_REQUIRED_FLAGS, XC_OPTIONAL_FLAGS)
        self.settings = settings

    @classmethod
    def from_options(cls,
                     xc: Xc,
                     symbol_set: tuple,
                     factor: int,
                     aexx: Optional[float] = 0.25,
                     hubbard_u: Optional[bool] = None,
                     ldauu: Optional[dict] = None,
                     ldaul: Optional[dict] = None,
                     ldaul_set_name: Optional[str] = "lda_gga_normal"):

        settings = \
            load_default_incar_settings(yaml_filename="xc_incar_set.yaml",
                                        required_flags=XC_REQUIRED_FLAGS,
                                        optional_flags=XC_OPTIONAL_FLAGS,
                                        key_name=str(xc))

        if not hubbard_u:
            hubbard_u = xc in LDA_OR_GGA

        if hubbard_u:
            settings.update({"LDAU": True, "LDAUTYPE": 2, "LDAUPRINT": 1})
            u_set = loadfn(SET_DIR / "u_parameter_set.yaml")
            ldauu_set = u_set["LDAUU"][ldaul_set_name]
            ldauu_set.update(ldauu)
            ldaul_set = u_set["LDAUL"][ldaul_set_name]
            ldaul_set.update(ldaul)

            settings["LDAUU"] = [ldauu_set[el] for el in symbol_set]
            settings["LDAUL"] = [ldaul_set[el] for el in symbol_set]

            # contains f-electrons
            if any([Element(el).Z > 56 for el in symbol_set]):
                settings["LMAXMIX"] = 6
            else:
                settings["LMAXMIX"] = 4

        if xc in HYBRID_FUNCTIONAL:
            settings["AEXX"] = aexx
            if factor > 1:
                settings["NKRED"] = factor

        return cls(settings=settings)


class XcTaskIncarSettings:

    def __init__(self, settings: dict):
        check_keys(settings, XC_TASK_REQUIRED_FLAGS, XC_TASK_OPTIONAL_FLAGS)
        self.settings = settings

    @classmethod
    def from_options(cls):
        settings = {}
        return cls(settings)


class CommonIncarSettings:
    def __init__(self, settings: dict):
        check_keys(settings, XC_REQUIRED_FLAGS, XC_OPTIONAL_FLAGS)
        self.settings = settings

    @classmethod
    def from_options(cls,
                     potcar: Potcar,
                     composition: Composition,
                     charge: int):
        """ """
        settings = {"NELM": 100, "LASPH": True, "LORBIT": 12, "LCHARG": False,
                    "SIGMA": 0.1}
        if charge:
            settings["NELCT"] = nelect(composition, potcar, charge)

        return cls(settings=settings)


