from copy import deepcopy
from typing import Optional

from monty.serialization import loadfn
from pymatgen import Structure, Composition, Element
from pymatgen.io.vasp import Potcar
from vise.config import BAND_GAP_CRITERION
from vise.input_set.vise_other_settings import SET_DIR, check_keys, logger, load_default_incar_settings, nbands, \
    nelect
from vise.input_set.xc import Xc, LDA_OR_GGA, HYBRID_FUNCTIONAL
from vise.input_set.task import LATTICE_RELAX_TASK, SPECTRA_TASK, Task

TASK_REQUIRED_FLAGS = {"LREAL", "ISPIN", "ISIF", "EDIFF", "IBRION", "ISMEAR",
                       "PREC"}
TASK_OPTIONAL_FLAGS = {"NSW", "EDIFFG", "POTIM", "ADDGRID", "KPAR", "ENCUT",
                       "NBANDS"}
TASK_FLAGS = TASK_REQUIRED_FLAGS | TASK_OPTIONAL_FLAGS
XC_REQUIRED_FLAGS = {"ALGO", "LWAVE"}
XC_OPTIONAL_FLAGS = {"LDAU", "LDAUTYPE", "LDAUPRINT", "LDAUU", "LDAUL",
                     "LMAXMIX", "NKRED", "LHFCALC", "PRECFOCK", "TIME",
                     "HFSCREEN", "AEXX", "NKRED"}
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
                     structure: Structure,
                     composition: Optional[Composition],
                     potcar: Potcar,
                     num_kpoints: int,
                     max_enmax: float,
                     is_magnetization: bool,
                     vbm_cbm: list,
                     npar_kpar: bool,
                     num_cores: list,
                     encut: Optional[float],
                     structure_opt_encut_factor: float):
        """
        See docstrings in the make_input method in ViseInputSet class
        """
        band_gap = vbm_cbm[1] - vbm_cbm[0] if vbm_cbm else None
        if band_gap:
            logger.info(f"Band gap: {round(band_gap, 3)} eV.")
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
            if hasattr(structure, "site_properties"):
                settings["MAGMOM"] = structure.site_properties["magmom"]

        if npar_kpar:
            settings["KPAR"] = 2

        if encut:
            settings["ENCUT"] = encut
        else:
            if task in LATTICE_RELAX_TASK:
                settings["ENCUT"] = \
                    round(max_enmax * structure_opt_encut_factor, 3)

        if task in SPECTRA_TASK:
            settings["NBANDS"] = nbands(composition, potcar)

        return cls(settings=settings)


class XcIncarSettings:

    def __init__(self, settings: dict):
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
                     ldaul_set_name: Optional[str] = None):

        settings = \
            load_default_incar_settings(yaml_filename="xc_incar_set.yaml",
                                        required_flags=XC_REQUIRED_FLAGS,
                                        optional_flags=XC_OPTIONAL_FLAGS,
                                        key_name=str(xc))

        if hubbard_u is None:
            hubbard_u = xc in LDA_OR_GGA
        ldaul_set_name = ldaul_set_name or "lda_gga_normal"
        ldauu = ldauu or {}
        ldaul = ldaul or {}

        if hubbard_u:
            u_set = loadfn(SET_DIR / "u_parameter_set.yaml")
            ldauu_set = u_set["LDAUU"][ldaul_set_name]
            ldauu_set.update(ldauu)
            ldauu = [ldauu_set.get(el, 0) for el in symbol_set]

            if sum(ldauu) > 0:
                settings["LDAUU"] = ldauu

                ldaul_set = u_set["LDAUL"][ldaul_set_name]
                ldaul_set.update(ldaul)
                settings["LDAUL"] = [ldaul_set.get(el, -1) for el in symbol_set]

                settings.update({"LDAU": True, "LDAUTYPE": 2, "LDAUPRINT": 1})
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
        check_keys(settings, COMMON_REQUIRED_FLAGS, COMMON_OPTIONAL_FLAGS)
        self.settings = settings

    @classmethod
    def from_options(cls,
                     potcar: Potcar,
                     composition: Composition,
                     charge: int):
        """ """
        settings = {"NELM": 100, "LASPH": True, "LORBIT": 12, "LCHARG": False,
                    "SIGMA": 0.1}
        # Structure charge is ignored.
        if charge:
            settings["NELECT"] = nelect(composition, potcar, charge)

        return cls(settings=settings)