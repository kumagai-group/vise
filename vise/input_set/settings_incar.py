# -*- coding: utf-8 -*-

from math import ceil
from pathlib import Path
import sys
from typing import Optional
from abc import abstractmethod, ABC
from collections import defaultdict

from monty.serialization import loadfn

from pymatgen import Structure, Composition, Element
from pymatgen.io.vasp import Potcar

from vise.config import BAND_GAP_CRITERION
from vise.input_set.settings_util import (
    load_default_incar_settings, check_keys, nelect, nbands, calc_npar_kpar)
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger

logger = get_logger(__name__)

SET_DIR = Path(__file__).parent / "datasets"

TASK_REQUIRED_FLAGS = {"LREAL", "ISPIN", "ISIF", "EDIFF", "IBRION", "ISMEAR",
                       "PREC"}
TASK_OPTIONAL_FLAGS = {"NSW", "EDIFFG", "POTIM", "ADDGRID", "KPAR", "NPAR",
                       "ENCUT", "NBANDS", "LEPSILON", "LCALCEPS", "LOPTICS",
                       "CSHIFT", "EMIN", "EMAX", "NEDOS", "MAGMOM"}
TASK_FLAGS = TASK_REQUIRED_FLAGS | TASK_OPTIONAL_FLAGS

XC_REQUIRED_FLAGS = {"ALGO", "LWAVE"}
XC_OPTIONAL_FLAGS = {"GGA", "METAGGA", "LDAU", "LDAUTYPE", "LDAUPRINT",
                     "LDAUU", "LDAUL", "LMAXMIX", "NKRED", "LHFCALC",
                     "PRECFOCK", "TIME", "HFSCREEN", "AEXX", "NKRED"}
XC_FLAGS = XC_REQUIRED_FLAGS | XC_OPTIONAL_FLAGS

XC_TASK_REQUIRED_FLAGS = set()
XC_TASK_OPTIONAL_FLAGS = set()
XC_TASK_FLAGS = XC_TASK_REQUIRED_FLAGS | XC_TASK_OPTIONAL_FLAGS

COMMON_REQUIRED_FLAGS = {"NELM", "LASPH", "LORBIT", "LCHARG", "SIGMA"}
COMMON_OPTIONAL_FLAGS = {"NELECT"}
COMMON_FLAGS = COMMON_REQUIRED_FLAGS | COMMON_OPTIONAL_FLAGS

incar_flag_list = sum(
    [list(i.keys()) for i in loadfn(SET_DIR / "incar_flags.yaml").values()], [])

ALL_FLAGS = set(incar_flag_list)

OTHER_FLAGS = ALL_FLAGS - (TASK_FLAGS | XC_FLAGS | XC_TASK_FLAGS | COMMON_FLAGS)


# incar_flags = defaultdict(list)

# task_required = set()
# task_optional = set()
# xc_required = set()
# xc_optional = set()
# common_required = set()
# common_optional = set()
# other_optional = set()

# for flag_type, flags in loadfn(SET_DIR / "incar_flags.yaml").items():
#     for flag, category in flags.items():
#         incar_flags[flag_type].append(flag)
#         if category is None:
#             other_optional.add(flag)
#         else:
#             try:
#                 getattr(sys.modules[__name__], category).add(flag)
#             except AttributeError:
#                 logger.error(f"{category} is not proper type.")
#                 raise


class IncarSettings(ABC):

    def __init__(self,
                 settings: dict):

        check_keys(settings, self.required_flags, self.optional_flags)
        self.settings = settings

    @property
    @abstractmethod
    def required_flags(self):
        pass

    @property
    @abstractmethod
    def optional_flags(self):
        pass


class TaskIncarSettings(IncarSettings):

    required_flags = TASK_REQUIRED_FLAGS
    optional_flags = TASK_OPTIONAL_FLAGS

    @classmethod
    def from_options(cls,
                     task: Task,
                     structure: Structure,
                     potcar: Potcar,
                     num_kpoints: int,
                     max_enmax: float,
                     is_magnetization: bool,
                     vbm_cbm: Optional[list],
                     npar_kpar: bool,
                     num_nodes: Optional[int],
                     encut: Optional[float],
                     structure_opt_encut_factor: float,
                     dos_step_size: float) -> "TaskIncarSettings":
        """Construct incar settings related to task with some options.

        Args: See ViseInputSet docstrings.

        Returns:
            TaskIncarSettings instance object.
        """
        band_gap = vbm_cbm[1] - vbm_cbm[0] if vbm_cbm else None
        if band_gap:
            logger.info(f"Band gap: {round(band_gap, 3)} eV.")
            is_band_gap = band_gap > BAND_GAP_CRITERION
        else:
            is_band_gap = None

        settings = \
            load_default_incar_settings(yaml_filename="task_incar_set.yaml",
                                        required_flags=TASK_REQUIRED_FLAGS,
                                        optional_flags=TASK_OPTIONAL_FLAGS,
                                        key_name=str(task))

        if task == Task.defect:
            settings["LREAL"] = "A"
        else:
            settings["LREAL"] = False

        # ISMEAR = -5 is acceptable when num_kpt >= 4 and the band gap exists.
        if is_band_gap is not True or num_kpoints < 4:
            settings["ISMEAR"] = 0

        if is_magnetization or task == Task.defect:
            settings["ISPIN"] = 2
            if hasattr(structure, "site_properties"):
                logger.error("Site properties are not inherited in the current "
                             "implementation.")

        if npar_kpar:
            kpar, npar = calc_npar_kpar(num_kpoints, num_nodes)
            settings["KPAR"] = kpar
            # now switch off NPAR
#            settings["NPAR"] = npar

        if not encut:
            if task.is_lattice_relax:
                encut = round(max_enmax * structure_opt_encut_factor, 3)
            else:
                encut = max_enmax
        settings["ENCUT"] = encut

        if task.is_plot_task:
            settings["NBANDS"] = nbands(structure.composition, potcar)

        if task.is_spectrum_task:
            if vbm_cbm:
                emin = ceil(vbm_cbm[0]) - 15 - dos_step_size
                emax = ceil(vbm_cbm[1]) + 15
            else:
                emin = -20 - dos_step_size
                emax = 20
            nedos = int(round((emax - emin) / dos_step_size, 0)) + 1
            settings.update({"EMIN": emin, "EMAX": emax, "NEDOS": nedos})

        if task == Task.dielectric_dfpt:
            settings["LEPSILON"] = True
        elif task == Task.dielectric_finite_field:
            settings["LCALCEPS"] = True
            settings["IBRION"] = 6
            settings["POTIM"] = 0.015
        elif task == Task.dielectric_function:
            settings["LOPTICS"] = True
            settings["CSHIFT"] = 0.01

        return cls(settings=settings)


class XcIncarSettings(IncarSettings):

    required_flags = XC_REQUIRED_FLAGS
    optional_flags = XC_OPTIONAL_FLAGS

    @classmethod
    def from_options(cls,
                     xc: Xc,
                     symbol_list: list,
                     factor: int,
                     aexx: Optional[float] = 0.25,
                     hubbard_u: Optional[bool] = None,
                     ldauu: Optional[dict] = None,
                     ldaul: Optional[dict] = None,
                     ldaul_set_name: Optional[str] = "default"
                     ) -> "XcIncarSettings":
        """Construct incar settings related to xc with some options.

        Args: See ViseInputSet docstrings

        Returns:
            XcIncarSettings instance object.
        """
        settings = \
            load_default_incar_settings(yaml_filename="xc_incar_set.yaml",
                                        required_flags=XC_REQUIRED_FLAGS,
                                        optional_flags=XC_OPTIONAL_FLAGS,
                                        key_name=str(xc))

        if ldauu:
            hubbard_u = True
        # By default Hubbard U is set for LDA or GGA.
        hubbard_u = xc.is_lda_or_gga if hubbard_u is None else hubbard_u
        ldauu = ldauu or {}
        ldaul = ldaul or {}

        if xc == Xc.pbesol:
            settings["GGA"] = "PS"
        elif xc == Xc.scan:
            settings["METAGGA"] = "SCAN"

        if hubbard_u:
            u_set = loadfn(SET_DIR / "u_parameter_set.yaml")
            ldauu_set = u_set["LDAUU"][ldaul_set_name]
            ldauu_set.update(ldauu)
            ldauu = [ldauu_set.get(el, 0) for el in symbol_list]

            if sum(ldauu) > 0:
                settings["LDAUU"] = ldauu
                settings.update({"LDAU": True, "LDAUTYPE": 2, "LDAUPRINT": 1})
                ldaul_set = u_set["LDAUL"][ldaul_set_name]
                ldaul_set.update(ldaul)
                settings["LDAUL"] = \
                    [ldaul_set.get(el, -1) for el in symbol_list]
                settings["LMAXMIX"] = \
                    6 if any([Element(el).Z > 56 for el in symbol_list]) else 4

        if xc.is_hybrid_functional:
            settings["AEXX"] = aexx
            if factor > 1:
                settings["NKRED"] = factor

        return cls(settings=settings)


class XcTaskIncarSettings(IncarSettings):

    required_flags = XC_TASK_REQUIRED_FLAGS
    optional_flags = XC_TASK_OPTIONAL_FLAGS

    @classmethod
    def from_options(cls) -> "XcTaskIncarSettings":
        """Construct incar settings related to task and xc with options.

        Note: When the GW related task and xc are implemented, this class will
              be used.

        Args: See ViseInputSet docstrings

        Returns:
            XcTaskIncarSettings instance object.
        """
        settings = {}
        return cls(settings)


class CommonIncarSettings(IncarSettings):

    required_flags = COMMON_REQUIRED_FLAGS
    optional_flags = COMMON_OPTIONAL_FLAGS

    @classmethod
    def from_options(cls,
                     potcar: Potcar,
                     composition: Composition,
                     charge: Optional[float] = None) -> "CommonIncarSettings":
        """Construct incar settings related to task and xc with options.

        Args: See ViseInputSet docstrings

        Returns:
            CommonIncarSettings instance object.
        """
        settings = {"NELM": 100,
                    "LASPH": True,
                    "LORBIT": 12,
                    "LCHARG": False,
                    "SIGMA": 0.1}
        # Structure charge is ignored.
        if charge:
            settings["NELECT"] = nelect(composition, potcar, charge)

        return cls(settings=settings)
