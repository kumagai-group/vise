# -*- coding: utf-8 -*-

from math import ceil
from pathlib import Path
from typing import Optional

from monty.serialization import loadfn

from pymatgen import Structure, Composition, Element
from pymatgen.io.vasp import Potcar

from vise.config import BAND_GAP_CRITERION
from vise.input_set.settings_util import (
    load_default_incar_settings, check_keys, nelect, nbands, calc_npar_kpar)
from vise.input_set.task import (
    LATTICE_RELAX_TASK, PLOT_TASK, SPECTRA_TASK, Task)
from vise.input_set.xc import Xc, LDA_OR_GGA, HYBRID_FUNCTIONAL
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

ALL_FLAGS = set(sum(loadfn(SET_DIR / "incar_flags.yaml").values(), []))

OTHER_FLAGS = ALL_FLAGS - (TASK_FLAGS | XC_FLAGS | XC_TASK_FLAGS | COMMON_FLAGS)


# TODO: Implement spin-orbit calc.

class TaskIncarSettings:

    def __init__(self,
                 settings: dict):

        check_keys(settings, TASK_REQUIRED_FLAGS, TASK_OPTIONAL_FLAGS)
        self.settings = settings

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

        settings["LREAL"] = "A" if task == Task.defect else False

        # -5 is acceptable for num_kpt >= 4 and is_band_gap = True
        if is_band_gap is not True or num_kpoints < 4:
            settings["ISMEAR"] = 0

        if task == Task.defect or is_magnetization:
            settings["ISPIN"] = 2
            if hasattr(structure, "site_properties"):
                logger.critical("MAGMOM is not inherited in the current "
                                "implementation.")
#                if "magmom" in structure.site_properties:
#                    settings["MAGMOM"] = structure.site_properties["magmom"]

        if npar_kpar:
            kpar, npar = calc_npar_kpar(num_kpoints, num_nodes)
            settings["KPAR"] = kpar
            # now switch off NPAR
#            settings["NPAR"] = npar

        if not encut:
            if task in LATTICE_RELAX_TASK:
                encut = round(max_enmax * structure_opt_encut_factor, 3)
            else:
                encut = max_enmax
        settings["ENCUT"] = encut

        if task in PLOT_TASK:
            settings["NBANDS"] = nbands(structure.composition, potcar)

        if task in SPECTRA_TASK:
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


class XcIncarSettings:

    def __init__(self, settings: dict):
        check_keys(settings, XC_REQUIRED_FLAGS, XC_OPTIONAL_FLAGS)
        self.settings = settings

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
        hubbard_u = xc in LDA_OR_GGA if hubbard_u is None else hubbard_u
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


class CommonIncarSettings:
    def __init__(self, settings: dict):
        check_keys(settings, COMMON_REQUIRED_FLAGS, COMMON_OPTIONAL_FLAGS)
        self.settings = settings

    @classmethod
    def from_options(cls,
                     potcar: Potcar,
                     composition: Composition,
                     charge: Optional[int] = None) -> "CommonIncarSettings":
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
