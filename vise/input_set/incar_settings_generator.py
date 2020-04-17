# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from math import ceil
from typing import Optional, Union, List, Dict

from pymatgen import Composition
from pymatgen.io.vasp.sets import Potcar

from vise.config import (
    DOS_STEP_SIZE, ENCUT_FACTOR_STR_OPT, DEFAULT_NUM_NODES,
    BAND_GAP_CRITERION)
from vise.input_set.datasets.dataset_util import num_bands, npar_kpar, LDAU
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger

logger = get_logger(__name__)


class IncarSettingsGenerator:
    def __init__(
            self,
            composition: Composition,
            symbol_list: list,
            num_kpts: int,
            num_kpt_multiplication_factor: int,
            potcar: Potcar,
            xc: Xc,
            task: Task,
            dos_step_size: float = DOS_STEP_SIZE,  # in eV
            charge: float = 0.0,  # in elementary charge as in vasp
            vbm_cbm: Optional[List[float]] = None,  # [vbm, cbm] in absolute eV
            aexx: float = 0.25,  # Exchange mixing parameter for hybrid func.
            ldauu: Optional[Dict[str, float]] = None,  # effective U, {"Ti": 4}
            ldaul: Optional[Dict[str, float]] = None,  # s:0, p:1, d:2 {"Ti": 2}
            set_hubbard_u: Optional[bool] = None,
            auto_npar_kpar: bool = True,  # Automatically set NPAR and KPAR.
            cutoff_energy: Optional[float] = None,  # ENCUT in eV.
            is_magnetization: bool = False,  # If the system is magnetic.
            num_nodes: int = DEFAULT_NUM_NODES,  # Used for KPAR and NPAR
            str_opt_encut_multi_factor: float = ENCUT_FACTOR_STR_OPT):

        self._composition = composition
        self._symbol_list = symbol_list
        self._num_kpt_multiplication_factor = num_kpt_multiplication_factor
        self._potcar = potcar
        self._num_kpts = num_kpts
        self._xc = xc
        self._task = task
        self._dos_step_size = dos_step_size
        self._charge = charge
        self._vbm_cbm = vbm_cbm
        self._aexx = aexx
        self._ldauu = ldauu
        self._ldaul = ldaul
        self._auto_npar_kpar = auto_npar_kpar
        self._cutoff_energy = cutoff_energy
        self._is_magnetization = is_magnetization
        self._num_nodes = num_nodes
        self._str_opt_encut_multi_factor = str_opt_encut_multi_factor

        # remove keys with None values
        self._settings = {
            "ALGO": self._xc.incar_algo,

            "PREC": self._task.incar_prec,
            "LREAL": self._task.incar_lreal,
            "EDIFF": self._task.incar_ediff,
            "ENCUT": self._encut,
            "ADDGRID": self._task.incar_addgrid_optional,
            "LASPH": True,
            "NELM": 100,

            "ISIF": self._task.incar_isif,
            "IBRION": self._task.incar_ibrion,
            "EDIFFG": self._task.incar_ediffg_optional,
            "NSW": self._task.incar_nsw,
            "POTIM": self._task.incar_potim_optional,

            "ISMEAR": self._ismear,
            "SIGMA": 0.1,

            "ISPIN": self._ispin,

            "LWAVE": self._xc.incar_lwave,
            "LCHARG": False,

            "NBANDS": self._nbands,
            "LORBIT": 12,

            "GGA": self._xc.incar_gga_optional,
            "METAGGA": self._xc.incar_metagga_optional,

            "NELECT": self._nelect,
        }

        if self._task.is_spectrum_task:
            self._set_spectrum_related_tags()
        if self._task.is_dielectric_task:
            self._set_dielectric_related_tags()
        if self._xc.is_hybrid_functional:
            self._set_hybrid_func_related_tags()
        if self._need_hubbard_u(set_hubbard_u):
            self._set_hubbard_u_related_tags()
        self._set_npar_kpar()
        self._remove_tags_with_none_values()

    def _need_hubbard_u(self, set_hubbard_u):
        if type(set_hubbard_u) == bool:
            return set_hubbard_u
        return (self._ldauu is not None) or self._xc.is_lda_or_gga

    @property
    def _encut(self) -> float:
        if self._cutoff_energy:
            return self._cutoff_energy
        max_enmax = max([p.enmax for p in self._potcar])

        if self._task.is_lattice_relaxed_task:
            max_enmax *= self._str_opt_encut_multi_factor
        return round(max_enmax, 3)

    @property
    def _ismear(self) -> Union[int, None]:
        # Tetrahedron method fails for irrep. NKPT<4 in vasp.
        if (is_band_gap(self._vbm_cbm)
                and self._num_kpts >= 4
                and self._task.is_spectrum_task):
            return -5
        return 0

    @property
    def _ispin(self) -> Optional[int]:
        if self._is_magnetization or (self._task.incar_ispin == 2):
            return 2
        return

    @property
    def _nbands(self) -> Union[int, None]:
        if self._task.is_plot_task:
            return num_bands(self._composition, self._potcar)
        return

    @property
    def _nelect(self) -> Union[float, None]:
        if self._charge:
            comp = self._composition.element_composition
            nelect = sum([comp[pt.element] * pt.ZVAL for pt in self._potcar])
            return nelect - self._charge
        return

    def _set_spectrum_related_tags(self) -> None:
        if self._vbm_cbm:
            emin = ceil(self._vbm_cbm[0]) - 15 - self._dos_step_size
            emax = ceil(self._vbm_cbm[1]) + 15
        else:
            emin = -20 - self._dos_step_size
            emax = 20
        self._settings["EMIN"] = emin
        self._settings["EMAX"] = emax
        num_steps = round((emax - emin) / self._dos_step_size)
        self._settings["NEDOS"] = num_steps + 1

    def _set_dielectric_related_tags(self) -> None:
        if self._task == Task.dielectric_dfpt:
            self._settings["LEPSILON"] = True
        elif self._task == Task.dielectric_finite_field:
            self._settings["LCALCEPS"] = True
            self._settings["IBRION"] = 6
            self._settings["POTIM"] = 0.015
        elif self._task == Task.dielectric_function:
            self._settings["LOPTICS"] = True
            self._settings["CSHIFT"] = 0.01

    def _set_hybrid_func_related_tags(self) -> None:
        self._settings["LHFCALC"] = True
        self._settings["PRECFOCK"] = "Fast"
        self._settings["TIME"] = 0.4
        self._settings["AEXX"] = self._aexx
        if self._num_kpt_multiplication_factor != 1:
            self._settings["NKRED"] = self._num_kpt_multiplication_factor
        if self._xc is Xc.hse:
            self._settings["HFSCREEN"] = 0.208

    def _set_hubbard_u_related_tags(self) -> None:
        ldau = LDAU(self._symbol_list, self._ldauu, self._ldaul)
        if ldau.is_ldau_needed:
            self._settings["LDAU"] = True
            self._settings["LDAUTYPE"] = 2
            self._settings["LMAXMIX"] = ldau.lmaxmix
            self._settings["LDAUPRINT"] = 1
            self._settings["LDAUU"] = ldau.ldauu
            self._settings["LDAUL"] = ldau.ldaul

    def _set_npar_kpar(self) -> None:
        if self._auto_npar_kpar:
            kpar, npar = npar_kpar(self._num_kpts, self._num_nodes)
            self._settings["KPAR"] = kpar
            # now switch off NPAR
            # self._settings["NPAR"] = npar

    def _remove_tags_with_none_values(self) -> None:
        for tag_name, value in self._settings.copy().items():
            if value is None:
                self._settings.pop(tag_name)

    @property
    def incar_settings(self) -> Dict[str, Union[str, bool, int, float]]:
        return self._settings


def is_band_gap(vbm_cbm: Optional[List[float]]) -> bool:
    band_gap = vbm_cbm[1] - vbm_cbm[0] if vbm_cbm else None
    if band_gap:
        logger.info(f"Band gap: {round(band_gap, 3)} eV.")
        return band_gap > BAND_GAP_CRITERION
    return False

