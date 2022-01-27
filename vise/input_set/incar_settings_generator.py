# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from math import ceil
from typing import Optional, Union, List, Dict

from pymatgen.core import Composition, Structure
from pymatgen.io.vasp.sets import Potcar

from vise.analyzer.band_edge_properties import is_band_gap
from vise.defaults import defaults
from vise.input_set.datasets.dataset_util import num_bands, npar_kpar, LDAU
from vise.input_set.fft_grids import vasp_grid
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger

logger = get_logger(__name__)


class IncarSettingsGenerator:
    def __init__(
            self,
            structure: Structure,
            symbol_list: list,
            num_kpts: int,
            num_kpt_multiplication_factor: int,
            potcar: Potcar,
            task: Task,
            xc: Xc,
            dos_step_size: float = defaults.dos_step_size,
            charge: float = 0.0,
            band_gap: Optional[float] = None,
            # [vbm, cbm] in absolute eV
            vbm_cbm: Optional[List[float]] = None,
            exchange_ratio: float = 0.25,
            set_hubbard_u: Optional[bool] = None,
            auto_npar_kpar: bool = True,
            cutoff_energy: Optional[float] = None,
            is_magnetization: bool = False,
            num_nodes_for_kpar: int = defaults.default_num_nodes,
            str_opt_encut_multi_factor: float = defaults.str_opt_encut_factor,
            multiples_for_grids: Optional[List[int]] = None):

        self._composition = structure.composition
        self._lattice = structure.lattice
        self._symbol_list = symbol_list
        self._num_kpt_multiplication_factor = num_kpt_multiplication_factor
        self._potcar = potcar
        self._num_kpts = num_kpts
        self._task = task
        self._xc = xc
        self._dos_step_size = dos_step_size
        self._charge = charge
        self._band_gap = band_gap
        self._vbm_cbm = vbm_cbm
        self._exchange_ratio = exchange_ratio
        self._auto_npar_kpar = auto_npar_kpar
        self._cutoff_energy = cutoff_energy
        self._is_magnetization = is_magnetization
        self._num_nodes_for_kpar = num_nodes_for_kpar
        self._str_opt_encut_multi_factor = str_opt_encut_multi_factor
        self._multiples_for_grids = multiples_for_grids

        self._incar_settings = {}
        self._set_incar_settings(set_hubbard_u)

    def _set_incar_settings(self, set_hubbard_u):
        self._set_default_settings()
        self._set_task_related_settings()
        self._set_xc_related_settings()
        self._set_options_related_settings()
        if self._task.is_spectrum_task:
            self._set_spectrum_related_settings()
        if self._task.is_dielectric:
            self._set_dielectric_related_settings()
        if self._xc.is_hybrid_functional:
            self._set_hybrid_func_related_settings()
        if self._need_hubbard_u(set_hubbard_u):
            self._set_hubbard_u_related_settings()
        if self._auto_npar_kpar:
            self._set_npar_kpar()
        self._remove_incar_settings_with_none_values()

    def _set_default_settings(self):
        self._incar_settings.update({
            "LASPH": True,
            "NELM": 100,
            "SIGMA": 0.1,
            "LCHARG": False,
        })

    def _set_task_related_settings(self):
        t = TaskIncarSettings(self._task)
        self._incar_settings.update({
            "PREC": t.prec,
            "LREAL": t.lreal,
            "EDIFF": t.ediff,
            "ADDGRID": t.addgrid_optional,
            "ISIF": t.isif,
            "IBRION": t.ibrion,
            "EDIFFG": t.ediffg_optional,
            "NSW": t.nsw,
            "POTIM": t.potim_optional,
            "LORBIT": t.lorbit,
        })

    def _set_xc_related_settings(self):
        x = XcIncarSettings(self._xc)
        self._incar_settings.update({
            "ALGO": x.algo,
            "LWAVE": x.lwave,
            "GGA": x.gga_optional,
            "METAGGA": x.metagga_optional,
        })

    def _set_options_related_settings(self):
        self._incar_settings.update({
            "ENCUT": self._encut,
            "ISMEAR": self._ismear,
            "ISPIN": self._ispin,
            "NBANDS": self._nbands,
            "NELECT": self._nelect,
        })
        if self._multiples_for_grids:
            grids = [vasp_grid(self._incar_settings["ENCUT"],
                               length,
                               self._incar_settings["PREC"])
                     for length in self._lattice.abc]
            m = self._multiples_for_grids
            ngx, ngy, ngz = (ceil(grids[i] / m[i]) * m[i] for i in range(3))
            self._incar_settings.update({"NGX": ngx})
            self._incar_settings.update({"NGY": ngy})
            self._incar_settings.update({"NGZ": ngz})

    def _need_hubbard_u(self, set_hubbard_u):
        if isinstance(set_hubbard_u, bool):
            return set_hubbard_u
        return False

    @property
    def _encut(self) -> float:
        if self._cutoff_energy:
            return self._cutoff_energy
        max_enmax = max([p.enmax for p in self._potcar])

        if self._task.is_lattice_relaxed:
            max_enmax *= self._str_opt_encut_multi_factor
        return round(max_enmax, 3)

    @property
    def _ismear(self) -> int:
        # Tetrahedron method fails for irrep. NKPT<4 in vasp.
        if self._task is Task.band:
            return 0
        elif is_band_gap(self._band_gap, self._vbm_cbm) and self._num_kpts >= 4:
            if self._task in (Task.dos, Task.dielectric_function):
                # tested -4 and -5 show the same results for spectra.
                return -4
            else:
                return -5
        return 0

    @property
    def _ispin(self) -> Optional[int]:
        if self._is_magnetization or self._task.need_spin:
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

    def _set_spectrum_related_settings(self) -> None:
        # Need one more step for vasp to remove weired huge value.
        if self._vbm_cbm:
            emin = ceil(self._vbm_cbm[0]) - 15 - self._dos_step_size
            emax = ceil(self._vbm_cbm[1]) + 15
        else:
            emin = -20 - self._dos_step_size
            emax = 20
        self._incar_settings["EMIN"] = emin
        self._incar_settings["EMAX"] = emax
        num_steps = round((emax - emin) / self._dos_step_size)
        self._incar_settings["NEDOS"] = num_steps + 1

    def _set_dielectric_related_settings(self) -> None:
        if self._task == Task.dielectric_dfpt:
            self._incar_settings["LEPSILON"] = True
        elif self._task == Task.dielectric_finite_field:
            self._incar_settings["LCALCEPS"] = True
#            self._incar_settings["IBRION"] = 6
            self._incar_settings["POTIM"] = 0.015
        elif self._task == Task.dielectric_function:
            self._incar_settings["LOPTICS"] = True
            # Since the real part of dielectric function depends on the vasp
            # version, we here adapt CSHIFT=0.0 and apply KK transformation
            # using vise.
            self._incar_settings["CSHIFT"] = 0.0

    def _set_hybrid_func_related_settings(self) -> None:
        self._incar_settings["LHFCALC"] = True
        self._incar_settings["PRECFOCK"] = "Fast"
        self._incar_settings["TIME"] = 0.4
        self._incar_settings["AEXX"] = self._exchange_ratio
        if self._num_kpt_multiplication_factor != 1:
            self._incar_settings["NKRED"] = self._num_kpt_multiplication_factor
        if self._xc is Xc.hse:
            self._incar_settings["HFSCREEN"] = 0.208

    def _set_hubbard_u_related_settings(self) -> None:
        ldau = LDAU(self._symbol_list)
        if ldau.is_ldau_needed:
            self._incar_settings["LDAU"] = True
            self._incar_settings["LDAUTYPE"] = 2
            self._incar_settings["LMAXMIX"] = ldau.lmaxmix
            self._incar_settings["LDAUPRINT"] = 1
            self._incar_settings["LDAUU"] = ldau.ldauu
            self._incar_settings["LDAUL"] = ldau.ldaul

    def _set_npar_kpar(self) -> None:
        kpar, npar = npar_kpar(self._num_kpts, self._num_nodes_for_kpar)
        if self._kpar_incompatible():
            return
        self._incar_settings["KPAR"] = kpar
        # now switch off NPAR
        # self._settings["NPAR"] = npar

    def _kpar_incompatible(self):
        if self._incar_settings.get("LELF", False):
            return True

    def _remove_incar_settings_with_none_values(self) -> None:
        for tag_name, value in self._incar_settings.copy().items():
            if value is None:
                self._incar_settings.pop(tag_name)

    @property
    def incar_settings(self) -> Dict[str, Union[str, bool, int, float]]:
        return self._incar_settings


class TaskIncarSettings:
    def __init__(self, task: Task):
        self._task = task

    @property
    def isif(self):
        if self._task.is_lattice_relaxed:
            return 3
        elif (self._task.is_atom_relaxed_lattice_fixed
              or self._task is Task.phonon_force):
            return 2
        elif self._task in (Task.band, Task.dos) or self._task.is_dielectric:
            return 0
        else:
            raise NotImplementedError

    # During dielectric_dfpt calculations, EDIFF is tightened automatically.
    @property
    def ediff(self):
        if self._task.is_tight_calc:
            return 1e-8
        elif self._task in (Task.structure_opt, Task.cluster_opt):
            return 1e-7
        elif self._task in (Task.dielectric_dfpt, Task.dielectric_finite_field):
            return 1e-6
        elif self._task in (Task.dielectric_function,
                            Task.band,
                            Task.dos,
                            Task.defect):
            return 1e-5
        elif self._task in (Task.structure_opt_rough,):
            return 1e-4
        else:
            raise NotImplementedError

    @property
    def ediffg_optional(self):
        if self._task.is_atom_relaxed:
            if self._task is Task.structure_opt_tight:
                return -0.001
            elif self._task in (Task.structure_opt, Task.cluster_opt):
                return -0.005
            elif self._task is Task.defect:
                return -0.03
            elif self._task is Task.structure_opt_rough:
                return -0.2
            else:
                raise NotImplementedError
        else:
            return

    @property
    def ibrion(self):
        return 8 if self._task is Task.dielectric_dfpt else 2

    @property
    def lorbit(self):
        if self._task is Task.phonon_force:
            return None
        return 10 if self._task is not Task.dos else 11

    @property
    def lreal(self):
        return "Auto" if self._task is Task.defect else False

    @property
    def prec(self):
        return "Accurate" if self._task.is_tight_calc else "Normal"

    @property
    def nsw(self):
        return 50 if self._task.is_atom_relaxed else 1

    @property
    def potim_optional(self):
        if self._task is Task.structure_opt_rough:
            return 0.1
        return

    @property
    def addgrid_optional(self):
        if self._task.is_tight_calc:
            return True
        return


class XcIncarSettings:
    def __init__(self, xc: Xc):
        self._xc = xc

    @property
    def algo(self):
        return "Damped" if self._xc.is_hybrid_functional else "Normal"

    @property
    def lwave(self):
        return True if self._xc.is_hybrid_functional else False

    @property
    def gga_optional(self):
        return "PS" if self._xc is Xc.pbesol else None

    @property
    def metagga_optional(self):
        return "SCAN" if self._xc is Xc.scan else None















