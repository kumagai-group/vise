# -*- coding: utf-8 -*-

from copy import deepcopy
from typing import Optional

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import (
    VaspInputSet)
from vise.core.config import (
    KPT_DENSITY, ENCUT_FACTOR_STR_OPT, ANGLE_TOL, SYMPREC)
from vise.input_set.sets import (
    Task, Xc, TaskStructureKpoints, XcTaskPotcar, TaskIncarSettings,
    XcIncarSettings, XcTaskIncarSettings, CommonIncarSettings)
from vise.util.logger import get_logger

logger = get_logger(__name__)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InputSet(VaspInputSet):
    def __init__(self,
                 structure: Structure,
                 orig_structure: Structure,
                 incar,
                 potcar,
                 kpoints,
                 task,
                 xc,
                 kpt_mode,
                 kpt_density,
                 band_ref_dist: float,
                 factor,
                 encut,
                 structure_opt_encut_factor: float,
                 is_magnetization: bool,
                 hubbard_u: Optional[bool],
                 ldauu: Optional[dict],
                 ldaul: Optional[dict],
                 **kwargs):
        pass

    @classmethod
    def make_input(cls,
                   structure: Structure,
                   prev_set: "InputSet" = None,
                   xc: Xc = Xc.pbe,
                   task: Task = Task.structure_opt,
                   files_to_transfer: Optional[dict] = None,
                   user_incar_settings: dict = None,
                   #  ---- TaskStructureKpoints
                   standardize_structure: bool = True,
                   sort_structure: bool = True,
                   is_magnetization: bool = False,
                   kpt_mode: str = "primitive_uniform",
                   kpt_density: float = KPT_DENSITY,
                   kpts_shift: list = None,
                   only_even: bool = False,
                   band_ref_dist: float = 0.03,
                   factor: int = None,
                   symprec: float = SYMPREC,
                   angle_tolerance: float = ANGLE_TOL,
                   #  ---- XcTaskPotcar
                   potcar_set_name: str = "normal",
                   override_potcar_set: dict = None,
                   #  ---- TaskIncarSettings
                   band_gap: float = None,
                   vbm_cbm: list = None,
                   npar_kpar: bool = True,
                   num_cores=None,
                   encut: float = None,
                   structure_opt_encut_factor: float = ENCUT_FACTOR_STR_OPT,
                   #  ---- XcIncarSettings
                   aexx: Optional[float] = 0.25,
                   hubbard_u: bool = True,
                   ldauu: Optional[dict] = None,
                   ldaul: Optional[dict] = None,
                   ldaul_set_name: Optional[dict] = None,
                   #  ---- CommonSettings
                   charge: int = 0):

        num_cores = num_cores or [36, 1]
        files_to_transfer = files_to_transfer or {}
        user_incar_settings = user_incar_settings or {}

        use_prev_xc = True if prev_set and prev_set.xc == xc else False
        use_prev_task = True if (prev_set and prev_set.task == task and
                                 prev_set.structure == structure) else False

        if use_prev_task:
            task_str_kpt = prev_set.task_str_kpt
        else:
            task_str_kpt = TaskStructureKpoints.from_options(
                task=task,
                original_structure=structure,
                standardize_structure=standardize_structure,
                sort_structure=sort_structure,
                is_magnetization=is_magnetization,
                kpt_mode=kpt_mode,
                kpt_density=kpt_density,
                kpt_shift=kpts_shift,
                only_even=only_even,
                band_ref_dist=band_ref_dist,
                factor=factor,
                symprec=symprec,
                angle_tolerance=angle_tolerance)

        orig_matrix = structure.lattice.matrix
        matrix = task_str_kpt.structure.lattice.matrix
        if not np.allclose(orig_matrix, matrix, atol=symprec):
            # The following files are useless when lattice is changed.
            for f in ["CHGCAR", "WAVECAR"]:
                files_to_transfer.pop(f, None)

        if use_prev_task and use_prev_xc:
            xc_task_potcar = prev_set.xc_task_potcar
        else:
            xc_task_potcar = XcTaskPotcar.from_options(
                xc=xc,
                task=task,
                symbol_set=task_str_kpt.structure.symbol_set,
                potcar_set_name=potcar_set_name,
                override_potcar_set=override_potcar_set)

        if use_prev_task:
            task_settings = prev_set.task_settings
        else:
            task_settings = TaskIncarSettings.from_options(
                task=task,
                composition=task_str_kpt.structure.composition,
                potcar=xc_task_potcar.potcar,
                num_kpoints=task_str_kpt.num_kpts,
                max_enmax=xc_task_potcar.max_enmax,
                is_magnetization=is_magnetization,
                band_gap=band_gap,
                vbm_cbm=vbm_cbm,
                npar_kpar=npar_kpar,
                num_cores=num_cores,
                encut=encut,
                structure_opt_encut_factor=structure_opt_encut_factor)

        if use_prev_xc:
            xc_settings = prev_set.xc_settings
        else:
            xc_settings = XcIncarSettings.from_options(
                xc=xc,
                symbol_set=task_str_kpt.structure.symbol_set,
                factor=factor,
                aexx=aexx,
                hubbard_u=hubbard_u,
                ldauu=ldauu,
                ldaul=ldaul,
                ldaul_set_name=ldaul_set_name)

        if use_prev_task and use_prev_xc:
            xc_task_settings = prev_set.xc_task_settings
        else:
            xc_task_settings = XcTaskIncarSettings()

        common_settings = CommonIncarSettings.from_options(
            potcar=xc_task_potcar.potcar,
            composition=task_str_kpt.structure.composition,
            charge=charge)



