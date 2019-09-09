# -*- coding: utf-8 -*-

from copy import deepcopy
from typing import Optional
from monty.io import zopen
import shutil
from pathlib import Path
import numpy as np
import os
import re
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import VaspInputSet, Poscar
from vise.core.config import (
    KPT_DENSITY, ENCUT_FACTOR_STR_OPT, ANGLE_TOL, SYMPREC)
from vise.input_set.incar import ViseIncar
from vise.input_set.sets import (
    Task, Xc, TaskStructureKpoints, XcTaskPotcar, TaskIncarSettings,
    XcIncarSettings, XcTaskIncarSettings, CommonIncarSettings, OTHER_FLAGS)
from vise.util.logger import get_logger

logger = get_logger(__name__)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InputSet(VaspInputSet):
    def __init__(self,
                 structure: Structure,
                 kpoints,
                 potcar,
                 incar_settings,
                 user_incar_settings,
                 files_to_transfer,
                 **kwargs):
        self.structure = structure
        self._kpoints = kpoints
        self._potcar = potcar
        self.incar_settings = incar_settings
        self.user_incar_settings = user_incar_settings
        self.files_to_transfer = files_to_transfer
        self.kwargs = kwargs

    @classmethod
    def make_input(cls,
                   structure: Structure,
                   prev_set: "InputSet" = None,
                   files_to_transfer: Optional[dict] = None,
                   user_incar_settings: Optional[dict] = None,
                   **kwargs) -> "InputSet":

        files_to_transfer = files_to_transfer or {}
        user_incar_settings = user_incar_settings or {}

        # First, set default.
        opts = {"xc": Xc.pbe,
                "task": Task.structure_opt,
                "sort_structure": True,
                "standardize_structure": True,
                "is_magnetization": False,
                "kpt_mode": "primitive_uniform",
                "kpt_density": KPT_DENSITY,
                "kpt_shift": None,
                "only_even": False,
                "band_ref_dist": 0.03,
                "factor": None,
                "symprec": SYMPREC,
                "angle_tolerance": ANGLE_TOL,
                "potcar_set_name": "normal",
                "override_potcar_set": None,
                "band_gap": None,
                "vbm_cbm": None,
                "npar_kpar": True,
                "num_cores": [36, 1],
                "encut": None,
                "structure_opt_encut_factor": ENCUT_FACTOR_STR_OPT,
                "aexx": 0.25,
                "hubbard_u": None,
                "ldauu": None,
                "ldaul": None,
                "ldaul_set_name": None,
                "charge": 0}

        # Second, override with previous condition
        if prev_set:
            opts.update(prev_set.kwargs)
        # Third, override with keyword arguments
        for k in kwargs:
            if k not in opts.keys():
                logger.warning(f"Keyword {k} is not adequate in vise InputSet.")
        opts.update(kwargs)

        task_str_kpt = TaskStructureKpoints.from_options(
            original_structure=structure, **opts)

        orig_matrix = structure.lattice.matrix
        matrix = task_str_kpt.structure.lattice.matrix
        if not np.allclose(orig_matrix, matrix, atol=opts["symprec"]):
            # The following files are useless when lattice is changed.
            pattern = \
                re.compile("|".join([r"CHGCAR$", r"WAVECAR$", r"WAVEDER"]))
            for f in files_to_transfer:
                if pattern.match(f):
                    files_to_transfer.pop(f, None)

        xc_task_potcar = XcTaskPotcar.from_options(
            symbol_set=task_str_kpt.structure.symbol_set, **opts)

        task_settings = TaskIncarSettings.from_options(
            composition=task_str_kpt.structure.composition,
            potcar=xc_task_potcar.potcar,
            num_kpoints=task_str_kpt.num_kpts,
            max_enmax=xc_task_potcar.max_enmax, **opts)

        xc_settings = XcIncarSettings.from_options(
            symbol_set=task_str_kpt.structure.symbol_set, **opts)

        xc_task_settings = XcTaskIncarSettings.from_options()

        common_settings = CommonIncarSettings.from_options(
            potcar=xc_task_potcar.potcar,
            composition=task_str_kpt.structure.composition, **opts)

        incar_settings = deepcopy(task_settings.settings)
        incar_settings.update(xc_settings.settings)
        incar_settings.update(xc_task_settings.settings)
        incar_settings.update(common_settings.settings)

        if prev_set:
            prev_other_settings = \
                {k: v for k, v in prev_set.incar_settings.items()
                 if k in OTHER_FLAGS}
            incar_settings.update(prev_other_settings)
            # Previous user_incar_settings is inherited, but current one is
            # prioritized.
            d = deepcopy(prev_set.user_incar_settings)
            user_incar_settings = d.update(user_incar_settings)

        incar_settings.update(user_incar_settings)

        return cls(structure=task_str_kpt.structure,
                   kpoints=task_str_kpt.kpoints,
                   potcar=xc_task_potcar.potcar,
                   incar_settings=incar_settings,
                   user_incar_settings=user_incar_settings,
                   files_to_transfer=files_to_transfer, **opts)

    @property
    def potcar(self):
        return self._potcar

    @property
    def kpoints(self):
        return self._kpoints

    @property
    def incar(self):
        incar = ViseIncar.from_dict(self.incar_settings)
        return incar

    @property
    def poscar(self):
        return Poscar(self.structure)

    def write_input(self,
                    output_dir: str,
                    make_dir_if_not_present: bool = True,
                    include_cif: bool = False) -> None:

        super().write_input(output_dir=output_dir,
                            make_dir_if_not_present=make_dir_if_not_present,
                            include_cif=include_cif)

        out_dir = Path(output_dir).absolute()
        for k, v in self.files_to_transfer.items():
            try:
                # full path is usually safer though it depends.
                filepath = Path(k).absolute()
                name = filepath.name
                with zopen(str(filepath), "rb") as fin, \
                        zopen(str(out_dir / name), "wb") as fout:
                    if v[0] == "c":
                        shutil.copyfileobj(fin, fout)
                    elif v[0] == "m":
                        shutil.move(fin, fout)
                    elif v[0] == "l":
                        os.symlink(fin, fout)

            except FileNotFoundError:
                logger.warning(f"{k} does not exist.")


