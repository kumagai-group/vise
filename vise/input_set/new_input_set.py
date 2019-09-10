# -*- coding: utf-8 -*-

import json
import os
import re
import shutil
from copy import deepcopy
from pathlib import Path
from typing import Optional

import numpy as np
from monty.io import zopen
from monty.json import MontyEncoder
from monty.serialization import loadfn
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import VaspInputSet, Poscar, Potcar
from pymatgen.io.vasp.sets import (
    get_vasprun_outcar, get_structure_from_prev_run)
from vise.analyzer.band_gap import band_gap_properties
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


__version__ = "0.0.1dev"


class InputSet(VaspInputSet):
    def __init__(self,
                 structure: Structure,
                 xc: Xc,
                 task: Task,
                 kpoints,
                 potcar,
                 incar_settings,
                 user_incar_settings,
                 files_to_transfer,
                 **kwargs):
        self.structure = structure
        self.xc = xc
        self.task = task
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

        # Don't forget to add keys below when add to opts.
        # Inherit those keys from prev_set if task is the same.
        task_keys = {"kpt_mode", "kpt_density", "kpt_shift", "only_even",
                     "band_ref_dist", "factor", "symprec", "angle_tolerance",
                     "charge"}
        # Inherit those keys from prev_set if xc is the same.
        xc_keys = {"potcar_set_name", "override_potcar_set", "band_gap",
                   "vbm_cbm", "aexx", "hubbard_u", "ldauu", "ldaul",
                   "ldaul_set_name"}
        # Inherit those keys from prev_set if xc and task are the same.
        # encut depends on the potcar which depends on xc.
        xc_task_keys = {"npar_kpar", "encut"}
        # Inherit those keys from prev_set anytime.
        common_keys = {"is_magnetization", "num_cores",
                       "structure_opt_encut_factor"}

        # Only for check.
        check = (set(opts.keys())
                 - (task_keys | xc_keys | xc_task_keys | common_keys))
        if check != {"standardize_structure", "sort_structure", "task", "xc"}:
            raise KeyError("Keys are not consistent.")

        # Second, override with previous condition
        if prev_set:
            key_set = common_keys
            if prev_set.task == opts["task"]:
                key_set.update(task_keys)
            if prev_set.xc == opts["xc"]:
                key_set.update(xc_keys)
            if prev_set.task == opts["task"] and prev_set.xc == opts["xc"]:
                key_set.update(xc_task_keys)
            for k in key_set:
                opts[k] = prev_set.kwargs[k]

        # Third, override with keyword arguments
        for k in kwargs:
            if "xc" in kwargs and type(kwargs["xc"]) == str:
                kwargs["xc"] = Xc.from_string(kwargs["xc"])
            if "task" in kwargs and type(kwargs["task"]) == str:
                kwargs["task"] = Xc.from_string(kwargs["task"])
            if k not in opts.keys():
                logger.warning(f"Keyword {k} is not adequate in vise InputSet.")
        opts.update(kwargs)

        task_str_kpt = TaskStructureKpoints.from_options(
            original_structure=structure, **opts)
        opts["factor"] = task_str_kpt.factor

        orig_matrix = structure.lattice.matrix
        matrix = task_str_kpt.structure.lattice.matrix
        structure_changed = \
            not np.allclose(orig_matrix, matrix, atol=opts["symprec"])

        if structure_changed and opts["charge"] != 0:
            raise ValueError("Structure is changed but charge is set.")

        if structure_changed:
            # The following files are useless when lattice is changed.
            pattern = \
                re.compile("|".join([r"CHGCAR$", r"WAVECAR$", r"WAVEDER"]))
            for f in files_to_transfer:
                if pattern.match(f):
                    files_to_transfer.pop(f, None)

        xc_task_potcar = XcTaskPotcar.from_options(
            symbol_set=task_str_kpt.structure.symbol_set, **opts)

        task_settings = TaskIncarSettings.from_options(
            structure=task_str_kpt.structure,
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

        # Note: user_incar_settings is not inherited from prev_set.
        #       This is a specification of vise.
        # user_incar_settings is the top prioritized.
        incar_settings.update(user_incar_settings)

        print(incar_settings)

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
        for key, value in self.files_to_transfer.items():
            try:
                # full path is usually safer for symbolic link.
                filepath = Path(key).absolute()
                name = filepath.name
                with zopen(filepath, "rb") as fin, \
                        zopen((out_dir / name), "wb") as fout:
                    if value == "c":
                        shutil.copyfileobj(fin, fout)
                    elif value == "m":
                        shutil.move(fin, fout)
                    elif value == "l":
                        os.symlink(fin, fout)

            except FileNotFoundError:
                logger.warning(f"{key} does not exist.")

    def as_dict(self, verbosity=2):
        # Xc and Task objects must be converted to string for to_json_file as
        # Enum is not compatible with MSONable.
        xc = str(self.xc)
        task = str(self.task)
        potcar = self.potcar.as_dict()

        d = {"@module":             self.__class__.__module__,
             "@class":              self.__class__.__name__,
             "structure":           self.structure,
             "xc":                  xc,
             "task":                task,
             "kpoints":             self.kpoints,
             "potcar":              potcar,
             "incar_settings":      self.incar_settings,
             "user_incar_settings": self.user_incar_settings,
             "files_to_transfer":   self.files_to_transfer,
             "kwargs":              self.kwargs}

        return d

    @classmethod
    def from_dict(cls, d):
        # Programmatic access to enumeration members in Enum class.
        xc = Xc.from_string(d["xc"])
        task = Task.from_string(d["task"])
        potcar = Potcar.from_dict(d["potcar"])

        structure = d["structure"]
        if isinstance(structure, dict):
            structure = Structure.from_dict(structure)

        return cls(structure=structure,
                   xc=xc,
                   task=task,
                   kpoints=d["kpoints"],
                   potcar=potcar,
                   incar_settings=d["incar_settings"],
                   user_incar_settings=d["user_incar_settings"],
                   files_to_transfer=d["files_to_transfer"],
                   **d["kwargs"])

    @classmethod
    def load_json(cls, filename: str = "vise.json"):
        return loadfn(filename)

    def to_json_file(self, filename: str = "vise.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def from_prev_calc(cls,
                       dirname,
                       json_filename: str = "vise.json",
                       parse_calc_results: bool = True,
                       sort_structure: bool = True,
                       standardize_structure: bool = False,
                       files_to_transfer: Optional[dict] = None,
                       user_incar_settings: Optional[dict] = None,
                       contcar_filename: str = "CONTCAR",
                       **kwargs):

        kwargs = kwargs or {}
        directory_path = Path(dirname)

        input_set = InputSet.load_json(json_filename)

        if parse_calc_results:
            vasprun, outcar = get_vasprun_outcar(dirname)
            structure = get_structure_from_prev_run(vasprun, outcar)

            gap_properties = band_gap_properties(vasprun)
            if gap_properties:
                band_gap, cbm, vbm = gap_properties
                kwargs["band_gap"] = band_gap
                kwargs["vbm_cbm"] = [vbm["energy"], cbm["energy"]]

            abs_max_mag = abs(max(structure.site_properties["magmom"], key=abs))
            if vasprun.is_spin and abs_max_mag > 0.05:
                kwargs["is_magnetization"] = True

        else:
            contcar = directory_path / contcar_filename
            poscar = directory_path / "POSCAR"
            if contcar.is_file and os.stat(contcar).st_size == 0:
                logger.info(f"{contcar} is parsed for structure.")
                structure = Structure.from_file(contcar)
            elif poscar.is_file and os.stat(poscar).st_size == 0:
                logger.warning(f"{poscar} is parsed for structure.")
                structure = Structure.from_file(poscar)
            else:
                raise FileNotFoundError("CONTCAR or POSCAR does not exist.")

        if files_to_transfer:
            abs_files_to_transfer = {}
            for filename, value in files_to_transfer.items():
                f = directory_path / filename
                if not f.is_file:
                    logger.warning(f"{f} does not exist.")
                elif os.stat(f).st_size == 0:
                    logger.warning(f"{f} is empty.")
                else:
                    if value[0].lower() not in ["l", "c", "m"]:
                        logger.warning(f"{value} option for {filename} is "
                                       f"invalid.")
                abs_files_to_transfer[str(f)] = value[0].lower()
            files_to_transfer = abs_files_to_transfer

        return cls.make_input(structure=structure,
                              prev_set=input_set,
                              files_to_transfer=files_to_transfer,
                              user_incar_settings=user_incar_settings,
                              sort_structure=sort_structure,
                              standardize_structure=standardize_structure)

