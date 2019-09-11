# -*- coding: utf-8 -*-

import json
import os
import re
import shutil
from copy import deepcopy
from pathlib import Path
from typing import Optional, Dict

import numpy as np
from monty.io import zopen
from monty.json import MontyEncoder
from monty.serialization import loadfn
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import VaspInputSet, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import (
    get_vasprun_outcar, get_structure_from_prev_run)

from vise.analyzer.band_gap import band_gap_properties
from vise.input_set.vise_incar import ViseIncar
from vise.input_set.vise_incar_settings import (
    TaskIncarSettings, XcIncarSettings, XcTaskIncarSettings,
    CommonIncarSettings, OTHER_FLAGS)
from vise.input_set.xc import Xc
from vise.input_set.task import Task
from vise.input_set.vise_other_settings import (
    TaskStructureKpoints, XcTaskPotcar)
from vise.config import (
    KPT_DENSITY, ENCUT_FACTOR_STR_OPT, ANGLE_TOL, SYMMETRY_TOLERANCE)
from vise.util.logger import get_logger

logger = get_logger(__name__)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


__version__ = "0.0.1dev"


class ViseInputSet(VaspInputSet):
    """

    Options:
        xc (Xc):
            Exchange-correlation (xc) defined in Xc.
        task (Task):
            Task defined in Task.
        sort_structure (bool):
            Whether to sort the elements using get_sorted_structure method
            of Structure class in pymatgen.
        standardize_structure (bool):
            Whether to convert the structure to a standardized primitive.
        kpt_mode (str):
            A string representing the k-point style that is used for the
            make_kpoints function. Now, "band", "primitive_uniform" and
            "manual_set" are supported. See make_kpoints docstring.
        kpt_density (float):
            K-point density in Angstrom along each direction used for
            determining k-point uniform mesh.
        kpt_shift (list):
            K-point shift in vasp definition.
        only_even (bool):
            Whether to ceil the numbers of kpoints to be even.
        band_ref_dist (float):
            A reference target distance between neighboring k-points in the
            path, in units of 1/Angstrom. The actual value will be as close
            as possible to this value, to have an integer number of points
            in each path. Note that seekpath default is 0.025.
        factor (int):
            Factor to be multiplied to the k-points along all the
            directions. This parameter should be distinguished from the
            kpt_density as it sets NKRED = factor for hybrid calculations:w
        symprec (float):
            Distance precision in Angstrom used for symmetry analysis.
        angle_tolerance (float):
            Angle precision in degrees used for symmetry analysis.
        charge (float):
            Charge state used for e.g., defect calculations.
        default_potcar (str / None):
            Default potcar name defined in potcar_set.yaml.
            By default, use "normal", but one needs to use GW potcars for GW
            calculation, which will be supported in the future.
        override_potcar_set (dict /None):
            Dict representation for overriding the potcar set, e.g.,
            {"Mg": "Mg_pv", "O": "O_h"}
        vbm_cbm (list):
            The absolute valence band maximum (vbm) and conduction band minimum
            (cbm) with [vbm, cbm] format in eV, which are used for determining
            ISMEAR and the spectra region of the dos and absorption spectra.
        aexx (float):
            Exchange mixing parameter, namely AEXX in INCAR for hybrid
            functional calculations.
        hubbard_u (bool):
            Whether to add Hubbard U parameter.
            By default (None), it is switched on for LDA and GGA calculations,
            while switched off for hybrid functional and GW approximations.
            Default U-parameters are determined in u_parameter_set.yaml.
        ldauu (dict):
            Set when users want to modify U parameters from default. E.g.,
            {"Ti": 4, "V": 3, "O": 5}
        ldaul (dict):
            Set when users want to modify orbitals from default.
            If ldauu add elements that are not set in u_parameter_set.yaml,
            ldaul also needs to be supplied for these elements. E.g.,
            {"Ti": 2, "V": 2, "O": 1}
        ldaul_set_name:
            LDAUU and LDAUL default set name defined in u_parameter_set.yaml.
        npar_kpar (bool):
            Whether to automatically set the NPAR and KPAR tags.
        encut (float):
            Cutoff energy in eV.
            This argument is useful when performing a set of calculations
            for models with different elements yet their energies to be
            compared, e.g., defect calculations with impurity.
        is_magnetization (bool):
            Whether to have magnetization. If exists, set
            ISPIN = 2 and time_reversal = False, the latter of which is
            necessary for band structure calculations.
        num_cores_per_node (int):
            Numbers of cores per node and nodes used to determine KPAR and NPAR
            INCAR setting.
        structure_opt_encut_factor (float):
            Encut multiplied factor for structure optimization, where encut
            needs to be increased to reduce the Pulay Stress errors.
    """
    GENERAL_OPTIONS = {"xc": Xc.pbe,
                       "task": Task.structure_opt,
                       "sort_structure": True,
                       "standardize_structure": True}

    TASK_OPTIONS = {"kpt_mode": "primitive_uniform",
                    "kpt_density": KPT_DENSITY,
                    "kpt_shift": None,
                    "only_even": False,
                    "band_ref_dist": 0.03,
                    "factor": None,
                    "symprec": SYMMETRY_TOLERANCE,
                    "angle_tolerance": ANGLE_TOL,
                    "charge": 0}

    XC_OPTIONS = {"potcar_set_name": "normal",
                  "override_potcar_set": None,
                  "band_gap": None,
                  "vbm_cbm": None,
                  "aexx": 0.25,
                  "hubbard_u": None,
                  "ldauu": None,
                  "ldaul": None,
                  "ldaul_set_name": None}

    XC_TASK_OPTIONS = {"npar_kpar": True,
                       "encut": None}

    COMMON_OPTIONS = {"is_magnetization": False,
                      "num_cores": [36, 1],
                      "structure_opt_encut_factor": ENCUT_FACTOR_STR_OPT}

    OPTIONS = {**GENERAL_OPTIONS, **TASK_OPTIONS, **XC_OPTIONS,
               **XC_TASK_OPTIONS, **COMMON_OPTIONS}

    def __init__(self,
                 structure: Structure,
                 xc: Xc,
                 task: Task,
                 kpoints: Kpoints,
                 potcar: Potcar,
                 incar_settings: dict,
                 user_incar_settings: dict,
                 files_to_transfer: dict,
                 **kwargs):
        """

        Args:
            structure (Structure): The Structure to create inputs for.
            xc (Xc): Exchange-correlation (xc) defined in Xc.
            task (Task): Task defined in Task.
            kpoints (Kpoints): Kpoints class object
            potcar (Potcar): Potcar class object
            incar_settings (dict): INCAR settings generated by some options
            user_incar_settings (dict): User INCAR settings.
            files_to_transfer (dict):
            kwargs (dict): OPTION arguments.
        """

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
                   prev_set: "ViseInputSet" = None,
                   files_to_transfer: Optional[dict] = None,
                   user_incar_settings: Optional[dict] = None,
                   **kwargs) -> "ViseInputSet":
        """

        Args:
            structure (Structure):
                The Structure to create inputs for.
            prev_set (ViseInputSet):
            files_to_transfer:
            user_incar_settings (dict): User INCAR settings.
            kwargs (dict): OPTION arguments.
        :return:
        """

        files_to_transfer = files_to_transfer or {}
        user_incar_settings = user_incar_settings or {}

        # First, set default.
        opts = deepcopy(cls.OPTIONS)

        # Second, override with previous condition
        if prev_set:
            key_set = set(cls.COMMON_OPTIONS.keys())
            if prev_set.task == opts["task"]:
                key_set.update(cls.TASK_OPTIONS.keys())
            if prev_set.xc == opts["xc"]:
                key_set.update(cls.XC_OPTIONS.keys())
            if prev_set.task == opts["task"] and prev_set.xc == opts["xc"]:
                key_set.update(cls.XC_TASK_OPTIONS)

            for k in key_set:
                opts[k] = prev_set.kwargs[k]

        # Third, override with keyword arguments
        for k in kwargs:
            if "xc" in kwargs and type(kwargs["xc"]) == str:
                kwargs["xc"] = Xc.from_string(kwargs["xc"])
            if "task" in kwargs and type(kwargs["task"]) == str:
                kwargs["task"] = Task.from_string(kwargs["task"])
            if k not in opts.keys():
                logger.warning(f"Keyword {k} is not adequate for ViseInputSet.")
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

        input_set = ViseInputSet.load_json(json_filename)

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

