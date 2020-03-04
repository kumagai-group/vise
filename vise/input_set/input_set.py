# -*- coding: utf-8 -*-

import json
import os
import re
import shutil
from copy import deepcopy
from importlib import import_module
from pathlib import Path
from typing import Optional, Union

import numpy as np
from monty.io import zopen
from monty.json import MontyEncoder
from monty.serialization import loadfn
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import VaspInputSet, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import (
    get_vasprun_outcar, get_structure_from_prev_run)
from vise.analyzer.band_gap import band_gap_properties
from vise.config import (
    DOS_STEP_SIZE, KPT_DENSITY, ENCUT_FACTOR_STR_OPT, ANGLE_TOL,
    SYMMETRY_TOLERANCE, BAND_REF_DIST, DEFAULT_NUM_NODES)
from vise.input_set.incar import ViseIncar
from vise.input_set.settings_incar import (
    TaskIncarSettings, XcIncarSettings, XcTaskIncarSettings,
    CommonIncarSettings)
from vise.input_set.settings_potcar import XcTaskPotcar
from vise.input_set.settings_structure_kpoints import TaskStructureKpoints
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger
from vise.util.structure_handler import get_symbol_list
from vise import __version__


logger = get_logger(__name__)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class ViseInputSet(VaspInputSet):
    """Vise version of VaspInputSet

    Use the kwargs for the options that control vasp input set as the dict
    can be easily updated. Here, a list of options are listed as below.

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
        num_nodes (int):
            Numbers of nodes used to determine KPAR and NPAR INCAR setting.
        structure_opt_encut_factor (float):
            ENCUT multiplied factor for structure optimization, where encut
            needs to be increased to reduce the Pulay Stress errors.
    """
    GENERAL_OPTIONS = {"sort_structure": True,
                       "standardize_structure": True}

    TASK_OPTIONS = {"kpt_mode": "primitive_uniform",
                    "kpt_density": KPT_DENSITY,
                    "kpt_shift": None,
                    "only_even": False,
                    "band_ref_dist": BAND_REF_DIST,
                    "factor": None,
                    "symprec": SYMMETRY_TOLERANCE,
                    "angle_tolerance": ANGLE_TOL,
                    "dos_step_size": DOS_STEP_SIZE,
                    "charge": 0}

    XC_OPTIONS = {"potcar_set_name": "normal",
                  "override_potcar_set": None,
                  "vbm_cbm": None,
                  "aexx": 0.25,
                  "hubbard_u": None,
                  "ldauu": None,
                  "ldaul": None,
                  "ldaul_set_name": "default"}

    XC_TASK_OPTIONS = {"npar_kpar": True,
                       "encut": None}

    COMMON_OPTIONS = {"is_magnetization": False,
                      "num_nodes": DEFAULT_NUM_NODES,
                      "structure_opt_encut_factor": ENCUT_FACTOR_STR_OPT}

    ALL_OPTIONS = {**GENERAL_OPTIONS, **TASK_OPTIONS, **XC_OPTIONS,
                   **XC_TASK_OPTIONS, **COMMON_OPTIONS}

    def __init__(self,
                 structure: Structure,
                 xc: Xc,
                 task: Task,
                 kpoints: Kpoints,
                 potcar: Potcar,
                 version: str,
                 incar_settings: dict,
                 files_to_transfer: dict,
                 **kwargs):
        """
        Args:
            structure (Structure): The Structure to create inputs for.
            xc (Xc): Exchange-correlation (xc) defined in Xc.
            task (Task): Task defined in Task.
            kpoints (Kpoints): Kpoints class object
            potcar (Potcar): Potcar class object
            version (str): Vise version
            incar_settings (dict): INCAR settings generated by some options
            files_to_transfer (dict):
            kwargs (dict): OPTION arguments.
        """

        self.structure = structure
        self.xc = xc
        self.task = task
        self._kpoints = kpoints
        self._potcar = potcar
        self._version = version
        self.incar_settings = incar_settings
        self.files_to_transfer = files_to_transfer
        self.kwargs = kwargs

    @classmethod
    def make_input(cls,
                   structure: Structure,
                   task: Union[str, Task] = "structure_opt",
                   xc: Union[str, Xc] = "pbe",
                   prev_set: "ViseInputSet" = None,
                   abs_files_to_transfer: Optional[dict] = None,
                   user_incar_settings: Optional[dict] = None,
                   **kwargs) -> "ViseInputSet":
        """Construct ViseInputSet from some options including xc and task.

        To make a simple but practical rule for inheriting the previous
        calculation condition and effectively use its results, we adopt
        fail-safe rule.

        When, the prev_set is set, we compare task and xc between current
        input and prev_set, and inherit some options depending on whether
        task and/or xc are common. For instance, when we calculate the band
        structure with the HSE06 hybrid functional, it would be better to
        generate the WAVECAR file with e.g., PBE functional to reduce the
        computational cost. Then, the initial task was done with task=Task.band
        and xc=Xc.pbe with user_incar_setting={"LWAVE": True}, and next task is
        performed with task=Task.band, xc=Xc.hs. Then, task is common, so
        input set related to TaskStructureKpoints and TaskIncarSettings are
        inherited. Note that, CommonIncarSettings options are always inherited.

        Other notes are as follows.

        Note1: Charge set to structure is ignored when determining NELECT. For
               this purpose, charge option needs to be used.
        Note2: When the structure is changed via find_spglib_standard_primitive,
               all the site properties are removed.
        Note3: When different version of ViseInputSet with different defaults,
               the consistency of input set is destroyed, so the same version
               must be used when options are inherited.
        Note4: Other INCAR flags than those defined by TaskIncarSettings,
               XcIncarSettings, XcTaskIncarSettings, and CommonIncarSettings are
               not inherited. When some of them need to be inherited, they
               should be added to COMMON_OPTIONAL_FLAGS.
        Note5: user_incar_settings is not inherited from prev_set. One needs
               to explicitly specify it, again, if needed.
        Note6: Since default args are used in main function, set
               "structure_opt" and "pbe" strings as defaults.

        Args:
            structure (Structure):
                The Structure to create inputs for.
            task (str/Task):
                Task defined in Task.
            xc (str/Xc):
                Exchange-correlation (xc) defined in Xc.
            prev_set (ViseInputSet):
                Previous ViseInputSet.
            abs_files_to_transfer (dict):
                Keys are file names with absolute paths to be transferred.
                Values mean the transfer modes, where there are three modes,
                "c": copy file, "m": move file, "l": make symbolic link
                e.g., {"/..../CHGCAR": "c", "/..../WAVECAR": "l", ..}
            user_incar_settings (dict):
                User INCAR settings.
                e.g., {"NSW": 100, "LWAVE": True, "LASPH": "A", ..}
            kwargs (dict):
                OPTION arguments.

        Return:
            ViseInputSet class object.
        """
        abs_files_to_transfer = abs_files_to_transfer or {}
        user_incar_settings = user_incar_settings or {}

        if isinstance(task, str):
            task = Task.from_string(task)
        if isinstance(xc, str):
            xc = Xc.from_string(xc)

        # First, set default.
        opts = deepcopy(cls.ALL_OPTIONS)

        # Second, override with previous condition
        if prev_set:
            key_set = set(cls.COMMON_OPTIONS.keys())
            if prev_set.task == task:
                key_set.update(cls.TASK_OPTIONS.keys())
            if prev_set.xc == xc:
                key_set.update(cls.XC_OPTIONS.keys())
            if prev_set.task == task and prev_set.xc == xc:
                key_set.update(cls.XC_TASK_OPTIONS)

            for k in key_set:
                opts[k] = prev_set.kwargs[k]

        # Third, override with keyword arguments
        for k in kwargs:
            if k not in opts.keys():
                logger.warning(f"Keyword {k} is not adequate for ViseInputSet.")
        opts.update(kwargs)

        task_str_kpt = TaskStructureKpoints.from_options(
            task=task,
            original_structure=structure,
            standardize_structure=opts["standardize_structure"],
            sort_structure=opts["sort_structure"],
            is_magnetization=opts["is_magnetization"],
            kpt_mode=opts["kpt_mode"],
            kpt_density=opts["kpt_density"],
            kpt_shift=opts["kpt_shift"],
            only_even=opts["only_even"],
            band_ref_dist=opts["band_ref_dist"],
            factor=opts["factor"],
            symprec=opts["symprec"],
            angle_tolerance=opts["angle_tolerance"])

        orig_matrix = structure.lattice.matrix
        matrix = task_str_kpt.structure.lattice.matrix
        structure_changed = \
            not np.allclose(orig_matrix, matrix, atol=opts["symprec"])

        if structure_changed and opts["charge"] != 0.0:
            raise ValueError("Structure is changed but charge is set.")

        if structure_changed:
            # The following files are useless when lattice is changed.
            pattern = \
                re.compile("|".join([r"CHGCAR$", r"WAVECAR$", r"WAVEDER"]))
            for f in abs_files_to_transfer:
                if pattern.match(f):
                    abs_files_to_transfer.pop(f, None)

        # unique_justseen https://docs.python.org/ja/3/library/itertools.html
        # ["H", "H", "O", "O", "H"] -> ['H', 'O', 'H']
        symbol_list = get_symbol_list(task_str_kpt.structure)

        xc_task_potcar = XcTaskPotcar.from_options(
            xc=xc,
            symbol_list=symbol_list,
            potcar_set_name=opts["potcar_set_name"],
            override_potcar_set=opts["override_potcar_set"])

        task_settings = TaskIncarSettings.from_options(
            task=task,
            structure=task_str_kpt.structure,
            potcar=xc_task_potcar.potcar,
            num_kpoints=task_str_kpt.num_kpts,
            max_enmax=xc_task_potcar.max_enmax,
            is_magnetization=opts["is_magnetization"],
            vbm_cbm=opts["vbm_cbm"],
            npar_kpar=opts["npar_kpar"],
            num_nodes=opts["num_nodes"],
            encut=opts["encut"],
            structure_opt_encut_factor=opts["structure_opt_encut_factor"],
            dos_step_size=opts["dos_step_size"])

        xc_settings = XcIncarSettings.from_options(
            xc=xc,
            symbol_list=symbol_list,
            factor=task_str_kpt.factor,
            aexx=opts["aexx"],
            hubbard_u=opts["hubbard_u"],
            ldauu=opts["ldauu"],
            ldaul=opts["ldaul"],
            ldaul_set_name=opts["ldaul_set_name"])

        xc_task_settings = XcTaskIncarSettings.from_options()

        common_settings = CommonIncarSettings.from_options(
            potcar=xc_task_potcar.potcar,
            composition=task_str_kpt.structure.composition,
            charge=opts["charge"])

        incar_settings = \
            {**task_settings.settings, **xc_settings.settings,
             **xc_task_settings.settings, **common_settings.settings}

        # user_incar_settings is the top prioritized.
        incar_settings.update(user_incar_settings)

        # TODO: tweak the unfavorable combination of the input set.
        # e.g., Avoiding ICHARG = 11 is a must for hybrid functional.

        return cls(structure=task_str_kpt.structure,
                   task=task,
                   xc=xc,
                   kpoints=task_str_kpt.kpoints,
                   potcar=xc_task_potcar.potcar,
                   version=__version__,
                   incar_settings=incar_settings,
                   files_to_transfer=abs_files_to_transfer,
                   **opts)

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
                    include_cif: bool = False,
                    to_json_file: bool = True,
                    json_filename: str = "vise.json") -> None:
        """Write vasp input files and handle transferred files.

        Args:
            See docstring of write_input of VaspInputSet

        Return:
            None
        """
        super().write_input(output_dir=output_dir,
                            make_dir_if_not_present=make_dir_if_not_present,
                            include_cif=include_cif)

        # Overwrite only POSCAR to increase the digit for frac coords.
        poscar_file = Path(output_dir) / "POSCAR"
        self.poscar.write_file(str(poscar_file), significant_figures=10)

        out_dir = Path(output_dir).absolute()
        for key, value in self.files_to_transfer.items():
            try:
                # full path is usually safer for symbolic link.
                filepath = Path(key).absolute()
                name = filepath.name
                if value == "c" or value == "m":
                    with zopen(filepath, "rb") as fin, \
                            zopen((out_dir / name), "wb") as fout:
                        if value == "c":
                            shutil.copyfileobj(fin, fout)
                        elif value == "m":
                            shutil.move(fin, fout)
                elif value == "l":
                    os.symlink(filepath, out_dir / name)

            except FileNotFoundError:
                logger.warning(f"{key} does not exist.")

        if to_json_file:
            self.to_json_file(json_filename)

    def as_dict(self, **kwargs):

        # Xc and Task objects must be converted to string for to_json_file as
        # Enum is not compatible with MSONable.
        d = {"@module":             self.__class__.__module__,
             "@class":              self.__class__.__name__,
             "structure":           self.structure,
             "xc":                  str(self.xc),
             "task":                str(self.task),
             "kpoints":             self.kpoints,
             "potcar":              self.potcar.as_dict(),
             "version":             self._version,
             "incar_settings":      self.incar_settings,
             "files_to_transfer":   self.files_to_transfer,
             "kwargs":              self.kwargs}

        return d

    def __eq__(self, other):
        try:
            return self.as_dict() == other.as_dict()
        except AttributeError:
            return False

    def __repr__(self):
        kwargs = "\n".join(
            [f"  {str(k):26s}: {str(v):26s}" for k, v in self.kwargs.items()])
        out = [f"Vise version: {self._version}",
               f"task: {self.task}",
               f"xc: {self.xc}",
               f"potcar: {self.potcar.symbols}",
               f"kwargs: ", f"{kwargs}"]

        return "\n".join(out)

    def __str__(self):
        return self.__repr__()

    @property
    def version(self):
        return self._version

    @classmethod
    def from_dict(cls, d):
        # Programmatic access to enumeration members in Enum class.
        structure = d["structure"]
        if isinstance(structure, dict):
            structure = Structure.from_dict(structure)

        return cls(structure=structure,
                   xc=Xc.from_string(d["xc"]),
                   task=Task.from_string(d["task"]),
                   kpoints=d["kpoints"],
                   potcar=Potcar.from_dict(d["potcar"]),
                   version=d["version"],
                   incar_settings=d["incar_settings"],
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
                       dirname: str,
                       task: Optional[Task] = None,
                       xc: Optional[Xc] = None,
                       json_filename: str = "vise.json",
                       parse_calc_results: bool = True,
                       parse_incar: bool = True,
                       sort_structure: bool = True,
                       standardize_structure: bool = False,
                       files_to_transfer: Optional[dict] = None,
                       user_incar_settings: Optional[dict] = None,
                       contcar_filename: str = "CONTCAR",
                       **kwargs) -> "ViseInputSet":
        """Constructor based on the previous calculations.

        If the task and/or xc are not set, the previous ones are assumed.

        Args:
            dirname (str):
                Directory name to be parsed.
            xc (Xc):
                Exchange-correlation (xc) defined in Xc.
            task (Task):
                Task defined in Task.
            json_filename (str):
                Json filename including ViseInputSet data.
            parse_calc_results (bool):
                Whether to parse the previous calculation results for band
                edges, magnetization.
            parse_incar (bool):
                Whether to parse the previous INCAR files, which could include
                other options that are not determined in vise.
            sort_structure (bool):
                Whether to sort the elements using get_sorted_structure method
                of Structure class in pymatgen.
            standardize_structure (bool):
                Whether to convert the structure to a standardized primitive.
            files_to_transfer (dict):
                Keys are file names to be transferred in the directory.
                Values mean the transfer modes, where there are three modes,
                "c": copy file, "m": move file, "l": make symbolic link
                e.g., {"CHGCAR": "c", "WAVECAR": "l", ..}
            user_incar_settings (dict):
                User INCAR settings.
                e.g., {"NSW": 100, "LWAVE": True, "LASPH": "A", ..}
            contcar_filename (str):
                CONTCAR type filename. Used if the user change the name.
            kwargs:
                Other OPTION arguments.

        Returns:
            ViseInputSet class object.
        """
        kwargs = kwargs or {}
        kwargs["sort_structure"] = sort_structure
        kwargs["standardize_structure"] = standardize_structure

        path = Path(dirname)

        input_set = ViseInputSet.load_json(path / json_filename)
        if input_set.version != __version__:
            logger.critical(f"The current vise version is {__version__}, "
                            f"while the previous version used in {path} is "
                            f"{input_set.version}. "
                            f"You must know what you're doing.")

        if parse_calc_results:
            vasprun, outcar = get_vasprun_outcar(dirname)
            structure = get_structure_from_prev_run(vasprun, outcar)

            gap_properties = band_gap_properties(vasprun, outcar)

            # For metall,
            # gap_properties = ({'energy': 0.0,
            #                    'direct': False,
            #                    'transition': None}, None, None)
            if gap_properties[2]:
                _, vbm, cbm = gap_properties
                kwargs["vbm_cbm"] = [vbm["energy"], cbm["energy"]]

            if "magmom" in structure.site_properties.keys():
                abs_max_mag = abs(max(structure.site_properties["magmom"],
                                      key=abs))
            else:
                abs_max_mag = 0
            if vasprun.is_spin and abs_max_mag > 0.05:
                kwargs["is_magnetization"] = True

        else:
            contcar = path / contcar_filename
            poscar = path / "POSCAR"
            if contcar.is_file and os.stat(contcar).st_size != 0:
                logger.info(f"{contcar} is parsed for structure.")
                structure = Structure.from_file(contcar)
            elif poscar.is_file and os.stat(poscar).st_size != 0:
                logger.warning(f"{poscar} is parsed for structure.")
                structure = Structure.from_file(poscar)
            else:
                raise FileNotFoundError("CONTCAR or POSCAR does not exist.")

        if parse_incar:
            input_set.incar_settings = ViseIncar.from_file(path / "INCAR")

        if files_to_transfer:
            abs_files_to_transfer = {}
            for filename, value in files_to_transfer.items():
                f = path / filename
                if not f.is_file():
                    logger.warning(f"{f} does not exist.")
                elif os.stat(f).st_size == 0:
                    logger.warning(f"{f} is empty.")
                elif value[0].lower() not in ["l", "c", "m"]:
                    logger.warning(f"{value} option for {filename} is invalid.")
                else:
                    abs_files_to_transfer[str(f)] = value[0].lower()
            files_to_transfer = abs_files_to_transfer

        return cls.make_input(structure=structure,
                              task=task or input_set.task,
                              xc=xc or input_set.xc,
                              prev_set=input_set,
                              abs_files_to_transfer=files_to_transfer,
                              user_incar_settings=user_incar_settings,
                              **kwargs)


