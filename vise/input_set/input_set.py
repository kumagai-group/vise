# -*- coding: utf-8 -*-

import json
import os
import re
import shutil
from copy import deepcopy
from pathlib import Path
from typing import Optional, List, Dict
from dataclasses import dataclass

import numpy as np
from monty.io import zopen
from monty.json import MontyEncoder
from monty.serialization import loadfn
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import Poscar, Potcar, Kpoints, VaspInput
from pymatgen.io.vasp.sets import (
    get_vasprun_outcar, get_structure_from_prev_run)
from pymatgen.util.serialization import pmg_serialize

from vise import __version__
from vise.analyzer.band_gap import band_gap_from_vasp
from vise.config import (
    DOS_STEP_SIZE, INSULATOR_KPT_DENSITY, ENCUT_FACTOR_STR_OPT, ANGLE_TOL,
    SYMMETRY_TOLERANCE, BAND_MESH_DISTANCE, DEFAULT_NUM_NODES,
    MAGNETIZATION_THRESHOLD)
from vise.input_set.incar import ViseIncar
from vise.input_set.settings_incar import (
    TaskIncarSettings, XcIncarSettings, XcTaskIncarSettings,
    CommonIncarSettings)
from vise.input_set.settings_potcar import XcTaskPotcar
from vise.input_set.make_kpoints import KpointsMode, MakeKpoints
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.util.logger import get_logger
from vise.util.structure_handler import symbol_list
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


@dataclass
class InputOptions:
    """Vise version of VaspInputSet

    Use the kwargs for the options that control vasp input set as the dict
    can be easily updated. Here, a list of options are listed as below.

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
    dos_step_size (float):
        Density of states step size.
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
    initial_structure: Structure

    xc: Xc = Xc.pbe
    task: Task = Task.structure_opt

    standardize_structure: bool = True

    kpt_mode: KpointsMode = KpointsMode.primitive
    kpt_density: float = INSULATOR_KPT_DENSITY
    kpt_shift: Optional[List[float]] = None
    manual_kpts: Optional[List[int]] = None
    only_even: bool = False
    band_ref_dist: float = BAND_MESH_DISTANCE
    factor: Optional[int] = None
    symprec: float = SYMMETRY_TOLERANCE
    angle_tolerance: float = ANGLE_TOL
    dos_step_size: float = DOS_STEP_SIZE
    charge: float = 0.0

    potcar_set_name: str = "normal"
    override_potcar_set: Optional[dict] = None
    vbm_cbm: Optional[List[float]] = None
    aexx: float = 0.25
    hubbard_u: Optional[bool] = None
    ldauu: Optional[Dict[str, float]] = None
    ldaul: Optional[Dict[str, float]] = None
    ldaul_set_name: str = "default"

    npar_kpar: bool = True
    encut: Optional[float] = None

    is_magnetization: bool = False
    num_nodes: int = DEFAULT_NUM_NODES
    structure_opt_encut_factor: float = ENCUT_FACTOR_STR_OPT

    # use dataclasses.InitVar
    # https://qiita.com/tag1216/items/13b032348c893667862a
    def parse_prev_calc(self):
        pass


class ViseInputSet:
    """Vise version of VaspInputSet"""
    # def __init__(self,
    #              structure: Structure,
    #              xc: Xc,
    #              task: Task,
    #              kpoints: Kpoints,
    #              potcar: Potcar,
    #              version: str,
    #              incar_settings: dict,
    #              files_to_transfer: Optional[FileTransfers],
    #              **kwargs):
    #     """
    #     Args:
    #         structure (Structure): The Structure to create inputs for.
    #         xc (Xc): Exchange-correlation (xc) defined in Xc.
    #         task (Task): Task defined in Task.
    #         kpoints (Kpoints): Kpoints class object
    #         potcar (Potcar): Potcar class object
    #         version (str): Vise version
    #         incar_settings (dict): INCAR settings generated by some options
    #         files_to_transfer (dict):
    #         kwargs (dict): OPTION arguments.
    #     """

        # self.structure = structure
        # self.xc = xc
        # self.task = task
        # self.kpoints = kpoints
        # self.potcar = potcar
        # self.version = version
        # self.incar_settings = incar_settings
        # self.files_to_transfer = files_to_transfer
        # self.kwargs = kwargs

    def __init__(self,
                 structure: Structure,
                 prev_calc_dir: Optional[str] = None,
                 user_incar_settings: Optional[dict] = None,
                 **options):
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

        Note1: Charge set in structure is ignored when determining NELECT. For
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
            files_to_transfer (dict):
                Keys are file names with absolute paths to be transferred.
                Values mean the transfer modes, where there are three modes,
                "c": copy file, "m": move file, "l": make symbolic link
                e.g., {"/..../CHGCAR": "c", "/..../WAVECAR": "l", ..}
            user_incar_settings (dict):
                User INCAR settings.
                e.g., {"NSW": 100, "LWAVE": True, "LASPH": "A", ..}
            options (dict):
                OPTION arguments.

        Return:
            ViseInputSet class object.
        """
        self.version = __version__
        self._create_input_options(options, structure)
        self.structure, self.kpoints, num_kpts, factor = self._make_str_kpt()
        self._check_lattice(structure)
        # self.files_to_transfer = deepcopy(files_to_transfer)
        # if lattice_changed:
        #     self.files_to_transfer.delete_file_transfers_w_keywords(
        #         ["CHGCAR$", "WAVECAR$", "WAVEDER"])

        xc_task_potcar = XcTaskPotcar.from_options(
            xc=self.opts.xc,
            symbol_list=symbol_list(self.structure),
            potcar_set_name=self.opts.potcar_set_name,
            override_potcar_set=self.opts.override_potcar_set)
        self.potcar = xc_task_potcar.potcar

        task_settings = TaskIncarSettings.from_options(
            task=self.opts.task,
            structure=self.structure,
            potcar=xc_task_potcar.potcar,
            num_kpoints=num_kpts,
            max_enmax=xc_task_potcar.max_enmax,
            is_magnetization=self.opts.is_magnetization,
            vbm_cbm=self.opts.vbm_cbm,
            npar_kpar=self.opts.npar_kpar,
            num_nodes=self.opts.num_nodes,
            encut=self.opts.encut,
            structure_opt_encut_factor=self.opts.structure_opt_encut_factor,
            dos_step_size=self.opts.dos_step_size)

        xc_settings = XcIncarSettings.from_options(
            xc=self.opts.xc,
            symbol_list=symbol_list(self.structure),
            factor=factor,
            aexx=self.opts.aexx,
            hubbard_u=self.opts.hubbard_u,
            ldauu=self.opts.ldauu,
            ldaul=self.opts.ldaul,
            ldaul_set_name=self.opts.ldaul_set_name)

        xc_task_settings = XcTaskIncarSettings.from_options()

        common_settings = CommonIncarSettings.from_options(
            potcar=xc_task_potcar.potcar,
            composition=self.structure.composition,
            charge=self.opts.charge)

        self.incar_settings = \
            {**task_settings.settings, **xc_settings.settings,
             **xc_task_settings.settings, **common_settings.settings}

        # user_incar_settings is top prioritized.
        self.incar_settings.update(user_incar_settings or {})

        # TODO: tweak the unfavorable combination of the input set.
        # e.g., Avoiding ICHARG = 11 is a must for hybrid functional.

    def _check_lattice(self, structure):
        orig_matrix = structure.lattice.matrix
        matrix = self.structure.lattice.matrix
        lattice_changed = \
            not np.allclose(orig_matrix, matrix, atol=self.opts.symprec)
        if lattice_changed and self.opts.charge != 0.0:
            raise ValueError("Structure is changed although charge is set.")

    # todo: move to InputOptions.
    def _make_str_kpt(self):
        sorted_structure = self._sort_species()
        structure = self._make_structure(sorted_structure)
        factor = self.opts.task.kpt_factor(self.opts.factor)
        make_kpoints = MakeKpoints(
            mode=self.opts.task.kpt_mode(self.opts.kpt_mode),
            structure=structure,
            kpt_density=self.opts.kpt_density,
            only_even=self.opts.task.only_even(self.opts.only_even),
            manual_kpts=self.opts.task.kpts(self.opts.manual_kpts),
            ref_distance=self.opts.band_ref_dist,
            kpt_shift=self.opts.task.kpt_shift(self.opts.kpt_shift),
            factor=factor,
            symprec=self.opts.symprec,
            angle_tolerance=self.opts.angle_tolerance,
            is_magnetization=self.opts.is_magnetization)
        make_kpoints.make_kpoints()

        return structure, make_kpoints.kpoints, make_kpoints.num_kpts, factor

    def _sort_species(self):
        structure = self.opts.initial_structure.get_sorted_structure()
        orig_symbol_list = symbol_list(self.opts.initial_structure)
        if symbol_list(structure) != orig_symbol_list:
            logger.warning(
                "CAUTION: The sequence of the species is changed. \n"
                f"Symbol set in the original structure {orig_symbol_list}\n"
                f"Symbol set in the generated structure "
                f"{symbol_list(structure)}")
        return structure

    def _make_structure(self, sorted_structure):
        symmetrizer = StructureSymmetrizer(sorted_structure,
                                           self.opts.symprec,
                                           self.opts.angle_tolerance)
        if self.opts.standardize_structure:
            if symmetrizer.is_lattice_changed:
                org_lat = self.opts.initial_structure.lattice.matrix
                primitive_lat = symmetrizer.primitive.lattice.matrix
                with np.printoptions(precision=3, suppress=True):
                    logger.warning(
                        "The structure is changed.\n"
                        f"Original lattice\n {org_lat} \n"
                        f"Generated lattice\n {primitive_lat} \n")
            return symmetrizer.primitive
        else:
            if symmetrizer.is_lattice_changed and \
                    self.opts.kpt_mode != "manual_set":
                raise ValueError(
                    "Although the given structure is not a primitive cell, "
                    "the k-point set is not set to manual_set.")
            else:
                logger.info("The structure is a standardized primitive.")
            return sorted_structure

    def _create_input_options(self, options, structure):
        input_options = InputOptions(initial_structure=structure, **options)
        input_options.parse_prev_calc()
        self.opts = input_options

    @property
    def incar(self):
        incar = ViseIncar.from_dict(self.incar_settings)
        return incar

    @property
    def poscar(self):
        return Poscar(self.structure)

    def vasp_input(self) -> VaspInput:
        return VaspInput(incar=self.incar, kpoints=self.kpoints,
                         poscar=self.poscar, potcar=self.potcar)

    def write_vasp_input_files(self, output_dir: str) -> None:
        self.vasp_input().write_input(output_dir, make_dir_if_not_present=True)
        # Overwrite POSCAR to increase the digit for frac coords.
        poscar_file = Path(output_dir) / "POSCAR"
        self.poscar.write_file(str(poscar_file), significant_figures=10)

    def create_input(self,
                     output_dir: str,
                     json_filename: str = "vise.json") -> None:
        self.write_vasp_input_files(output_dir)
        if self.files_to_transfer:
            self.files_to_transfer.transfer(Path(output_dir))
#        self.to_json_file(Path(output_dir) / json_filename)

    @pmg_serialize
    def as_dict(self, **kwargs):
        # Xc and Task objects must be converted to string for to_json_file as
        # Enum is not compatible with MSONable.
        d = {"structure":           self.structure,
             "kpoints":             self.kpoints,
             "potcar":              self.potcar.as_dict(),
             "version":             self.version,
             "incar_settings":      self.incar_settings,
             "files_to_transfer":   self.files_to_transfer,
             "kwargs":              self.opts}

        return d

    def __eq__(self, other):
        try:
            return self.as_dict() == other.as_dict()
        except AttributeError:
            return False

    # def __repr__(self):
    #     kwargs = "\n".join(
    #         [f"  {str(k):26s}: {str(v):26s}" for k, v in self.kwargs.items()])
    #     out = [f"Vise version: {self.version}",
    #            f"task: {self.task}",
    #            f"xc: {self.xc}",
    #            f"potcar: {self.potcar.symbols}",
    #            f"kwargs: ", f"{kwargs}"]

        # return "\n".join(out)

    # def __str__(self):
    #     return self.__repr__()

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

    # @classmethod
    # def load_json(cls, filename: str = "vise.json") -> "ViseInputSet":
    #     return loadfn(filename)

    def to_json_file(self, filename: str = "vise.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    # # TODO: to class
    # @classmethod
    # def from_prev_calc(cls,
    #                    dirname: str,
    #                    task: Optional[Task] = None,
    #                    xc: Optional[Xc] = None,
    #                    json_filename: str = "vise.json",
    #                    parse_calc_results: bool = True,
    #                    parse_incar: bool = True,
    #                    sort_structure: bool = True,
    #                    standardize_structure: bool = False,
    #                    files_to_transfer_dict: Optional[dict] = None,
    #                    user_incar_settings: Optional[dict] = None,
    #                    contcar_filename: str = "CONTCAR",
    #                    **kwargs) -> "ViseInputSet":
    #     """Constructor based on the previous calculations.
    #
    #     If the task and/or xc are not set, the previous ones are assumed.
    #
    #     Args:
    #         dirname (str):
    #             Directory name to be parsed.
    #         xc (Xc):
    #             Exchange-correlation (xc) defined in Xc.
    #         task (Task):
    #             Task defined in Task.
    #         json_filename (str):
    #             Json filename including ViseInputSet data.
    #         parse_calc_results (bool):
    #             Whether to parse the previous calculation results for band
    #             edges, magnetization.
    #         parse_incar (bool):
    #             Whether to parse the previous INCAR files, which could include
    #             other options that are not determined in vise.
    #         sort_structure (bool):
    #             Whether to sort the elements using get_sorted_structure method
    #             of Structure class in pymatgen.
    #         standardize_structure (bool):
    #             Whether to convert the structure to a standardized primitive.
    #         files_to_transfer_dict (dict):
    #             Keys are file names to be transferred in the directory.
    #             Values mean the transfer modes, where there are three modes,
    #             "c": copy file, "m": move file, "l": make symbolic link
    #             e.g., {"CHGCAR": "c", "WAVECAR": "l", ..}
    #         user_incar_settings (dict):
    #             User INCAR settings.
    #             e.g., {"NSW": 100, "LWAVE": True, "LASPH": "A", ..}
    #         contcar_filename (str):
    #             CONTCAR type filename. Used if the user change the name.
    #         kwargs:
    #             Other OPTION arguments.
    #
    #     Returns:
    #         ViseInputSet class object.
    #     """
    #     kwargs = deepcopy(kwargs) if kwargs else {}
    #     kwargs["sort_structure"] = sort_structure
    #     kwargs["standardize_structure"] = standardize_structure
    #
    #     path = Path(dirname)
    #
    #     input_set = loadfn(path / json_filename)
    #     cls.validate_vise_version_consistency(input_set.version, path)
    #
    #     if parse_calc_results:
    #         vasprun, outcar = get_vasprun_outcar(dirname)
    #         structure = get_structure_from_prev_run(vasprun, outcar)
    #         gap_properties = band_gap_from_vasp(vasprun, outcar)
    #
    #         if gap_properties[1]:
    #             _, vbm, cbm = gap_properties
    #             kwargs["vbm_cbm"] = [vbm["energy"], cbm["energy"]]
    #
    #         site_magmom = structure.site_properties.get("magmom", [0.0])
    #         kwargs["is_magnetization"] = \
    #             abs(max(site_magmom, key=abs)) > MAGNETIZATION_THRESHOLD
    #
    #     else:
    #         contcar = path / contcar_filename
    #         if contcar.is_file and contcar.stat().st_size:
    #             logger.info(f"{contcar} is parsed for structure.")
    #             structure = Structure.from_file(contcar)
    #         else:
    #             raise FileNotFoundError(f"{contcar} does not exist.")
    #
    #     if parse_incar:
    #         input_set.incar_settings = ViseIncar.from_file(path / "INCAR")
    #
    #     files_to_transfer = \
    #         FileTransfers.from_dict(files_to_transfer_dict, path)
    #
    #     return cls.make_input(structure=structure,
    #                           task=task or input_set.task,
    #                           xc=xc or input_set.xc,
    #                           prev_set=input_set,
    #                           files_to_transfer=files_to_transfer,
    #                           user_incar_settings=user_incar_settings,
    #                           **kwargs)
    #
    # @staticmethod
    # def validate_vise_version_consistency(version: str, path: Path) -> None:
    #     if version != __version__:
    #         logger.critical(f"The current vise version is {__version__}, "
    #                         f"while the previous version used in {path} is "
    #                         f"{version}. You must know what you're doing.")
    #
    #
    #
    #
