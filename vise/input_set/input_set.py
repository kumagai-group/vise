# -*- coding: utf-8 -*-

import logging
import os
import re
import warnings
from copy import deepcopy
from enum import Enum, unique
from os.path import join, isfile, getsize
from typing import Optional

import numpy as np
from monty.serialization import loadfn
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Potcar, Kpoints, Poscar
from pymatgen.io.vasp.sets import (
    get_vasprun_outcar, DictSet, get_structure_from_prev_run)
from vise.core.config import (
    KPT_DENSITY, ENCUT_FACTOR_STR_OPT, ANGLE_TOL, SYMPREC)
from vise.input_set.incar import ObaIncar
from vise.input_set.kpoints import make_kpoints, num_irreducible_kpoints
from vise.util.logger import get_logger
from vise.util.structure_handler import find_spglib_standard_primitive

logger = get_logger(__name__)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_DIR = os.path.join(MODULE_DIR, "sets")


def _load_oba_yaml_config(set_list_name: str,
                          potcar_yaml: str,
                          override_potcar_set: dict = None):
    """Load the yaml setting files for config and POTCAR list.

    Args:
        set_list_name (str): Set yaml file name
        potcar_yaml (str): Potcar list yaml file name
        override_potcar_set (dict): User specifying POTCAR set

    Return:
          config (dict):
    """
    config = loadfn(os.path.join(CONFIG_DIR, "%s.yaml" % set_list_name))
    potcar_list_file = os.path.join(CONFIG_DIR, "%s.yaml" % potcar_yaml)
    config["POTCAR"] = loadfn(potcar_list_file)
    if override_potcar_set:
        config["POTCAR"].update(override_potcar_set)
    return config


@unique
class Task(Enum):
    """ Supported tasks """
    single_point = "single_point"
    structure_opt = "structure_opt"
    phonon_disp = "phonon_disp"
    defect = "defect"
    band = "band"
    dos = "dos"
    dielectric = "dielectric"
    dielectric_function = "dielectric_function"
    gw_pre_calc1 = "gw_pre_calc1"
    gw_pre_calc2 = "gw_pre_calc2"

    def __str__(self):
        return self.name

    @classmethod
    def from_string(cls, s):
        for m in Task:
            if m.value == s:
                return m
        raise AttributeError("Task: " + str(s) + " is not proper.\n" +
                             "Supported Task:\n" + cls.name_list())

    @classmethod
    def name_list(cls):
        return ', '.join([e.value for e in cls])


@unique
class Xc(Enum):
    """ Supported exchange-correlation treatment. """
    pbe = "pbe"
    pbesol = "pbesol"
    lda = "lda"
    scan = "scan"
    pbe0 = "pbe0"
    hse = "hse"
    gw0 = "gw0"

    def __str__(self):
        return self.name

    @classmethod
    def from_string(cls, s):

        for m in Xc:
            if m.value == s:
                return m
        if s == "perdew-zunger81":
            return Xc.lda
        raise AttributeError("Xc: " + str(s) + " is not proper.\n" +
                             "Supported Xc:\n" + cls.name_list())

    @classmethod
    def name_list(cls):
        return ', '.join([e.value for e in cls])

    @property
    def require_wavefunctions(self):
        return True if self in BEYOND_DFT else False


LDA_OR_GGA = (Xc.pbe, Xc.pbesol, Xc.lda)
DFT_FUNCTIONAL = (Xc.pbe, Xc.pbesol, Xc.lda, Xc.scan)
HYBRID_FUNCTIONAL = (Xc.pbe0, Xc.hse)
BEYOND_GGA = (Xc.pbe0, Xc.hse, Xc.scan)
GW = (Xc.gw0, )
BEYOND_DFT = HYBRID_FUNCTIONAL + GW


class ObaSet(DictSet):
    """ Implementation of DictSet tuned for the research in Oba group.

    Special cares are as follows.
    1. The magnetization is assumed to be non-magnetic or ferromagnetic. Since
       the other magnetic configurations lower the host symmetry, one can use
       only kpt_mode="manual_set" and mostly kpts_shift=[0, 0, 0].
    2. Many options can be set in the **kwargs. See also the document of the
       DictSet and VaspInputSet, which are superclasses of the ObaSet.

    Below, the rules written in the sets file in pymatgen DictSet.

    The above are recommendations. The following are UNBREAKABLE rules:
    1. All input sets must take in a structure or list of structures as the
       first argument.
    2. user_incar_settings and user_kpoints_settings are absolute. Any new sets
       you implement must obey this. If a user wants to override your settings,
       you assume he knows what he is doing. Do not magically override user
       supplied settings. You can issue a warning if you think the user is
       wrong.
    3. All input sets must save all supplied args and kwargs as instance
       variables. E.g., self.my_arg = my_arg and self.kwargs = kwargs in
       the __init__. This ensures the as_dict and from_dict work correctly.
    """

    def __init__(self,
                 structure: Structure,
                 orig_structure: Structure,
                 config_dict: dict,
                 incar_settings: dict,
                 potcar: Potcar,
                 potcar_functional: str,
                 kpoints: Kpoints,
                 xc: str,
                 task: str,
                 is_magnetization: bool,
                 sg: int,
                 kpt_mode: str,
                 kpt_density: float,
                 band_ref_dist: float,
                 factor: int,
                 rough: bool,
                 cluster: bool,
                 only_even: bool,
                 use_structure_charge: bool,
                 files_to_move: Optional[dict] = None,
                 files_to_link: Optional[dict] = None,
                 files_to_transfer: Optional[dict] = None,
                 **kwargs):

        """
        Args:
            structure (Structure):
                The Structure to create Poscar.
            orig_structure (Structure):
                The original Structure.
            config_dict (dict):
            incar_settings (dict):
            potcar (Potcar):
            potcar_functional (str):
            kpoints (Kpoints):
            xc (str):
                A string representing xc defined in Xc.
            task (str):
                A string representing task defined in Task.
            is_magnetization (bool):
                Whether to have magnetization.
            sg (int):
                Space group number.
            kpt_mode (str):
                A string representing the k-point style that is used for the
                "primitive_uniform", and "manual_set" exist.
                See make_kpoints docstring.
            kpt_density (float):
                K-point density in Angstrom along each direction used for
                k-point uniform mesh.
            band_ref_dist (float):
                A reference target distance between neighboring k-points in the
                path, in units of 1/Angstrom. The actual value will be as close
                as possible to this value, to have an integer number of points
                in each path. Note that seekpath default is 0.025.
            factor (int):
                Factor to be multiplied to the k-points along all the
                directions. This parameter should be distinguished from the
                kpt_density as it can be compatible with NKRED = factor.
            rough (bool):
                Whether to ease the convergence criterion.
                Set EDIFF = 1e-5, EDIFFG = -0.04, NKRED = 2 in INCAR.
            cluster (bool):
                Whether the calculation is for a cluster including a molecule.
            only_even (bool):
                Whether to ceil the kpoints to be even numbers.
            use_structure_charge (bool):
            files_to_move (dict):
                Files to be moved to the directly where files are to be written.
            files_to_link (dict):
                Files to be linked via symbolic link with full path name to the
                directly where files are to be written.
            files_to_transfer (dict):
                Files to be copied to the directly where files are to be written.
            kwargs (dict):
        """
        self.incar_settings = incar_settings
        super().__init__(deepcopy(structure), config_dict,
                         files_to_transfer=files_to_transfer,
                         user_incar_settings=incar_settings,
                         sort_structure=False,
                         potcar_functional=potcar_functional,
                         use_structure_charge=use_structure_charge, **kwargs)

        # Only when structure is not changed, site_properties["magmom"] exists.
        self.orig_structure = orig_structure
        self._potcar = potcar
        self.xc = xc
        self.task = task
        self.is_magnetization = is_magnetization
        self.sg = sg
        self.kpt_mode = kpt_mode
        self.kpt_density = kpt_density
        self._kpoints = kpoints
        self.band_ref_dist = band_ref_dist
        self.factor = factor
        self.rough = rough
        self.cluster = cluster
        self.only_even = only_even
        self.files_to_move = files_to_move or {}
        self.files_to_link = files_to_link or {}

    @classmethod
    def make_input(cls,
                   structure: Structure,
                   incar_from_prev_calc: bool = False,
                   standardize_structure: bool = True,
                   xc: str = "pbe",
                   task: str = "structure_opt",
                   is_cluster: bool = False,
                   rough: bool = False,
                   factor: int = None,
                   encut_factor_str_opt: float = ENCUT_FACTOR_STR_OPT,
                   charge: float = None,
                   encut: float = None,
                   vbm_cbm: list = None,
                   band_gap: float = None,
                   is_magnetization: bool = False,
                   kpt_density: float = KPT_DENSITY,
                   sort_structure: bool = True,
                   weak_incar_settings: dict = None,
                   only_even: bool = True,
                   hubbard_u: bool = True,
                   ldauu: dict = None,
                   ldaul: dict = None,
                   ionic_contribution: bool = True,
                   kpt_mode: str = "primitive_uniform",
                   band_ref_dist: float = 0.03,
                   kpts_shift: list = None,
                   npar_kpar: bool = True,
                   num_cores_per_node: int = 36,
                   num_nodes: int = 1,
                   files_to_move: bool = None,
                   files_to_link: bool = None,
                   files_to_transfer: bool = None,
                   default_potcar: str = None,
                   override_potcar_set: dict = None,
                   symprec: float = SYMPREC,
                   angle_tolerance: float = ANGLE_TOL,
                   **kwargs):
        """
         --------------------------------------------------------------------
        POSCAR: structure, standardize_structure
        POTCAR: xc, user_potfar_yaml
         INCAR: incar_from_prev_calc, xc, task, rough, charge, encut, vbm_cbm,
                band_gap hubbard_u, ionic_contribution, num_cores_per_node,
                num_nodes, auto_kpar_npar, is_magnetization
        KPOINTS: task, kpt_mode, generate_kpoints, is_magnetization
        --------------------------------------------------------------------
        Args:
            structure (Structure):
                The Structure to create Poscar.
            incar_from_prev_calc (bool):
                Whether from previous calc or not.
            standardize_structure (bool):
                Whether to convert the structure to a standardized primitive.
            xc (str):
                A string representing exchange-correlation (xc) defined in Xc.
            task (str):
                A string representing task defined in Task.
            is_cluster (bool):
                Whether the calculation is for a cluster including a molecule.
            rough (bool):
                Whether to ease the convergence criterion.
                Set EDIFF = 1e-5, EDIFFG = -0.04, NKRED = 2 in INCAR.
            factor (int):
                Factor to be multiplied to the k-points along all the
                directions. This parameter should be distinguished from the
                kpt_density as it can be compatible with NKRED = factor.
            encut_factor_str_opt (float):
            charge (float):
                Charge state used for e.g., defect calculations.
            encut (float):
                Cutoff energy in eV.
                This argument is useful when performing a set of calculations
                for models with different elements yet their energies to be
                compared, e.g., defect calculations with impurity.
                For static calculations, 400 eV would be appropriate.
            vbm_cbm (list):
                The valence band maximum (vbm) and conduction band minimum (cbm)
                with [vbm, cbm] format. Used for determining the dos and
                absorption spectra region.
            band_gap (float):
                band gap. If it is larger than 0.2 eV, set ISMEAR = -5.
                It should be determined from previous calculations.
                If None, it is determined from vbm_cbm if given.
            is_magnetization (bool):
                Whether to have magnetization. If exists, set
                ISPIN = 2 and time_reversal = False, the latter of which is
                necessary for band structure calculations.
            kpt_density (float):
                K-point density in Angstrom along each direction used for
                k-point uniform mesh.
            sort_structure (bool):
                Whether to sort the elements using get_sorted_structure method
                of Structure class in pymatgen.
            weak_incar_settings (dict):
                Incar settings with the weakest hierarchy, meaning the
                flags are overwritten by other settings.
                This can be used for inheriting the modification by custodian.
            only_even (bool):
                Whether to ceil the kpoints to be even numbers.
            hubbard_u (bool):
                Whether to add Hubbard U parameter for lda/gga calculations.
                By default, it is switched off for hybrid functional and GW
                approximations. If one wants to add it, use user_incar_setting.
                U-parameters are also prefixed in element_specific_parameters.py
            ldauu (dict):
                Set when users want to modify U parameters from default.
                LDAUL also needs to be supplied.
                {"Ti": 4, "V": 3, "O": 5}
            ldaul (dict):
                {"Ti": 2, "V": 2, "O": 1}
            ionic_contribution (bool):
                Whether to calculate ionic contribution for dielectric constant.
                IBRION = 8 (DFPT w/ symmetry) is set.
            kpt_mode (str):
                A string representing the k-point style that is used for the
                make_kpoints function. See make_kpoints docstring.
                "band":
                  Kpoints with the band path will be returned based on the
                  seekpath program. The space group is analyzed and primitive
                  unitcell that must be used for the band structure calculation is
                  returned as well.
                "primitive_uniform":
                  Kpoints with uniform k-point sampling. The k-point sampling mesh
                  and centering are determined based on the standardized primitive
                  unitcell. Structure is also changed if not primitive.
                "manual_set":
                  Kpoints with uniform k-point sampling. The k-point sampling mesh
                  and centering are determined based on the given lattice. Note
                  that only when the angles are 90 degrees, the centering is
                  shifted along the perpendicular direction.
                  This mode is useful when calculating the supercells.
            band_ref_dist (float):
                This is used for band structure calculations.
                A reference target distance between neighboring k-points in the
                path, in units of 1/Angstrom. The actual value will be as close
                as possible to this value, to have an integer number of points
                in each path. Note that seekpath default is 0.025.
            kpts_shift (list):
                Kpoint shift in vasp definition.
            npar_kpar (bool):
                Whether to automatically set the NPAR and KPAR tags.
            num_cores_per_node (int):
                Number of cores per node used for automatic KPAR/NPAR setting.
            num_nodes (int):
                Number of nodes used for automatic KPAR/NPAR setting.
            files_to_move (dict):
                Files to be moved to the directory where files are to be written.
            files_to_link (dict):
                Files to be linked via symbolic link with full path name to the
                directory where files are to be written.
            files_to_transfer (dict):
                Files to be copied to the directory where files are to be
                written.
            default_potcar (str / None):
                Default potcar yaml-type file name without suffix.
                By default, use "default_POTCAR_list" but for GW approximations
                use "default_GW_POTCAR_list"
            override_potcar_set (dict /None):
                {"Mg": "Mg_pv"}
            kwargs (dict):
                keyword arguments. Some variables are quite useful as __init__
                and incar methods of DictSet superclass are called.
                Some examples are shown:
                    user_incar_settings (dict):
                    user_kpoints_settings (dict):
                    user_potcar_settings (dict):
                    files_to_transfer (dict):
        """
        structure = deepcopy(structure)

        kwargs = kwargs or {}
        files_to_move = files_to_move or {}
        files_to_link = files_to_link or {}
        files_to_transfer = files_to_transfer or {}

        if isinstance(task, str):
            task = Task.from_string(task)
        elif not isinstance(task, Task):
            raise TypeError("task needs to be str or Task object")

        xc = Xc.from_string(xc)

        if not default_potcar:
            if xc in GW or task in (Task.gw_pre_calc1, Task.gw_pre_calc2):
                default_potcar = "default_GW_POTCAR_list"
            else:
                default_potcar = "default_POTCAR_list"

        # read config_dict and override some values.
        config_dict = _load_oba_yaml_config("ObaRelaxSet", default_potcar,
                                            override_potcar_set)
        if ldauu:
            config_dict["INCAR"]["LDAUU"].update(ldauu)
        if ldaul:
            config_dict["INCAR"]["LDAUL"].update(ldaul)

        # Reset only INCAR. Be careful if the INCAR flags depend on the
        # POTCAR files, e.g., NBANDS.
        if incar_from_prev_calc is True:
            config_dict["INCAR"] = {}

        # ---- Handling kpt_mode depending on the task--------------------------
        if task == Task.defect:
            logger.warning("For defect task, kpt_mode is set to 'manual_set'.")
            kpt_mode = "manual_set"

        # ---- Compatibility check. ---------------
        # kpt_mode is overwritten if task is set to band and vice versa.
        if task is Task.band and kpt_mode != "band":
            logging.warning("For band task, kpt_mode is set to 'manual_set'.")
            kpt_mode = "band"

        if kpt_mode == "band" and task != Task.band:
            logging.warning("For band task, kpt_mode is set to 'manual_set'.")
            task = Task.band

        if band_gap:
            logger.info(f"Band gap : {band_gap} ")

        potcar_functional = "LDA" if xc == Xc.lda else "PBE"

        # ---- Handling of weak_incar_settings--------------------------------
        if weak_incar_settings:
            config_dict["INCAR"].update(weak_incar_settings)

        # remove the hubbard U parameter for hybrid functionals
        if not (xc in (Xc.lda, Xc.pbe, Xc.pbesol) and hubbard_u):
            for u_tag in ("LDAU", "LDAUU", "LDAUL", "LDAUTYPE", "LDAUPRINT"):
                config_dict["INCAR"].pop(u_tag, None)

        config_dict["INCAR"].pop("ICHARG", None)
        # ------------------------------------------

        # - charge
        if structure.charge:
            if charge:
                if charge != structure.charge:
                    raise ValueError(f"structure's charge {structure.charge} "
                                     f"and kwargs' charge {charge} are "
                                     f"simultaneously specified but different.")
            use_structure_charge = True
        elif charge:
            use_structure_charge = True
            structure._charge = charge
        else:
            use_structure_charge = False

        # structure
        original_structure = structure

        # Sorting the structure makes reading make_input difficult, and so
        # should be controlled only here.
        # VaspInputSet is also possible to sort elements but it should be
        # avoided for confusion.
        if sort_structure:
            structure = structure.get_sorted_structure()

        # For band structure calculations, structure is determined later.
        if task is not Task.band:

            if task is Task.defect:
                standardize_structure = False

            # Note that when applying find_spglib_standard_primitive, site
            # properties are removed.
            primitive_structure, is_structure_changed = \
                find_spglib_standard_primitive(
                    structure, symprec, angle_tolerance)

            if standardize_structure:
                structure = primitive_structure
            else:
                if is_structure_changed:
                    if kpt_mode != "manual_set":
                        logger.warning(
                            "Standardizaion is set to False and the given "
                            "structure is not a primitive cell. Thus, the "
                            "kpoint set is switched to manual_set.")
                        kpt_mode = "manual_set"
                else:
                    logger.info("The structure is a standardized primitive.")

        # Potcar
        potcar_symbols = []
        elements = structure.symbol_set
        if "user_potcar_settings" in kwargs:
            potcar_settings = kwargs["user_potcar_settings"]
        else:
            potcar_settings = config_dict["POTCAR"]

        if isinstance(potcar_settings[elements[-1]], dict):
            for el in elements:
                potcar_symbols.append(potcar_settings[el]['symbol']
                                      if el in potcar_settings else el)
        else:
            for el in elements:
                potcar_symbols.append(potcar_settings.get(el, el))

        potcar = Potcar(potcar_symbols, functional=potcar_functional)

        # if factor is None, set depending on the task.
        incar_settings, factor = \
            cls.make_incar_setting(structure, xc, task, rough, encut,
                                   vbm_cbm, potcar, band_gap, only_even,
                                   ionic_contribution, factor, encut_factor_str_opt)

        # - kpoints shift
        # To sample the band edges, Gamma-centered mesh is forced for dos and
        # dielectric function.
        # For the gw calcs, Gamma center is a must due to vasp implementation.
        # For tetrahedron method, Gamma-centered mesh is forced.

        if kpts_shift is None:
            if task in (Task.dos, Task.dielectric_function, Task.gw_pre_calc1,
                        Task.gw_pre_calc2) or \
                    incar_settings.get("ISMEAR", None) == -5:
                kpts_shift = (0, 0, 0)
        else:
            try:
                kpts_shift = tuple(kpts_shift)
            except:
                raise ValueError("The givenkpts_shift {} is not adequate."
                                 .format(kpts_shift))

        # is_cluster flag controls some INCAR tags.
        if is_cluster:
            incar_settings.update({"ISIF": 2})
            incar_settings.pop("ENCUT", None)
            incar_settings.pop("NKRED", None)
            kpts_shift = (0, 0, 0)
            kpt_density = 0.000001
            only_even = False

        # Refactor here.
        if incar_settings.get("ISIF", None) not in [3, 4, 5, 6] \
                and config_dict["INCAR"].get("ISIF", None) not in [3, 4, 5, 6]:
            incar_settings.pop("ENCUT", None)

        # The manual setting has the utmost priority.
        if "user_incar_settings" in kwargs:
            incar_settings.update(kwargs["user_incar_settings"])

        # Avoiding ICHARG = 11 is a must for hybrid functional.
        if incar_settings.get("LHFCALC", None) is True or \
                config_dict["INCAR"].get("LHFCALC", None) is True:
            incar_settings.pop("ICHARG", None)

        kwargs.pop("user_incar_settings", None)

        # - KPOINTS construction
        user_kpoints_settings = kwargs.get("user_kpoints_settings", {})
        logger.info("Magnetization: {}".format(is_magnetization))
        # Settings constructed by super().__init__ will be overwritten.
        # If a Kpoints object is supplied, it will be used.
        if not isinstance(user_kpoints_settings, Kpoints):
            kpt_settings = config_dict["KPOINTS"] or {}
            if user_kpoints_settings:
                kpt_settings.update(user_kpoints_settings)
            # Split of KPOINTS file for band structure calc is not supported,
            # so num_split_kpoints = 1.
            if "kpts_shift" in kpt_settings:
                kpts_shift = kpt_settings.pop("kpts_shift", None)
            kpoints, structure, sg = \
                make_kpoints(mode=kpt_mode,
                             structure=structure,
                             kpts_density=kpt_density,
                             only_even=only_even,
                             num_split_kpoints=1,
                             ref_distance=band_ref_dist,
                             kpts_shift=kpts_shift,
                             factor=factor,
                             symprec=SYMPREC,
                             angle_tolerance=ANGLE_TOL,
                             is_magnetization=is_magnetization, **kpt_settings)
        else:
            kpoints = user_kpoints_settings
            sg = None

        # calc num_kpoints
        num_kpoints = num_irreducible_kpoints(structure, kpoints,
                                              SYMPREC, ANGLE_TOL)
        kpoints.comment += ", # kpts: {}".format(num_kpoints)

        # Tetrahedron method is valid when num_kpoints >= 4
        if num_kpoints < 4 and incar_settings.get("ISMEAR", None) == -5:
            incar_settings.pop("ISMEAR")

        # - kpar and npar
        # If KPAR or NPAR exists in INCAR, this flow is switched off.
        # This routine should be followed after constructing KPOINTS to check
        # the irreducible number of kpoints.
        if npar_kpar and not (("KPAR" or "NPAR") in incar_settings):
            kpar, npar, num_kpoints = \
                calc_npar_kpar(num_kpoints, num_cores_per_node, num_nodes)
            incar_settings["KPAR"] = kpar

            # NPAR is now switched off as it causes some errors.
            # # Phonon calculations are not compatible with NPAR tag.
            # if incar_settings.get("IBRION", None) in [5, 6, 7, 8] or \
            #         (incar_settings.get("NBANDS", None) and
            #          xc in HYBRID_FUNCTIONAL):
            #     pass
            # else:
            #     incar_settings["NPAR"] = npar

        if cls.structure_changed(structure, original_structure):
            # The following files are useless when the lattice is changed.
            for f in ["CHGCAR", "WAVECAR", "WAVEDER"]:
                files_to_move.pop(f, None)
                files_to_link.pop(f, None)
                files_to_transfer.pop(f, None)

        return cls(structure, original_structure, config_dict, incar_settings,
                   potcar, potcar_functional, kpoints, xc,
                   task, is_magnetization, sg, kpt_mode, kpt_density,
                   band_ref_dist, factor, rough, is_cluster, only_even,
                   use_structure_charge, files_to_move, files_to_link,
                   files_to_transfer, **kwargs)

    # This is needed to overwrite the VaspInputSet potcar.
    @property
    def potcar(self):
        return self._potcar

    @property
    def kpoints(self):
        return self._kpoints

    @property
    def incar(self):
        # Do a lot of things with super().incar. See incar property in DictSet.
        settings = dict(self._config_dict["INCAR"])
        for k, v in self.user_incar_settings.items():
            if v is None:
                try:
                    del settings[k]
                except KeyError:
                    settings[k] = v
            else:
                settings[k] = v
        structure = self.structure
        incar = ObaIncar()
        comp = structure.composition
        elements = sorted([el for el in comp.elements if comp[el] > 0],
                          key=lambda e: e.X)
        most_electroneg = elements[-1].symbol
        poscar = Poscar(structure)
        hubbard_u = settings.get("LDAU", False)

        for k, v in settings.items():
            if k == "MAGMOM":
                mag = []
                for site in structure:
                    # YK: Add this
                    if isinstance(v, list):
                        mag = v
                    elif hasattr(site, 'magmom'):
                        mag.append(site.magmom)
                    elif hasattr(site.specie, 'spin'):
                        mag.append(site.specie.spin)
                    elif str(site.specie) in v:
                        mag.append(v.get(str(site.specie)))
                    else:
                        mag.append(v.get(site.specie.symbol, 0.6))
                incar[k] = mag
            elif k in ('LDAUU', 'LDAUJ', 'LDAUL'):
                if hubbard_u:
                    if hasattr(structure[0], k.lower()):
                        m = dict([(site.specie.symbol, getattr(site, k.lower()))
                                  for site in structure])
                        incar[k] = [m[sym] for sym in poscar.site_symbols]
                    # lookup specific LDAU if specified for most_electroneg atom
                    elif most_electroneg in v.keys() and \
                            isinstance(v[most_electroneg], dict):
                        incar[k] = [v[most_electroneg].get(sym, 0)
                                    for sym in poscar.site_symbols]
                    # else, use fallback LDAU value if it exists
                    else:
                        incar[k] = [v.get(sym, 0)
                                    if isinstance(v.get(sym, 0), (float, int))
                                    else 0 for sym in poscar.site_symbols]
            elif k.startswith("EDIFF") and k != "EDIFFG":
                if "EDIFF" not in settings and k == "EDIFF_PER_ATOM":
                    incar["EDIFF"] = float(v) * structure.num_sites
                else:
                    incar["EDIFF"] = float(settings["EDIFF"])
            else:
                incar[k] = v

        has_u = hubbard_u and sum(incar['LDAUU']) > 0
        if has_u:
            # modify LMAXMIX if LSDA+U and you have d or f electrons
            # note that if the user explicitly sets LMAXMIX in settings it will
            # override this logic.
            if 'LMAXMIX' not in settings.keys():
                # contains f-electrons
                if any([el.Z > 56 for el in structure.composition]):
                    incar['LMAXMIX'] = 6
                # contains d-electrons
                elif any([el.Z > 20 for el in structure.composition]):
                    incar['LMAXMIX'] = 4
        else:
            for key in list(incar.keys()):
                if key.startswith('LDAU'):
                    del incar[key]

        if self.constrain_total_magmom:
            nupdown = sum([mag if abs(mag) > 0.6 else 0
                           for mag in incar['MAGMOM']])
            incar['NUPDOWN'] = nupdown

        if self.use_structure_charge:
            incar["NELECT"] = self.nelect

        if np.product(self.kpoints.kpts) < 4 and incar.get("ISMEAR", 0) == -5:
            incar["ISMEAR"] = 0

        return incar

    @property
    def poscar(self):
        # The structure change should be checked here as it can happen at the
        # standardization and seekpath.
        if self.structure.symbol_set != self.orig_structure.symbol_set:
            logger.warning("CAUTION: The sequence of the species is changed.")
            logger.warning("Symbol set in the original structure")
            logger.warning(self.orig_structure.symbol_set)
            logger.warning("Symbol set in the generated structure")
            logger.warning(self.structure.symbol_set)

        if self.structure_changed(self.structure, self.orig_structure):
            logger.warning("CAUTION: The structure is changed.")
            logger.warning("Lattice of the original structure.")
            logger.warning(self.orig_structure.lattice)
            logger.warning("Lattice of the generated structure.")
            logger.warning(self.structure.lattice)

        return Poscar(self.structure)

    @staticmethod
    def structure_changed(structure1, structure2, symprec=SYMPREC):
        return not np.allclose(structure1.lattice.matrix,
                               structure2.lattice.matrix, atol=symprec)

    def write_input(self, output_dir, make_dir_if_not_present=True,
                    include_cif=False):
        """
        """

        super(DictSet, self).write_input(
            output_dir=output_dir,
            make_dir_if_not_present=make_dir_if_not_present,
            include_cif=include_cif)
        from monty.io import zopen
        import shutil
        for k, v in self.files_to_transfer.items():
            try:
                with zopen(v, "rb") as fin, zopen(os.path.join(output_dir, k),
                                                  "wb") as fout:
                    shutil.copyfileobj(fin, fout)
            except FileNotFoundError:
                warnings.warn("{} does not exist.".format(v))

        for k, v in self.files_to_move.items():
            try:
                shutil.move(v, os.path.join(output_dir, k))
            except FileNotFoundError:
                warnings.warn("{} does not exist.".format(v))

        for k, v in self.files_to_link.items():
            try:
                # full path is usually safer though it depends.
                os.symlink(os.path.realpath(v), os.path.join(output_dir, k))
            except FileNotFoundError:
                logger.warning("{} does not exist.".format(v))

    @classmethod
    def from_prev_calc(cls,
                       dirname,
                       parse_calc_results=True,
                       parse_magnetization=True,
                       standardize_structure=False,
                       sort_structure=True,
                       parse_potcar=True,
                       parse_incar=False,
                       parse_kpoints=False,
                       copied_file_names=None,
                       **kwargs):
        """
        Args:
            dirname (str):
                Directory name to be parsed.
            parse_calc_results (bool):
                Whether to parse the calculation results.
                Obtain the band gap and magnetization at all the atomic sites.
                The Structure is also constructed from vasprun.xml.
            parse_magnetization (bool):
                Whether to parse the calculation results.
            standardize_structure (bool):
            sort_structure (bool):
            parse_potcar (bool):
            parse_incar (bool):
            parse_kpoints (bool / str):
                True/False/"density"
            copied_file_names (dict):
                Key values show if files are
                moved (first letter="M")
                linked (first letter="L")
                copied (first letter="C")
                e.g., {"CHGCAR": "Moved", "WAVECAR": "L", "WAVEDER": "Co"}
            kwargs (dict):
        """
        kwargs = kwargs or {}
        kwargs["standardize_structure"] = standardize_structure
        kwargs["sort_structure"] = sort_structure

        if parse_calc_results:
            vasprun, outcar = get_vasprun_outcar(dirname)
            # When sym_prec=None, the primitive finder does not work.
            # see get_structure_from_prev_run for details.
            structure = get_structure_from_prev_run(vasprun, outcar)
            band_gap, cbm, vbm = vasprun.eigenvalue_band_properties[0:3]
            kwargs.update({"band_gap": band_gap})
            if band_gap > 0.1:
                kwargs.update({"vbm_cbm": [vbm, cbm]})

            if parse_magnetization and vasprun.is_spin:
                # magmom
                max_mag = max(structure.site_properties["magmom"])
                # max_mag = max([abs(v) for i in range(len(magnetization))
                #                for v in magnetization[i].values()])
                if "user_incar_settings" not in kwargs:
                    kwargs["user_incar_settings"] = {}
                # The criterion of the magnetic system is controlled below.
                # INCAR flags may be overwritten if INCAR is parsed.
                if max_mag > 0.05:
                    kwargs["user_incar_settings"].update({"ISPIN": 2})
                    kwargs["is_magnetization"] = True
                else:
                    kwargs["user_incar_settings"].pop("ISPIN", None)
                    kwargs["user_incar_settings"].pop("MAGMOM", None)
                    kwargs["is_magnetization"] = False
        else:
            if isfile(join(dirname, "CONTCAR")) and \
                    getsize(join(dirname, "CONTCAR")) > 0:
                structure = Structure.from_file(join(dirname, "CONTCAR"))
            elif isfile(join(dirname, "POSCAR-finish")) and \
                    getsize(join(dirname, "POSCAR-finish")) > 0:
                structure = Structure.from_file(join(dirname, "POSCAR-finish"))
            elif isfile(join(dirname, "POSCAR")) and \
                    getsize(join(dirname, "POSCAR")) > 0:
                structure = Structure.from_file(join(dirname, "POSCAR"))
            # This file name is only used in Oba group.
            else:
                raise OSError("CONTCAR or POSCAR does not exist.")

        if parse_potcar:
            try:
                potcar = Potcar.from_file(join(dirname, "POTCAR"))
                potcar_dict = {p.element: p.symbol for p in potcar}
                # Manually set potcar and xc are always prioritized.
                if "user_potcar_settings" not in kwargs:
                    kwargs["user_potcar_settings"] = potcar_dict

                if "xc" not in kwargs:
                    kwargs["xc"] = potcar.functional.lower()
            except:
                raise FileNotFoundError("POTCAR files is not found.")

        if parse_incar:
            kwargs["incar_from_prev_calc"] = True
            try:
                incar = ObaIncar.from_file(join(dirname, "INCAR"))
                weak_incar_settings = incar.as_dict()
                # Manually set potcar and xc are always prioritized.
                if "weak_incar_settings" in kwargs:
                    weak_incar_settings.update(kwargs["user_incar_settings"])
                for k in ('LDAUU', 'LDAUJ', 'LDAUL'):
                    if k in weak_incar_settings:
                        d = {}
                        if structure.symbol_set:
                            for i, e in enumerate(structure.symbol_set):
                                d[e] = weak_incar_settings[k][i]
                            weak_incar_settings[k] = d
                        else:
                            raise ValueError("Symbol set is required to parse "
                                             "LDAU tag in INCAR.")

                    # If structure doesn't have magmom property, set based on
                    # MAGMOM in INCAR.
                if 'MAGMOM' in weak_incar_settings and \
                        not hasattr(structure[0], 'magmom'):
                    structure.add_site_property('magmom',
                                                weak_incar_settings['MAGMOM'])
                else:
                    weak_incar_settings.pop("MAGMOM", None)

                if weak_incar_settings.get("ISPIN", None) == 2:
                    kwargs["is_magnetization"] = True

                kwargs["weak_incar_settings"] = weak_incar_settings

            except:
                raise FileNotFoundError("INCAR files is not found.")

        if parse_kpoints is True:
            try:
                kpoints = Kpoints.from_file(join(dirname, "KPOINTS"))
                kwargs["user_kpoints_settings"] = kpoints
            except:
                raise FileNotFoundError("KPOINTS file is not found.")
        # TODO: add kpt-density parser.
#        elif parse_kpoints == "density":

        kwargs["files_to_move"] = {}
        kwargs["files_to_link"] = {}
        kwargs["files_to_transfer"] = {}

        if copied_file_names:
            for filename in copied_file_names.keys():
                f = os.path.join(dirname, filename)
                if not isfile(f):
                    logger.warning("{} does not exist.".format(f))
                elif getsize(f) == 0:
                    logger.warning("{} is empty.".format(f))
                else:
                    if copied_file_names[filename][0].lower() == "m":
                        kwargs["files_to_move"][filename] = f
                    elif copied_file_names[filename][0].lower() == "l":
                        kwargs["files_to_link"][filename] = f
                    elif copied_file_names[filename][0].lower() == "c":
                        kwargs["files_to_transfer"][filename] = f
                    else:
                        logger.warning("{} option for {} file is invalid.".
                                      format(copied_file_names[filename][0],
                                             filename))

        return cls.make_input(structure, **kwargs)

    @classmethod
    def from_prev_calc_gw(cls, dirname, wavecar_transfer_style="Move",
                          waveder_transfer_style="Move", **kwargs):

        files = {"WAVECAR": wavecar_transfer_style,
                 "WAVEDER": waveder_transfer_style}

        return cls.from_prev_calc(dirname=dirname,
                                  parse_calc_results=True,
                                  parse_magnetization=True,
                                  standardize_structure=False,
                                  sort_structure=False,
                                  parse_potcar=True,
                                  parse_incar=True,
                                  parse_kpoints=True,
                                  copied_file_names=files,
                                  task="single_point",
                                  xc="gw0",
                                  **kwargs)

    @classmethod
    def from_prev_calc_gw_pre2(cls, dirname, chgcar_transfer_style="Move",
                               wavecar_transfer_style="Move", **kwargs):

        vasprun, _ = get_vasprun_outcar(dirname)

        if "LHFCALC" in vasprun.incar and vasprun.incar["LHFCALC"] is True:
            files = {"WAVECAR": wavecar_transfer_style}
        else:
            files = {"CHGCAR": chgcar_transfer_style}

        # need to calculate the number of the minimum G-vectors (plane waves)
        outcar = os.path.join(dirname, "OUTCAR")
        with open(outcar) as f:
            for line in f:
                if re.search("maximum and minimum number of plane-waves", line):
                    minimum_num_plane_waves = line.split()[10]
                    break
            if not minimum_num_plane_waves:
                raise ValueError("Number of plane wwaves is unavailable.")

            nbands = int(int(minimum_num_plane_waves) / 36) * 36

            if "user_incar_settings" not in kwargs:
                user_incar_settings = {"NBANDS": nbands}
            else:
                kwargs["user_incar_settings"]["NBANDS"] = nbands

        return cls.from_prev_calc(dirname=dirname,
                                  parse_calc_results=True,
                                  parse_magnetization=True,
                                  standardize_structure=False,
                                  sort_structure=False,
                                  parse_potcar=True,
                                  parse_incar=True,
                                  parse_kpoints=True,
                                  copied_file_names=files,
                                  user_incar_settings=user_incar_settings,
                                  task="gw_pre_calc2",
                                  **kwargs)


def calc_npar_kpar(num_kpoints, num_cores_per_node, num_nodes):
    """
    :param structure:
    :param kpoints:
    :param num_cores_per_node:
    :param num_nodes:
    :param symprec:
    :return:
    """

    kpar_set = {
        1:        [1, 1, 1],
        2:        [2, 2, 2],
        3:        [3, 3, 3],
        4:        [2, 4, 4],
        5:        [1, 2, 2],
        6:        [3, 6, 6],
        7:        [2, 2, 4],
        8:        [4, 8, 8],
        9:        [3, 3, 3],
        10:       [2, 2, 2],
        11:       [2, 2, 2],
        12:       [3, 4, 4],
        13:       [2, 2, 2],
        14:       [2, 2, 2],
        15:       [3, 2, 2],
        16:       [4, 4, 8],
        17:       [2, 4, 4],
        18:       [3, 6, 4],
        19:       [2, 4, 4],
        20:       [2, 4, 4],
        21:       [3, 3, 3],
        22:       [2, 4, 4],
        23:       [2, 4, 4],
        24:       [3, 6, 8],
        25:       [2, 4, 4],
        26:       [2, 4, 4],
        27:       [3, 4, 4],
        28:       [2, 4, 4],
        29:       [2, 4, 4],
        30:       [3, 6, 4],
        31:       [3, 4, 4],
        32:       [3, 4, 8],
        33:       [3, 4, 4],
        34:       [3, 4, 4],
        35:       [3, 4, 4],
        36:       [4, 8, 16],
        48:       [4, 8, 16],
        60:       [4, 8, 16],
        72:       [4, 8, 16],
        "others": [4, 8, 16]}

    if num_kpoints in kpar_set:
        num_kpt_key = num_kpoints
    else:
        num_kpt_key = "others"

    if num_nodes == 2:
        kpar = kpar_set[num_kpt_key][1]
    elif num_nodes == 4:
        kpar = kpar_set[num_kpt_key][2]
    else:
        kpar = kpar_set[num_kpt_key][0]

    if kpar == 1:
        npar = 2
    else:
        npar = 1

    return kpar, npar, num_kpoints

