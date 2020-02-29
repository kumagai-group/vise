# -*- coding: utf-8 -*-

from pathlib import Path
from typing import Optional

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Kpoints
from vise.input_set.make_kpoints import MakeKpoints
from vise.input_set.task import Task, SPECTRA_TASK
from vise.util.logger import get_logger
from vise.util.structure_handler import find_spglib_primitive, get_symbol_list

logger = get_logger(__name__)
SET_DIR = Path(__file__).parent / "datasets"


class TaskStructureKpoints:

    def __init__(self,
                 structure: Structure,
                 kpoints: Kpoints,
                 is_structure_changed: bool,
                 sg: int,
                 num_kpts: int,
                 factor: int):
        self.structure = structure
        self.kpoints = kpoints
        self.is_structure_changed = is_structure_changed
        self.sg = sg
        self.num_kpts = num_kpts
        self.factor = factor

    @classmethod
    def from_options(cls,
                     task: Task,
                     original_structure: Structure,
                     standardize_structure: bool,
                     sort_structure: bool,
                     is_magnetization: bool,
                     kpt_mode: str,
                     kpt_density: float,
                     kpt_shift: list,
                     only_even: bool,
                     band_ref_dist: float,
                     factor: Optional[int],
                     symprec: float,
                     angle_tolerance: float):
        """Construct Structure and Kpoints from task and some options.

        When task and kpt_mode are not consistent e.g., task=Task.band,
        kpt_mode="manual", task is prioritized.

        Args: See ViseInputSet docstrings

        Returns:
            TaskStructureKpoints instance object
        """
        if sort_structure:
            structure = original_structure.get_sorted_structure()
            symbol_list = get_symbol_list(structure)
            orig_symbol_list = get_symbol_list(original_structure)
            if symbol_list != orig_symbol_list:
                logger.warning(
                    "CAUTION: The sequence of the species is changed. \n"
                    f"Symbol set in the original structure {symbol_list} \n"
                    f"Symbol set in the generated structure {orig_symbol_list}")
        else:
            structure = original_structure.copy()

        is_structure_changed = False
        manual_kpts = None
        if task == Task.defect:
            kpt_mode = "manual_set"
        elif task == Task.cluster_opt:
            kpt_mode = "manual_set"
            manual_kpts = [1, 1, 1]
            only_even = False
            kpt_shift = [0, 0, 0]
        elif task == Task.band:
            kpt_mode = "band"
        elif task == Task.phonon_force:
            kpt_mode = "manual_set"
            if kpt_shift != [0, 0, 0]:
                logger.warning("For phonon force calculations, Gamma centering "
                               "is forced for k-point sampling.")
            kpt_shift = [0, 0, 0]
        else:
            primitive_structure, is_structure_changed = \
                find_spglib_primitive(structure, symprec, angle_tolerance)

            if standardize_structure:
                org = original_structure.lattice.matrix
                primitive = primitive_structure.lattice.matrix
                if is_structure_changed:
                    with np.printoptions(precision=3, suppress=True):
                        logger.warning(
                            "CAUTION: The structure is changed.\n"
                            f"Original lattice\n {org} \n"
                            f"Generated lattice\n {primitive} \n")
                structure = primitive_structure
            else:
                if is_structure_changed and kpt_mode != "manual_set":
                    logger.warning(
                        "Standardizaion is set to False and the given "
                        "structure is not a primitive cell. Thus, the "
                        "kpoint set is switched to manual_set.")
                    kpt_mode = "manual_set"
                else:
                    logger.info("The structure is a standardized primitive.")

        # Gamma-center mesh is a must for GW calculations due to vasp
        # implementation and tetrahedron method, while is a strong recommend for
        # dos and dielectric function to sample the band edges.
        if task in SPECTRA_TASK and kpt_shift != [0, 0, 0]:
            logger.warning("Gamma centering is forced for k-point sampling.")
            kpt_shift = [0, 0, 0]

        if factor is None:
            if task == Task.dielectric_function:
                factor = 3
            elif task in (Task.dos, Task.dielectric_dfpt,
                          Task.dielectric_finite_field):
                factor = 2
            else:
                factor = 1

        kpoints = MakeKpoints(mode=kpt_mode,
                              structure=structure,
                              kpt_density=kpt_density,
                              only_even=only_even,
                              manual_kpts=manual_kpts,
                              ref_distance=band_ref_dist,
                              kpt_shift=kpt_shift,
                              factor=factor,
                              symprec=symprec,
                              angle_tolerance=angle_tolerance,
                              is_magnetization=is_magnetization)
        kpoints.make_kpoints()

        return cls(structure=structure,
                   kpoints=kpoints.kpoints,
                   is_structure_changed=kpoints.is_structure_changed,
                   sg=kpoints.sg,
                   num_kpts=kpoints.num_kpts,
                   factor=factor)


