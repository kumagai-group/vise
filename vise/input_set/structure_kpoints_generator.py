# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from math import ceil
from typing import Optional, List, Union

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import Kpoints
from vise.defaults import defaults
from vise.input_set.kpoints_mode import KpointsMode
from vise.input_set.task import Task
from vise.util.bravais_lattice import BravaisLattice
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


class StructureKpointsGenerator:
    def __init__(
            self,
            initial_structure: Structure,
            task: Task,
            kpt_density: float,  # in Å
            kpt_mode: Union[KpointsMode, str, None] = None,  # None for default
            gamma_centered: Optional[bool] = None,  # Vasp definition
            only_even_num_kpts: bool = False,  # If ceil kpt numbers to be even.
            num_kpt_factor: Optional[int] = None,  # Set NKRED to this as well.
            band_ref_dist: float = defaults.band_mesh_distance,
            symprec: float = defaults.symmetry_length_tolerance,
            angle_tolerance: float = defaults.symmetry_angle_tolerance,
            is_magnetization: bool = False):  # Whether the system is magnetic.

        if kpt_mode:
            kpt_mode = KpointsMode.from_string(str(kpt_mode))

        self._initial_structure = initial_structure.copy()
        self._task = task
        self._kpt_density = kpt_density
        self._is_magnetization = is_magnetization
        self._num_kpt_factor = num_kpt_factor or self._task.default_kpt_factor
        if self._num_kpt_factor != 1:
            logger.info(f"kpoint factor is set to {self._num_kpt_factor}")
        self._symmetrizer = \
            StructureSymmetrizer(structure=self._initial_structure,
                                 symprec=symprec,
                                 angle_tolerance=angle_tolerance,
                                 band_mesh_distance=band_ref_dist,
                                 time_reversal=(not self._is_magnetization))
        # overwrite options fixed by other options
        self._gamma_centered = task.requisite_gamma_centered or gamma_centered
        self._adjust_only_even_num_kpts(only_even_num_kpts)
        self._adjust_kpt_mode(kpt_mode)

    def _adjust_only_even_num_kpts(self, option):
        task_requisite = self._task.requisite_only_even_num_kpts
        if isinstance(task_requisite, bool):
            result = task_requisite
        else:
            result = option
        self._only_even_num_kpts = result

    def _adjust_kpt_mode(self, option):
        self._kpt_mode = option or self._task.default_kpt_mode

    def generate_input(self):
        self._make_structure()
        self._make_kpoints()

    def _make_structure(self):
        if self._kpt_mode == KpointsMode.uniform:
            self._structure = self._initial_structure.copy()
        elif self._kpt_mode == KpointsMode.primitive:
            self._structure = self._symmetrizer.primitive
        elif self._kpt_mode == KpointsMode.band:
            self._structure = self._symmetrizer.band_primitive
        else:
            raise NotImplementedError

    def _make_kpoints(self):
        self._set_num_kpt_list()
        self._set_kpt_shift()
        self._set_kpoints()

    def _set_num_kpt_list(self) -> None:
        if self._task.requisite_num_kpt_list:
            self._num_kpt_list = self._task.requisite_num_kpt_list
            return

        kpt_list = []
        for reciprocal_lattice_length in self._reciprocal_lat_abc:
            a_kpt_num = self._kpt_density * reciprocal_lattice_length

            if self._only_even_num_kpts:
                rounded_up_kpt_num = ceil(a_kpt_num / 2) * 2
            else:
                rounded_up_kpt_num = ceil(a_kpt_num)

            kpt_list.append(rounded_up_kpt_num * self._num_kpt_factor)
        self._num_kpt_list = kpt_list

    @property
    def _reciprocal_lat_abc(self) -> List[float]:
        recipro_lat_abc = self._structure.lattice.reciprocal_lattice.abc
        self.bravais = BravaisLattice.from_sg_num(self._symmetrizer.sg_number)
        if self._kpt_mode.band_or_primitive and self.bravais.need_same_num_kpt:
            logger.warning("To keep the space group symmetry, the number "
                           "of k-points along three directions must be the "
                           "same for oI and tI Bravais lattice.")
            geometric_mean_abc = pow(float(np.prod(recipro_lat_abc)), 1 / 3)
            result = (geometric_mean_abc,) * 3
        else:
            result = recipro_lat_abc

        return result

    def _set_kpt_shift(self):
        if self._gamma_centered:
            self._kpt_shift = [0.0, 0.0, 0.0]
            return

        if self._kpt_mode is KpointsMode.uniform:
            kpt_shift = []
            angles = self._structure.lattice.angles
            for i in range(3):
                normal_to_plane = max([abs(angles[i - 2] - 90),
                                       abs(angles[i - 1] - 90)]) < 1e-5
                kpt_shift.append(0.5 if normal_to_plane else 0.0)
        else:
            kpt_shift = self.bravais.kpt_centering

        even_num_kpt_pos = [self._num_kpt_list[i] % 2 == 0 for i in range(3)]
        self._kpt_shift = [s * a for s, a in zip(kpt_shift, even_num_kpt_pos)]

    def _set_kpoints(self):
        # symmetrizer must be recreated since the fractional coordinates are
        # different in different lattices.
        # Symmetrized structure doesn't need the symprec and angle_tolerance
        symmetrizer = StructureSymmetrizer(self._structure)
        irreducible_kpoints = symmetrizer.irreducible_kpoints(
            self._num_kpt_list, self._kpt_shift)

        self.comment = ""
        if self._kpt_mode is KpointsMode.band:
            self._add_weighted_kpts(irreducible_kpoints)
            self._add_band_path_kpts()
            self._num_kpts = len(self._kpoints.kpts)
        else:
            self._kpoints = Kpoints(comment=self.comment,
                                    kpts=(self._num_kpt_list,),
                                    kpts_shift=self._kpt_shift)
            self._num_kpts = len(irreducible_kpoints)

            self.kpoints.comment += f"Num irrep kpoints: {self._num_kpts}"

    def _add_weighted_kpts(self, irreducible_kpoints):
        kpts = [k[0] for k in irreducible_kpoints]
        weights = [k[1] for k in irreducible_kpoints]
        labels = [None] * len(irreducible_kpoints)
        self._kpoints = Kpoints(comment=self.comment,
                                style=Kpoints.supported_modes.Reciprocal,
                                num_kpts=len(kpts),
                                kpts=kpts,
                                kpts_weights=weights,
                                labels=labels)

    def _add_band_path_kpts(self):
        k_path = self._symmetrizer.seekpath_data["explicit_kpoints_rel"]
        k_labels = self._symmetrizer.seekpath_data["explicit_kpoints_labels"]
        formula = self._structure.composition.reduced_formula
        sg = self._symmetrizer.sg_number
        self.kpoints.comment += \
            f"k-path added by seekpath. Formula: {formula} SG: {sg} "
        self._kpoints.num_kpts += len(k_path)
        self._kpoints.kpts += list(k_path)
        self._kpoints.labels += k_labels
        self._kpoints.kpts_weights += [0] * len(k_path)

    @property
    def structure(self):
        return self._structure

    @property
    def kpoints(self):
        return self._kpoints

    @property
    def num_kpts(self):
        return self._num_kpts

    @property
    def num_kpt_factor(self):
        return self._num_kpt_factor

