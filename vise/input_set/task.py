# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import Optional, List

from monty.json import MSONable
from vise.input_set.kpoints_mode import KpointsMode
from vise.util.enum import ExtendedEnum


class Task(MSONable, ExtendedEnum):
    structure_opt = "structure_opt"
    structure_opt_rough = "structure_opt_rough"
    structure_opt_tight = "structure_opt_tight"
    cluster_opt = "cluster_opt"
    phonon_force = "phonon_force"
    defect = "defect"
    band = "band"
    dos = "dos"
    dielectric_dfpt = "dielectric_dfpt"
    dielectric_finite_field = "dielectric_finite_field"
    dielectric_function = "dielectric_function"

    @property
    def is_lattice_relaxed(self) -> bool:
        return "structure_opt" in self.name

    @property
    def is_atom_relaxed_lattice_fixed(self) -> bool:
        return self in (self.cluster_opt, self.defect)

    @property
    def is_atom_relaxed(self) -> bool:
        return (self.is_lattice_relaxed
                or self.is_atom_relaxed_lattice_fixed)

    @property
    def is_dielectric(self) -> bool:
        return "dielectric" in self.name

    @property
    def is_tight_calc(self) -> bool:
        return self in (self.structure_opt_tight, self.phonon_force)

    @property
    def is_plot_task(self) -> bool:
        return self in (self.band, self.dos, self.dielectric_function)

    @property
    def is_spectrum_task(self) -> bool:
        return self in (self.dos, self.dielectric_function)

    @property
    def default_kpt_factor(self) -> Optional[int]:
        if self is self.dielectric_function:
            return 3
        elif self in (self.dos,
                      self.dielectric_dfpt,
                      self.dielectric_finite_field):
            return 2
        elif self.condition_kpt_factor_is_one:
            return 1
        else:
            raise NotImplementedError

    @property
    def condition_kpt_factor_is_one(self):
        if self in (self.band,
                    self.cluster_opt,
                    self.phonon_force,
                    self.defect) or \
                self.is_lattice_relaxed:
            return True
        else:
            return False

    @property
    def default_kpt_mode(self) -> KpointsMode:
        if self == self.band:
            return KpointsMode.band
        elif self in (self.defect, self.cluster_opt, self.phonon_force):
            return KpointsMode.uniform
        elif self.condition_kpoints_mode_is_primitive:
            return KpointsMode.primitive
        else:
            raise NotImplementedError

    @property
    def condition_kpoints_mode_is_primitive(self):
        if self is self.dos or \
                self.is_dielectric or \
                self.is_lattice_relaxed:
            return True
        else:
            return False

    @property
    def requisite_num_kpt_list(self) -> Optional[List[int]]:
        if self == self.cluster_opt:
            return [1, 1, 1]
        return

    @property
    def requisite_only_even_num_kpts(self) -> Optional[bool]:
        if self == self.cluster_opt:
            return False
        return

    # Gamma-center mesh is a must for GW calculations due to vasp
    # implementation and tetrahedron method, while is a strong recommend for
    # dos and dielectric function to sample the band edges.
    @property
    def requisite_gamma_centered(self) -> Optional[bool]:
        if self in (self.cluster_opt, self.phonon_force) \
                or self.is_spectrum_task:
            return True
        return

    @property
    def need_spin(self):
        return True if self is self.defect else False

    @property
    def fine_to_inherit_structure_from_prev(self):
        # For phonon_calculations, the POSCAR must no be inherited.
        return self is not self.phonon_force
