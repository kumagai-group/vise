# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import Optional, List

from vise.util.enum import ExtendedEnum
from vise.input_set.make_kpoints import KpointsMode


class Task(ExtendedEnum):
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
    def is_lattice_relaxed_task(self) -> bool:
        return "structure_opt" in self.name

    @property
    def is_tight_calc(self) -> bool:
        return self in (self.structure_opt_tight, self.phonon_force)

    @property
    def is_atom_relaxed_task(self) -> bool:
        return (self.is_lattice_relaxed_task or
                self in (self.cluster_opt, self.defect))

    @property
    def is_dielectric_task(self) -> bool:
        return "dielectric" in self.name

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
        elif (self in (self.band, self.cluster_opt,
                       self.phonon_force, self.defect)
              or self.is_lattice_relaxed_task):
            return 1
        else:
            raise NotImplementedError

    @property
    def default_kpt_mode(self) -> KpointsMode:
        if self == self.band:
            return KpointsMode.band
        elif self in (self.defect, self.cluster_opt, self.phonon_force):
            return KpointsMode.uniform
        elif (self is self.dos
              or self.is_dielectric_task
              or self.is_lattice_relaxed_task):
            return KpointsMode.primitive
        else:
            raise NotImplementedError

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
    def requisite_kpt_shift(self) -> Optional[List[float]]:
        if self in (self.cluster_opt, self.phonon_force) \
                or self.is_spectrum_task:
            return [0, 0, 0]
        return

    @property
    def incar_isif(self):
        if self.is_lattice_relaxed_task:
            return 3
        elif self.is_atom_relaxed_task:
            return 2
        elif self in (self.band, self.dos) or self.is_dielectric_task:
            return 0
        else:
            raise NotImplementedError

    # During dielectric_dfpt calculations, EDIFF is tightened automatically.
    @property
    def incar_ediff(self):
        if self.is_tight_calc:
            return 1e-8
        elif self in (self.structure_opt, self.cluster_opt):
            return 1e-7
        elif self in (self.dielectric_dfpt, self.dielectric_finite_field):
            return 1e-6
        elif self in (self.dielectric_function,):
            return 1e-5
        elif self in (self.structure_opt_rough,):
            return 1e-4
        else:
            raise NotImplementedError

    @property
    def incar_ediffg_optional(self):
        if self.is_atom_relaxed_task:
            if self is self.structure_opt_tight:
                return -0.001
            elif self in (self.structure_opt, self.cluster_opt):
                return -0.005
            elif self is self.defect:
                return -0.03
            elif self is self.structure_opt_rough:
                return -0.2
            else:
                raise NotImplementedError
        else:
            return

    @property
    def incar_ibrion(self):
        return 8 if self is self.dielectric_dfpt else 2

    @property
    def incar_lreal(self):
        return "Auto" if self is Task.defect else False

    @property
    def incar_ispin(self):
        return 2 if self is self.defect else 1

    @property
    def incar_prec(self):
        return "Accurate" if self.is_tight_calc else "Normal"

    @property
    def incar_nsw(self):
        return 50 if self.is_atom_relaxed_task else 0

    @property
    def incar_potim_optional(self):
        if self is self.structure_opt_rough:
            return 0.1
        return

    @property
    def incar_addgrid_optional(self):
        if self.is_tight_calc:
            return 0.1
        return















