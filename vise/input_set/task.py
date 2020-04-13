# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import Optional, List

from vise.util.enum import ExtendedEnum
from vise.input_set.make_kpoints import KpointsMode


class Task(ExtendedEnum):
    """ Supported tasks """
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

    def kpt_factor(self, factor: int) -> int:
        if factor:
            return factor
        elif self is self.dielectric_function:
            return 3
        elif self in (self.dos,
                      self.dielectric_dfpt,
                      self.dielectric_finite_field):
            return 2

        return 1

    @property
    def is_lattice_relax(self) -> bool:
        return "structure_opt" in self.name

    @property
    def is_plot_task(self) -> bool:
        return self in (self.band, self.dos, self.dielectric_function)

    @property
    def is_spectrum_task(self) -> bool:
        return self in (self.dos, self.dielectric_function)

    def kpt_mode(self, kpt_mode: Optional[KpointsMode]) -> KpointsMode:
        if kpt_mode:
            return kpt_mode
        elif self == self.band:
            return KpointsMode.band
        elif self in (self.defect, self.cluster_opt, self.phonon_force):
            return KpointsMode.manual_set
        else:
            return KpointsMode.uniform

    def kpts(self, kpts: Optional[List[int]]) -> List[int]:
        if self == self.cluster_opt:
            return [1, 1, 1]
        return kpts

    def only_even(self, only_even: bool) -> bool:
        if self == self.cluster_opt:
            return True
        return only_even

    # Gamma-center mesh is a must for GW calculations due to vasp
    # implementation and tetrahedron method, while is a strong recommend for
    # dos and dielectric function to sample the band edges.

    def kpt_shift(self, kpt_shift: Optional[List[float]]) -> List[float]:
        if self in (self.cluster_opt, self.phonon_force) \
                or self.is_spectrum_task:
            return [0, 0, 0]
        return kpt_shift

