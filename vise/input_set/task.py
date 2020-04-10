# -*- coding: utf-8 -*-

from vise.util.enum import ExtendedEnum


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

    @property
    def kpt_factor(self):
        if self is self.dielectric_function:
            return 3
        elif self in (self.dos,
                      self.dielectric_dfpt,
                      self.dielectric_finite_field):
            return 2

        return 1

    @property
    def is_lattice_relax(self):
        return "structure_opt" in self.name

    @property
    def is_plot_task(self):
        return self in (self.band, self.dos, self.dielectric_function)

    @property
    def is_spectrum_task(self):
        return self in (self.dos, self.dielectric_function)

