# -*- coding: utf-8 -*-

from enum import unique, Enum


@unique
class Task(Enum):
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

    def __str__(self):
        return self.name

    @classmethod
    def from_string(cls, s):
        for m in Task:
            if m.value == s:
                return m
        raise AttributeError(f"Task: {s} is not proper.\n"
                             f"Supported Task:\n {cls.name_list()}")

    @classmethod
    def name_list(cls):
        return ', '.join([e.value for e in cls])


LATTICE_RELAX_TASK = (Task.structure_opt, Task.structure_opt_rough,
                      Task.structure_opt_tight)
PLOT_TASK = (Task.band, Task.dos, Task.dielectric_function)
SPECTRA_TASK = (Task.dos, Task.dielectric_function)

