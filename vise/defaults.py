# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from abc import ABC

from monty.design_patterns import singleton
from pathlib import Path
from vise.input_set.datasets.potcar_set import PotcarSet
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.user_settings import UserSettings
from vise.util.logger import get_logger

logger = get_logger(__name__)


class DefaultsBase(ABC):
    def set_user_settings(self, yaml_filename):
        user_settings = UserSettings(yaml_filename=yaml_filename)
        self.yaml_files = user_settings.yaml_files_from_root_dir
        self.user_settings = user_settings.user_settings

        for k, v in self.user_settings.items():
            if hasattr(self, k):
                self.__setattr__("_" + k, v)
            else:
                logger.warning(f"{k} key in {yaml_filename} is not valid.")


@singleton
class Defaults(DefaultsBase):
    def __init__(self):
        self._symmetry_length_tolerance = 0.01
        self._symmetry_angle_tolerance = 5.0
        self._dos_step_size = 0.01
        self._kpoint_density = 5.0
        self._insulator_kpoint_density = 2.5
        self._defect_kpoint_density = 1.8
        self._band_mesh_distance = 0.025
        self._str_opt_encut_factor = 1.3
        self._band_gap_criterion = 0.2  # in eV
        self._integer_criterion = 0.1
        self._default_num_nodes = 1
        self._task = str(Task.structure_opt)
        self._xc = str(Xc.pbe)
        self._options = {}
        self._user_incar_settings = {}
        self._ldauu = {}
        self._ldaul = {}
        self._outcar = "OUTCAR"
        self._contcar = "CONTCAR"
        self._vasprun = "vasprun.xml"
        self._procar = "PROCAR"
        self._overridden_potcar = []
        self._potcar_set = str(PotcarSet.normal)

        self.set_user_settings(yaml_filename="vise.yaml")

    @property
    def symmetry_length_tolerance(self):
        return self._symmetry_length_tolerance

    @property
    def symmetry_angle_tolerance(self):
        return self._symmetry_angle_tolerance

    @property
    def dos_step_size(self):
        return self._dos_step_size

    @property
    def kpoint_density(self):
        return self._kpoint_density

    @property
    def insulator_kpoint_density(self):
        return self._insulator_kpoint_density

    @property
    def defect_kpoint_density(self):
        return self._defect_kpoint_density

    @property
    def band_mesh_distance(self):
        return self._band_mesh_distance

    @property
    def str_opt_encut_factor(self):
        return self._str_opt_encut_factor

    @property
    def band_gap_criterion(self):
        return self._band_gap_criterion

    @property
    def integer_criterion(self):
        return self._integer_criterion

    @property
    def default_num_nodes(self):
        return self._default_num_nodes

    @property
    def task(self):
        return Task(self._task)

    @property
    def xc(self):
        return Xc(self._xc)

    @property
    def options(self):
        return self._options

    @property
    def user_incar_settings(self):
        return self._user_incar_settings

    @property
    def ldauu(self):
        return self._ldauu

    @property
    def ldaul(self):
        return self._ldaul

    @property
    def outcar(self):
        return Path(self._outcar)

    @property
    def contcar(self):
        return Path(self._contcar)

    @property
    def vasprun(self):
        return Path(self._vasprun)

    @property
    def procar(self):
        return Path(self._procar)

    @property
    def overridden_potcar(self):
        return self._overridden_potcar

    @property
    def potcar_set(self):
        return PotcarSet(self._potcar_set)


defaults = Defaults()
