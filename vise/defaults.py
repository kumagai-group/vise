# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from monty.design_patterns import singleton

from vise.input_set.datasets.dataset_util import PotcarSet
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.user_settings import UserSettings


@singleton
class Defaults:
    def __init__(self):
        self.symmetry_length_tolerance = 0.01
        self.symmetry_angle_tolerance = 5.0
        self.dos_step_size = 0.01
        self.kpoint_density = 5.0  # for insulator 2.5 is okay.
        self.band_mesh_distance = 0.025
        self.str_opt_encut_factor = 1.3
        self.band_gap_criterion = 0.2  # in eV
        self.integer_criterion = 0.1
        self.default_num_nodes = 1
        self.task = Task.structure_opt
        self.xc = Xc.pbe
        self.options = {}
        self.user_incar_settings = {}
        self.ldauu = {}
        self.ldaul = {}
        self.outcar = "OUTCAR"
        self.contcar = "CONTCAR"
        self.vasprun = "vasprun.xml"
        self.procar = "PROCAR"
        self.overridden_potcar = {}
        self.potcar_set_name = str(PotcarSet.normal)

        user_settings = UserSettings(yaml_filename="vise.yaml")
        self.yaml_files = user_settings.yaml_files_from_root_dir
        self.user_settings = user_settings.user_settings

        for k, v in self.user_settings.items():
            if hasattr(self, k):
                self.__setattr__(k, v)


defaults = Defaults()


