# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os
from pathlib import Path
from monty.design_patterns import singleton

from vise.cli.main_tools import get_user_settings
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.input_set.datasets.dataset_util import PotcarSet

SYMMETRY_TOLERANCE = 0.01
ANGLE_TOL = 5

DOS_STEP_SIZE = 0.01

INSULATOR_KPT_DENSITY = 2.5
METAL_KPT_FACTOR = 5.0  # Multiplied factor used for metallic systems.

BAND_MESH_DISTANCE = 0.025

ENCUT_FACTOR_STR_OPT = 1.3  # This value * ENMAX is used for structure opt calc.

BAND_GAP_CRITERION = 0.2  # Criterion in eV to determine if a band gap exists.

DEFAULT_NUM_NODES = 1

MAGNETIZATION_THRESHOLD = 0.05

MAIN_SETTINGS = {}
# The following keys are set by vise.yaml
main_setting_keys = ["xc",
                     "vise_opts",
                     "user_incar_settings",
                     "ldauu",
                     "ldaul",
                     "outcar",
                     "contcar",
                     "vasprun",
                     "procar",
                     "potcar_dict",
                     "potcar_set_name",
                     "removed_files",
                     "left_files"]

var_keys = [i for i in globals().keys() if i[0].isupper()]
all_keys = var_keys + main_setting_keys

user_settings, VISE_YAML_FILES = \
    get_user_settings(yaml_filename="vise.yaml",
                      setting_keys=all_keys,
                      home_hidden_directory=".vise")


@singleton
class Defaults:
    def __init__(self):
        self.symmetry_length_tolerance = 0.01
        self.symmetry_angle_tolerance = 5.0
        self.dos_step_size = 0.01
        self.kpoint_density = 5.0  # for insulator 2.5 is okay.
        self.band_mesh_distance = 0.025
        self.encut_factor_str_opt = 1.3
        self.band_gap_criterion = 0.2  # in eV
        self.magnetization_criterion = 0.05
        self.default_num_nodes = 1
        self.task = str(Task.structure_opt)
        self.xc = str(Xc.pbe)
        self.vise_opts = {}
        self.user_incar_settings = {}
        self.ldauu = {}
        self.ldaul = {}
        self.outcar = "OUTCAR"
        self.contcar = "CONTCAR"
        self.vasprun = "vasprun.xml"
        self.procar = "PROCAR"
        self.potcar_dict = {}
        self.potcar_set_name = str(PotcarSet.normal)

        a_user_settings, self.vise_yaml_files = \
            get_user_settings(yaml_filename="vise.yaml",
                              setting_keys=all_keys,
                              home_hidden_directory=".vise")

        # for k, v in user_settings.items():
        #     if hasattr(self, k):
        #         setattr(self, k) = v


defaults = Defaults()


class UserSettings:
    def __init__(self, yaml_filename: str):
        self.cwd = Path.cwd()
        self.yaml_filename = yaml_filename
        self.yaml_file_list = self._yaml_filenames

    @property
    def _yaml_filenames(self):
        result = []

        dirname = self.cwd
        while True:
            target_yaml_filename = dirname / self.yaml_filename
            if target_yaml_filename.exists():
                result.append(target_yaml_filename)

            if dirname == Path("/"):
                break
            else:
                dirname = dirname.parent

        return list(reversed(result))

#    def parse_yaml_files(self):

