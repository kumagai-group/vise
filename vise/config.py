# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.cli.main_tools import get_user_settings

SYMMETRY_TOLERANCE = 0.01
ANGLE_TOL = 5

DOS_STEP_SIZE = 0.01

INSULATOR_KPT_DENSITY = 2.5
METAL_KPT_FACTOR = 5.0  # Multiplied factor used for metallic systems.

BAND_MESH_DISTANCE = 0.025

ENCUT_FACTOR_STR_OPT = 1.3  # This value * ENMAX is used for structure opt calc.

BAND_GAP_CRITERION = 0.2  # Criterion in eV to determine if a band gap exists.

DEFAULT_NUM_NODES = 1

# chempotdiag
ROOM_TEMPERATURE = 298.15
REFERENCE_PRESSURE = 1e5

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
                     "potcar_set",
                     "potcar_set_name",
                     "max_relax_num",
                     "removed_files",
                     "left_files"]

var_keys = [i for i in globals().keys() if i[0].isupper()]
all_keys = var_keys + main_setting_keys

user_settings, VISE_YAML_FILES = \
    get_user_settings(yaml_filename="vise.yaml",
                      setting_keys=all_keys,
                      home_hidden_directory=".vise")

for k, v in user_settings.items():
    if k in globals():
        globals()[k] = v
    else:
        MAIN_SETTINGS[k] = v



