# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from vise.cli.main_tools import get_user_settings

SYMMETRY_TOLERANCE = 0.01
ANGLE_TOL = 5

DOS_STEP_SIZE = 0.01

KPT_DENSITY = 2.5
KPT_INIT_DENSITY = 2.5  # Initial k-point density used for k-point convergence
KPT_FACTOR = 1.2  # Multiplied factor used for incrementing the k-point density.

BAND_MESH_DISTANCE = 0.025

ENCUT_FACTOR_STR_OPT = 1.3  # This times ENMAX is used for structure opt calc.

BAND_GAP_CRITERION = 0.2  # Criterion in eV to determine if band gap exists.

DEFAULT_NUM_NODES = 1

TIMEOUT = 518400

# chempotdiag
ROOM_TEMPERATURE = 298.15
REFERENCE_PRESSURE = 1e5


MAIN_SETTINGS = {}
# The following keys are set by vise.yaml
main_setting_keys = ["vasp_cmd",
                     "xc",
                     "vise_opts",
                     "user_incar_settings",
                     "ldauu",
                     "ldaul",
                     "outcar",
                     "convergence_criterion",  # kpt convergence
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



