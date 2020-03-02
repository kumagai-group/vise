# -*- coding: utf-8 -*-

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


SYMMETRY_TOLERANCE = 0.01
ANGLE_TOL = 5

DOS_STEP_SIZE = 0.01

KPT_DENSITY = 2.5

KPT_INIT_DENSITY = 2.5  # Initial k-point density used for k-point convergence
KPT_FACTOR = 1.2  # Multiplied factor used for incrementing the k-point density.

BAND_REF_DIST = 0.025

ENCUT_FACTOR_STR_OPT = 1.3  # This times ENMAX is used for structure opt calc.

BAND_GAP_CRITERION = 0.2  # Criterion in eV to determine if band gap exists.

DEFAULT_NUM_NODES = 1

TIMEOUT = 518400

# chempotdiag
ROOM_TEMPERATURE = 298.15
REFERENCE_PRESSURE = 1e+5
MOLECULE_SUFFIX = "molecule_chempotdiag"
