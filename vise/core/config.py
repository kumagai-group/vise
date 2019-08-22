# -*- coding: utf-8 -*-

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


# Following defaults determine the condition of automatic defect calculations.
# electronegativity difference for antisites and substitutional impurities
ELECTRONEGATIVITY_DIFFERENCE = 1.0
# Maximum displacement displacement_distance
DISPLACEMENT_DISTANCE = 0.2
# Cutoff radius in which atoms are considered as neighbors of a defect, and
# perturbed when DISPLACEMENT_DISTANCE is set.
CUTOFF_RADIUS = 3.2
SYMPREC = 0.01
# The following must be used after structure optimization anytime.
DEFECT_SYMMETRY_TOLERANCE = 0.1

ANGLE_TOL = 5

KPT_DENSITY = 2.5
ENCUT_FACTOR_STR_OPT = 1.3

COLOR = [
    "xkcd:blue",
    "xkcd:brown",
    "xkcd:crimson",
    "xkcd:darkgreen",
    "xkcd:gold",
    "xkcd:magenta",
    "xkcd:orange",
    "xkcd:darkblue",
    "xkcd:navy",
    "xkcd:red",
    "xkcd:olive",
    "xkcd:black",
    "xkcd:indigo"
]