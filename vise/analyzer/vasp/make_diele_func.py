# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.dielectric_function import DieleFuncData
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties

import numpy as np


def make_diele_func(vasprun: Vasprun, outcar: Outcar):
    energies, real, imag = vasprun.dielectric_data["density"]
    real, imag = np.array(real), np.array(imag)
    band_gap = VaspBandEdgeProperties(vasprun, outcar).band_gap
    return DieleFuncData(energies, real.tolist(), imag.tolist(), band_gap)


