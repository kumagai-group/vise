# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np

from pymatgen.electronic_structure.boltztrap2 import VasprunBSLoader, \
    BztInterpolator, BztTransportProperties
from pymatgen.io.vasp import Vasprun
from vise.analyzer.effective_mass import EffectiveMass


def make_effective_mass(vasprun, temps, concentrations):
    vl = VasprunBSLoader(vasprun)
    bi = BztInterpolator(vl)
    btp = BztTransportProperties(bi, temp_r=np.array(temps))
    btp.compute_properties_doping(concentrations)

    return EffectiveMass(p=btp.Effective_mass_doping["p"],
                         n=btp.Effective_mass_doping["n"],
                         temps=temps,
                         concentrations=concentrations)