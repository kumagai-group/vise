# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.dielectric_function import DieleFuncData, \
    kramers_kronig_trans
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties

import numpy as np


def make_diele_func(vasprun: Vasprun,
                    outcar: Outcar,
                    use_vasp_real: bool = True,
                    ita: float = 0.01) -> DieleFuncData:

    energies, real, imag = vasprun.dielectric_data["density"]
    imag = np.array(imag)
    if use_vasp_real:
        real = np.array(real)
    else:
        # When CSHIFT = 0.0, the first component becomes 99999.0
        # the following lines are copied from the vasprun.xml
        # <dielectricfunction>
        # <imag>
        # <array>
        # <dimension dim=“1”>gridpoints</dimension>
        # <field>energy</field>
        # <field>xx</field>
        # <field>yy</field>
        # <field>zz</field>
        # <field>xy</field>
        # <field>yz</field>
        # <field>zx</field>
        # <set>
        # <r>     0.0000     0.0000     0.0000    -0.0000 99999.0000 99999.0000 99999.0000 </r>
        imag[0] = 0.0
        real = kramers_kronig_trans(imag, energies, ita)
    band_gap = VaspBandEdgeProperties(vasprun, outcar).band_gap
    return DieleFuncData(energies, real.tolist(), imag.tolist(), band_gap)


