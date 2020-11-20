# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.dielectric_function import AbsorptionCoeffPlotlyPlotter
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.analyzer.vasp.make_diele_func import make_diele_func, \
    make_shifted_diele_func

import numpy as np


def test_make_diele_func(test_data_files):
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun_gamma.xml")
    o = Outcar(test_data_files / "MgSe_absorption_OUTCAR_gamma")
    actual = make_diele_func(v, o)


#    print(actual.diele_func_imag)
#    print(actual.band_gap)

    print(VaspBandEdgeProperties(v, o))
    assert actual.energies[1] == 0.0407


def test_make_diele_func_2(test_data_files):
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun_gamma.xml")
    o = Outcar(test_data_files / "MgSe_absorption_OUTCAR_gamma")
    actual = make_shifted_diele_func(make_diele_func(v, o),
                                     original_band_gap=1.997,
                                     shift=1.0)
    assert isinstance(actual.diele_func_imag[0], list)
    assert isinstance(actual.diele_func_imag[0][0], float)
#    print(type(actual.diele_func_imag[0]))
#    print(type(actual.diele_func_imag[0][0]))
    #    print(actual.diele_func_real)
    # print(actual.diele_func_imag[1000])
    # plotter = AbsorptionCoeffPlotlyPlotter(actual)
    # fig = plotter.create_figure()
    # fig.show()

