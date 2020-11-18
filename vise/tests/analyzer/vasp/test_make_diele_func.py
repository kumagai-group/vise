# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.dielectric_function import AbsorptionCoeffPlotlyPlotter
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.analyzer.vasp.make_diele_func import make_diele_func, \
    make_shifted_diele_func


def test_make_diele_func(test_data_files):
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun.xml")
    o = Outcar(test_data_files / "MgSe_absorption_OUTCAR")
    actual = make_diele_func(v, o)
    print(VaspBandEdgeProperties(v, o))
    assert actual.energies[1] == 0.0407


def test_make_diele_func_2(test_data_files):
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun.xml")
    o = Outcar(test_data_files / "MgSe_absorption_OUTCAR")
    actual = make_shifted_diele_func(make_diele_func(v, o),
                                     original_band_gap=1.997,
                                     shift=1.0)
    print(len(actual.energies))
    print(len(actual.diele_func_imag))
    print(len(actual.diele_func_real))
    # print(actual.diele_func_imag[1000])
    # plotter = AbsorptionCoeffPlotlyPlotter(actual)
    # fig = plotter.create_figure()
    # fig.show()

