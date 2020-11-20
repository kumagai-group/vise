# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest
from monty.serialization import loadfn
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.dielectric_function import DieleFuncData, \
    AbsorptionCoeffMplPlotter, AbsorptionCoeffPlotlyPlotter, eV_to_inv_cm

import numpy as np
from math import pi, sqrt

from vise.analyzer.vasp.make_diele_func import make_diele_func
from vise.tests.conftest import test_data_files, assert_msonable
from vise.util.dash_helper import show_png


try:
    import psutil
    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


@pytest.fixture
def diele_func_data():
    array_real = [1, 2, 3, 0, 0, 0]
    array_imag = [4, 5, 6, 0, 0, 0]
    return DieleFuncData(energies=list(np.linspace(0.0, 10.0, num=11)),
                         diele_func_real=[array_real]*11,
                         diele_func_imag=[array_imag]*11,
                         band_gap=1.0)


def test_msonable(diele_func_data):
    assert_msonable(diele_func_data)


def test_json_file_mixin(diele_func_data, tmpdir):
    tmpdir.chdir()
    diele_func_data.to_json_file()
    actual = loadfn("diele_func_data.json")
    assert actual.diele_func_real == diele_func_data.diele_func_real


def test_ave_absorption_coeff(diele_func_data):
    actual = diele_func_data.ave_absorption_coeff
    expected = [2 * sqrt(2) * pi * sqrt(sqrt(2 ** 2 + 5 ** 2) - 2)
                * i * eV_to_inv_cm for i in range(0, 11)]
    assert actual == expected


@pytest.fixture
def actual_diele_func_data(test_data_files):
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun.xml")
    o = Outcar(test_data_files / "MgSe_absorption_OUTCAR")
    return make_diele_func(v, o)


def test_target_coeff_e_from_band_gap(actual_diele_func_data):
    actual = actual_diele_func_data.target_coeff_e_from_band_gap
    np.testing.assert_almost_equal(actual, 0.4014000000000002)


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_absorption_coeff_plotly_plotter(actual_diele_func_data):
    plotter = AbsorptionCoeffPlotlyPlotter(actual_diele_func_data,
                                           materials=["GaAs"])
    #                                        align_gap=False)
    fig = plotter.create_figure()
    show_png(fig)


def test_absorption_coeff_plotly_plotter_wo_alignment(actual_diele_func_data):
    plotter = AbsorptionCoeffPlotlyPlotter(actual_diele_func_data,
                                           materials=["GaAs"])
    fig = plotter.create_figure()
    fig.show()
#    show_png(fig)


def test_absorption_coeff_mpl_plotter(actual_diele_func_data):
    plotter = AbsorptionCoeffMplPlotter(actual_diele_func_data,
                                        materials=["GaAs", "Si"])
#                                        align_gap=False)
    plotter.construct_plot()
    plotter.plt.show()
