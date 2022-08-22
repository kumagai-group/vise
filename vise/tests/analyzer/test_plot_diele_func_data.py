# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import copy

import pytest
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.plot_diele_func_data import AbsorptionCoeffPlotlyPlotter, \
    AbsorptionCoeffMplPlotter, TensorDirection

from vise.analyzer.vasp.make_diele_func import make_diele_func
from vise.util.dash_helper import show_png


try:
    import psutil
    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


def test_tensor_direction():
    tensors = [[0.0, 1.0, 2.0, 3.0, 4.0, 5.0]]
    actual = TensorDirection.average.val(tensors)
    assert actual == [1.0]

    actual = TensorDirection.xx.val(tensors)
    assert actual == [0.0]

    actual = TensorDirection.xy.val(tensors)
    assert actual == [3.0]


@pytest.fixture
def actual_diele_func_data(test_data_files):
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun.xml")
    o = Outcar(test_data_files / "MgSe_absorption_OUTCAR")
    return make_diele_func(v, o)


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_absorption_coeff_plotly_plotter(actual_diele_func_data):
    data = copy(actual_diele_func_data)
    data.band_gap = None
    plotter = AbsorptionCoeffPlotlyPlotter(data)
    fig = plotter.create_figure(materials=["GaAs"])
    show_png(fig)


def test_absorption_coeff_mpl_plotter(actual_diele_func_data):
    plotter = AbsorptionCoeffMplPlotter(actual_diele_func_data)
    #                                        align_gap=False)
    plotter.construct_plot(materials=["GaAs", "Si"])
    plotter.plt.show()
