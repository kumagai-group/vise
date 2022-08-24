# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import copy

import pytest
from monty.serialization import loadfn
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.plot_diele_func_data import DieleFuncPlotlyPlotter, \
    DieleFuncMplPlotter, TensorDirection, DieleFuncPlotType

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


def test_diele_func_plot_type():
    actual = DieleFuncPlotType.absorption_coeff.y_axis_label("plotly")
    expected = "Absorption coefficient (cm<sup>-1</sup>)"
    assert actual == expected
    actual = DieleFuncPlotType.absorption_coeff.y_axis_label("matplotlib")
    expected = "Absorption coefficient (cm$^{-1}$)"
    assert actual == expected


@pytest.fixture
def actual_diele_func_data(test_data_files):
    return loadfn(test_data_files / "MgSe_diele_func_data.json")


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_absorption_coeff_plotly_plotter(actual_diele_func_data):
    data = copy(actual_diele_func_data)
    data.band_gap = None
    plotter = DieleFuncPlotlyPlotter(data)
    plotter.create_figure()
    show_png(plotter.fig)


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_diele_func_plotly_plotter(actual_diele_func_data):
    plotter = DieleFuncPlotlyPlotter(actual_diele_func_data)
    plotter.create_figure(
        directions=[TensorDirection.xx, TensorDirection.zz],
        plot_type=DieleFuncPlotType.diele_func)
    show_png(plotter.fig)


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_refraction_plotly_plotter(actual_diele_func_data):
    plotter = DieleFuncPlotlyPlotter(actual_diele_func_data)
    plotter.create_figure(
        directions=[TensorDirection.yy],
        plot_type=DieleFuncPlotType.refraction)
    show_png(plotter.fig)


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_reflectivity_plotly_plotter(actual_diele_func_data):
    plotter = DieleFuncPlotlyPlotter(actual_diele_func_data)
    plotter.create_figure(
        directions=[TensorDirection.yy],
        plot_type=DieleFuncPlotType.reflectivity)
    show_png(plotter.fig)


def test_absorption_coeff_mpl_plotter(actual_diele_func_data):
    plotter = DieleFuncMplPlotter(actual_diele_func_data)
    plotter.construct_plot()
    plotter.plt.show()


def test_refraction_mpl_plotter(actual_diele_func_data):
    plotter = DieleFuncMplPlotter(actual_diele_func_data)
    plotter.construct_plot(plot_type=DieleFuncPlotType.refraction)
    plotter.plt.show()
