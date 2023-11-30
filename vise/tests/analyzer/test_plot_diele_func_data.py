# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import copy

import pytest
from monty.serialization import loadfn
from vise.analyzer.plot_diele_func_data import (DieleFuncMplPlotter,
                                                DieleFuncPlotType)


def test_diele_func_plot_type():
    actual = DieleFuncPlotType.absorption_coeff.y_axis_label("plotly")
    expected = "Absorption coefficient (cm<sup>-1</sup>)"
    assert actual == expected
    actual = DieleFuncPlotType.absorption_coeff.y_axis_label("matplotlib")
    expected = "Absorption coefficient (cm$^{-1}$)"
    assert actual == expected


@pytest.fixture
def actual_diele_func_data(test_data_files):
    result = loadfn(test_data_files / "MgSe_diele_func_data.json")
    result.title = "test diele"
    return result


def test_absorption_coeff_mpl_plotter(actual_diele_func_data):
    plotter = DieleFuncMplPlotter(actual_diele_func_data)
    plotter.construct_plot()
    plotter.plt.show()


def test_refraction_mpl_plotter(actual_diele_func_data):
    plotter = DieleFuncMplPlotter(actual_diele_func_data)
    plotter.construct_plot(plot_type=DieleFuncPlotType.refraction)
    plotter.plt.show()
