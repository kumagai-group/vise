# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.plot_absorption_coeff import AbsorptionCoeffPlotlyPlotter, \
    AbsorptionCoeffMplPlotter

from vise.analyzer.vasp.make_diele_func import make_diele_func
from vise.util.dash_helper import show_png


try:
    import psutil
    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


@pytest.fixture
def actual_diele_func_data(test_data_files):
    v = Vasprun(test_data_files / "MgSe_absorption_vasprun.xml")
    o = Outcar(test_data_files / "MgSe_absorption_OUTCAR")
    return make_diele_func(v, o)


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_absorption_coeff_plotly_plotter(actual_diele_func_data):
    plotter = AbsorptionCoeffPlotlyPlotter(actual_diele_func_data,
                                           materials=["GaAs"])
    fig = plotter.create_figure()
    show_png(fig)


def test_absorption_coeff_mpl_plotter(actual_diele_func_data):
    plotter = AbsorptionCoeffMplPlotter(actual_diele_func_data,
                                        materials=["GaAs", "Si"])
    #                                        align_gap=False)
    plotter.construct_plot()
    plotter.plt.show()
