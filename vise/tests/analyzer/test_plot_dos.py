# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from collections import OrderedDict

from pymatgen import Spin

from vise.analyzer.plot_dos import OrbitalDos, DosInfo, DosPlotter, DosMplDefaults

"""
TODO:
+ Add labels

+ Allow to set vbm, cbm, efermi
+ Shift energy zero to vbm or efermi.


DONE:
+ Show one spectrum.
+ Show spin down as well.
+ Add zero at y=0.
+ Show two data.
+ Allow to set y-scale for total
+ Allow to set y-scale for pdos
+ Allow to set x-scale  
"""


@pytest.fixture
def energies():
    return [i for i in range(-5, 6)]


@pytest.fixture
def total_up():
    return [0.0] * 2 + [2] * 2 + [4] * 4 + [2] * 2 + [0]


@pytest.fixture
def total_down():
    return [0.0] * 4 + [6] * 6 + [0]


@pytest.fixture
def h_s_up():
    return [0.0] * 2 + [2] * 8 + [0]


@pytest.fixture
def h_s_down():
    return [0.0] * 4 + [1] * 6 + [0]


@pytest.fixture
def h_p_up():
    return [0.0] * 2 + [3] * 8 + [0]


@pytest.fixture
def h_p_down():
    return [0.0] * 3 + [2] * 7 + [0]


@pytest.fixture
def doses(total_up, total_down, h_s_up, h_s_down, h_p_up, h_p_down):

    total_dos = [OrbitalDos("total", [total_up, total_down])]

    h_s = OrbitalDos("H-s", [h_s_up, h_s_down])
    h_p = OrbitalDos("H-p", [h_p_up, h_p_down])
    h_dos = [h_s, h_p]

    return [total_dos, h_dos]


@pytest.fixture
def dos_info(doses, energies):
    return DosInfo(energies=energies, doses=doses)


@pytest.fixture
def xlim():
    return [-6, 6]


@pytest.fixture
def ylim_set():
    return [[-10, 10], [-5, 5]]


@pytest.fixture
def dos_info_manual_axis(doses, energies, xlim, ylim_set):
    return DosInfo(energies=energies, doses=doses, xlim=xlim, ylim_set=ylim_set)


def test_orbital_dos(total_up, total_down):
    orbital_dos = OrbitalDos("total", [total_up, total_down])
    assert orbital_dos.max_each_dos() == 6.0


def test_doses(dos_info, energies, total_up, total_down):
    assert dos_info.energies == energies
    assert dos_info.doses[0][0].name == "total"
    assert dos_info.doses[0][0].dos == [total_up, total_down]
    assert dos_info.xlim == [-10, 10]
    assert dos_info.ylim_set == [[-6.6, 6.6], [-3.3, 3.3]]


def test_doses_manual_axis(dos_info_manual_axis, xlim, ylim_set):
    assert dos_info_manual_axis.xlim == xlim
    assert dos_info_manual_axis.ylim_set == ylim_set


def colors():
    return DosMplDefaults().colors


def test_dos_mpl_defaults():
    dos_defaults = DosMplDefaults()
    assert dos_defaults.colors == ['#36454f', '#E15759', '#4E79A7', '#F28E2B',
                                   '#76B7B2']
    assert dos_defaults.linewidth == 1.0
    assert dos_defaults.title_font_size == 15
    assert dos_defaults.label_font_size == 15
    assert dos_defaults.vline == \
           {"linewidth": 0.75, "color": "black", "linestyle": "-."}

    dos_defaults = DosMplDefaults(colors=["black"],
                                  linewidth=2.0,
                                  band_edge_line_width=1.0,
                                  band_edge_line_color="blue",
                                  band_edge_line_style=":",
                                  title_font_size=30,
                                  label_font_size=40)

    assert dos_defaults.colors == ["black"]
    assert dos_defaults.linewidth == 2.0
    assert dos_defaults.title_font_size == 30
    assert dos_defaults.label_font_size == 40
    assert dos_defaults.vline == \
           {"linewidth": 1.0, "color": "blue", "linestyle": ":"}


@pytest.fixture
def mock_plt_1st_ax(dos_info_manual_axis, mocker):
    mock_plt = mocker.patch("vise.analyzer.plot_dos.plt", auto_spec=True)
    mock_1st_ax = mocker.MagicMock()
    mock_2nd_ax = mocker.MagicMock()
    other_axs_len = (len(dos_info_manual_axis.doses) - 2)

    mock_axs = [mock_1st_ax, mock_2nd_ax] + [mocker.MagicMock()] * other_axs_len

    mock_plt.subplots.return_value = (None, mock_axs)
    plotter = DosPlotter(dos_info=dos_info_manual_axis)
    plotter.construct_plot()

    return mock_plt, mock_1st_ax, mock_2nd_ax


def test_plot_dos(energies, total_up, total_down, h_s_up, h_p_down, dos_info,
                  mock_plt_1st_ax):
    num_axs = len(dos_info.doses)
    mock_plt, mock_1st_ax, mock_2nd_ax = mock_plt_1st_ax
    mock_plt.subplots.assert_called_once_with(num_axs, 1, sharex=True)

    reversed_total_down = [dos * -1 for dos in total_down]
    reversed_h_p_down = [dos * -1 for dos in h_p_down]

    args = DosMplDefaults().dos_line(0)
    args["label"] = "total"
    mock_1st_ax.plot.assert_any_call(energies, total_up, **args)

    args = DosMplDefaults().dos_line(0)
    mock_1st_ax.plot.assert_any_call(energies, reversed_total_down, **args)

    mock_1st_ax.axhline.assert_called_once_with(0, linestyle=":", color="black")

    args = DosMplDefaults().dos_line(0)
    args["label"] = "H-s"
    mock_2nd_ax.plot.assert_any_call(energies, h_s_up, **args)

    args = DosMplDefaults().dos_line(1)
    mock_2nd_ax.plot.assert_any_call(energies, reversed_h_p_down, **args)


def test_set_figure_legend(mock_plt_1st_ax):
    _, mock_1st_ax, mock_2nd_ax = mock_plt_1st_ax
    mock_1st_ax.legend.assert_called_with(loc="best", markerscale=0.1)
    mock_2nd_ax.legend.assert_called_with(loc="best", markerscale=0.1)


def test_set_x_and_y_range(mock_plt_1st_ax):
    _, mock_1st_ax, _ = mock_plt_1st_ax
    mock_1st_ax.set_ylim.assert_called_with([-10, 10])
    mock_1st_ax.set_xlim.assert_called_with([-6, 6])


def test_actual_plot(dos_info_manual_axis):
    plotter = DosPlotter(dos_info_manual_axis)
    plotter.construct_plot()
    plotter.plt.show()
