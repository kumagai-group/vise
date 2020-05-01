# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from collections import OrderedDict

from pymatgen import Spin

from vise.analyzer.plot_dos import DosPlotter

"""
TODO:
+ Use +- same scale for total with the ylim max is 1.1 * max of dos
+ Use +- same scale for pdos with the ylim max is 1.1 * max of dos
+ Use consistent scale for pdos.

+ Allow to set x-scale  
+ Allow to set y-scale for total
+ Allow to set y-scale for pdos
+ Allow to set vbm, cbm, efermi
+ Shift energy zero to vbm or efermi.
+ Add labels


DONE:
+ Show one spectrum.
+ Show spin down as well.
+ Add zero at y=0.
+ Show two data.
"""


@pytest.fixture
def energies():
    return [i for i in range(-5, 6)]


@pytest.fixture
def total_up():
    return [0] * 2 + [4] * 8 + [0]


@pytest.fixture
def total_down():
    return [0] * 4 + [6] * 6 + [0]


@pytest.fixture
def h_s_up():
    return [0] * 2 + [2] * 8 + [0]


@pytest.fixture
def h_s_down():
    return [0] * 4 + [1] * 6 + [0]


@pytest.fixture
def h_p_up():
    return [0] * 2 + [3] * 8 + [0]


@pytest.fixture
def h_p_down():
    return [0] * 4 + [1] * 6 + [0]


@pytest.fixture
def doses(total_up, total_down, h_s_up, h_s_down, h_p_up, h_p_down):

    total_dos = {Spin.up: total_up, Spin.down: total_down}

    h_s = {Spin.up: h_s_up, Spin.down: h_s_down}
    h_p = {Spin.up: h_p_up, Spin.down: h_p_down}

    total_dos_dict = OrderedDict({"total": total_dos})
    h_dos_dict = OrderedDict({"H-s": h_s, "H-p": h_p})

    return [total_dos_dict, h_dos_dict]


@pytest.fixture
def mock_plt_1st_ax(energies, doses, mocker):
    mock_plt = mocker.patch("vise.analyzer.plot_dos.plt", auto_spec=True)
    mock_1st_ax = mocker.MagicMock()
    mock_2nd_ax = mocker.MagicMock()
    other_mocks = [mocker.MagicMock()] * (len(doses) - 2)
    mock_axs = [mock_1st_ax, mock_2nd_ax] + other_mocks

    mock_plt.subplots.return_value = (None, mock_axs)
    plotter = DosPlotter(energies, doses)
    plotter.construct_plot()

    return mock_plt, mock_1st_ax, mock_2nd_ax


def test_plot_dos(energies, total_up, total_down, doses, mock_plt_1st_ax):
    num_axs = len(doses)
    mock_plt, mock_1st_ax, _ = mock_plt_1st_ax

    mock_plt.subplots.assert_called_once_with(num_axs, 1, sharex=True)

    total_dos_sign_reversed = [dos * -1 for dos in total_down]

    mock_1st_ax.plot.assert_any_call(energies, total_up)
    mock_1st_ax.plot.assert_any_call(energies, total_dos_sign_reversed)
    mock_1st_ax.axhline.assert_called_once_with(0)


def test_set_y_range(distances, mock_plt_axis):
    _, mock_1st_ax, _ = mock_plt_1st_ax
    mock_1st_ax.xlim.assert_called_with(-6.6, 6.6)


def test_actual_plot(energies, doses):
    plotter = DosPlotter(energies, doses)
    plotter.construct_plot()
    plotter.plt.show()
