# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from copy import deepcopy
from unittest.mock import MagicMock

import numpy as np
import pytest

from vise.analyzer.plot_band import (
    BandPlotter, BandInfo, BandEdge, XTicks, BandMplSettings,
    BandPlotInfo, ViseBandInfoError)
from vise.tests.helpers.assertion import assert_msonable
from vise.util.matplotlib import float_to_int_formatter

"""
TODO:
+ Not use Spin

+ Raise error when 5 bands are drawn
+ Set band gap arrow(s).
+ Check multiple bands with the same distances and x_ticks?
+ To csv
+ Brillouin zone
+ Add spin info to band edge 

DONE:
"""

distances = [[0, 1, 2, 3, 4, 5, 6], [6, 7, 8, 9]]
x_ticks = XTicks(["A", "B", "C$\\mid$D", "E"], [0.0, 5.0, 6.0, 9.0])
y_range = [-10, 10]
title = "Title"
base_energy = -1.0
# [by branch][by spin][by band idx][by k-path idx]
band_energies = [[[[-3.0, -2, -1, -1, -1, -2, -3], [7.0, 6, 5, 4, 3, 2, 3]]],
                 [[[-4.0, -5, -6, -7], [5, 2, 7, 8]]]]
shifted_band_energies = [(np.array(e) + 1.0).tolist() for e in band_energies]

band_edge = BandEdge(vbm=-1, cbm=2, vbm_distances=[2, 3, 4], cbm_distances=[5, 7])
fermi_level = 1.5

band_info_fermi = BandInfo(band_energies=band_energies, fermi_level=fermi_level)
band_info_edge = BandInfo(band_energies=band_energies, band_edge=band_edge)

colors = ['#E15759', '#4E79A7', '#F28E2B', '#76B7B2']


@pytest.fixture
def band_info():
    return BandInfo(band_energies=deepcopy(band_energies), band_edge=band_edge,
                    fermi_level=fermi_level)


@pytest.fixture
def band_info_set(band_info):
    return [band_info]


@pytest.fixture
def band_plot_info(band_info_set):
    return BandPlotInfo(band_info_set, distances, x_ticks, title)


def test_band_plot_info_msonable(band_info, band_plot_info):
    assert_msonable(band_edge)
    assert_msonable(band_info)
    assert_msonable(band_plot_info)


@pytest.fixture
def mock_band_plt_list(mocker, band_plot_info):
    mock_plt = mocker.patch("vise.analyzer.plot_band.plt", auto_spec=True)
    mock_axis = MagicMock()
    mock_plt.gca.return_value = mock_axis
    plotter = BandPlotter(band_plot_info, y_range)
    plotter.construct_plot()
    return mock_plt, mock_axis


def test_band_info_slide_energies(band_info):
    band_info.slide_energies(base_energy=base_energy)
    expected_band_edge = BandEdge(vbm=0, cbm=3,
                                  vbm_distances=[2, 3, 4], cbm_distances=[5, 7])

    assert band_info.band_energies == shifted_band_energies
    assert band_info.band_edge == expected_band_edge
    assert band_info.fermi_level == fermi_level - base_energy


def test_raise_error_when_both_band_edge_fermi_level_absent():
    # [spin][branch idx][band idx][k-path idx]
    with pytest.raises(ViseBandInfoError):
        BandInfo(band_energies=band_energies)


def test_slide_energies_when_band_edge_is_none():
    band_info_fermi.slide_energies(base_energy=base_energy)
    assert band_info_fermi.band_edge is None
    assert band_info_fermi.fermi_level == fermi_level - base_energy


def test_slide_energies_when_fermi_is_none():
    band_info_edge.slide_energies(base_energy=base_energy)
    assert band_info_edge.band_edge.vbm == band_edge.vbm - base_energy
    assert band_info_edge.fermi_level is None


def test_band_mpl_defaults():
    band_defaults = BandMplSettings()

    assert band_defaults.colors == colors
    assert next(band_defaults.linewidth) == 1.0
    assert band_defaults.circle_size == 70
    assert band_defaults.circle_colors == ["pink", "green"]
    assert band_defaults.title_font_size == 15
    assert band_defaults.label_font_size == 15

    band_defaults = BandMplSettings(colors=["black"],
                                    linewidth=[2.0],
                                    circle_size=200,
                                    circle_colors=["blue"],
                                    title_font_size=30,
                                    label_font_size=40)

    assert band_defaults.colors == ["black"]
    assert next(band_defaults.linewidth) == 2.0
    assert band_defaults.circle_size == 200
    assert band_defaults.circle_colors == ["blue"]
    assert band_defaults.title_font_size == 30
    assert band_defaults.label_font_size == 40


def test_add_band_structures(mock_band_plt_list):
    mock_plt, _ = mock_band_plt_list
    linewidth = BandMplSettings().linewidth
    # 1st branch, 1st band
    args = {"color": colors[0], "linewidth": next(linewidth)}
    mock_plt.plot.assert_any_call(distances[0], shifted_band_energies[0][0][0], **args)

    # 2nd branch, 1st band
    args.pop("label", None)
    mock_plt.plot.assert_any_call(distances[1], shifted_band_energies[1][0][0], **args)

    # 1st branch, 2nd band
    mock_plt.plot.assert_any_call(distances[0], shifted_band_energies[0][0][1], **args)

    # 2nd branch, 2nd band
    mock_plt.plot.assert_any_call(distances[1], shifted_band_energies[1][0][1], **args)


def test_band_plot_info_add(band_plot_info, band_info_set):
    band_plot_info_2 = BandPlotInfo(band_info_set,
                                    distances,
                                    x_ticks,
                                    "a")
    added = band_plot_info + band_plot_info_2
    assert added.band_info_set[0] == band_info_set[0]
    assert added.band_info_set[1] == band_info_set[0]
    assert added.distances_by_branch == distances
    assert added.x_ticks == x_ticks
    assert added.title == title  # title is set to the original one.


def test_add_band_edge_circles(mock_band_plt_list, band_info_set):
    mock_plt, _ = mock_band_plt_list
    edge = band_info_set[0].band_edge
    defaults = BandMplSettings()

    mock_plt.scatter.assert_any_call(edge.vbm_distances[0], 0,
                                     **defaults.circle(0))
    mock_plt.scatter.assert_any_call(edge.cbm_distances[0], 3.0,
                                     **defaults.circle(0))
    mock_plt.scatter.assert_any_call(edge.cbm_distances[1], 3.0,
                                     **defaults.circle(0))
    mock_plt.axhline.assert_any_call(y=0.0, **defaults.hline)
    mock_plt.axhline.assert_any_call(y=3.0, **defaults.hline)


def test_figure_legends(mock_band_plt_list):
    mock_plt, _ = mock_band_plt_list
#    mock_plt.legend.assert_called_once_with(loc="lower right")
    mock_plt.legend.assert_not_called()


def test_set_x_range(mock_band_plt_list):
    mock_plt, _ = mock_band_plt_list
    mock_plt.xlim.assert_called_with(distances[0][0], distances[-1][-1])


def test_set_y_range(mock_band_plt_list):
    mock_plt, _ = mock_band_plt_list
    mock_plt.ylim.assert_called_with(y_range[0], y_range[1])


def test_set_labels(mock_band_plt_list):
    mock_plt, _ = mock_band_plt_list
    defaults = BandMplSettings()
    mock_plt.xlabel.assert_called_with("Wave vector",
                                       size=defaults.label_font_size)
    mock_plt.ylabel.assert_called_with("Energy (eV)",
                                       size=defaults.label_font_size)


def test_set_x_tics(mock_band_plt_list):
    mock_plt, mock_axis = mock_band_plt_list

    mock_axis.set_xticks.assert_called_once_with(x_ticks.distances)
    mock_axis.set_xticklabels.assert_called_once_with(x_ticks.labels)
    mock_plt.axvline.assert_any_call(x=x_ticks.distances[2], linestyle="-")
    mock_plt.axvline.assert_any_call(x=x_ticks.distances[1], linestyle="--")


def test_set_title(mock_band_plt_list):
    mock_plt, _ = mock_band_plt_list
    defaults = BandMplSettings()
    mock_plt.title.assert_called_once_with(title, size=defaults.title_font_size)


def test_set_float_to_int_formatter(mock_band_plt_list):
    _, mock_axis = mock_band_plt_list
    print(mock_axis)
    mock_axis.yaxis.set_major_formatter.assert_called_once_with(
        float_to_int_formatter)


def test_plot_tight_layout(mock_band_plt_list):
    mock_plt, _ = mock_band_plt_list
    mock_plt.tight_layout.assert_called_once_with()


#@pytest.mark.skip()
def test_draw_bands(band_plot_info):
    band_plotter = BandPlotter(band_plot_info, y_range)
    band_plotter.construct_plot()
    band_plotter.plt.show()


@pytest.mark.parametrize("ref_energy,subtracted_energy",
                         [(1.0, 1.0), (None, band_edge.vbm)])
def test_reference_energy(ref_energy, subtracted_energy, mocker, band_plot_info):
    mock_plt = mocker.patch("vise.analyzer.plot_band.plt", auto_spec=True)
    mock_axis = MagicMock()
    mock_plt.gca.return_value = mock_axis

    plotter = BandPlotter(band_plot_info, y_range, ref_energy)
    plotter.construct_plot()

    band_defaults = BandMplSettings()
    mock_plt.axhline.assert_called_with(y=band_edge.cbm - subtracted_energy,
                                        **band_defaults.hline)


@pytest.fixture
def two_band_set():
    first_branch_energies = [[[-3.0, -2, -1, -1, -1, -2, -3],
                               [7.0, 6, 5, 4, 3, 2, 3]],
                              [[-2.0, -1, 0, 0, 0, -1, -2],
                               [6.0, 5, 4, 3, 2, 1, 2]]]
    second_branch_energies = [[[-4, -5, -6, -7],
                               [5, 2, 7, 8]],
                              [[-3, -4, -5, -6],
                               [4, 1, 6, 7]]]
    first_band_energies = [first_branch_energies, second_branch_energies]

    band_edge_1 = BandEdge(vbm=0, cbm=1,
                           vbm_distances=[2, 3, 4], cbm_distances=[5, 7])
    band_info_1 = BandInfo(band_energies=first_band_energies,
                           band_edge=band_edge_1,
                           fermi_level=None)

    first_branch_energies = [[[-5, -5, -5, -5, -4, -5, -5],
                              [ 5,  5,  5,  5,  4,  5,  5]],
                             [[-6, -6, -6, -6, -5, -6, -6],
                              [6, 6, 6, 6, 5, 6, 6]]]
    second_branch_energies = [[[-5, -5, -5, -5],
                               [5, 5, 5, 5]],
                              [[-6, -6, -6, -6],
                               [6, 6, 6, 6]]]
    second_band_energies = [first_branch_energies, second_branch_energies]

    band_edge_2 = BandEdge(vbm=-4, cbm=4,
                           vbm_distances=[4], cbm_distances=[4])
    band_info_2 = BandInfo(band_energies=second_band_energies,
                           band_edge=band_edge_2,
                           fermi_level=None)

    return [band_info_1, band_info_2]


@pytest.fixture
def mock_two_bands_plt_axis(two_band_set, mocker):
    mock_plt = mocker.patch("vise.analyzer.plot_band.plt", auto_spec=True)
    mock_axis = MagicMock()
    mock_plt.gca.return_value = mock_axis
    band_plot_info = BandPlotInfo(two_band_set, distances, x_ticks, title)
    plotter = BandPlotter(band_plot_info, y_range)
    plotter.construct_plot()
    return mock_plt, mock_axis


def test_check_color_generator_with_last_call(
        two_band_set, mock_two_bands_plt_axis):
    mock_plt, _ = mock_two_bands_plt_axis
    mock_plt.plot.assert_called_with(
        distances[1],
        two_band_set[-1].band_energies[-1][-1][-1],
        color=colors[1],
        linewidth=1.0,
        linestyle=":")


def test_draw_two_bands(two_band_set):
    band_plot_info = BandPlotInfo(two_band_set, distances, x_ticks, title)
    band_plotter = BandPlotter(band_plot_info, y_range)
    band_plotter.construct_plot()
    band_plotter.plt.show()


