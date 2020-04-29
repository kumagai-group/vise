# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from unittest.mock import MagicMock

import pytest
from pymatgen import Spin
from vise.analyzer.plot_band import (
    BandPlotter, BandInfo,
    BandEdge, XTicks, BandPlotterDefaults)
from vise.util.matplotlib import float_to_int_formatter

"""
TODO:

+ Raise error when 5 bands are drawn
+ Set color of vbm, cbm, efermi lines
+ Set band gap arrow(s).
+ Enable multiple bands (distances and x_ticks must be the same)
+ To csv
+ Brillouin zone
+ Add spin info to band edge 

DONE:
+ Draw a non-magnetic band with two bands
+ Add x-label and y-label
+ Add x-axis and y-axis
+ Change x-label to unicode.
+ Add vbm and cbm circles
+ Add line when two x_ticks conflict
+ Enable to set y-range
+ Enable spin polarized band
+ Add Fermi level
+ Set horizontal lines (vbm, cbm)
+ Add dashed line for the positions with labels
+ Use set_major_formatter
+ Create a class that holds default values for drawing band structures.
+ Set tuning parameters with BandPlotterDefaults class
+ Create VaspBandPlotter inheriting BandPlotter
+ Move a|b part to VaspBandPlotter
+ From a vasprun.xml and a KPOINTS, generate minimum BandPlotter data.
+ Set title size
+ Set vbm to zero of energy.
+ Generate cyclic color for band structures.
+ Draw band edge lines.
+ Add figure legend
"""


@pytest.fixture
def distances():
    return [[0, 1, 2, 3, 4, 5, 6], [6, 7, 8, 9]]


@pytest.fixture
def x_ticks():
    labels = ["A", "B", "C$\\mid$D", "E"]
    distances = [0.0, 5.0, 6.0, 9.0]
    return XTicks(labels, distances)


@pytest.fixture
def y_range():
    return [-10, 10]


@pytest.fixture
def title():
    return "Title"


@pytest.fixture
def reference_energy():
    return 1.0


@pytest.fixture
def band_info():
    # [spin][branch idx][band idx][k-path idx]
    band_energies = {Spin.up: [[[-3.0, -2, -1, -1, -1, -2, -3],
                                [7.0, 6, 5, 4, 3, 2, 3]],
                               [[-4, -5, -6, -7],
                                [5, 2, 7, 8]]]}

    band_edge = BandEdge(vbm=-1, cbm=2,
                         vbm_distances=[2, 3, 4], cbm_distances=[5, 7])
    fermi_level = 1.5
    return BandInfo(band_energies=band_energies,
                    band_edge=band_edge,
                    fermi_level=fermi_level)


def test_band_info_slide_energies(band_info, reference_energy):
    band_info.slide_energies(reference_energy=reference_energy)
    expected_band_energies = {Spin.up: [[[-4.0, -3, -2, -2, -2, -3, -4],
                                          [6.0, 5, 4, 3, 2, 1, 2]],
                                         [[-5, -6, -7, -8],
                                          [4, 1, 6, 7]]]}
    expected_band_edge = BandEdge(vbm=-2, cbm=1,
                                  vbm_distances=[2, 3, 4], cbm_distances=[5, 7])
    expected_fermi_level = 0.5

    assert band_info.band_energies == expected_band_energies
    assert band_info.band_edge == expected_band_edge
    assert band_info.fermi_level == expected_fermi_level


@pytest.fixture
def band_info_set(band_info):
    return [band_info]


@pytest.fixture
def mock_plt_axis(band_info_set, distances, x_ticks, y_range, title, mocker):
    mock_plt = mocker.patch("vise.analyzer.plot_band.plt", auto_spec=True)
    mock_axis = MagicMock()
    mock_plt.gca.return_value = mock_axis
    plotter = BandPlotter(band_info_set, distances, x_ticks, y_range, title)
    plotter.construct_plot()
    return mock_plt, mock_axis


@pytest.fixture
def colors():
    return ['#E15759', '#4E79A7', '#F28E2B', '#76B7B2']


def test_defaults_class(colors):
    band_defaults = BandPlotterDefaults()

    assert band_defaults.band_colors == colors
    assert band_defaults.band_linewidth == 1.0
    assert band_defaults.band_edge_circle_size == 100
    assert band_defaults.title_font_size == 15
    assert band_defaults.label_font_size == 15

    generator = band_defaults.band_args
    assert next(generator) == {"color": colors[0], "linewidth": 1.0}
    assert next(generator) == {"color": colors[1], "linewidth": 1.0}
    assert next(generator) == {"color": colors[2], "linewidth": 1.0}
    assert next(generator) == {"color": colors[3], "linewidth": 1.0}

    band_defaults = BandPlotterDefaults(band_colors=["black"],
                                        band_linewidth=2.0,
                                        band_edge_circle_size=200,
                                        title_font_size=30,
                                        label_font_size=40)

    assert band_defaults.band_colors == ["black"]
    assert band_defaults.band_linewidth == 2.0
    assert band_defaults.band_edge_circle_size == 200
    assert band_defaults.title_font_size == 30
    assert band_defaults.label_font_size == 40


def test_set_band_structures(band_info_set, distances, colors, mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    mock_plt.plot.assert_any_call(distances[0],
                                  band_info_set[0].band_energies[Spin.up][0][0],
                                  color=colors[0], linewidth=1.0,
                                  label="1th")
    mock_plt.plot.assert_any_call(distances[1],
                                  band_info_set[0].band_energies[Spin.up][1][0],
                                  color=colors[0], linewidth=1.0)
    mock_plt.plot.assert_any_call(distances[0],
                                  band_info_set[0].band_energies[Spin.up][0][1],
                                  color=colors[0], linewidth=1.0)
    mock_plt.plot.assert_any_call(distances[1],
                                  band_info_set[0].band_energies[Spin.up][1][1],
                                  color=colors[0], linewidth=1.0)


def test_set_band_edge_circles(band_info_set, mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    band_edge = band_info_set[0].band_edge
    band_defaults = BandPlotterDefaults()

    mock_plt.scatter.assert_any_call(band_edge.vbm_distances[0], band_edge.vbm,
                                     **band_defaults.band_edge_circles_args)
    mock_plt.scatter.assert_any_call(band_edge.cbm_distances[0], band_edge.cbm,
                                     **band_defaults.band_edge_circles_args)
    mock_plt.scatter.assert_any_call(band_edge.cbm_distances[1], band_edge.cbm,
                                     **band_defaults.band_edge_circles_args)
    mock_plt.axhline.assert_any_call(y=band_edge.vbm)
    mock_plt.axhline.assert_any_call(y=band_edge.cbm)


def test_figure_legends(mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    mock_plt.legend.assert_called_once_with(loc="lower right")


def test_set_fermi_level(band_info_set, mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    mock_plt.axhline.assert_called_with(y=band_info_set[0].fermi_level)


def test_set_x_range(distances, mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    mock_plt.xlim.assert_called_with(distances[0][0], distances[-1][-1])


def test_set_y_range(y_range, mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    mock_plt.ylim.assert_called_with(y_range[0], y_range[1])


def test_set_axes(mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    defaults = BandPlotterDefaults()
    mock_plt.xlabel.assert_called_with("Wave vector",
                                       size=defaults.label_font_size)
    mock_plt.ylabel.assert_called_with("Energy (eV)",
                                       size=defaults.label_font_size)


def test_set_x_tics(x_ticks, mock_plt_axis):
    mock_plt, mock_axis = mock_plt_axis

    mock_axis.set_xticks.assert_called_once_with(x_ticks.distances)
    mock_axis.set_xticklabels.assert_called_once_with(x_ticks.labels)
    mock_plt.axvline.assert_any_call(x=x_ticks.distances[2], linestyle="-")
    mock_plt.axvline.assert_any_call(x=x_ticks.distances[1], linestyle="--")


def test_set_title(title, mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    defaults = BandPlotterDefaults()
    mock_plt.title.assert_called_once_with(title, size=defaults.title_font_size)


def test_set_float_to_int_formatter(mock_plt_axis):
    _, mock_axis = mock_plt_axis
    mock_axis.yaxis.set_major_formatter.assert_called_once_with(
        float_to_int_formatter)


def test_plot_tight_layout(mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    mock_plt.tight_layout.assert_called_once_with()


#@pytest.mark.skip()
def test_draw_bands(band_info_set, distances, x_ticks, y_range, title):
    band_plotter = \
        BandPlotter(band_info_set, distances, x_ticks, y_range, title)
    band_plotter.construct_plot()
    band_plotter.plt.show()


@pytest.mark.parametrize("reference_energy,subtracted_energy",
                         [(1.0, 1.0), (None, -1.0)])
def test_reference_energy(band_info_set, distances, x_ticks, y_range,
                          title, reference_energy, subtracted_energy, mocker):
    original_fermi_level = band_info_set[0].fermi_level
    shifted_fermi_level = original_fermi_level - subtracted_energy

    mock_plt = mocker.patch("vise.analyzer.plot_band.plt", auto_spec=True)
    mock_axis = MagicMock()
    mock_plt.gca.return_value = mock_axis

    plotter = BandPlotter(band_info_set, distances, x_ticks, y_range, title,
                          reference_energy)
    plotter.construct_plot()

    mock_plt.axhline.assert_called_with(y=shifted_fermi_level)