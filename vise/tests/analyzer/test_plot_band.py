# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from vise.analyzer.plot_band import (
    greek_to_unicode, italic_to_roman, pairwise, BandPlotter, BandInfo, BandEdge, XTics
)
from unittest.mock import MagicMock

from pymatgen import Spin
from vise.util.matplotlib import float_to_int_formatter

"""
TODO:
+ Set tuning parameters

+ Enable multiple bands
+ Create VaspBandPlotter inheriting BandPlotter
+ Move a|b part to VaspBandPlotter

DONE:
+ Draw a non-magnetic band with a single band
+ Draw a non-magnetic band with two bands
+ Add x-label and y-label
+ Add x-axis and y-axis
+ Change x-label to Unicodes.
+ Add vbm and cbm circles
+ Add line when two x_tics conflict
+ Enable to set y-range
+ Enable spin polarized band
+ Add Fermi level
+ Set horizontal lines (vbm, cbm)
+ Add dashed line for the positions with labels
+ Use set_major_formatter
"""


@pytest.fixture
def distance():
    return [i for i in range(10)]


@pytest.fixture
def x_tics():
    labels = ["A", "B", "C|D", "E"]
    distances = [0.0, 5.0, 6.0, 9.0]
    return XTics(labels, distances)


@pytest.fixture
def y_range():
    return [-10, 10]


@pytest.fixture
def band_set():
    band_energies = {Spin.up: [[-3.0, -2, -1,  0, -1, -2, -3, -4, -5, -6],
                               [ 7.0,  6,  5,  4,  3,  2,  3,  2,  3,  4]]}
    band_edge = BandEdge(vbm=0, cbm=2, vbm_distances=[3], cbm_distances=[5, 7])
    fermi_level = 1.5
    return [BandInfo(band_energies=band_energies,
                     band_edge=band_edge,
                     fermi_level=fermi_level)]


@pytest.fixture
def mock_plt(band_set, distance, x_tics, y_range, mocker):
    mock_plt = mocker.patch("vise.analyzer.plot_band.plt", auto_spec=True)
    plotter = BandPlotter(band_set, distance, x_tics, y_range)
    plotter.construct_plot()
    return mock_plt


def test_set_band_structures(band_set, distance, mock_plt):
    mock_plt.plot.assert_any_call(distance,
                                  band_set[0].band_energies[Spin.up][0],
                                  **band_set[0].band_structure_args,
                                  )
    mock_plt.plot.assert_any_call(distance,
                                  band_set[0].band_energies[Spin.up][1],
                                  **band_set[0].band_structure_args,
                                  )


def test_set_band_edge_circles(band_set, mock_plt):
    mock_plt.scatter.assert_any_call(band_set[0].band_edge.vbm_distances[0],
                                     band_set[0].band_edge.vbm)
    mock_plt.scatter.assert_any_call(band_set[0].band_edge.cbm_distances[0],
                                     band_set[0].band_edge.cbm)
    mock_plt.scatter.assert_any_call(band_set[0].band_edge.cbm_distances[1],
                                     band_set[0].band_edge.cbm)
    mock_plt.axhline.assert_any_call(y=band_set[0].band_edge.vbm)
    mock_plt.axhline.assert_any_call(y=band_set[0].band_edge.cbm)


def test_set_fermi_level(band_set, mock_plt):
    mock_plt.axhline.assert_called_with(y=band_set[0].fermi_level)


def test_set_x_range(distance, mock_plt):
    mock_plt.xlim.assert_called_with(distance[0], distance[-1])


def test_set_y_range(y_range, mock_plt):
    mock_plt.ylim.assert_called_with(y_range[0], y_range[1])


def test_set_axes(mock_plt):
    mock_plt.xlabel.assert_called_with("Wave vector")
    mock_plt.ylabel.assert_called_with("Energy (eV)")


@pytest.fixture
def mocker_axis(band_set, distance, x_tics, y_range, mocker):
    mock_plt = mocker.patch("vise.analyzer.plot_band.plt", auto_spec=True)
    mock_axis = MagicMock()
    mock_plt.gca.return_value = mock_axis
    plotter = BandPlotter(band_set, distance, x_tics, y_range)
    plotter.construct_plot()
    return mock_plt, mock_axis


def test_set_x_tics(x_tics, mocker_axis):
    mock_plt, mock_axis = mocker_axis

    mock_axis.set_xticks.assert_called_once_with(x_tics.distances)
    mock_axis.set_xticklabels.assert_called_once_with(x_tics.labels)
    mock_plt.axvline.assert_any_call(x=x_tics.distances[2])
    mock_plt.axvline.assert_any_call(x=x_tics.distances[1], linestyle="--")


def test_set_float_to_int_formatter(mocker_axis):
    mock_plt, mock_axis = mocker_axis
    mock_axis.xaxis.set_major_formatter.assert_called_once_with(
        float_to_int_formatter)
    mock_axis.yaxis.set_major_formatter.assert_called_once_with(
        float_to_int_formatter)


def test_plot_tight_layout(mock_plt):
    mock_plt.tight_layout.assert_called_once_with()


def test_draw_bands(band_set, distance, x_tics, y_range):
    band_plotter = BandPlotter(band_set, distance, x_tics, y_range)
    band_plotter.construct_plot()
    band_plotter.plt.show()





def test_greek_to_unicode():
    assert greek_to_unicode("GAMMA") == "Γ"
    assert greek_to_unicode("SIGMA") == "Σ"
    assert greek_to_unicode("DELTA") == "Δ"


def test_italic_to_roman():
    assert italic_to_roman("A_0") == "${\\rm A}_0$"






def test_pairwise():
    assert list(pairwise([0, 1, 2, 3])) == [(0, 1), (1, 2), (2, 3)]




# @pytest.fixture
# def xtics_with_overlap():
#     labels = ["A", "B", "C", "D", "E", "F"]
#     distances = [0.0, 5.0, 5.0, 7.0, 7.0, 9.0]
#     return Xtics(labels, distances)


# def test_assert_incorrect_distances():
#     labels = ["A", "B", "C"]
#     distances = [0.0, 6.0, 5.0]
#     with pytest.raises(AssertionError):
#         Xtics(labels, distances)


# def test_xtics_with_overlap(xtics_with_overlap):
#     assert xtics_with_overlap.labels == ["A", "B|C", "D|E", "F"]
#     assert xtics_with_overlap.distances == [0.0, 5.0, 7.0, 9.0]
#     assert xtics_with_overlap.vertical_line_positions == [5.0, 7.0]


