# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path
from unittest.mock import MagicMock

import pytest
from pymatgen import Spin
from pymatgen.io.vasp import Vasprun
from vise.analyzer.plot_band import (
    greek_to_unicode, italic_to_roman, BandPlotter, BandInfo,
    BandEdge, XTicks, BandPlotterDefaults,
    VaspBandPlotter)
from vise.util.matplotlib import float_to_int_formatter

"""
TODO:
+ Enable multiple bands (distances and x_ticks must be the same)
+ Generate cyclic color and linestyle for band structures.
+ Consider if the structure is parsed or not.
+ Set title

DONE:
+ Draw a non-magnetic band with a single band
+ Draw a non-magnetic band with two bands
+ Add x-label and y-label
+ Add x-axis and y-axis
+ Change x-label to Unicodes.
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
def band_info_set():
    # band_energies = [{Spin.up: [[-3.0, -2, -1,  0, -1, -2, -3, -4, -5, -6],
    #                            [ 7.0,  6,  5,  4,  3,  2,  3,  2,  3,  4]]}]
    band_energies = {Spin.up: [[[-3.0, -2, -1,  0, -1, -2, -3],
                                [ 7.0,  6,  5,  4,  3,  2,  3]],
                               [[-4, -5, -6, -7],
                                [ 5,  2,  7,  8]]]}

    band_edge = BandEdge(vbm=0, cbm=2, vbm_distances=[3], cbm_distances=[5, 7])
    fermi_level = 1.5
    return [BandInfo(band_energies=band_energies,
                     band_edge=band_edge,
                     fermi_level=fermi_level)]


@pytest.fixture
def mock_plt_axis(band_info_set, distances, x_ticks, y_range, mocker):
    mock_plt = mocker.patch("vise.analyzer.plot_band.plt", auto_spec=True)
    mock_axis = MagicMock()
    mock_plt.gca.return_value = mock_axis
    plotter = BandPlotter(band_info_set, distances, x_ticks, y_range)
    plotter.construct_plot()
    return mock_plt, mock_axis


def test_defaults_class():
    band_defaults = BandPlotterDefaults()

    assert band_defaults.band_linestyle == ["-"]
    assert band_defaults.band_colors == ["red", "blue"]
    assert band_defaults.band_linewidth == 1.0
    assert band_defaults.band_edge_circle_size == 100

    band_defaults = BandPlotterDefaults(band_linestyle=["--"],
                                        band_colors=["black"],
                                        band_linewidth=2.0,
                                        band_edge_circle_size=200)

    assert band_defaults.band_linestyle == ["--"]
    assert band_defaults.band_colors == ["black"]
    assert band_defaults.band_linewidth == 2.0
    assert band_defaults.band_edge_circle_size == 200


def test_set_band_structures(band_info_set, distances, mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    band_defaults = BandPlotterDefaults()
    mock_plt.plot.assert_any_call(distances[0],
                                  band_info_set[0].band_energies[Spin.up][0][0],
                                  **band_defaults.band_args)
    mock_plt.plot.assert_any_call(distances[1],
                                  band_info_set[0].band_energies[Spin.up][1][0],
                                  **band_defaults.band_args)
    mock_plt.plot.assert_any_call(distances[0],
                                  band_info_set[0].band_energies[Spin.up][0][1],
                                  **band_defaults.band_args)
    mock_plt.plot.assert_any_call(distances[1],
                                  band_info_set[0].band_energies[Spin.up][1][1],
                                  **band_defaults.band_args)


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
    mock_plt.xlabel.assert_called_with("Wave vector")
    mock_plt.ylabel.assert_called_with("Energy (eV)")


def test_set_x_tics(x_ticks, mock_plt_axis):
    mock_plt, mock_axis = mock_plt_axis

    mock_axis.set_xticks.assert_called_once_with(x_ticks.distances)
    mock_axis.set_xticklabels.assert_called_once_with(x_ticks.labels)
    mock_plt.axvline.assert_any_call(x=x_ticks.distances[2], linestyle="-")
    mock_plt.axvline.assert_any_call(x=x_ticks.distances[1], linestyle="--")


def test_set_float_to_int_formatter(mock_plt_axis):
    _, mock_axis = mock_plt_axis
    mock_axis.yaxis.set_major_formatter.assert_called_once_with(
        float_to_int_formatter)


def test_plot_tight_layout(mock_plt_axis):
    mock_plt, _ = mock_plt_axis
    mock_plt.tight_layout.assert_called_once_with()


#@pytest.mark.skip()
def test_draw_bands(band_info_set, distances, x_ticks, y_range):
    band_plotter = BandPlotter(band_info_set, distances, x_ticks, y_range)
    band_plotter.construct_plot()
    band_plotter.plt.show()


def test_greek_to_unicode():
    assert greek_to_unicode("GAMMA") == "Γ"
    assert greek_to_unicode("SIGMA") == "Σ"
    assert greek_to_unicode("DELTA") == "Δ"


def test_italic_to_roman():
    assert italic_to_roman(r"S$\mid$$S_0$") == "S$\mid$${\\rm S}_0$"


test_data = [(True,  None),
             (False, BandEdge(vbm=-100, cbm=100,
                              vbm_distances=[0], cbm_distances=[1, 2]))]


@pytest.mark.parametrize("is_metal,expected", test_data)
def test_vasp_band_plotter(is_metal, expected, mocker):
    stub_vasprun = mocker.MagicMock(spec=Vasprun)
    mock_bs = mocker.MagicMock()
    mock_bs.efermi = 10
    mock_bs.is_metal.return_value = is_metal

    stub_vasprun.get_band_structure.return_value = mock_bs

    mock_bsp = mocker.patch("vise.analyzer.plot_band.BSPlotter", auto_spec=True)
    energy = [{Spin.up: [[0.1]]}]
    distances = [[0.0, 0.1, 0.2]]
    labels = ["A", "$A_0$", "GAMMA"]
    label_distances = [0.0, 0.1, 0.2]
    plot_data = {"ticks": {"label": labels, "distance": label_distances},
                 "energy": energy,
                 "distances": distances,
                 "vbm": [[0, -100]],
                 "cbm": [[1, 100], [2, 100]]}
    mock_bsp().bs_plot_data.return_value = plot_data

    plotter = VaspBandPlotter(stub_vasprun, "KPOINTS")

    assert plotter.band_info_set == [BandInfo(band_energies={Spin.up: [[[0.1]]]},
                                              band_edge=expected,
                                              fermi_level=10)]
    assert plotter.distances_by_branch == distances
    assert plotter.x_ticks == XTicks(labels=["A", "${\\rm A}_0$", "Γ"],
                                     distances=label_distances)
    assert plotter.y_range == [-10, 10]


def test_vasp_band_plotter_with_actual_files(test_data_files: Path):
    vasprun = Vasprun(str(test_data_files / "KO2_band_vasprun.xml"))
    plotter = VaspBandPlotter(vasprun, str(test_data_files / "KO2_band_KPOINTS"))
    plotter.construct_plot()
    plotter.plt.show()

