# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

import pytest
from pymatgen.core import Composition
from pymatgen.io.vasp import Vasprun

from vise.analyzer.plot_band import BandEdgeForPlot, XTicks, BandMplSettings
from vise.analyzer.plot_brillouin_zone import BZPlotlyPlotter
from vise.analyzer.vasp.plot_band import greek_to_unicode, italic_to_roman, \
    BandPlotInfoFromVasp
from vise.analyzer.plot_band import BandMplPlotter
import numpy as np
from vise.util.dash_helper import show_png

try:
    import psutil
    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


def test_greek_to_unicode():
    assert greek_to_unicode("GAMMA") == "Γ"
    assert greek_to_unicode("SIGMA") == "Σ"
    assert greek_to_unicode("DELTA") == "Δ"


def test_italic_to_roman():
    assert italic_to_roman(r"S$\mid$$S_0$") == "S$\mid$${\\rm S}_0$"


test_data = [(True, None),
             (False, BandEdgeForPlot(vbm=-100, cbm=100,
                                     vbm_distances=[0], cbm_distances=[1, 2]))]


@pytest.mark.parametrize("is_metal,expected_band_edge", test_data)
def test_vasp_band_plotter(is_metal, expected_band_edge, mocker):
    mock_bs = mocker.MagicMock()
    mock_bs.efermi = 10
    mock_bs.is_metal.return_value = is_metal

    stub_vasprun = mocker.MagicMock()
    stub_vasprun.final_structure.composition = Composition("MgO2")
    stub_vasprun.get_band_structure.return_value = mock_bs

    energy = {"1": [np.array([[0.1], [0.2], [0.3]])]}
    distances = [np.array([0.0, 0.1, 0.2])]
    labels = ["A", "$A_0$", "GAMMA"]
    label_distances = [0.0, 0.1, 0.2]
    plot_data = {"energy": energy,
                 "distances": distances,
                 "vbm": [[0, -100]],
                 "cbm": [[1, 100], [2, 100]]}

    mock_bsp = mocker.patch("vise.analyzer.vasp.plot_band.BSPlotter",
                            auto_spec=True)
    mock_bsp.return_value.bs_plot_data.return_value = plot_data
    mock_bsp.return_value.get_ticks_old.return_value = \
        {"label": labels, "distance": [0.0, 0.1, 0.2]}

    band_plot_info = BandPlotInfoFromVasp(stub_vasprun, "KPOINTS")
    plot_info = band_plot_info.make_band_plot_info()

    expected_x_ticks = XTicks(labels=["A", "${\\rm A}_0$", "Γ"],
                              distances=label_distances)

    assert plot_info.band_energy_infos["1"].band_energies == \
           [[[[0.1], [0.2], [0.3]]]]
    assert plot_info.band_energy_infos["1"].band_edge == expected_band_edge
    assert plot_info.distances_by_branch == [[0.0, 0.1, 0.2]]
    assert plot_info.x_ticks == expected_x_ticks
    assert plot_info.title == "MgO$_{2}$"


def test_draw_band_plotter_with_actual_vasp_files(test_data_files: Path):
    vasprun_file = str(test_data_files / "KO2_band_vasprun.xml")
    kpoints_file = str(test_data_files / "KO2_band_KPOINTS")
    vasprun = Vasprun(vasprun_file)
    plot_info = BandPlotInfoFromVasp(vasprun, kpoints_file).make_band_plot_info()
    plotter = BandMplPlotter(plot_info, [-10, 10])
    plotter.construct_plot()
    plotter.plt.show()


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_bz_plotter_with_actual_vasp_files(test_data_files: Path):
    vasprun_file = str(test_data_files / "H_band_vasprun.xml")
    kpoints_file = str(test_data_files / "H_band_KPOINTS")
    vasprun = Vasprun(vasprun_file)
    band_plot_info = BandPlotInfoFromVasp(vasprun, kpoints_file)
    bz_plot_info = band_plot_info.make_bz_plot_info()
    fig = BZPlotlyPlotter(bz_plot_info).create_figure()
    # fig.show()
    show_png(fig)


# def test_irrep_actual_vasp_files(test_data_files: Path):


def test_draw_two_bands(test_data_files: Path):
    mpl_settings = BandMplSettings(linewidth=[0.3, 1.0],
                                   circle_size=50,
                                   show_legend=True)

    vasprun_file = str(test_data_files / "CdAs2O6-vasprun1.xml")
    vasprun2_file = str(test_data_files / "CdAs2O6-vasprun2.xml")
    kpoints_file = str(test_data_files / "CdAs2O6-KPOINTS")
    vasprun = Vasprun(vasprun_file)
    vasprun2 = Vasprun(vasprun2_file)
    plot_info = BandPlotInfoFromVasp(
        vasprun, kpoints_file, first_band_plot_name="pbe",
        vasprun2=vasprun2, second_band_plot_name="dd hybrid",
        energy_window=[-20.0, 20.0]).make_band_plot_info()

    plotter = BandMplPlotter(plot_info, [-10, 10], mpl_defaults=mpl_settings)
    plotter.construct_plot()
    plotter.plt.show()


def test_energy_window(mocker):
    mock_bs = mocker.MagicMock()
    mock_bs.efermi = 0
    mock_bs.is_metal.return_value = True

    stub_vasprun = mocker.MagicMock()
    stub_vasprun.final_structure.composition = Composition("MgO2")
    stub_vasprun.get_band_structure.return_value = mock_bs

    energy = {"1": [np.array([[-0.2, -0.1, -0.3, -0.1],
                              [-0.2, -0.1, -0.3, 0.1],
                              [1.2, 1.3, 1.1, 1.2]])]}
    distances = [[0.0, 1.0]]
    labels = ["A", "GAMMA"]
    label_distances = [0.0, 1.0]
    plot_data = {"energy": energy,
                 "distances": distances,
                 "vbm": None,
                 "cbm": None}

    mock_bsp = mocker.patch("vise.analyzer.vasp.plot_band.BSPlotter",
                            auto_spec=True)
    mock_bsp.return_value.bs_plot_data.return_value = plot_data
    mock_bsp.return_value.get_ticks_old.return_value = \
        {"label": labels, "distance": label_distances}

    plot_info = BandPlotInfoFromVasp(
        stub_vasprun, "KPOINTS", energy_window=[0.0, 1.0]).make_band_plot_info()

    assert (plot_info.band_energy_infos["1"].band_energies
            == [[[[-0.2, -0.1, -0.3, 0.1]]]])

