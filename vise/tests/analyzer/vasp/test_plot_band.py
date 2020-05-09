# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

import pytest
from pymatgen import Composition
from pymatgen.io.vasp import Vasprun

from vise.analyzer.plot_band import BandEdge, XTicks
from vise.analyzer.vasp.plot_band import greek_to_unicode, italic_to_roman, \
    VaspBandPlotInfo
from vise.analyzer.plot_band import BandPlotter


def test_greek_to_unicode():
    assert greek_to_unicode("GAMMA") == "Γ"
    assert greek_to_unicode("SIGMA") == "Σ"
    assert greek_to_unicode("DELTA") == "Δ"


def test_italic_to_roman():
    assert italic_to_roman(r"S$\mid$$S_0$") == "S$\mid$${\\rm S}_0$"


test_data = [(True, None),
             (False, BandEdge(vbm=-100, cbm=100,
                              vbm_distances=[0], cbm_distances=[1, 2]))]


@pytest.mark.parametrize("is_metal,expected_band_edge", test_data)
def test_vasp_band_plotter(is_metal, expected_band_edge, mocker):
    mock_bs = mocker.MagicMock()
    mock_bs.efermi = 10
    mock_bs.is_metal.return_value = is_metal

    stub_vasprun = mocker.MagicMock()
    stub_vasprun.final_structure.composition = Composition("MgO2")
    stub_vasprun.get_band_structure.return_value = mock_bs

    energy = [{"1": [[0.1]]}]
    distances = [[0.0, 0.1, 0.2]]
    labels = ["A", "$A_0$", "GAMMA"]
    label_distances = [0.0, 0.1, 0.2]
    plot_data = {"ticks": {"label": labels, "distance": label_distances},
                 "energy": energy,
                 "distances": distances,
                 "vbm": [[0, -100]],
                 "cbm": [[1, 100], [2, 100]]}

    mock_bsp = mocker.patch("vise.analyzer.vasp.plot_band.BSPlotter", auto_spec=True)
    mock_bsp().bs_plot_data.return_value = plot_data

    plot_info = VaspBandPlotInfo(stub_vasprun, "KPOINTS")

    expected_x_ticks = XTicks(labels=["A", "${\\rm A}_0$", "Γ"],
                              distances=label_distances)

    assert plot_info.band_info_set[0].band_energies == [[[[0.1]]]]
    assert plot_info.distances_by_branch == distances
    assert plot_info.x_ticks == expected_x_ticks
    assert plot_info.title == "MgO$_{2}$"


def test_draw_band_plotter_with_actual_vasp_files(test_data_files: Path):
    vasprun_file = str(test_data_files / "KO2_band_vasprun.xml")
    kpoints_file = str(test_data_files / "KO2_band_KPOINTS")
    vasprun = Vasprun(vasprun_file)
    plot_info = VaspBandPlotInfo(vasprun, kpoints_file)
    plotter = BandPlotter(plot_info, [-10, 10])
    plotter.construct_plot()
    plotter.plt.show()