# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from monty.serialization import loadfn
from vise.analyzer.plot_band import BandPlotInfo, Irrep, Irreps
from vise.analyzer.plot_band_dos import BandDosPlotlyPlotter, \
    plotly_sanitize_label
from vise.util.dash_helper import show_png

try:
    import psutil

    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


def test_plotly_sanitize_label():
    assert plotly_sanitize_label(r"A$\mid$B") == "A|B"
    assert plotly_sanitize_label("A_0") == "A<sub>0</sub>"
    assert plotly_sanitize_label("A_i1") == "A<sub>i1</sub>"




@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_no_data():
    band_dos_component = BandDosPlotlyPlotter()
#    band_dos_component.fig.show()
    show_png(band_dos_component.fig)


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_plotly_band_dos_actual_files(test_data_files):
    band_plot_data: BandPlotInfo = loadfn(
        str(test_data_files / "MgSe_band_plot_info.json"))
    dos_plot_data = loadfn(str(test_data_files / "MgSe_dos_plot_data.json"))
    irrep = Irrep(frac_coords=[0.0, 0.0, 0.0], symbols=["Γ1+"], energies=[2.96],
                  degeneracies=[1])
    irreps = Irreps(sg_num=225, irreps={"Γ": irrep})
    band_plot_data.band_energy_infos["1"].irreps = irreps
    band_dos_component = BandDosPlotlyPlotter(dos_plot_data, band_plot_data)
    band_dos_component.fig.show()
    # show_png(band_dos_component.fig)
