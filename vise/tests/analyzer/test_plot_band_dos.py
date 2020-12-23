# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from monty.serialization import loadfn
from vise.analyzer.plot_band_dos import BandDosPlotlyPlotter
from vise.util.dash_helper import show_png

try:
    import psutil
    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test(test_data_files):
    band_plot_data = loadfn(str(test_data_files / "KAlSi3O8_band_plot_info.json"))
    dos_plot_data = loadfn(str(test_data_files / "KAlSi3O8_dos_plot_data.json"))
    band_dos_component = BandDosPlotlyPlotter(dos_plot_data, band_plot_data)
#    band_dos_component.fig.show()
    show_png(band_dos_component.fig)
