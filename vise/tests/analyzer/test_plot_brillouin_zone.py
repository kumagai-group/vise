# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest
from vise.analyzer.plot_brillouin_zone import BZPlotInfo, BZPlotlyPlotter
from vise.tests.conftest import assert_msonable


@pytest.fixture
def bz_plot_info():
    f1 = [0.5, 0.5, 0.5]
    f2 = [-0.5, 0.5, 0.5]
    f3 = [0.5, -0.5, 0.5]
    f4 = [0.5, 0.5, -0.5]
    f5 = [0.5, -0.5, -0.5]
    f6 = [-0.5, 0.5, -0.5]
    f7 = [-0.5, -0.5, 0.5]
    f8 = [-0.5, -0.5, -0.5]

    faces = [[f1, f3, f5, f4], [f2, f6, f8, f7], [f1, f2, f6, f4],
             [f3, f5, f8, f7], [f1, f2, f7, f3], [f4, f5, f8, f6]]
    return BZPlotInfo(faces=faces,
                      labels={"Γ": [0.0, 0.0, 0.0],
                              "X": [0.5, 0.0, 0.0],
                              "M": [0.5, 0.5, 0.0]})


def test_bz_plot_info(bz_plot_info):
    assert_msonable(bz_plot_info)


def test_bz_plot_plotly(bz_plot_info):
    fig = BZPlotlyPlotter(bz_plot_info).create_figure()
    fig.show()
