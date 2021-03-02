# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest
from vise.analyzer.plot_brillouin_zone import BZPlotInfo, BZPlotlyPlotter, \
    pairwise
from vise.tests.helpers.assertion import assert_msonable
from vise.util.dash_helper import show_png

try:
    import psutil
    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


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
                      band_paths=[
                          [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
                          [[0.5, 0.0, 0.0], [0.5, 0.5, 0.0]],
                          [[0.5, 0.5, 0.0], [0.0, 0.0, 0.0]]],
                      rec_lat_vec=[[1, 1, 0], [1, 0, 1], [0, 1, 1]],
                      labels={"Γ": {"cart":[0.0, 0.0, 0.0], "frac":[0.0, 0.0, 0.0]},
                              "X": {"cart":[0.5, 0.0, 0.0], "frac":[0.5, 0.0, 0.0]},
                              "M": {"cart":[0.5, 0.5, 0.0], "frac":[0.5, 0.5, 0.0]}})


def test_pairwise():
    actual = list(pairwise([[-0.43, -0.74, -0.643], [-0.43, -0.74, 0.643],
                            [-0.858, 0.0, 0.643], [-0.858, 0.0, -0.643]]))
    expected = [([-0.43, -0.74, -0.643], [-0.43, -0.74, 0.643]),
                ([-0.43, -0.74, 0.643], [-0.858, 0.0, 0.643]),
                ([-0.858, 0.0, 0.643], [-0.858, 0.0, -0.643]),
                ([-0.858, 0.0, -0.643], [-0.43, -0.74, -0.643])]
    assert actual == expected


def test_bz_plot_info(bz_plot_info):
    assert_msonable(bz_plot_info)


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_bz_plot_plotly(bz_plot_info):
    fig = BZPlotlyPlotter(bz_plot_info).create_figure()
    show_png(fig)
