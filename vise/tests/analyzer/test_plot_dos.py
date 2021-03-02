# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
import pytest

from vise.analyzer.dos_data import DosBySpinEnergy, DosPlotData
from vise.analyzer.plot_dos import DosPlotter, DosMplSettings
from vise.tests.helpers.assertion import assert_msonable
from vise.util.matplotlib import float_to_int_formatter

"""
TODO:
+ Create PDos class composed of orbital-decomposed pdos
+ Create from FullDos instance with total and atom- and orbital-decomposed pdos
+ Add labels
+ Allow to set vbm, cbm, efermi
+ Shift energy zero to vbm or efermi.
+ Create DosInfo from vasp output.

DONE:
"""

relative_energies = [i for i in range(-5, 6)]

total_up = [0.0] * 2 + [2] * 2 + [4] * 4 + [2] * 2 + [0]
total_down = [0.0] * 4 + [6] * 6 + [0]
h_s_up = [0.0] * 2 + [2] * 8 + [0]
h_s_down = [0.0] * 4 + [1] * 6 + [0]
h_p_up = [0.0] * 2 + [3] * 8 + [0]
h_p_down = [0.0] * 3 + [2] * 7 + [0]

xlim = [-6, 6]
ylim_set = [[-10, 10], [-5, 5]]
colors = DosMplSettings().colors


doses = [[DosBySpinEnergy("", [total_up, total_down])],
         [DosBySpinEnergy("s", np.array([h_s_up, h_s_down])),
          DosBySpinEnergy("p", np.array([h_p_up, h_p_down]))]]

dos_plot_data = DosPlotData(relative_energies=relative_energies,
                            doses=doses,
                            names=["total", "H"],
                            energy_range=xlim,
                            dos_ranges=ylim_set,
                            energy_lines=[0.0, 1.0])
dos_data_len = len(dos_plot_data.doses)


def test_dos_plot_data_msonable():
    assert_msonable(dos_plot_data)


@pytest.fixture
def mock_plt_list(mocker):
    mock_plt = mocker.patch("vise.analyzer.plot_dos.plt", auto_spec=True)
    mock_1st_ax = mocker.MagicMock()
    mock_2nd_ax = mocker.MagicMock()
    other_axs_len = dos_data_len - 2

    mock_axs = [mock_1st_ax, mock_2nd_ax] + [mocker.MagicMock()] * other_axs_len

    mock_plt.subplots.return_value = (None, mock_axs)
    plotter = DosPlotter(dos_plot_data=dos_plot_data)
    plotter.construct_plot()

    return mock_plt, mock_1st_ax, mock_2nd_ax


def test_dos_mpl_settings_defaults():
    mpl_settings = DosMplSettings()
    assert mpl_settings.colors == ['#36454f', '#E15759', '#4E79A7', '#F28E2B']
    assert mpl_settings.linewidth == 1.0
    assert mpl_settings.title_font_size == 15
    assert mpl_settings.label_font_size == 12
    assert mpl_settings.vline == {"linewidth": 0.75, "color": "black",
                                  "linestyle": "-."}


def test_dos_mpl_settings_modify():
    mpl_settings = DosMplSettings(colors=["black"],
                                  linewidth=2.0,
                                  band_edge_line_width=1.0,
                                  band_edge_line_color="blue",
                                  band_edge_line_style=":",
                                  title_font_size=30,
                                  label_font_size=40)

    assert mpl_settings.colors == ["black"]
    assert mpl_settings.linewidth == 2.0
    assert mpl_settings.title_font_size == 30
    assert mpl_settings.label_font_size == 40
    assert mpl_settings.vline == \
           {"linewidth": 1.0, "color": "blue", "linestyle": ":"}


def test_plot_dos(mock_plt_list):
    mock_plt, mock_1st_ax, mock_2nd_ax = mock_plt_list
    mock_plt.subplots.assert_called_once_with(dos_data_len, 1, sharex=True,
                                              gridspec_kw={'hspace': 0.1})

    reversed_total_down = [dos * -1 for dos in total_down]
    reversed_h_p_down = [dos * -1 for dos in h_p_down]

    args = DosMplSettings().dos_line(0)
    args["label"] = "total"
    mock_1st_ax.plot.assert_any_call(relative_energies, total_up, **args)

    args = DosMplSettings().dos_line(0)
    mock_1st_ax.plot.assert_any_call(
        relative_energies, reversed_total_down, **args)

    mock_1st_ax.axhline.assert_called_once_with(0, linestyle=":", color="black")

    args = DosMplSettings().dos_line(0)
    args["label"] = "H-s"
    mock_2nd_ax.plot.assert_any_call(relative_energies, h_s_up, **args)

    args = DosMplSettings().dos_line(1)
    mock_2nd_ax.plot.assert_any_call(
        relative_energies, reversed_h_p_down, **args)


def test_axs_is_list_when_single_dos_passed():
    single_dos = [[DosBySpinEnergy("", [total_up, total_down])]]
    dos_info = DosPlotData(relative_energies=relative_energies,
                           doses=single_dos, names=["total"], energy_range=xlim,
                           dos_ranges=ylim_set, energy_lines=[0.0, 1.0])
    plotter = DosPlotter(dos_plot_data=dos_info)
    assert isinstance(plotter._axs, list)
    assert len(plotter._axs) == 1


def test_set_figure_legend(mock_plt_list):
    _, mock_1st_ax, mock_2nd_ax = mock_plt_list
    mock_1st_ax.legend.assert_called_with(bbox_to_anchor=(0.9, 1),
                                          loc='upper left',
                                          borderaxespad=0,
                                          markerscale=0.1)
    mock_2nd_ax.legend.assert_called_with(bbox_to_anchor=(0.9, 1),
                                          loc='upper left',
                                          borderaxespad=0,
                                          markerscale=0.1)


def test_set_x_and_y_range(mock_plt_list):
    _, mock_1st_ax, _ = mock_plt_list
    mock_1st_ax.set_ylim.assert_called_with([-10, 10])
    mock_1st_ax.set_xlim.assert_called_with([-6, 6])


def test_set_vlines(mock_plt_list):
    _, mock_1st_ax, _ = mock_plt_list
    mpl_settings = DosMplSettings()
    mock_1st_ax.axvline(x=0.0, **mpl_settings.vline)


def test_set_labels(mock_plt_list):
    mock_plt, mock_1st_ax, mock_2nd_ax = mock_plt_list
    defaults = DosMplSettings()
    mock_plt.xlabel.assert_called_with("Energy (eV)",
                                       size=defaults.label_font_size)
    mock_1st_ax.set_ylabel.assert_called_with("Dos (1/eV)",
                                              size=defaults.label_font_size)


def test_set_float_to_int_formatter(mock_plt_list):
    _, mock_1st_ax, mock_2nd_ax = mock_plt_list
    mock_1st_ax.xaxis.set_major_formatter.assert_called_once_with(
        float_to_int_formatter)
    mock_1st_ax.yaxis.set_major_formatter.assert_called_once_with(
        float_to_int_formatter)
    mock_2nd_ax.yaxis.set_major_formatter.assert_called_once_with(
        float_to_int_formatter)


def test_actual_plot():
    plotter = DosPlotter(dos_plot_data)
    plotter.construct_plot()
    plotter.plt.show()


# def test_plotly_actual_plot():
#     pploter = PlotlyDosPlotter(dos_plot_data)
#     pploter.show()


# # 1. imports of your dash app
# import dash
# import dash_html_components as html


# # 2. give each testcase a tcid, and pass the fixture
# # as a function argument, less boilerplate
# def test_bsly001_falsy_child(dash_duo):

    # # 3. define your app inside the test function
    # app = dash.Dash(__name__)
    # app.layout = html.Div(id="nully-wrapper", children=0)

    # # 4. host the app locally in a thread, all dash server configs could be
    # # passed after the first app argument
    # dash_duo.start_server(app)

    # # 5. use wait_for_* if your target element is the result of a callback,
    # # keep in mind even the initial rendering can trigger callbacks
    # dash_duo.wait_for_text_to_equal("#nully-wrapper", "0", timeout=4)

    # # 6. use this form if its present is expected at the action point
    # assert dash_duo.find_element("#nully-wrapper").text == "0"

    # # 7. to make the checkpoint more readable, you can describe the
    # # acceptance criterion as an assert message after the comma.
    # assert dash_duo.get_logs() == [], "browser console should contain no error"

    # # 8. visual testing with percy snapshot
    # dash_duo.percy_snapshot("bsly001-layout")