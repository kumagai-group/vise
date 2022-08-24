# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from math import log10
from typing import List

import numpy as np
from matplotlib import pyplot as plt
from plotly import graph_objects as go
from vise.analyzer.dielectric_function import DieleFuncData
from vise.analyzer.dielectric_function_data.exp_dielectric_func import \
    ExpDieleFunc
from vise.util.enum import ExtendedEnum
from vise.util.matplotlib import float_to_int_formatter


class TensorDirection(ExtendedEnum):
    average = "average"
    xx = "xx"
    yy = "yy"
    zz = "zz"
    xy = "xy"
    yz = "yz"
    xz = "xz"

    def val(self, tensors: List[List[float]]) -> List[float]:
        if self is self.average:
            return [(t[0] + t[1] + t[2]) / 3 for t in tensors]
        else:
            idx = ["xx", "yy", "zz", "xy", "yz", "xz"].index(self.name)
            return [tensor[idx] for tensor in tensors]


class DieleFuncPlotType(ExtendedEnum):
    diele_func = "Dielectric function"
    absorption_coeff = "Absorption coefficient"
    refraction = "Refraction"
    reflectivity = "Reflectivity"

    def y_axis_label(self, plot_engine):
        if self is self.absorption_coeff:
            if plot_engine == "plotly":
                sup_before, sup_after = "<sup>", "</sup>"
            elif plot_engine == "matplotlib":
                sup_before, sup_after = "$^{", "}$"
            else:
                raise ValueError("The plotter engine is not adequate.")
            return f"{self} (cm{sup_before}-1{sup_after})"

        return str(self)

    @property
    def y_axis_default_min(self):
        if self is self.absorption_coeff:
            return 10**3, 10**7
        return None

    def tensors(self, data: DieleFuncData):
        if self is self.diele_func:
            return data.diele_func_real, data.diele_func_imag
        elif self is self.absorption_coeff:
            return data.absorption_coeff,
        elif self is self.refraction:
            return data.refractive_idx_real, data.refractive_idx_imag
        elif self is self.reflectivity:
            return data.reflectivity,


def auto_y_range(plot_type, values):
    _max = max(values)
    if plot_type is DieleFuncPlotType.absorption_coeff:
        y_min, y_max = 10 ** 3, _max * 2
    else:
        _min = min(values)
        y_min, y_max = 1.05 * _min - .05 * _max, 1.05 * _max - .05 * _min
    return y_min, y_max


class DieleFuncPlotter:

    def __init__(self, diele_func_data: DieleFuncData,
                 energy_range: List[float] = None):

        self.diele_func_data = diele_func_data
        self.energies = diele_func_data.energies
        self.band_gap = diele_func_data.band_gap
        self.energy_range = energy_range or [0, 10]
        self._x_axis_title = "Energy (eV)"

class DieleFuncPlotlyPlotter(DieleFuncPlotter):

    def create_figure(self,
                      directions=(TensorDirection.average,),
                      plot_type=DieleFuncPlotType.absorption_coeff,
                      y_range: List[float] = None):
        fig = go.Figure()
        y_axis = plot_type.y_axis_label("plotly")
        fig.update_layout(xaxis_title=self._x_axis_title,
                          yaxis_title=y_axis,
                          font_size=25, width=800, height=700)

        _idx = np.where(np.array(self.energies) > self.energy_range[1])[0][0]
        values = []
        a = plot_type.tensors(self.diele_func_data)
        for b, tensors in zip(["real", "imag"], a):
            for direction in directions:
                tensor = direction.val(tensors)
                values.extend(tensor[:_idx])
                if len(a) == 1:
                    name = direction.name
                else:
                    name = f"{b}_{direction.name}"
                fig.add_trace(go.Scatter(x=self.energies,
                                         y=tensor,
                                         line=dict(width=2.5),
                                         name=name))

        y_min, y_max = y_range or auto_y_range(plot_type, values)

        if plot_type is DieleFuncPlotType.absorption_coeff:
            kwargs = {"type": "log", "range": [log10(y_min), log10(y_max)]}
        else:
            kwargs = {"range": [y_min, y_max]}

        if self.band_gap:
            fig.add_trace(go.Scatter(x=[self.band_gap, self.band_gap],
                                     y=[y_min, y_max],
                                     line=dict(width=2, dash="dash"),
                                     line_color="black",
                                     name="Band gap"))

        fig.update_xaxes(range=self.energy_range, tickfont_size=20)
        fig.update_yaxes(tickfont_size=20,
                         showexponent="all", exponentformat='power',
                         **kwargs)
        return fig


class DieleFuncMplPlotter(DieleFuncPlotter):
    def construct_plot(self,
                       directions=(TensorDirection.average,),
                       plot_type=DieleFuncPlotType.absorption_coeff,
                       y_range=None):
        self.plt = plt
        self._add_coeffs(directions, plot_type, y_range)
        self._add_band_gap()
        # self._set_figure_legend()
        self._set_x_range()
        self._set_labels()
        self._set_formatter()
        self.plt.tight_layout()

    def _add_coeffs(self, directions, plot_type, y_range):
        _idx = np.where(np.array(self.energies) > self.energy_range[1])[0][0]
        values = []
        a = plot_type.tensors(self.diele_func_data)
        for b, tensors in zip(["real", "imag"], a):
            for direction in directions:
                tensor = direction.val(tensors)
                values.extend(tensor[:_idx])
                if len(a) == 1:
                    name = direction.name
                else:
                    name = f"{b}_{direction.name}"
                self.plt.semilogy(self.energies, tensor, label=name)

        y_min, y_max = y_range or auto_y_range(plot_type, values)
        self.plt.gca().set_ylim(ymin=y_min, ymax=y_max)

    def _add_band_gap(self):
        self.plt.axvline(x=self.band_gap, linestyle="dashed", color="black",
                         linewidth=1)

    def _set_figure_legend(self):
        self.plt.legend(loc="lower right")

    def _set_x_range(self):
        self.plt.xlim(self.energy_range[0], self.energy_range[1])

    def _set_labels(self):
        self.plt.xlabel(self._x_axis_title)
        y_axis = DieleFuncPlotType.absorption_coeff.y_axis_label("matplotlib")
        self.plt.ylabel(y_axis)

    def _set_formatter(self):
        axis = self.plt.gca()
        axis.xaxis.set_major_formatter(float_to_int_formatter)