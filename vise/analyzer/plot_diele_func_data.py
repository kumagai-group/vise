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
    diele_real = "diele_real"
    diele_imag = "diele_imag"
    absorption_coeff = "absorption_coeff"
    refraction_real = "refraction_real"
    refraction_imag = "refraction_imag"
    reflection = "reflection"


class DieleFuncPlotter:

    def __init__(self,
                 diele_func_data: DieleFuncData,
                 energy_range: List[float] = None,
                 yranges: List[float] = None):

        self.diele_func_data = diele_func_data
        self.energies = diele_func_data.energies
        self.band_gap = diele_func_data.band_gap
        self.energy_range = energy_range or [0, 10]

        self._xaxis_title = "Energy (eV)"

        if yranges:
            self.ymin, self.ymax = yranges
        else:
            self.ymin, self.ymax = 10**3, 10**7

        self.plt = plt


class DieleFuncPlotlyPlotter(DieleFuncPlotter):

    _yaxis_title = "Absorption coefficient (cm <sup>-1</sup>)"

    def create_figure(self,
                      directions=(TensorDirection.average,),
                      materials: List[str] = None):
        fig = go.Figure()
        fig.update_layout(xaxis_title=self._xaxis_title,
                          yaxis_title=self._yaxis_title,
                          font_size=25, width=800, height=700)

        for direction in directions:
            abs_coeff = direction.val(self.diele_func_data.absorption_coeff)

            fig.add_trace(go.Scatter(x=self.energies,
                                     y=abs_coeff,
                                     line=dict(width=2.5),
                                     name=direction.name))
        if self.band_gap:
            fig.add_trace(go.Scatter(x=[self.band_gap, self.band_gap],
                                     y=[self.ymin, self.ymax],
                                     line=dict(width=2, dash="dash"),
                                     line_color="black",
                                     name="Band gap"))

        if materials:
            for material in materials:
                exp = ExpDieleFunc(material)
                energies = exp.energies
                fig.add_trace(go.Scatter(x=energies, y=exp.absorption_coeff,
                                         line=dict(width=1, dash="dashdot"),
                                         name=material))
                fig.add_trace(go.Scatter(x=[exp.band_gap, exp.band_gap],
                                         y=[self.ymin, self.ymax],
                                         line=dict(width=1, dash="dash"),
                                         showlegend=False,
                                         line_color="black",
                                         name=f"{material} band gap"))

        fig.update_xaxes(range=self.energy_range, tickfont_size=20)
        fig.update_yaxes(type="log",
                         range=[log10(self.ymin), log10(self.ymax)],
                         tickfont_size=20,
                         showexponent="all", exponentformat='power'
                         )
        return fig


class DieleFuncMplPlotter(DieleFuncPlotter):
    _yaxis_title = "Absorption coefficient (cm$^{-1}$)"

    def construct_plot(self, directions=(TensorDirection.average,),
                       materials: List[str] = None):
        self._add_coeffs(directions)
        self._add_band_gap()
        if materials:
            self._add_materials(materials)
        # self._set_figure_legend()
        self._set_x_range()
        self._set_y_range()
        self._set_labels()
        self._set_formatter()
        self.plt.tight_layout()

    def _add_coeffs(self, directions):
        for direction in directions:
            abs_coeff = direction.val(self.diele_func_data.absorption_coeff)
            self.plt.semilogy(self.energies, abs_coeff)

    def _add_band_gap(self):
        self.plt.axvline(x=self.band_gap, linestyle="dashed", color="black",
                         linewidth=1)

    def _add_materials(self, materials):
        for material in materials:
            exp = ExpDieleFunc(material)
            energies = exp.energies
            self.plt.axvline(x=exp.band_gap,
                             linestyle="dashed",
                             color="black",
                             linewidth=1)

            self.plt.semilogy(energies,
                              list(exp.absorption_coeff),
                              label=material,
                              linestyle="dashdot",
                              )

    def _set_figure_legend(self):
        self.plt.legend(loc="lower right")

    def _set_x_range(self):
        self.plt.xlim(0, 10)

    def _set_y_range(self):
        self.plt.gca().set_ylim(ymin=self.ymin, ymax=self.ymax)

    def _set_labels(self):
        self.plt.xlabel(self._xaxis_title)
        self.plt.ylabel(self._yaxis_title)

    def _set_formatter(self):
        axis = self.plt.gca()
        axis.xaxis.set_major_formatter(float_to_int_formatter)