# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from dataclasses import dataclass
from typing import List, Dict, Optional

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class NamedSingleDos:
    name: str
    dos: List[List[float]]

    def max_dos(self):
        return np.max(np.array(self.dos))


class DosPlotInfo:

    def __init__(self,
                 energies: List[float],
                 doses: List[List[NamedSingleDos]],
                 xlim: Optional[List[float]] = None,
                 ylim_set: Optional[List[List[float]]] = None):
        self.energies = energies
        self.doses = doses  # [ax][orbital]
        self.xlim = xlim or [-10, 10]
        self.ylim_set = ylim_set or [[-y, y] for y in self.max_y_ranges()]

    def max_y_ranges(self, multi=1.1, round_digit=2):
        return [round(i * multi, round_digit) for i in self.max_modulus_dos()]

    def max_modulus_dos(self):
        result = []
        for dos in self.doses:
            result.append(np.max([orb_dos.max_dos() for orb_dos in dos]))
        return result


class DosMplDefaults:
    def __init__(self,
                 colors: Optional[List[str]] = None,
                 linewidth: float = 1.0,
                 band_edge_line_width: float = 0.75,
                 band_edge_line_color: str = "black",
                 band_edge_line_style: str = "-.",
                 title_font_size: int = 15,
                 label_font_size: int = 15,
                 ):
        self.colors = colors or ['#36454f', '#E15759', '#4E79A7', '#F28E2B',
                                 '#76B7B2']
        self.linewidth = linewidth

        self.vline = {"linewidth": band_edge_line_width,
                      "color": band_edge_line_color,
                      "linestyle": band_edge_line_style}

        self.title_font_size = title_font_size
        self.label_font_size = label_font_size

    def dos_line(self, index):
        return {"color": self.colors[index],
                "linewidth": self.linewidth}


class DosPlotter:
    def __init__(self,
                 dos_info: DosPlotInfo,
                 mpl_defaults: Optional[DosMplDefaults] = DosMplDefaults()):
        self._dos_info = dos_info
        self.mpl_defaults = mpl_defaults

        self.plt = plt
        num_axs = len(self._dos_info.doses)
        _, self._axs = self.plt.subplots(num_axs, 1, sharex=True)
        if num_axs == 1:
            self._axs = [self._axs]

    def construct_plot(self):
        self._axs[0].set_xlim(self._dos_info.xlim)
        for i in range(len(self._dos_info.doses)):
            self._add_ax(i)

    def _add_ax(self, i):
        self._add_dos(i)
        self._set_y_range(i)
        self._set_dos_zero_line(i)

    def _add_dos(self, i):
        for j, dos_each_name in enumerate(self._dos_info.doses[i]):
            for k, dos_each_spin in enumerate(dos_each_name.dos):
                sign = - (2 * k) + 1
                dos_for_plot = [d * sign for d in dos_each_spin]
                args = self.mpl_defaults.dos_line(j)
                if k == 0:
                    args["label"] = self._dos_info.doses[i][j].name
                self._axs[i].plot(self._dos_info.energies, dos_for_plot, **args)

        self._axs[i].legend(loc="best", markerscale=0.1)

    def _set_dos_zero_line(self, i):
        self._axs[i].axhline(0, linestyle=":", color="black")

    def _set_y_range(self, i):
        self._axs[i].set_ylim(self._dos_info.ylim_set[i])

