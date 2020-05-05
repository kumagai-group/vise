# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import List, Optional

import matplotlib.pyplot as plt
from vise.analyzer.dos_info import DosPlotData


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
                 dos_data: DosPlotData,
                 mpl_defaults: Optional[DosMplDefaults] = DosMplDefaults()):

        self._dos_info = dos_data
        self.mpl_defaults = mpl_defaults

        self.plt = plt
        num_axs = len(self._dos_info.doses)
        fig, self._axs = self.plt.subplots(num_axs, 1, sharex=True)
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
        for j, by_name_dos in enumerate(self._dos_info.doses[i]):
            for k, by_spin_dos in enumerate(by_name_dos.dos):
                sign = - (2 * k) + 1
                dos_for_plot = [d * sign for d in by_spin_dos]
                args = self.mpl_defaults.dos_line(j)
                if k == 0:
                    args["label"] = self._dos_info.doses[i][j].name
                self._axs[i].plot(self._dos_info.relative_energies,
                                  dos_for_plot, **args)

        self._axs[i].legend(loc="best", markerscale=0.1)

    def _set_dos_zero_line(self, i):
        self._axs[i].axhline(0, linestyle=":", color="black")

    def _set_y_range(self, i):
        self._axs[i].set_ylim(self._dos_info.ylim_set[i])

