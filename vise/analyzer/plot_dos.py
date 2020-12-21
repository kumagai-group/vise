# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import List, Optional

import matplotlib.pyplot as plt

from vise.analyzer.dos_data import DosPlotData
from vise.util.matplotlib import float_to_int_formatter


class DosMplSettings:
    def __init__(self,
                 colors: Optional[List[str]] = None,
                 linewidth: float = 1.0,
                 band_edge_line_width: float = 0.75,
                 band_edge_line_color: str = "black",
                 band_edge_line_style: str = "-.",
                 title_font_size: int = 15,
                 label_font_size: int = 12,
                 ):
        self.colors = colors or ['#36454f', '#E15759', '#4E79A7', '#F28E2B']
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
                 dos_plot_data: DosPlotData,
                 show_legend: bool = True,
                 mpl_defaults: Optional[DosMplSettings] = DosMplSettings()):

        self._dos_plot_data = dos_plot_data
        self._show_legend = show_legend
        self.mpl_defaults = mpl_defaults
        self.plt = plt

        num_axs = len(self._dos_plot_data.doses)
        fig, self._axs = self.plt.subplots(num_axs, 1,
                                           sharex=True,
                                           gridspec_kw={'hspace': 0.1})
        if num_axs == 1:
            self._axs = [self._axs]

    def construct_plot(self):
        self._axs[0].set_xlim(self._dos_plot_data.energy_range)
        self._axs[0].xaxis.set_major_formatter(float_to_int_formatter)
        for i in range(len(self._dos_plot_data.doses)):
            self._add_ax(i)
        self._set_x_labels()
        middle_idx = int(i / 2)
        self._set_ylabel(middle_idx)

    def _add_ax(self, i):
        self._add_dos(i)
        self._set_y_range(i)
        self._set_dos_zero_line(i)
        self._set_vline(i)
        self._set_formatter(i)

    def _add_dos(self, i):
        for j, by_name_dos in enumerate(self._dos_plot_data.doses[i]):
            for k, by_spin_dos in enumerate(by_name_dos.dos):
                sign = - (2 * k) + 1  # up: 1, down: -1
                dos_for_plot = [d * sign for d in by_spin_dos]
                args = self.mpl_defaults.dos_line(j)
                if k == 0 and self._show_legend:
                    if self._dos_plot_data.doses[i][j].name:
                        args["label"] = "-".join(
                            [self._dos_plot_data.names[i],
                             self._dos_plot_data.doses[i][j].name])
                    else:
                        args["label"] = self._dos_plot_data.names[i]
                self._axs[i].plot(self._dos_plot_data.relative_energies,
                                  dos_for_plot, **args)

        if self._show_legend:
            self._axs[i].legend(bbox_to_anchor=(0.9, 1), loc='upper left',
                                borderaxespad=0, markerscale=0.1)

    def _set_ylabel(self, i):
        self._axs[i].set_ylabel("Dos (1/eV)",
                                size=self.mpl_defaults.label_font_size)

    def _set_dos_zero_line(self, i):
        self._axs[i].axhline(0, linestyle=":", color="black")

    def _set_y_range(self, i):
        self._axs[i].set_ylim(self._dos_plot_data.dos_ranges[i])

    def _set_vline(self, i):
        for line_energy in self._dos_plot_data.energy_lines:
            self._axs[i].axvline(x=line_energy, **self.mpl_defaults.vline)

    def _set_formatter(self, i):
        self._axs[i].yaxis.set_major_formatter(float_to_int_formatter)

    def _set_x_labels(self):
        self.plt.xlabel("Energy (eV)", size=self.mpl_defaults.label_font_size)

