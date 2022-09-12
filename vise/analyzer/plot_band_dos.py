# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import re
from typing import Optional, Tuple

from plotly.subplots import make_subplots
from vise.analyzer.dos_data import DosPlotData
from vise.analyzer.plot_band import BandPlotInfo, BandEdgeForPlot, \
    BandEnergyInfo, get_base_energy, slide_band_energies
from plotly import graph_objects as go


def plotly_sanitize_label(label: str) -> str:
    str_replace = {"$\mid$": "|", r"${\rm ": "", "}": "", "$": ""}
    result = re.sub(r'(\w*)_(\w*)', '\\1<sub>\\2</sub>', label)
    for key in str_replace.keys():
        if key in result:
            result = result.replace(key, str_replace[key])
    return result


class BandDosPlotlySettings:
    def __init__(self,
                 band_colors: Optional[Tuple[str, str]] = ("black", "green"),
                 band_opacity: Optional[Tuple[float, float]] = (0.7, 0.2),
                 band_line_widths: Optional[Tuple[float, float]] = (2.0, 3.0),
                 layout_width: Optional[int] = 1600,
                 layout_height: Optional[int] = 700,
                 yaxis_range: Optional[Tuple[float, float]] = (-5.0, 10.0),

                 vbm_color: Optional[str] = "blue",
                 cbm_color: Optional[str] = "red",

                 h_line_colors: Optional[Tuple[str, str]] =
                 ("royalblue", "green"),
                 h_line_width: Optional[int] = 3,

                 subplot_font_size: Optional[int] = 24,
                 tickfont_size: Optional[int] = 20,

                 edge_circle_mode: Optional[str] = "markers",
                 edge_circle_size: Optional[int] = 20
                 ):
        self.band_colors = band_colors
        self.band_opacity = band_opacity
        self.band_line_widths = band_line_widths
        self.layout_width = layout_width
        self.layout_height = layout_height
        self.yaxis_range = yaxis_range
        self.vbm_color = vbm_color
        self.cbm_color = cbm_color
        self.h_line_colors = h_line_colors
        self.h_line_width = h_line_width
        self.subplot_font_size = subplot_font_size
        self.tickfont_size = tickfont_size
        self.band_edge_circles = {"mode": edge_circle_mode,
                                  "marker_size": edge_circle_size,
                                  "showlegend": False}


class BandDosPlotlyPlotter:

    def __init__(self,
                 dos_plot_data: DosPlotData = None,
                 band_plot_info: BandPlotInfo = None,
                 plotly_defaults: Optional[BandDosPlotlySettings]
                 = BandDosPlotlySettings()):

        self.dos_plot_data = dos_plot_data
        self.band_plot_info = band_plot_info
        if band_plot_info:
            base_energy = get_base_energy(self.band_plot_info.band_energy_infos)
            slide_band_energies(
                self.band_plot_info.band_energy_infos, base_energy)

        self.plotly_defaults = plotly_defaults

        self._create_fig_w_subplots(dos_plot_data)
        self._set_subplot_font_size()
        self._set_y_axis()
        self._add_band()
        self._add_dos()
        self.fig.update_layout(width=plotly_defaults.layout_width,
                               height=plotly_defaults.layout_height)

    def _set_subplot_font_size(self):
        for i in self.fig['layout']['annotations']:
            i['font'] = dict(size=self.plotly_defaults.subplot_font_size)

    def _set_y_axis(self):
        self.fig["layout"]["yaxis"]["range"] = self.plotly_defaults.yaxis_range
        self.fig.update_yaxes(tickfont_size=self.plotly_defaults.tickfont_size)

    def _create_fig_w_subplots(self, dos_plot_data):
        if dos_plot_data:
            num_dos_panels = len(dos_plot_data.doses)
            names = [""] + [plotly_sanitize_label(name)
                            for name in dos_plot_data.names]
            column_width = [0.25] + [0.75 / num_dos_panels] * num_dos_panels
        else:
            num_dos_panels = 1
            names = [""]
            column_width = [0.4, 0.6]

        self.fig = make_subplots(rows=1,
                                 cols=1 + num_dos_panels,
                                 column_widths=column_width,
                                 shared_yaxes=True,
                                 horizontal_spacing=0.005,
                                 subplot_titles=names,
                                 y_title="Energy (eV)")
        if self.dos_plot_data:
            self.fig.update_xaxes(title_text="DOS (1/eV/unit cell)",
                                  row=1, col=2)

    def _add_band(self):
        if self.band_plot_info is None:
            # add empty trace when band_plot_info is None.
            self.fig.add_trace(go.Scatter(x=[], y=[], name="band"),
                               row=1, col=1)
            return

        self._last_kpt_distance = self.band_plot_info.x_ticks.distances[-1]
        self.fig["layout"]["xaxis"]["range"] = [0.0, self._last_kpt_distance]

        for (subtitle, band_info), width, color, opacity, h_line_color in \
                zip(self.band_plot_info.band_energy_infos.items(),
                    self.plotly_defaults.band_line_widths,
                    self.plotly_defaults.band_colors,
                    self.plotly_defaults.band_opacity,
                    self.plotly_defaults.h_line_colors):
            self._add_bands(band_info, width, color, opacity)
            if band_info.band_edge:
                self._add_band_edges(band_info.band_edge, opacity, h_line_color)
            else:
                self._add_horizontal_lines([band_info.fermi_level],
                                           h_line_color, "dash")

        self._add_special_points_and_border_lines()
        self._add_irrep_points()

    def _add_band_edges(self, band_edge: BandEdgeForPlot, opacity,
                        h_line_color):
        self._add_horizontal_lines(
            [band_edge.vbm, band_edge.cbm], h_line_color, "dash")
        self._add_band_edge_circles(band_edge.vbm_distances,
                                    opacity=opacity, energy=band_edge.vbm,
                                    color=self.plotly_defaults.vbm_color)
        self._add_band_edge_circles(band_edge.cbm_distances,
                                    opacity=opacity, energy=band_edge.cbm,
                                    color=self.plotly_defaults.cbm_color)

    def _add_band_edge_circles(self, distances, opacity, energy, color):
        self.fig.add_trace(
            go.Scatter(x=distances, y=[round(energy, 2)] * len(distances),
                       line_color=color, fillcolor=color, opacity=opacity,
                       **self.plotly_defaults.band_edge_circles),
            row=1, col=1),

    def _add_horizontal_lines(self, line_heights, color, dash):
        for height in line_heights:
            self.fig.add_shape(
                dict(type="line", x0=0, x1=self._last_kpt_distance,
                     y0=height, y1=height,
                     line=dict(color=color, dash=dash,
                               width=self.plotly_defaults.h_line_width)),
                row=1, col=1)

    def _add_special_points_and_border_lines(self):
        new_labels = [plotly_sanitize_label(label)
                      for label in self.band_plot_info.x_ticks.labels]
        distances = self.band_plot_info.x_ticks.distances
        self.fig.update_xaxes(tickvals=distances,
                              ticktext=new_labels,
                              tickfont_size=self.plotly_defaults.tickfont_size,
                              row=1, col=1)
        for label, distance in zip(new_labels, distances):
            if "|" in label:
                self.fig.add_shape(dict(type="line",
                                        x0=distance, x1=distance, y0=-20, y1=20,
                                        line=dict(color="black",
                                                  width=2)), row=1, col=1)

    def _add_irrep_points(self):
        for band_e_info in self.band_plot_info.band_energy_infos.values():
            if band_e_info.irreps is None:
                continue
            distances = band_e_info.irreps.get_distances(self.band_plot_info.x_ticks)
            for irrep, d in zip(band_e_info.irreps().values(), distances):
                for symbol, energy, degeneracy in irrep.irrep_info_set:
                    hover = f"{energy:.2f} eV {symbol} " \
                            f"({degeneracy})<extra></extra>"
                    self.fig.add_trace(
                        go.Scatter(x=d, y=[energy] * len(d), mode='markers',
                                   hovertemplate=hover,
                                   showlegend=False,
                                   line_color="pink"),
                        row=1, col=1),

    def _add_bands(self, band_info: BandEnergyInfo, width, color, opacity=1.0):
        try:
            assert len(band_info.band_energies[0]) == 1
        except AssertionError:
            print("Currently, spin-polarization is not allowed.")

        for branch_idx, (distances, eigvals) in enumerate(zip(
                self.band_plot_info.distances_by_branch,
                band_info.band_energies)):
            for eigvals_along_kpath in eigvals[0]:  # only spin-up
                x_wo_irrep, y_wo_irrep = [], []
                for x, y in zip(distances, eigvals_along_kpath):
                    x_wo_irrep.append(x)
                    y_wo_irrep.append(y)

                self.fig.add_trace(
                    go.Scatter(x=x_wo_irrep,
                               y=y_wo_irrep,
                               hoverinfo="skip",
                               line_color=color,
                               showlegend=False,
                               name="band",
                               mode="lines",
                               opacity=opacity,
                               line={"width": width}), row=1, col=1),

    def _add_dos(self):
        if self.dos_plot_data is None:
            return
        total_dos = self.dos_plot_data.doses[0][0].dos[0]
        energies = self.dos_plot_data.relative_energies
        cbm = min([e for d, e in zip(total_dos, energies)
                   if e > 0 and d > 1e-5])

        colors = ["#1f77b4",  # muted blue
                  "#ff7f0e",  # safety orange
                  "#2ca02c",  # cooked asparagus green
                  "#9467bd",  # muted purple
                  ]
        for i, (doses, y_lim) in enumerate(zip(self.dos_plot_data.doses,
                                               self.dos_plot_data.dos_ranges),
                                           2):
            for dos, color in zip(doses, colors):
                self.fig.add_trace(
                    go.Scatter(y=self.dos_plot_data.relative_energies,
                               x=dos.dos[0],
                               name=f"<i>{dos.name}</i>",
                               line_color=color,
                               line={"width": 3},
                               legendgroup="total" if i == 2 else "pdos",
                               showlegend=True if i == 3 else False,
                               ),
                    row=1, col=i
                )
            for energy in [0.0, cbm]:
                self.fig.add_shape(dict(type="line",
                                        x0=y_lim[0], x1=y_lim[1],
                                        y0=energy, y1=energy,
                                        line=dict(color="royalblue",
                                                  width=3,
                                                  dash="dash")), row=1, col=i)
            self.fig.update_xaxes(
                range=y_lim, row=1, col=i,
                tickfont_size=self.plotly_defaults.tickfont_size)
        self.fig.update_layout(legend=dict(yanchor="top", y=0.99,
                                           xanchor="right", x=0.99,
                                           font_size=20))
