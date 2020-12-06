# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
from crystal_toolkit.core.mpcomponent import MPComponent
from crystal_toolkit.helpers.layouts import Columns, Column
from dash_core_components import Loading
from plotly.subplots import make_subplots
from vise.analyzer.dos_data import DosPlotData
from vise.analyzer.plot_band import BandPlotInfo


class BandDosComponent(MPComponent):
    def __init__(self,
                 dos_plot_data: DosPlotData,
                 band_plot_info: BandPlotInfo,
                 band_plot_info_2: BandPlotInfo = None,
                 id: str = None,
                 *args, **kwargs):
        super().__init__(id=id, *args, **kwargs)
        self.dos_plot_data = dos_plot_data
        self.band_plot_info = band_plot_info
        self.band_plot_info_2 = band_plot_info_2
#        self.create_store("dos_data", initial_data=dos_data)

    @property
    def _sub_layouts(self):

        # pdos_range = self.get_slider_input(
        #     kwarg_label="pdos_range",
        #     state={},
        #     label="PDOS max range",
        #     min=0.1,
        #     max=
        # )

        if self.dos_plot_data:
            num_dos = len(self.dos_plot_data.doses)
            names = ["band"] + self.dos_plot_data.names
            column_width = [0.25] + [0.75 / num_dos] * num_dos
        else:
            num_dos = 1
            names = ["band"]
            column_width = [0.4, 0.6]

        fig = make_subplots(rows=1,
                            cols=1 + num_dos,
                            column_widths=column_width,
                            shared_yaxes=True,
                            horizontal_spacing=0.01,
                            subplot_titles=names,
                            font_size=25,
                            x_title="DOS (1/eV/unit cell)",
                            y_title="Energy (eV)")

        if self.dos_plot_data:
            fig["layout"]["yaxis"]["range"] = self.dos_plot_data.energy_range
        else:
            fig["layout"]["yaxis"]["range"] = [-5, 10]

        if self.band_plot_info:
            self.add_band(fig)
        else:
            fig.add_trace(
                go.Scatter(x=[], y=[], name="band"), row=1, col=1)

        fig.update_layout(width=2000, height=700)

        if self.dos_plot_data:
            self.add_dos(fig)

        # Main plot
        graph = Loading([dcc.Graph(figure=fig)])

        return {"graph": graph}

    def add_band(self, fig):
        last_path = self.band_plot_info.x_ticks.distances[-1]
        fig["layout"]["xaxis"]["range"] = [0.0, last_path]

        band_edge = self.band_plot_info.band_info_set[0].band_edge
        h_lines_2 = None

        if band_edge:
            base_energy = band_edge.vbm
            cbm = band_edge.cbm - base_energy
            h_lines = [0.0, cbm]
            # Add vbm and cbm.
            fig.add_trace(
                go.Scatter(x=band_edge.vbm_distances,
                           y=[0.0] * len(band_edge.vbm_distances),
                           line_color="blue",
                           fillcolor="blue",
                           mode="markers",
                           marker_size=20,
                           opacity=0.7,
                           showlegend=False),
                row=1, col=1),
            fig.add_trace(
                go.Scatter(x=band_edge.cbm_distances,
                           y=[cbm] * len(band_edge.cbm_distances),
                           line_color="red",
                           fillcolor="red",
                           mode="markers",
                           marker_size=20,
                           opacity=0.7,
                           showlegend=False),
                row=1, col=1),

            if self.band_plot_info_2:
                band_edge_2 = self.band_plot_info_2.band_info_set[0].band_edge
                vbm = band_edge_2.vbm - base_energy
                cbm = band_edge_2.cbm - base_energy
                h_lines_2 = [vbm, cbm]
                # Add vbm and cbm.
                fig.add_trace(
                    go.Scatter(x=band_edge_2.vbm_distances,
                               y=[vbm] * len(band_edge_2.vbm_distances),
                               line_color="blue",
                               fillcolor="blue",
                               mode="markers",
                               marker_size=20,
                               opacity=0.5,
                               showlegend=False),
                    row=1, col=1),
                fig.add_trace(
                    go.Scatter(x=band_edge_2.cbm_distances,
                               y=[cbm] * len(band_edge_2.cbm_distances),
                               line_color="red",
                               fillcolor="red",
                               mode="markers",
                               marker_size=20,
                               opacity=0.5,
                               showlegend=False),
                    row=1, col=1),

        else:
            base_energy = self.band_plot_info.band_info_set[0].fermi_level
            h_lines = [0.0]

        # Add vbm and cbm lines.
        for energy in h_lines:
            fig.add_shape(dict(type="line",
                               x0=0, x1=last_path, y0=energy, y1=energy,
                               line=dict(color="royalblue", width=3, dash="dash")),
                          row=1, col=1)

        if h_lines_2:
            for energy in h_lines:
                fig.add_shape(dict(type="line",
                                   x0=0, x1=last_path, y0=energy, y1=energy,
                                   line=dict(color="royalblue", width=3, dash="dot")),
                              row=1, col=1)


        # Add x labels
        str_replace = {
            "$\mid$": "|",
            "_1": "₁",
            "_2": "₂",
            "_3": "₃",
            "_4": "₄",
            r"${\rm ": "",
            "}": "",
            "$": "",
        }
        labels = self.band_plot_info.x_ticks.labels
        distances = self.band_plot_info.x_ticks.distances
        new_labels = []
        for label in labels:
            for key in str_replace.keys():
                if key in label:
                    label = label.replace(key, str_replace[key])
            new_labels.append(label)

        fig.update_xaxes(tickvals=self.band_plot_info.x_ticks.distances,
                         ticktext=new_labels,
                         tickfont_size=16,
                         row=1, col=1)

        # Add k-path border lines.
        for label, distance in zip(new_labels, distances):
            if "|" in label:
                fig.add_shape(dict(type="line",
                                   x0=distance, x1=distance, y0=-20, y1=20,
                                   line=dict(color="black",
                                             width=2)), row=1, col=1)

        for d, e in zip(self.band_plot_info.distances_by_branch,
                        self.band_plot_info.band_info_set[0].band_energies):
            for i in e[0]:
                fig.add_trace(
                    go.Scatter(x=d,
                               y=[j - base_energy for j in i],
                               hoverinfo="skip",
                               line_color="black",
                               showlegend=False,
                               name="band",
                               mode="lines",
                               line={"width": 2.0}), row=1, col=1),

        if self.band_plot_info_2:
            for d, e in zip(self.band_plot_info.distances_by_branch,
                            self.band_plot_info_2.band_info_set[0].band_energies):
                for i in e[0]:
                    fig.add_trace(
                        go.Scatter(x=d,
                                   y=[j - base_energy for j in i],
                                   hoverinfo="skip",
                                   line_color="green",
                                   showlegend=False,
                                   name="band",
                                   mode="lines",
                                   opacity=0.5,
                                   line={"width": 3.0}), row=1, col=1),

    def add_dos(self, fig):
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
                                               self.dos_plot_data.dos_ranges), 2):
            for dos, color in zip(doses, colors):
                fig.add_trace(
                    go.Scatter(y=self.dos_plot_data.relative_energies,
                               x=dos.dos[0],
                               name=dos.name,
                               line_color=color,
                               line={"width": 3},
                               legendgroup="total" if i == 2 else "pdos",
                               showlegend=True if i == 3 else False),
                    row=1, col=i
                )
            for energy in [0.0, cbm]:
                fig.add_shape(dict(type="line",
                                   x0=y_lim[0], x1=y_lim[1],
                                   y0=energy, y1=energy,
                                   line=dict(color="royalblue",
                                             width=3,
                                             dash="dash")), row=1, col=i)
            fig.update_xaxes(range=y_lim, row=1, col=i)
        fig.update_layout(legend=dict(yanchor="top", y=0.99,
                                      xanchor="right", x=0.99))

    def layout(self):
        return html.Div([Columns(Column(self._sub_layouts["graph"]))])

    def generate_callbacks(self, app, cache):
        pass
    #     @app.callback(
    #         Output(self.id("dos-plot"), "figure"),
    #         [
    #             Input(self.get_kwarg_id("grouping"), "value"),
    #         ],
    #     )
    #     def update_graph(grouping):

            # plot_data: DosPlotData = self.dos_plot_data.dos_plot_data({"O": [10, 11]})

            # fig = make_subplots(rows=1, cols=2, shared_yaxes=True,
            #                     horizontal_spacing=0.01)

            # fig.add_trace(
            #     go.Scatter(y=plot_data.relative_energies,
            #                x=plot_data.doses[0][0].dos[0],
            #                name="total", line={"width": 4}),
            #     row=1, col=1
            # )

            # fig.add_trace(
            #     go.Scatter(y=plot_data.relative_energies,
            #                x=plot_data.doses[1][0].dos[0],
            #                name="s", line={"width": 4}),
            #     row=1, col=2
            # )

            # fig.add_trace(
            #     go.Scatter(y=plot_data.relative_energies,
            #                x=plot_data.doses[1][1].dos[0],
            #                name="p", line={"width": 4}),
            #     row=1, col=2
            # )

            # fig.add_trace(
            #     go.Scatter(y=plot_data.relative_energies,
            #                x=plot_data.doses[1][2].dos[0],
            #                name="d", line={"width": 4}),
            #     row=1, col=2
            # )

            # return fig

