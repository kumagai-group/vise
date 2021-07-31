# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import itertools
from dataclasses import dataclass
from math import sqrt
from typing import List, Dict

from monty.json import MSONable

import plotly.graph_objects as go
from vise.analyzer.plot_band_dos import plotly_sanitize_label
from vise.util.plotly_util import sort_coords, make_triangles
import numpy as np


@dataclass
class BZPlotInfo(MSONable):
    faces: List[List[List[float]]]
    # {"X": {"cart": [0.5, 0, 0], "frac": [0.7514, 0, 0]"}}
    labels: Dict[str, Dict[str, List[float]]]
    band_paths: List[List[List[float]]] = None
    rec_lat_vec: List[List[float]] = None


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itertools.tee(iterable)
    next(b)
    return itertools.zip_longest(a, b, fillvalue=iterable[0])


class BZPlotlyPlotter:
    def __init__(self, bz_plot_info: BZPlotInfo):
        self._bz_plot_info = bz_plot_info

    def create_figure(self):
        data = []
        for face_vertices in self._bz_plot_info.faces:
            vertex_coords = sort_coords(np.array(face_vertices))
            # data.append(go.Mesh3d(**make_triangles(vertex_coords),
            #                       alphahull=-1, opacity=0.1, hoverinfo='skip',
            #                       color="blue"))
            for s, t in pairwise(vertex_coords):
                data.append(go.Scatter3d(x=[s[0], t[0]],
                                         y=[s[1], t[1]],
                                         z=[s[2], t[2]], mode="lines",
                                         hoverinfo="none",
                                         marker_color="black"))

        for label, coords in self._bz_plot_info.labels.items():
            c_coords, f_coords = coords["cart"], coords["frac"]
            data.append(go.Scatter3d(x=[c_coords[0]], y=[c_coords[1]], z=[c_coords[2]],
                                     text=[plotly_sanitize_label(label)],
                                     mode="markers+text",
                                     marker=dict(color="orange", size=6),
                                     meta=f_coords,
                                     hovertemplate="(%{meta[0]:.2f} %{meta[1]:.2f} %{meta[2]:.2f})",
                                     textposition="middle center",
                                     textfont=dict(size=32)))

        for i, l in self._bz_plot_info.band_paths:
            data.append(go.Scatter3d(x=[i[0], l[0]],
                                     y=[i[1], l[1]],
                                     z=[i[2], l[2]], mode="lines",
                                     opacity=0.5,
                                     hoverinfo="none",
                                     line_width=8,
                                     marker_color="purple"))

        _max = np.amax(np.array(sum(self._bz_plot_info.faces, [])))
        c_max = _max * 1.3
        c_cone = _max * 0.3

        data.append(go.Cone(x=[c_max], y=[0], z=[0], u=[c_cone], v=[0], w=[0],
                            colorscale=[[0, 'rgb(200,200,200)'],
                                        [1, 'rgb(200,200,200)']],
                            showscale=False, hovertemplate="<b>x</b>"))
        data.append(go.Cone(x=[0], y=[c_max], z=[0], u=[0], v=[c_cone], w=[0],
                            colorscale=[[0, 'rgb(200,200,200)'],
                                        [1, 'rgb(200,200,200)']],
                            showscale=False, hovertemplate="<b>y</b>"))
        data.append(go.Cone(x=[0], y=[0], z=[c_max], u=[0], v=[0], w=[c_cone],
                            colorscale=[[0, 'rgb(200,200,200)'],
                                        [1, 'rgb(200,200,200)']],
                            showscale=False, hovertemplate="<b>z</b>"))
        data.append(go.Scatter3d(x=[0, c_max], y=[0, 0], z=[0, 0], mode="lines",
                                 hoverinfo="none",
                                 marker_color="black"))
        data.append(go.Scatter3d(x=[0, 0], y=[0, c_max], z=[0, 0], mode="lines",
                                 hoverinfo="none",
                                 marker_color="black"))
        data.append(go.Scatter3d(x=[0, 0], y=[0, 0], z=[0, c_max], mode="lines",
                                 hoverinfo="none",
                                 marker_color="black"))

        for i, direct in zip(self._bz_plot_info.rec_lat_vec, ["kx", "ky", "kz"]):
            nn = sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2)
            kx, ky, kz = np.array(i) * c_max / nn * 1.15
            norm_kx, norm_ky, norm_kz = np.array(i) * c_cone / nn
            data.append(go.Cone(x=[kx], y=[ky], z=[kz],
                                u=[norm_kx], v=[norm_ky], w=[norm_kz],
                                colorscale=[[0, 'rgb(100,100,100)'],
                                            [1, 'rgb(100,100,100)']],
                                showscale=False, hovertemplate=f"<b>{direct}</b>"))
            data.append(go.Scatter3d(x=[0, kx], y=[0, ky], z=[0, kz],
                                     mode="lines",
                                     hoverinfo="none",
                                     marker_color="black"))

        range_max = c_max * 1.4
        fig = go.Figure(data=data)

        fig.update_layout(
            title_font_size=30,
            font_size=24,
            width=900, height=900,
            showlegend=False,
            scene=dict(xaxis=dict(showspikes=False, range=[-range_max, range_max],
                                  showticklabels=False, visible=False),
                       yaxis=dict(showspikes=False, range=[-range_max, range_max],
                                  showticklabels=False, visible=False),
                       zaxis=dict(showspikes=False, range=[-range_max, range_max],
                                  showticklabels=False, visible=False))
        )

        return fig


