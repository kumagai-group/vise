# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Dict

from monty.json import MSONable

import plotly.graph_objects as go
from vise.analyzer.plot_band_dos import plotly_sanitize_label
from vise.util.plotly_util import sort_coords, make_triangles
import numpy as np


@dataclass
class BZPlotInfo(MSONable):
    faces: List[List[List[float]]]
    labels: Dict[str, List[float]]


class BZPlotlyPlotter:
    def __init__(self, bz_plot_info: BZPlotInfo):
        self._bz_plot_info = bz_plot_info

    def create_figure(self):
        data = []
        for face_vertices in self._bz_plot_info.faces:
            vertex_coords = sort_coords(np.array(face_vertices))
            data.append(go.Mesh3d(**make_triangles(vertex_coords),
                                  alphahull=-1, opacity=0.3, hoverinfo='skip',
                                  color="blue"))

        for label, coords in self._bz_plot_info.labels.items():
            data.append(go.Scatter3d(x=[coords[0]], y=[coords[1]], z=[coords[2]],
                                     text=[plotly_sanitize_label(label)],
                                     mode="markers+text",
                                     marker=dict(color="orange", size=5),
                                     textposition="middle center",
                                     textfont=dict(size=24),
                                     hoverinfo='skip'))

        fig = go.Figure(data=data)
        fig.update_layout(
            title=f"Brillouin Zone",
            title_font_size=30,
            font_size=24,
            width=700, height=700,
            showlegend=False)

        return fig


