# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
import dash
from pydefect.analyzer.band_edge_states import BandEdgeEigenvalues
from pydefect.analyzer.defect_energy import DefectEnergy
from pydefect.analyzer.defect_energy_plotter import DefectEnergyPlotlyPlotter
from pydefect.analyzer.eigenvalue_plotter import EigenvaluePlotlyPlotter
from pydefect.corrections.efnv_correction.efnv_correction import PotentialSite, \
    ExtendedFnvCorrection
from pydefect.corrections.efnv_correction.site_potential_plotter import \
    SitePotentialPlotlyPlotter
from pymatgen import loadfn, Structure, Lattice
from pymatgen.util.string import latexify
from vise.analyzer.dash_components.structure_component import StructureComponent
from vise.analyzer.dos_data import DosData
import dash_html_components as html

import crystal_toolkit.components as ctc

app = dash.Dash()
ctc.register_app(app)

structure = Structure(lattice=Lattice.cubic(10), species=["H"], coords=[[0.0]*3])

component = StructureComponent(structure)
my_layout = html.Div(
    [
        component.layout(),
    ]
)

# app.layout = ctc.crystal_toolkit_layout(my_layout)
ctc.register_crystal_toolkit(app=app, layout=my_layout, cache=None)


if __name__ == "__main__":
    app.run_server(debug=True, port=8051)

