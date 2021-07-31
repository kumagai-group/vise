# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

from dash import Dash
from monty.serialization import loadfn
from vise.analyzer.dash_components.band_dos_dash import BandDosComponent
import dash_html_components as html

import crystal_toolkit.components as ctc

#
#def test_plotly_actual_plot(test_files):
#     dos_data = loadfn(str(test_files / "KAlSi3O8_dos_data.json"))
#     pploter = PlotlyDosPlotter(dos_data)
#     fig = pploter.show({"O": [10, 11]})
#     fig.show()


app = Dash(suppress_callback_exceptions=True)
ctc.register_app(app)

test_files = Path(__file__).parent / ".." / ".." / "tests" / "test_data_files"
band_plot_data = loadfn(str(test_files / "KAlSi3O8_band_plot_info.json"))
dos_plot_data = loadfn(str(test_files / "KAlSi3O8_dos_plot_data.json"))
band_dos_component = BandDosComponent(dos_plot_data, band_plot_data)
#dos_component = BandDosComponent(None, band_plot_info)
#dos_component = BandDosComponent(dos_plot_data, band_plot_info, band_plot_info_2)

# example layout to demonstrate capabilities of component
my_layout = html.Div(
    [
        html.H1("BandDosComponent Example"),
        band_dos_component.layout(),
    ]
)
# {'K_i1': [0, 1],
#  'Al_i1': [2, 3],
#  'Si_i1': [4, 5],
#  'Si_i2': [6, 7],
#  'Si_i3': [8, 9],
#  'O_i1': [10, 11],
#  'O_i2': [12, 13],
#  'O_i3': [14, 15],
#  'O_i4': [16, 17],
#  'O_i5': [18, 19],
#  'O_i6': [20, 21],
#  'O_i7': [22, 23],
#  'O_i8': [24, 25]}
    # wrap your app.layout with crystal_toolkit_layout()
    # to ensure all necessary components are loaded into layout
# app.layout = ctc.crystal_toolkit_layout(my_layout)
ctc.register_crystal_toolkit(app=app, layout=my_layout, cache=None)

    # allow app to be run using "python structure.py"
    # in production, deploy behind gunicorn or similar
    # see Dash documentation for more information
if __name__ == "__main__":
    app.run_server(debug=True, port=8051)


"""
TODO
- Enable to select inequivalent sites shown
- Enable to change the max range of dos

DONE
- Show vbm and cbm
- Show both pbesol and dd-hybrid.
- Change size of dos (smaller width)
- Set max range of default DOS.
- Add lines 
- Add band figure
- Add x and y labels
- title
- show s, p, d orbitals
- show multiple figures for K, Al, Si, and O
"""