# coding: utf-8
#  Copyright (c) 2020 Kumagai group.
import argparse
import sys
from fractions import Fraction
from pathlib import Path

import crystal_toolkit.components as ctc
import numpy as np
from crystal_toolkit.components import StructureMoleculeComponent
from crystal_toolkit.helpers.layouts import *
from dash import Dash
from pymatgen.core import Structure, Lattice, Molecule
from pymatgen.ext.matproj import MPRester

from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import \
    SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import \
    AllCoordinationGeometries
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import \
    LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import \
    LightStructureEnvironments
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.structure_matcher import StructureMatcher
from vise.analyzer.dash_components.band_dos_dash import BandDosComponent
from vise.analyzer.dash_components.structure_component import StructureComponent
from vise.analyzer.dos_data import DosPlotData
from vise.analyzer.plot_band import BandPlotInfo
from vise.defaults import defaults
from vise.util.structure_symmetrizer import StructureSymmetrizer


def create_ctk(struct_component,
               crystal_symmetry,
               sites,
               band_dos_component):
    box_size = "30vmin"
    _layout = Container(
        [
            Section(
                [
                    Columns(
                        [
                            Column(
                                [struct_component.title_layout()]
                            )
                        ]
                    ),
                    Columns(
                        [
                            Column(
                                [
                                    Box(
                                        struct_component.layout(size="100%"),
                                        style={
                                            "width": box_size,
                                            "height": box_size,
                                            "minWidth": "300px",
                                            "minHeight": "300px",
                                            "maxWidth": "600px",
                                            "maxHeight": "600px",
                                            "overflow": "hidden",
                                            "padding": "0.25rem",
                                            "marginBottom": "0.5rem",
                                        },
                                    ),
                                    html.Div(
                                        [
                                            html.Div(
                                                struct_component.legend_layout(),
                                                style={"float": "left"},
                                            ),
                                        ],
                                        style={
                                            "width": box_size,
                                            "minWidth": "300px",
                                            "marginBottom": "40px",
                                        },
                                    ),
                                ],
                                narrow=True,
                            ),
                            Column(
                                [
                                    crystal_symmetry,
                                ],
                                style={"width": box_size, "max-width": box_size},
                            ),
                        ],
                        desktop_only=False,
                        centered=False,
                    ),
                    Columns(Column([sites])),
                ]
            ),
            Section([band_dos_component.layout()]
            ),
        ]
    )

    return _layout


def symmetry_layout(structure: Structure):
    data = dict()
    symmetrizer = StructureSymmetrizer(structure)
    data["Space Group"] = f"{symmetrizer.spglib_sym_data['international']} " \
                          f"({symmetrizer.sg_number})"
    data["Point Group"] = symmetrizer.point_group

    datalist = get_data_list(data)

    return Columns([Column([H4("Crystal symmetry"), datalist])])


def mpid_and_link(symmetrizer: StructureSymmetrizer):

    reduced_formula = symmetrizer.structure.composition.reduced_formula
    criteria = {"reduced_cell_formula": reduced_formula,
                "spacegroup.number": symmetrizer.sg_number}
    properties = ["task_id", "structure"]
    with MPRester() as m:
        for doc in m.query(criteria, properties):
            sm = StructureMatcher()
            if sm.fit(doc["structure"], symmetrizer.structure):
                mpid = doc["task_id"]
                break
        else:
            return H4("None")

    return dcc.Link(f'mp-id {mpid}',
                    href=f'https://materialsproject.org/materials/{mpid}/',
                    style={'font-weight': 'bold', "font-size": "20px"})


def pretty_frac_format(x):
    x = x % 1
    fraction = Fraction(x).limit_denominator(8)
    if np.allclose(x, 1):
        x_str = "0"
    elif not np.allclose(x, float(fraction)):
        x = np.around(x, decimals=3)
        x_str = f"{x:.3g}"
    else:
        x_str = str(fraction)
    return x_str


def site_layout(symmetrizer: StructureSymmetrizer):
    columns = []
    distance_cutoff = 2.0
    angle_cutoff = 0.3
    cnn = CrystalNN()
    nn_info = cnn.get_all_nn_info(structure=symmetrizer.structure)

    for name, site in symmetrizer.sites.items():
        wyckoff_contents = []
        data = dict()
        data["Site Symmetry"] = site.site_symmetry
        data["Wyckoff Letter"] = site.wyckoff_letter
        wyckoff_contents.append(html.Label(f"{name}", className="mpc-label"))
        site_data = []
        repr_idx = site.equivalent_atoms[0]
        for idx in site.equivalent_atoms:
            coords = symmetrizer.structure[idx].frac_coords
            site_data += [(pretty_frac_format(coords[0]),
                           pretty_frac_format(coords[1]),
                           pretty_frac_format(coords[2]))]

        # lgf = LocalGeometryFinder()
        # lgf.setup_structure(structure=symmetrizer.structure)

        # se = lgf.compute_structure_environments(
        #     maximum_distance_factor=distance_cutoff + 0.01)
        # strategy = SimplestChemenvStrategy(
        #     distance_cutoff=distance_cutoff, angle_cutoff=angle_cutoff
        # )
        # lse = LightStructureEnvironments.from_structure_environments(
        #     strategy=strategy, structure_environments=se)

    # represent the local environment as a molecule

        mol = Molecule.from_sites(
            [symmetrizer.structure[repr_idx]] + [i["site"] for i in nn_info[repr_idx]]
        )
        mol = mol.get_centered_molecule()
        mg = MoleculeGraph.with_empty_graph(molecule=mol)
        for i in range(1, len(mol)):
            mg.add_edge(0, i)

        view = html.Div(
            [
                StructureMoleculeComponent(
                    struct_or_mol=mg,
                    disable_callbacks=True,
                    id=f"{symmetrizer.structure.composition.reduced_formula}_site_{repr_idx}",
                    scene_settings={"enableZoom": False, "defaultZoom": 0.6},
                )._sub_layouts["struct"]
            ],
            style={"width": "300px", "height": "300px"},
        )

        # env = lse.coordination_environments[repr_idx]
        # all_ce = AllCoordinationGeometries()
        # co = all_ce.get_geometry_from_mp_symbol(env[0]["ce_symbol"])
        # name = co.name

        data.update(
            {
                # "Environment": name,
                # "IUPAC Symbol": co.IUPAC_symbol_str,
                "Structure": view,
            }
        )

        wyckoff_contents.append(get_data_list(data))
        wyckoff_contents.append(Reveal(get_table(site_data),
                                       title=f"Positions ({len(site_data)})",
                                       id=f"{name}-positions"))
        columns.append(Column(html.Div(wyckoff_contents)))
    _layout = Columns(columns)

    return Columns([Column([H4("Sites"), html.Div(_layout)])])


def create_app(structure: Structure, dos_plot_data, band_plot_data):
    app = Dash(suppress_callback_exceptions=True,
               assets_folder=SETTINGS.ASSETS_PATH)

    structure_component = StructureComponent(structure)
    comp = structure.composition.reduced_formula
    band_dos_component = BandDosComponent(dos_plot_data, band_plot_data, id=f"band_dos_{comp}",)
    symmetrizer = StructureSymmetrizer(structure)
    delayed_layout = create_ctk(structure_component,
                                symmetry_layout(structure),
                                site_layout(symmetrizer),
                                band_dos_component)
    ctc.register_crystal_toolkit(layout=delayed_layout, app=app, cache=None)

    return app


def parse_args(args):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-b", "--band", type=Path)
    parser.add_argument("-d", "--dos", type=Path)
    parser.add_argument("-p", "--port", type=int)
    return parser.parse_args(args)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    structure = Structure.from_file("POSCAR")
    dos_plot_data = loadfn(args.dos / "dos_plot_data.json")
    band_plot_data = loadfn(args.band / "band_plot_info.json")
    my_app = create_app(structure, dos_plot_data, band_plot_data)
    my_app.run_server(debug=True, port=args.port)
    # symmetrizer = StructureSymmetrizer(structure)
    #
    # app = Dash(suppress_callback_exceptions=True,
    #            assets_folder=SETTINGS.ASSETS_PATH)
    # ctc.register_crystal_toolkit(layout=html.Div([mpid_and_link(symmetrizer)]),
    #                              app=app, cache=None)
    # app.run_server(debug=True, port=args.port)

