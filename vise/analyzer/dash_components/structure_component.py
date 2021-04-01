# coding: utf-8
#  Copyright (c) 2020 Kumagai group.
import re
from collections import OrderedDict
from itertools import combinations_with_replacement, chain
from typing import Tuple

import numpy as np
from crystal_toolkit.core.legend import Legend
from crystal_toolkit.core.mpcomponent import MPComponent
from crystal_toolkit.core.scene import Scene
from crystal_toolkit.helpers.layouts import *
from crystal_toolkit.settings import SETTINGS
from dash_mp_components import CrystalToolkitScene
from pymatgen.analysis.graphs import StructureGraph, MoleculeGraph
from pymatgen.analysis.local_env import NearNeighbors
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.core.structure import Structure, Molecule

DEFAULTS = {
    "bonding_strategy": "CrystalNN",
    "radius_strategy": "specified_or_average_ionic",
    "draw_image_atoms": True,
    "bonded_sites_outside_unit_cell": False,
    "hide_incomplete_bonds": True,
    "show_compass": True,
    "unit_cell_choice": "input",
}


class StructureComponent(MPComponent):
    """
    A component to display pymatgen Structure objects.
    """

    available_bonding_strategies = \
        {subcls.__name__: subcls for subcls in NearNeighbors.__subclasses__()}

    default_scene_settings = {
        "extractAxis": True,
        # For visual diff testing, we change the renderer
        # to SVG since this WebGL support is more difficult
        # in headless browsers / CI.
        "renderer": "svg" if SETTINGS.TEST_MODE else "webgl",
    }

    default_title = "Structure"

    def __init__(
            self,
            structure: Optional[Structure] = None,
            id: str = None,
            scene_additions: Optional[Scene] = None,
            bonding_strategy: str = DEFAULTS["bonding_strategy"],
            bonding_strategy_kwargs: Optional[dict] = None,
            color_scale: Optional[str] = None,
            radius_strategy: str = DEFAULTS["radius_strategy"],
            unit_cell_choice: str = DEFAULTS["unit_cell_choice"],
            draw_image_atoms: bool = DEFAULTS["draw_image_atoms"],
            bonded_sites_outside_unit_cell: bool = DEFAULTS[
                "bonded_sites_outside_unit_cell"
            ],
            hide_incomplete_bonds: bool = DEFAULTS["hide_incomplete_bonds"],
            show_compass: bool = DEFAULTS["show_compass"],
            scene_settings: Optional[Dict] = None,
            **kwargs,
    ):

        super().__init__(id=id, default_data=structure, **kwargs)

        self.initial_scene_settings = self.default_scene_settings.copy()
        if scene_settings:
            self.initial_scene_settings.update(scene_settings)

        self.create_store("scene_settings",
                          initial_data=self.initial_scene_settings)

        # unit cell choice and bonding algorithms need to come from a settings
        # object (in a dcc.Store) guaranteed to be present in layout, rather
        # than from the controls themselves -- since these are optional and
        # may not be present in the layout
        self.create_store(
            "graph_generation_options",
            initial_data={"bonding_strategy": bonding_strategy,
                          "bonding_strategy_kwargs": bonding_strategy_kwargs,
                          "unit_cell_choice": unit_cell_choice,
            },
        )

        self.create_store(
            "display_options",
            initial_data={
                "color_scale"                   : color_scale,
                "radius_strategy"               : radius_strategy,
                "draw_image_atoms"              : draw_image_atoms,
                "bonded_sites_outside_unit_cell": bonded_sites_outside_unit_cell,
                "hide_incomplete_bonds"         : hide_incomplete_bonds,
                "show_compass"                  : show_compass,
            },
        )

        if scene_additions:
            initial_scene_additions = Scene(name="scene_additions",
                                            contents=scene_additions).to_json()
        else:
            initial_scene_additions = None
        self.create_store("scene_additions",
                          initial_data=initial_scene_additions)

        if structure:
            # graph is cached explicitly, this isn't necessary but is an
            # optimization so that graph is only re-generated if bonding
            # algorithm changes
            graph = self._preprocess_input_to_graph(
                structure,
                bonding_strategy=bonding_strategy,
                bonding_strategy_kwargs=bonding_strategy_kwargs,
            )
            scene, legend = self.get_scene_and_legend(
                graph,
                scene_additions=self.initial_data["scene_additions"],
                **self.initial_data["display_options"],
            )
            if hasattr(structure, "lattice"):
                self._lattice = structure.lattice
        else:
            # component could be initialized without a structure, in which case
            # an empty scene should be displayed
            graph = None
            scene, legend = self.get_scene_and_legend(
                None,
                scene_additions=self.initial_data["scene_additions"],
                **self.initial_data["display_options"])

        self.create_store("legend_data", initial_data=legend)
        self.create_store("graph", initial_data=graph)

        # this is used by a Simple3DScene component, not a dcc.Store
        self._initial_data["scene"] = scene

        # hide axes inset for molecules
        if isinstance(structure, Molecule) or isinstance(
                structure, MoleculeGraph
        ):
            self.scene_kwargs = {"axisView": "HIDDEN"}
        else:
            self.scene_kwargs = {}

    def _make_legend(self, legend):

        if not legend:
            return html.Div(id=self.id("legend"))

        def get_font_color(hex_code):
            # ensures contrasting font color for background color
            c = tuple(int(hex_code[1:][i: i + 2], 16) for i in (0, 2, 4))
            if 1 - (c[0] * 0.299 + c[1] * 0.587 + c[2] * 0.114) / 255 < 0.5:
                font_color = "#000000"
            else:
                font_color = "#ffffff"
            return font_color

        try:
            formula = Composition.from_dict(
                legend["composition"]).reduced_formula
        except:
            # TODO: fix legend for Dummy Specie compositions
            formula = "Unknown"

        legend_colors = OrderedDict(
            sorted(list(legend["colors"].items()),
                   key=lambda x: formula.find(x[1]))
        )

        legend_elements = [
            Button(
                html.Span(
                    name, className="icon",
                    style={"color": get_font_color(color)}
                ),
                kind="static",
                style={"backgroundColor": color},
            )
            for color, name in legend_colors.items()
        ]

        return Field(
            [Control(el, style={"marginRight": "0.2rem"}) for el in
             legend_elements],
            id=self.id("legend"),
            grouped=True,
        )

    def _make_title(self, legend):

        if not legend or (not legend.get("composition", None)):
            return H1(self.default_title, id=self.id("title"))

        composition = legend["composition"]
        if isinstance(composition, dict):

            try:
                composition = Composition.from_dict(composition)

                # strip DummySpecie if present (TODO: should be method in pymatgen)
                composition = Composition(
                    {
                        el: amt
                        for el, amt in composition.items()
                        if not isinstance(el, DummySpecie)
                    }
                )
                composition = composition.get_reduced_composition_and_factor()[
                    0]
                formula = composition.reduced_formula
                formula_parts = re.findall(r"[^\d_]+|\d+", formula)
                formula_components = [
                    html.Sub(part.strip())
                    if part.isnumeric()
                    else html.Span(part.strip())
                    for part in formula_parts
                ]
            except:
                formula_components = list(map(str, composition.keys()))

        return H1(
            formula_components, id=self.id("title"),
            style={"display": "inline-block"}
        )

    @staticmethod
    def _make_bonding_algorithm_custom_cuffoff_data(graph):
        if not graph:
            return [{"A": None, "B": None, "A—B": None}]
        struct_or_mol = StructureComponent._get_structure(graph)
        # can't use type_of_specie because it doesn't work with disordered structures
        species = set(
            map(
                str,
                chain.from_iterable(
                    [list(c.keys()) for c in struct_or_mol.species_and_occu]
                ),
            )
        )
        rows = [
            {"A": combination[0], "B": combination[1], "A—B": 0}
            for combination in combinations_with_replacement(species, 2)
        ]
        return rows

    @property
    def _sub_layouts(self):

        struct_layout = html.Div(
            CrystalToolkitScene(
                id=self.id("scene"),
                data=self.initial_data["scene"],
                settings=self.initial_scene_settings,
                sceneSize="100%",
                **self.scene_kwargs,
            ),
            style={
                "width"   : "100%",
                "height"  : "100%",
                "overflow": "hidden",
                "margin"  : "0 auto",
            },
        )

        title_layout = html.Div(
            self._make_title(self._initial_data["legend_data"]),
            id=self.id("title_container"),
        )

        legend_layout = html.Div(
            self._make_legend(self._initial_data["legend_data"]),
            id=self.id("legend_container"),
        )

        return {"struct": struct_layout,
                "title": title_layout,
                "legend": legend_layout}

    def layout(self, size: str = "500px") -> html.Div:
        """
        :param size: a CSS string specifying width/height of Div
        :return: A html.Div containing the 3D structure or molecule
        """
        return html.Div(
            self._sub_layouts["struct"], style={"width": size, "height": size}
        )

    @staticmethod
    def _preprocess_input_to_graph(
            input: Union[Structure, StructureGraph],
            bonding_strategy: str = DEFAULTS["bonding_strategy"],
            bonding_strategy_kwargs: Optional[Dict] = None,
    ) -> Union[StructureGraph, MoleculeGraph]:

        # ensure fractional co-ordinates are normalized to be in [0,1)
        # (this is actually not guaranteed by Structure)
        try:
            input = input.as_dict(verbosity=0)
        except TypeError:
            # TODO: remove this, necessary for Slab(?), some structure subclasses don't have verbosity
            input = input.as_dict()
        for site in input["sites"]:
            site["abc"] = np.mod(site["abc"], 1)
        input = Structure.from_dict(input)

        if not input.is_ordered:
            # calculating bonds in disordered structures is currently very flaky
            bonding_strategy = "CutOffDictNN"

        if (bonding_strategy
                not in StructureComponent.available_bonding_strategies.keys()):
            raise ValueError(
                "Bonding strategy not supported. Please supply a name "
                "of a NearNeighbor subclass, choose from: {}".format(
                    ", ".join(
                        StructureComponent.available_bonding_strategies.keys()
                    )
                )
            )
        else:
            bonding_strategy_kwargs = bonding_strategy_kwargs or {}
            if bonding_strategy == "CutOffDictNN":
                if "cut_off_dict" in bonding_strategy_kwargs:
                    # TODO: remove this hack by making args properly JSON serializable
                    bonding_strategy_kwargs["cut_off_dict"] = {
                        (x[0], x[1]): x[2]
                        for x in bonding_strategy_kwargs["cut_off_dict"]
                    }
            bonding_strategy = StructureComponent.available_bonding_strategies[bonding_strategy](**bonding_strategy_kwargs)
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    graph = StructureGraph.with_local_env_strategy(input, bonding_strategy)
            except:
                # for some reason computing bonds failed, so let's not have any bonds(!)
                graph = StructureGraph.with_empty_graph(input)

        return graph

    def generate_callbacks(self, app, cache):
        pass

    @staticmethod
    def _get_structure(graph: Union[StructureGraph, Structure]) -> Union[Structure, Molecule]:
        if isinstance(graph, StructureGraph):
            return graph.structure
        elif isinstance(graph, Structure):
            return graph
        else:
            raise ValueError

    @staticmethod
    def get_scene_and_legend(
            graph: Optional[Union[StructureGraph, MoleculeGraph]],
            color_scale=None,
            radius_strategy=DEFAULTS["radius_strategy"],
            draw_image_atoms=DEFAULTS["draw_image_atoms"],
            bonded_sites_outside_unit_cell=DEFAULTS[
                "bonded_sites_outside_unit_cell"],
            hide_incomplete_bonds=DEFAULTS["hide_incomplete_bonds"],
            explicitly_calculate_polyhedra_hull=False,
            scene_additions=None,
            show_compass=DEFAULTS["show_compass"],
    ) -> Tuple[Scene, Dict[str, str]]:

        scene = Scene(name="AceStructureMoleculeComponentScene")

        if graph is None:
            return scene, {}

        structure = StructureComponent._get_structure(graph)

        # TODO: add radius_scale
        legend = Legend(
            structure,
            color_scheme="VESTA",
            radius_scheme=radius_strategy,
            cmap_range=color_scale,
        )

        scene = graph.get_scene(
                draw_image_atoms=draw_image_atoms,
                bonded_sites_outside_unit_cell=bonded_sites_outside_unit_cell,
                hide_incomplete_edges=hide_incomplete_bonds,
                explicitly_calculate_polyhedra_hull=explicitly_calculate_polyhedra_hull,
                legend=legend,
            )

        scene.name = "StructureComponentScene"

        if hasattr(structure, "lattice"):
            axes = structure.lattice._axes_from_lattice()
            axes.visible = show_compass
            scene.contents.append(axes)

        scene = scene.to_json()
        if scene_additions:
            # TODO: need a Scene.from_json() to make this work
            # raise NotImplementedError
            scene["contents"].append(scene_additions)

        return scene, legend.get_legend()

    def title_layout(self):
        """
        :return: A layout including the composition of the structure/molecule as a title.
        """
        return self._sub_layouts["title"]

    def legend_layout(self):
        """
        :return: A layout including a legend for the structure/molecule.
        """
        return self._sub_layouts["legend"]
