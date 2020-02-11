#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.
import argparse
from collections import OrderedDict
from copy import deepcopy
import os
import string
from typing import Optional, List, Dict, Tuple, Union, Any

import ruamel.yaml
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from scipy.spatial import HalfspaceIntersection
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition, CompositionError
from pymatgen.analysis.phase_diagram import PhaseDiagram
from vise.chempotdiag.compound import (
    Compound, DummyCompoundForDiagram, CompoundsList, ElemOrderType)
from vise.chempotdiag.vertex import (
    Vertex, VertexOnBoundary, VerticesList)

molecule_directory = os.path.dirname(__file__) + "/molecules"


# TODO: Need to write test


class ChemPotDiag:
    """ Object for chemical potentials diagram. """

    def __init__(self,
                 element_free_energy: dict,
                 stable_compounds: CompoundsList,
                 unstable_compounds: CompoundsList,
                 vertices: VerticesList,
                 compounds_to_vertex_list: List[List[int]],
                 vertex_to_compounds_list: List[List[int]],
                 temperature: float = 0,
                 pressure: Optional[Dict[str, float]] = None):
        """

        Args:
            element_free_energy (dict):
            stable_compounds (CompoundsList):
            unstable_compounds (CompoundsList):
            vertices (VerticesList):
            compounds_to_vertex_list (List[List[int]]):
            vertex_to_compounds_list (List[List[int]]):
        """
        self.element_free_energy = element_free_energy
        self.stable_compounds = stable_compounds
        self.unstable_compounds = unstable_compounds
        self.vertices = vertices
        self.compounds_to_vertex_list = compounds_to_vertex_list
        self.vertex_to_compounds_list = vertex_to_compounds_list
        self.temperature = temperature
        self.pressure = pressure

    @classmethod
    def from_compound_list(cls,
                           compound_list: CompoundsList,
                           temperature: Optional[float] = None,
                           pressure: Optional[Dict[str, float]] = None
                           ) -> "ChemPotDiag":
        """Create a object of ChemPot.

        Args:
            compound_list (CompoundsList):
                List of considered compounds
            temperature (float):
            pressure (Dict[str, float]):
        """
        compound_list = deepcopy(compound_list)

        # Alphabetically sort elements' names
        elements = sorted(compound_list.elements)
        dim = len(elements)

        # Standardize energies
        compound_list, element_energy = \
            compound_list.standard_energy_array(temperature, pressure)

        #  unary system
        if dim == 1:
            min_energy = min(compound_list, key=lambda x: x.energy).energy
            stable_compounds \
                = CompoundsList([comp for comp in compound_list
                                 if comp.energy == min_energy])
            unstable_compounds \
                = CompoundsList([comp for comp in compound_list
                                 if comp.energy != min_energy])
            vertices = VerticesList([])
            for comp in compound_list:
                d = {list(compound_list.elements)[0]: comp.energy}
                vertex = Vertex(None, d)
                vertices.append(vertex)

            compounds_to_vertex_list = [[i] for i in range(len(compound_list))]
            vertex_to_compounds_list = [[i] for i in range(len(compound_list))]
            return cls(element_energy, stable_compounds, unstable_compounds,
                       vertices,
                       compounds_to_vertex_list, vertex_to_compounds_list)

        # set boundary range
        intersections = np.zeros((len(compound_list), dim))
        free_energies = compound_list.free_energies(temperature, pressure)
        for i, (comp_dat, fe) in enumerate(zip(compound_list, free_energies)):
            c = comp_dat.composition_vector(elements)
            for j in range(dim):
                if c[j] == 0:  # This does not related to search_vertices_range.
                    intersections[i][j] = float("inf")
                else:
                    intersections[i][j] = fe / c[j]

        boundary_rate = 1.1  # can be arbitrary number larger than 1.0
        search_vertices_range = np.min(intersections) * boundary_rate

#        elements = compound_list.elements
        for e in compound_list.elements:
            energy = search_vertices_range
            boundary = DummyCompoundForDiagram.construct_boundary(e, -energy)
            compound_list.append(boundary)
            free_energies.append(-energy)

        # Calculate HalfspaceIntersection
        eps = 1e-2
        interior_point = np.array([search_vertices_range + eps] * dim)

        halfspaces = []
        for i, comp_dat in enumerate(compound_list):
            comp = comp_dat.composition_vector(compound_list.elements)
            halfspaces.append(np.append(comp, -free_energies[i]))

        halfspaces = np.array(halfspaces)
        hi = HalfspaceIntersection(halfspaces, interior_point)
        # facets_by_halfspace
        n = len(hi.halfspaces)
        facets_by_halfspace = []
        for i in range(n):
            indices = []
            for j, facet in enumerate(hi.dual_facets):
                if i in facet:
                    indices.append(j)
            facets_by_halfspace.append(indices)

        # classify all compounds into dummy, stable, and unstable.
        stable_compounds = CompoundsList([])
        unstable_compounds = CompoundsList([])
        stable_original_index_list = []
        unstable_original_index_list = []
        which_vertex_on_boundary = set()
        for i, compound in enumerate(compound_list):
            if isinstance(compound, DummyCompoundForDiagram):
                which_vertex_on_boundary |= set(facets_by_halfspace[i])
            elif facets_by_halfspace[i]:
                stable_compounds.append(compound_list[i])
                stable_original_index_list.append(i)
            else:
                unstable_compounds.append(compound_list[i])
                unstable_original_index_list.append(i)

        # make vertices_array
        a = [item for item in hi.intersections.flatten() if abs(item - search_vertices_range) > 1e-6]
        draw_criterion = min(a)
        draw_vertices = deepcopy(hi.intersections)
        draw_range = draw_criterion * 1.1

        # sometimes boundary range too small like -100
        # due to an unstable substance, the value is changed to draw_range
        vertices = VerticesList([])
        for i in range(len(draw_vertices)):
            if i in which_vertex_on_boundary:

                is_element_boundary = [abs(v - search_vertices_range) < 1e-6 for v in draw_vertices[i]]
                for j, flag in enumerate(is_element_boundary):
                    if flag:
                        draw_vertices[i][j] = draw_range
                    vertex = VertexOnBoundary(None,
                                              {e: en for e, en in zip(elements, draw_vertices[i])},
                                              draw_range,
                                              draw_criterion
                    )
            else:
                vertex = Vertex(
                    None, {e: en for e, en in zip(elements, draw_vertices[i])})
            vertices.append(vertex)

        # make compounds_to_vertex_list, vertex_to_compounds_list
        compounds_to_vertex_list = [l for l in facets_by_halfspace if l]
        vertex_to_compounds_list = []
        for i, l in enumerate(hi.dual_facets):
            stable_index_list = [stable_original_index_list.index(j) for j in l if j in stable_original_index_list]
            vertex_to_compounds_list.append(stable_index_list)

        return cls(element_energy, stable_compounds, unstable_compounds,
                   vertices,
                   compounds_to_vertex_list, vertex_to_compounds_list,
                   temperature=temperature, pressure=pressure)

    @classmethod
    def from_file(cls,
                  filename: str,
                  temperature: Optional[float] = None,
                  pressure: Optional[dict] = None) -> "ChemPotDiag":
        """
        Args:
            filename (str):
            temperature (float): (K)
            pressure (dict): e.g. {"O2": 1e+3} (Pa)

        Returns:
        """
        compounds_array = CompoundsList.from_file(filename)
        return cls.from_compound_list(compounds_array, temperature, pressure)

    @classmethod
    def from_vasp_calculations_files(
            cls,
            poscar_paths: List[str],
            output_paths: List[str],
            fmt: str = "outcar",
            temperature: Optional[float] = None,
            pressure: Optional[Dict[str, float]] = None,
            energy_shift_dict: Union[Dict[str, float], None] = None
            ) -> "ChemPotDiag":
        """
        Args:
            poscar_paths (list of str):
            output_paths (list of str):
            fmt (str): "outcar" or "vasprun".
            temperature (float): (K)
            pressure (dict): (Pa)
            energy_shift_dict (dict): unit: (eV), e.g. {N2_molecule/OUTCAR: 1}
        Returns:
            (ChemPotDiag) ChemPotDiag object from vasp files.
        """
        energy_shift_dict = energy_shift_dict if energy_shift_dict else {}
        compounds_list = \
            CompoundsList.from_vasp_files(poscar_paths,
                                          output_paths,
                                          fmt=fmt,
                                          energy_shift_dict=
                                                       energy_shift_dict)
        return cls.from_compound_list(compounds_list,
                                      temperature=temperature,
                                      pressure=pressure)

    @classmethod
    def from_vasp_and_materials_project(
            cls,
            vasp_target_poscar: str,
            vasp_target_output: str,
            vasp_element_poscar: List[str],
            vasp_element_output: List[str],
            fmt: str = "outcar",
            temperature: Optional[float] = None,
            pressure: Union[str, float, None] = None,
            energy_shift_dict: Union[Dict[str, float], None] = None
            ) -> "ChemPotDiag":
        """

        Args:
            vasp_target_poscar (str): path
            vasp_target_output (str): path
            vasp_element_poscar (list): path
            vasp_element_output (list): path
            fmt (str):
            temperature (float):
            pressure (dict):
            energy_shift_dict (dict):

        Returns:
            (CompoundsList) CompoundsList object from materials project.

        """
        compounds_list = CompoundsList.\
            from_vasp_and_materials_project(vasp_target_poscar,
                                            vasp_target_output,
                                            vasp_element_poscar,
                                            vasp_element_output,
                                            fmt=fmt,
                                            energy_shift_dict=energy_shift_dict)
        return cls.from_compound_list(compounds_list, temperature, pressure)

    @property
    def dim(self) -> int:
        """ (int) Number of considered atoms. """
        return self.stable_compounds.dim

    @property
    def elements(self) -> List[Element]:
        """Considered elements. Order of list is used as order of coordinates"""
        return self.stable_compounds.elements

    @property
    def all_compounds(self) -> CompoundsList:
        """ (CompoundsList) CompoundList """
        return CompoundsList(self.stable_compounds + self.unstable_compounds)

    @property
    def equilibrium_points(self) -> VerticesList:
        """
            (VerticesList) Vertices, excluding vertices on drawing boundary
                           and only physically meaningful points are included.
        """
        return VerticesList(
            [v for v in self.vertices if not isinstance(v, VertexOnBoundary)])

    def get_neighbor_vertices(self,
                              compound: Union[int, str, Compound],
                              elements: Optional[ElemOrderType] = None
                              ) -> VerticesList:
        """ Search equilibrium_points on remarked compound.

        Args:
            compound (int/str/Compound):
                Index of compound or compound name or Compound object.
            elements (None/[Element or str]):
                Order of elements.

        Returns (VerticesList):
            Vertices on input compound.
        """
        index = None
        if isinstance(compound, int):
            index = compound
        elif isinstance(compound, str):
            matched = self.stable_compounds.get_indices_and_compounds(compound)
            if matched is None:
                raise ValueError(f"No such compounds in diagram {compound}")
            if len(matched) >= 2:
                raise ValueError(f"More than one compounds matched {compound}.")
            if matched:
                index = matched[0][0]
        elif isinstance(compound, Compound):
            index = self.stable_compounds.index(compound)
        else:
            raise TypeError("Input compound must be int or str or Compound.")

        if index is None:
            raise ValueError(f"{compound} did not found in stable_compounds.")

        to_return = VerticesList(
            [self.vertices[i] for i in self.compounds_to_vertex_list[index]])
        if self.dim == 3:
            to_return = to_return.sorted_to_loop_in_3d(elements)

        return to_return

    def get_neighbor_compounds(self,
                               vertex: Union[int, str, Vertex]
                               ) -> CompoundsList:
        """ Search equilibrium_points on remarked vertex.

        Args:
            vertex (int/str/Vertex)/
                Input compound (index or label or vertex).
        Returns (CompoundsList):
            Compounds on the vertex.
        """
        if isinstance(vertex, int):
            index = vertex
        elif isinstance(vertex, str):
            matched = self.vertices.get_indices_and_vertices(vertex)
            if len(matched) >= 2:
                raise ValueError(f"More than one equilibrium_points "
                                 f"matched the name {vertex}.")
            if matched:
                index = matched[0][0]
        elif isinstance(vertex, Vertex):
            index = self.vertices.index(vertex)
        else:
            raise TypeError("Input vertex must be int or str or Vertex.")

        if index is None:
            raise ValueError(
                f"{vertex} did not found in equilibrium_points.")

        return CompoundsList([self.stable_compounds[i]
                              for i in self.vertex_to_compounds_list[index]])

    def get_neighbor_vertices_as_dict(self,
                                      remarked_compound: Union[int,
                                                               str,
                                                               Compound],
                                      elements,
                                      **kwargs) -> Dict[str, Any]:
        d = {"compound": remarked_compound}
        # standard_energy = {str(el): float(en) for el, en
        #                    in zip(self.elements, self.element_energy)}
        standard_energy = self.element_energy
        d["standard_energy"] = {k: v for k, v in standard_energy.items()}
        vertices_list =\
            VerticesList(self.get_neighbor_vertices(remarked_compound))
        if self.dim == 3:
            vertices_list.sorted_to_loop_in_3d(elements)
        vertices_list.set_alphabetical_label()
        for i, v in enumerate(vertices_list):
            if not isinstance(v, VertexOnBoundary):
                if v.label is not None:
                    d[v.label] = v.coords
                else:
                    d[f"vertex {i}"] = v.coords

        for k, v in kwargs.items():
            d[k] = v
        return d

    def dump_vertices_yaml(self,
                           file_path: str,
                           remarked_compound: str,
                           elements: ElemOrderType,
                           **kwargs) -> None:
        """Dumps coordination of vertex, compound, and standard_energy.

        Labels of vertices must be one capital alphabet.

        Args:
            file_path(str):
            remarked_compound(str):
            elements([Element or str]): Order of elements
            **kwargs(dict):
                other property, like {"supercell_comment" : "foo,var"}

        Returns:
            None
        """
        d = self.get_neighbor_vertices_as_dict(remarked_compound,
                                               elements,
                                               **kwargs)

        def elem_key_to_str(dict_):
            if isinstance(dict_, dict):
                d_temp = {}
                for k, v in dict_.items():
                    k_ = str(k) if isinstance(k, Element) else k
                    if isinstance(v, np.ndarray):
                        v_ = list(float(val) for val in v)
                    elif isinstance(v, np.float):
                        v_ = float(v)
                    elif isinstance(v, dict):
                        v_ = elem_key_to_str(v)
                    else:
                        v_ = v
                    d_temp[k_] = v_
                return d_temp
            else:
                return dict_

        d = elem_key_to_str(d)
        d["temperature"] = self.temperature
        d["pressure"] = self.pressure
        f = f"{file_path}/vertices.yaml" \
            if os.path.isdir(file_path) else file_path
        with open(f, 'w') as fw:
            ruamel.yaml.dump(d, fw)

    @staticmethod
    def load_vertices_yaml(file_name: str
                           ) -> Tuple[VerticesList, Dict[Element, float]]:
        """Read yaml file return VerticesList and standard_energy.

        Labels of vertices must be one capital alphabet.

        Args:
            file_name(str):

        Returns:
            tuple of (VerticesList, standard_energy)
        """
        with open(file_name, 'r') as fr:
            yaml_data = ruamel.yaml.safe_load(fr)
        vl = VerticesList()
        name_list = list(string.ascii_uppercase)
        for name in name_list:
            if name in yaml_data.keys():
                if not isinstance(yaml_data[name], dict):
                    raise TypeError(f"Failed to read yaml_data[{name}] as dict")
                vertex_dict = OrderedDict(yaml_data[name])
                vertex = Vertex(name, vertex_dict)
                vl.append(vertex)
        return \
            vl, {Element(e): v for e, v in yaml_data["standard_energy"].items()}

    def draw_diagram(self,
                     title: str = None,
                     save_file_name: Optional[str] = None,
                     remarked_compound: Optional[str] = None,
                     with_label: bool = True,
                     draw_range: Optional[float] = None,
                     elements: Optional[ElemOrderType] = None):
        """ Draw chemical potential diagram.

        Args:
            title (str):
                Title of diagram.
            save_file_name (None/str):
                If you will save diagram as image file, specify name of the file
                by this arg.
            remarked_compound (None/str):
                Vertices on the specified compound will be labeled.
            with_label (bool):
                Whether if you will display names of compounds.
            draw_range (None/float):
                If none, range will be determined automatically.
            elements (None/[Element or str]):
                Order of elements.
        """
        elements = elements or sorted([str(e) for e in self.elements])

        if self.dim >= 4:
            # TODO: fix some parameter and plot in 3D or 2D?
            raise NotImplementedError("4D or more data can not be drawn.")

        if self.dim != 1:
            if not draw_range:
                draw_range = self.vertices.boundary_range_limit * 1.1
            self.vertices.set_boundary_range(draw_range)

        #  1D, 2D, and 3D dimension. More than 4D has not yet implemented.
        if self.dim == 1:
            ax = plt.figure().add_subplot(111)
            x = np.zeros(len(self.stable_compounds)
                         + len(self.unstable_compounds))
            y = [cd.energy
                 for cd in self.stable_compounds + self.unstable_compounds]
            ax.scatter(x, y)
            ax.set_ylabel(f"Chemical potential of {elements[0]}")
            ax.set_xlim(-1, 1)
            y_max = np.max(y)
            ax.set_ylim(-y_max * 0.1, y_max * 1.2)
            for cd in self.all_compounds:
                x_shift = -0.2
                ax.text(x_shift,
                        cd.energy,
                        cd.name,
                        size='smaller')
            if title:
                ax.set_title(title)
            else:
                ax.set_title(f"Chemical potential diagram of {elements[0]}")

        elif self.dim == 2:
            ax = plt.figure().add_subplot(111)
            num_line = len(self.stable_compounds)
            for i, compound in enumerate(self.stable_compounds):
                vertices \
                    = VerticesList([self.vertices[j]
                                    for j in self.compounds_to_vertex_list[i]])
                vertices_coords = [v.coords_vector(elements) for v in vertices]
                x = [v[0] for v in vertices_coords]
                y = [v[1] for v in vertices_coords]
                color = "black"
                plt.plot(x, y, color=color)
                mean = np.mean(vertices_coords, axis=0)
                x_shift = 0
                y_shift = 0
                if with_label:
                    ax.text(mean[0] + x_shift,
                            mean[1] + y_shift,
                            compound.name,
                            size='smaller',
                            zorder=num_line + i)
                    if compound.name is None:
                        compound_name = \
                            compound.composition.reduced_formula
                    else:
                        compound_name = compound.name

                    remarked_compound_name = None
                    if remarked_compound:
                        remarked_compound_name = \
                            Composition(remarked_compound).reduced_formula

                    if compound_name == remarked_compound_name:
                        vertices.set_alphabetical_label()
                        for j, v in enumerate(vertices):
                            if v.label:
                                ax.text(v.coords_vector(elements)[0] + x_shift,
                                        v.coords_vector(elements)[1] + y_shift,
                                        v.label,
                                        size="smaller",
                                        zorder=2 * num_line + j,
                                        color="red",
                                        weight="bold")
            if title:
                ax.set_title(title)
            else:
                ax.set_title(f"Chemical potential diagram of "
                             f"{elements[0]} and {elements[1]}")
            ax.set_xlabel(f"Chemical potential of {elements[0]}")
            ax.set_ylabel(f"Chemical potential of {elements[1]}")
            ax.set_xlim(draw_range, 0)
            ax.set_ylim(draw_range, 0)

        elif self.dim == 3:
            ax = plt.figure().add_subplot(111, projection='3d')
            num_plane = len(self.stable_compounds)
            for i, compound in enumerate(self.stable_compounds):
                sorted_vertices = self.get_neighbor_vertices(i, elements)
                sorted_vertices_coords = \
                    [v.coords_vector(elements) for v in sorted_vertices]
                face = Poly3DCollection([sorted_vertices_coords])
                color = [(c + 3) / 4 for c in compound.composition_vector(elements)]
                face.set_color(color)
                face.set_edgecolor("black")
                ax.add_collection3d(face)
                mean = np.mean(sorted_vertices_coords, axis=0)
                if with_label:
                    ax.text(mean[0], mean[1], mean[2],
                            compound.name,
                            size='smaller',
                            zorder=num_plane+i,
                            ha='center',
                            va='center')
                    try:
                        compound_name = \
                            Composition(compound.name).reduced_formula
                    except CompositionError:
                        compound_name = compound.name

                    remarked_compound_name = None
                    if remarked_compound:
                        remarked_compound_name = \
                            Composition(remarked_compound).reduced_formula

                    if compound_name == remarked_compound_name:
                        sorted_vertices.set_alphabetical_label()
                        # vertices_coords.set_alphabetical_label()
                        for j, v in enumerate(sorted_vertices):
                            if v.label:
                                ax.text(v.coords_vector(elements)[0],
                                        v.coords_vector(elements)[1],
                                        v.coords_vector(elements)[2],
                                        v.label,
                                        size="smaller",
                                        zorder=2 * num_plane + j,
                                        weight="bold",
                                        color="red")
            if title:
                ax.set_title(title)
            else:
                ax.set_title(f'Chemical potential diagram of {elements[0]}, '
                             f'{elements[1]}, and {elements[2]}')
            ax.set_xlabel(f'Chemical potential of {elements[0]} (eV)')
            ax.set_ylabel(f'Chemical potential of {elements[1]} (eV)')
            ax.set_zlabel(f'Chemical potential of {elements[2]} (eV)')
            ax.set_xlim3d(draw_range, 0)
            ax.set_ylim3d(0, draw_range)
            ax.set_zlim3d(draw_range, 0)

        ax.grid(color='b', alpha=0.2, linestyle='dashed', linewidth=0.5)
        plt.tight_layout()

        if save_file_name:
            plt.savefig(save_file_name)
        else:
            plt.show()


