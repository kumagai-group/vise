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
from pymatgen.util.string import latexify
from pymatgen.analysis.phase_diagram import PDEntry

from vise.chempotdiag.compound import (
    Compound, DummyCompoundForDiagram, CompoundsList, ElemOrderType)
from vise.chempotdiag.vertex import (
    Vertex, VertexOnBoundary, VerticesList)


class ChemPotDiag:
    """chemical potential diagrams (CPD).

   .. attribute: elements (List[Elements]):

        List of pmg Element involved.

   .. attribute: el_refs (Dict[Element, PDEntry]):

        List of PDEntry for each Element.

   .. attribute: dim (int):

        Dimension of chemical potential diagram, i.e., number of elements.

   .. attribute: vertices (List[Dict[Element, float]]):

       Vertices in the CPD surrounded by compounds.

   .. attribute: comp_facets (Dict[Composition, List[int]] ):

        List of vertex indices comprising each Composition.

   .. attribute: target_composition (Composition):

       Target Composition.

   .. attribute: target_comp_chempot:

       Chemical potentials at the vertices of the target Composition.

   .. attribute: target_comp_abs_chempot:

       Absolute chemical potentials at the vertices of the target Composition.
    """

    def __init__(self,
                 pd: PhaseDiagram,
                 target_composition: Composition = None,
                 ):
        """

        Args:
            pd (PhaseDiagram):
                Pmg PhaseDiagram object.
            target_composition (Composition):
                Target composition.
        """
        self.elements: List[Element] = sorted(pd.elements)
        self.el_refs: Dict[Element, PDEntry] = pd.el_refs
        self.dim: int = len(self.elements)

        self.vertices: List[Dict[Element, float]] = []
        adjoin_comps = []

        for f in pd.facets:
            comp_list = [pd.qhull_entries[i].composition for i in f]
            energy_list = \
                [pd.get_form_energy_per_atom(pd.qhull_entries[i]) for i in f]
            atomic_frac = [[c.get_atomic_fraction(e) for e in self.elements]
                           for c in comp_list]
            chempot = dict(zip(self.elements,
                               np.linalg.solve(atomic_frac, energy_list)))

            self.vertices.append(
                {e: round(chempot[e], 5) for e in self.elements})
            adjoin_comps.append([pd.qhull_entries[i].composition for i in f])

        self.comp_facets: Dict[Composition, List[int]] = {}
        for e in pd.qhull_entries:
            self.comp_facets[e.composition] = \
                [i for i, cs in enumerate(adjoin_comps) if e.composition in cs]

        self.target_composition: Composition = target_composition
        self.target_comp_chempot: Dict[str, Dict[Element, float]] = {}
        self.target_comp_abs_chempot: Dict[str, Dict[Element, float]] = {}

        if target_composition:
            alphabet_list = list(string.ascii_uppercase)
            for i, f in enumerate(self.comp_facets[target_composition]):
                self.target_comp_chempot[alphabet_list[i]] = self.vertices[f]
                self.target_comp_abs_chempot[alphabet_list[i]] = \
                    {e: self.el_refs[e].energy_per_atom + self.vertices[f][e]
                     for e in self.elements}

    def draw_diagram(self,
                     title: str = None,
                     filename: Optional[str] = None,
                     draw_range: Optional[float] = None,
                     elements: List[Element] = None):
        """ Draw chemical potential diagram.

        Args:
            title (str):
                Title of the diagram.
            filename (None/str):
                Specify this if one will save diagram as image file.
            draw_range (None/float):
                Lower limit for the diagram.
                If none, range will be determined from the lowest chemical
                potential.
            elements (None/[Element or str]):
                Specify this when changing order of elements in the diagram.
        """
        elements = elements or self.elements
        draw_range = draw_range or min(sum([list(i.values())
                                            for i in self.vertices], [])) * 1.1

        if self.dim == 2:
            ax = self._plot_2d(draw_range, elements)
        elif self.dim == 3:
            ax = self._plot_3d(draw_range, elements)
        else:
            raise NotImplementedError(f"Drawing {self.dim} is impossible.")

        if title:
            ax.set_title(title)
        else:
            ele = "-".join([str(e) for e in elements])
            ax.set_title(f"Chemical potential diagram: {ele}")

        ax.set_xlabel(f"Chemical potential of {elements[0]}")
        ax.set_ylabel(f"Chemical potential of {elements[1]}")
        if self.dim == 2:
            ax.set_xlim(draw_range, 0)
            ax.set_ylim(draw_range, 0)
            plt.gca().set_aspect('equal', adjustable='box')
        else:
            ax.set_zlabel(f"Chemical potential of {elements[2]}")
            ax.set_xlim3d(draw_range, 0)
            ax.set_ylim3d(0, draw_range)
            ax.set_zlim3d(draw_range, 0)

        ax.grid(color='b', alpha=0.2, linestyle='dashed', linewidth=0.5)
        plt.tight_layout()

        if filename:
            plt.savefig(filename)
        else:
            plt.show()

    def _plot_2d(self, draw_range, elements):
        ax = plt.figure().add_subplot(111)

        for i, (comp, vertex_indices) in enumerate(self.comp_facets.items()):
            coords = np.array([[self.vertices[i][e]
                                for i in vertex_indices] for e in elements])
            if len(vertex_indices) == 1:
                if comp.elements[0] == elements[0]:
                    np.append(coords, [[0.0], [draw_range]], axis=1)
                else:
                    np.append(coords, [[draw_range], [0.0]], axis=1)

            color = "black"
            if self.target_composition and comp == self.target_composition:
                color = "blue"

            plt.plot(coords[0], coords[1], zorder=1, color=color)
            label = latexify(comp.get_reduced_formula_and_factor()[0])
            ax.text(np.average(coords[0]), np.average(coords[1]), label,
                    size='smaller', zorder=2, color=color)

        return ax

    def _plot_3d(self, draw_range, elements):
        ax = plt.figure().add_subplot(111, projection='3d')

        for i, (comp, vertex_indices) in enumerate(self.comp_facets.items()):
            coords = np.array([[self.vertices[i][e]
                                for i in vertex_indices] for e in elements])

            if len(comp.elements) == 1:
                inds = []
                for ind, e in enumerate(elements):
                    if e not in comp.elements:
                        inds.append(ind)

                c1 = coords.copy()
                c2 = coords.copy()
                c1[inds[0], :] = draw_range
                c2[inds[1], :] = draw_range
                c3 = np.array([[0.0], [0.0], [0.0]])
                c3[inds[0], 0] = draw_range
                c3[inds[1], 0] = draw_range
                coords = np.concatenate((coords, c1, c2, c3), axis=1)

            elif len(comp.elements) == 2:
                for ind, e in enumerate(elements):
                    if e not in comp.elements:
                        continue
                c = coords.copy()
                c[ind, :] = draw_range
                coords = np.append(coords, c, axis=1)

            sorted_t_coords = sort_coords(np.transpose(coords))
            face = Poly3DCollection([sorted_t_coords])

            color = [(comp.get_atomic_fraction(e) + 3) / 4 for e in elements]
            face.set_color(color)
            face.set_edgecolor("black")
            ax.add_collection3d(face)

            center = np.average(coords, axis=1)
            label = latexify(comp.get_reduced_formula_and_factor()[0])
            color = "black"
            if self.target_composition and comp == self.target_composition:
                color = "blue"

            ax.text(center[0], center[1], center[2], label, size='smaller',
                    zorder=100, ha='center', va='center', color=color)

        for label, cp in self.target_comp_chempot.items():
            ax.text(cp[elements[0]], cp[elements[1]], cp[elements[2]], label,
                    zorder=200, color="blue", weight="bold")

        return ax


def sort_coords(coords: np.ndarray):
    """Sort coordinates based on the angle with first coord from the center.

    Args:
        coords (np.ndarray):
            Coordinates to be sorted
            np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]
        sort indices of vertex to loop (for drawing diagram)

    Returns:
        np.ndarray for sorted coordinates.
    """
    if len(coords[0]) != 3:
        raise ValueError("Only valid for 3D vector")

    center = np.average(coords, axis=0)
    rel_coords = coords - center
    n = np.cross(rel_coords[0], rel_coords[1])
    if abs(np.linalg.norm(n)) < 1e-8:  # unfortunately, get parallel vectors
        n = np.cross(rel_coords[0], rel_coords[2])
    n = n / np.linalg.norm(n)

    def angle_between_v0(index: int) -> float:
        """
        Args:
            index (int): index of vertices_coords.
        Returns (float):
            angle between from_mean[index] and from_mean[0]
        """
        v0 = rel_coords[0] / np.linalg.norm(rel_coords[0])
        v = rel_coords[index] / np.linalg.norm(rel_coords[index])
        det = np.linalg.det(np.concatenate(([v0], [v], [n]), axis=0))
        angle = np.arctan2(np.clip(np.dot(v0, v), -1.0, 1.0), det)
        return angle

    indices = [i for i in range(len(coords))]
    indices.sort(key=angle_between_v0)

    return coords[indices]
