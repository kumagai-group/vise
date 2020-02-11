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

from vise.chempotdiag.compound import (
    Compound, DummyCompoundForDiagram, CompoundsList, ElemOrderType)
from vise.chempotdiag.vertex import (
    Vertex, VertexOnBoundary, VerticesList)


class ChemPotDiag:
    """chemical potential diagrams (CPD).

   .. attribute: vertices:

       Vertices in CPD surrounded by compounds.

   .. attribute: remarked_vertices:

       Vertices in CPD equilibrating with a given composition.
       When the given composition is unstable, the nearest chemical potential
       is set.

   .. attribute: compound_instability

        Show how the given compound is unstable, set zero if it is stable.

   .. attribute: boundary_vertices:

       Vertices representing a border in CPD.

   .. attribute: facets:

       Facets representing compounds.


    """

    def __init__(self,
                 pd: PhaseDiagram,
                 composition: Composition = None,
                 equilibrating_phases: List[Composition] = None,
                 ):
        """

        Args:
            phase_diagram:
            composition:
        """
        self.elements = pd.elements
        self.el_refs = pd.el_refs
        self.dim = len(self.elements)

        self.vertices = []
        adjoin_comps = []

        self.absolute_chempot = {}
        for f in pd.facets:
            complist = [pd.qhull_entries[i].composition for i in f]
            energylist = \
                [pd.get_form_energy_per_atom(pd.qhull_entries[i]) for i in f]
            m = [[c.get_atomic_fraction(e) for e in self.elements] for c in
                 complist]
            chempot = dict(zip(self.elements, np.linalg.solve(m, energylist)))

            self.vertices.append({e: chempot[e] for e in self.elements})
            adjoin_comps.append([pd.qhull_entries[i].composition for i in f])

        self.comp_facets = {}
        for e in pd.qhull_entries:
            self.comp_facets[e.composition] = \
                [i for i, cs in enumerate(adjoin_comps) if e.composition in cs]

    def absolute_chem_pot(self, vertex: int):
        return {e: self.el_refs[e].energy_per_atom + self.vertices[vertex][e]
                for e in self.elements}

#        self.all_facets = None
#        self.remarked_vertices = None
#        self.compound_instability = None
#        self.facets = None

    def draw_diagram(self,
                     title: str = None,
                     save_file_name: Optional[str] = None,
                     remarked_comp: Optional[str] = None,
                     draw_range: Optional[float] = None,
                     elements: List[Element] = None):
        """ Draw chemical potential diagram.

        Args:
            title (str):
                Title of diagram.
            save_file_name (None/str):
                If you will save diagram as image file, specify name of the file
                by this arg.
            remarked_comp (None/str):
                Vertices on the specified compound will be labeled.
            draw_range (None/float):
                If none, range will be determined automatically.
            elements (None/[Element or str]):
                Order of elements.
        """
        elements = elements or sorted(self.elements)

        draw_range = draw_range or min(sum([list(i.values())
                                            for i in self.vertices], [])) * 1.1

        #  1D, 2D, and 3D dimension. More than 4D has not yet implemented.
        if self.dim == 2:
            ax = self._plot_2d(draw_range, elements, remarked_comp)
        # elif self.dim == 3:
        #     ax = self._plot_3d(draw_range, elements, remarked_compound, with_label)
        else:
            raise NotImplementedError(f"Drawing {self.dim} is impossible.")

        if title:
            ax.set_title(title)
        else:
            ele = "-".join([str(e) for e in self.elements])
            ax.set_title(f"Chemical potential diagram {ele}")
        ax.set_xlabel(f"Chemical potential of {elements[0]}")
        ax.set_ylabel(f"Chemical potential of {elements[1]}")
        if self.dim == 2:
            ax.set_xlim(draw_range, 0)
            ax.set_ylim(draw_range, 0)
            plt.gca().set_aspect('equal', adjustable='box')
        else:
            ax.set_zlabel(f"Chemical potential of {elements[0]}")
            ax.set_xlim3d(draw_range, 0)
            ax.set_ylim3d(0, draw_range)
            ax.set_zlim3d(draw_range, 0)

        ax.grid(color='b', alpha=0.2, linestyle='dashed', linewidth=0.5)
        plt.tight_layout()

        if save_file_name:
            plt.savefig(save_file_name)
        else:
            plt.show()

    def _plot_2d(self, draw_range, elements, remarked_compound):
        ax = plt.figure().add_subplot(111)

        for i, (comp, vertex_indices) in enumerate(self.comp_facets.items()):
            x = [self.vertices[i][elements[0]] for i in vertex_indices]
            y = [self.vertices[i][elements[1]] for i in vertex_indices]

            if len(vertex_indices) == 1:
                if comp.elements[0] == elements[0]:
                    y.append(draw_range)
                    x.append(0.0)
                else:
                    y.append(0.0)
                    x.append(draw_range)

            if remarked_compound and comp == remarked_compound:
                color = "blue"
            else:
                color = "black"

            plt.plot(x, y, zorder=1, color=color)
            label = latexify(comp.get_reduced_formula_and_factor()[0])
            ax.text(sum(x) / 2, sum(y) / 2, label, size='smaller', zorder=2,
                    color=color)

        return ax

    def _plot_3d(self, draw_range, elements, remarked_compound):
        ax = plt.figure().add_subplot(111, projection='3d')

        for i, (comp, vertex_indices) in enumerate(self.comp_facets.items()):
            coords = np.array([[self.vertices[i][elements[x]]
                                for i in vertex_indices] for x in range(3)])

            if len(comp.elements) == 1:
                idx = elements.index(comp.elements[0])
                coords[idx].append([0.0, 0.0, 0.0])
                coords[idx-1].append([draw_range, draw_range, 0.0])
                coords[idx-2].append([draw_range, 0.0, draw_range])

            elif len(comp.elements) == 2:
                for i in range(3):
                    if elements[i] not in comp.elements:
                        idx = i
                        continue
                c = coords.copy()
                c[idx, :] = draw_range
                np.append(coords, c, axis=1)

            sorted_coords = sort_coords(coords)
            face = Poly3DCollection([sorted_coords])
#            color = [(c + 3) / 4 for c in compound.composition_vector(elements)]
#            face.set_color(color)
            face.set_edgecolor("black")
            ax.add_collection3d(face)
            # mean = np.mean(sorted_vertices_coords, axis=0)
            # ax.text(mean[0], mean[1], mean[2],
            #         compound.name,
            #         size='smaller',
            #         zorder=1,
            #         ha='center',
            #         va='center')
            # try:
            #     compound_name = \
            #         Composition(compound.name).reduced_formula
            # except CompositionError:
            #     compound_name = compound.name

            # remarked_compound_name = None
            # if remarked_compound:
            #     remarked_compound_name = \
            #         Composition(remarked_compound).reduced_formula

            # if compound_name == remarked_compound_name:
            #     sorted_vertices.set_alphabetical_label()
            #     # vertices_coords.set_alphabetical_label()
            #     for j, v in enumerate(sorted_vertices):
            #         if v.label:
            #             ax.text(v.coords_vector(elements)[0],
            #                     v.coords_vector(elements)[1],
            #                     v.coords_vector(elements)[2],
            #                     v.label,
            #                     size="smaller",
            #                     zorder=2,
            #                     weight="bold",
            #                     color="red")
        return ax


def sort_coords(coords):
    """
    Args:
        sort indices of vertex to loop (for drawing diagram)
    """
    if len(coords) != 3:
        raise ValueError("Only valid for 3D vector")

    t_coords = np.transpose(coords)
    mean = np.average(coords)
    rel_coords = t_coords - mean
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
        v0 = rel_coords[0]
        v = rel_coords[index]
        v0 = v0 / np.linalg.norm(v0)
        v = v / np.linalg.norm(v)
        if np.linalg.norm(v - v0) < 1e-10:
            return 0
        dot = np.dot(v0, v)
        det = v0[0]*v[1]*n[2] + v[0]*n[1]*v0[2] + n[0]*v0[1]*v[2] - \
            v0[2]*v[1]*n[0] - v[2]*n[1]*v0[0] - n[2]*v0[1]*v[0]
        angle = np.arctan2(det, dot)
        # For easy debug
        angle = np.rad2deg(angle)
        return angle

    indices = [i for i in range(len(vertices_coords))]
    indices.sort(key=angle_between_v0)
    return VerticesList([self[i] for i in indices])