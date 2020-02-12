# -*- coding: utf-8 -*-
#  Copyright (c) Oba-group
#  Distributed under the terms of the MIT License.

import json
import string
from typing import Optional, List, Dict, Tuple, Union, Any

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from monty.json import MSONable, MontyEncoder
from monty.serialization import loadfn

import numpy as np

from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.util.string import latexify


# TODO: Add CompoundChemPotDiag


class ChemPotDiag(MSONable):
    """chemical potential diagrams (CPD)"""
    def __init__(self,
                 elements: List[str],
                 el_refs: Dict[str, PDEntry],
                 dim: int,
                 vertices: List[Dict[str, float]],
                 comp_facets: Dict[str, List[int]],
                 target_composition: str,
                 target_comp_chempot: Dict[str, Dict[str, float]]):
        """
        Args:
            elements (List[str]):
                List of pmg Element involved.
            el_refs (Dict[str, PDEntry]):
                List of PDEntry for each Element.
            dim (int):
                Dimension of chemical potential diagram, or number of elements.
            vertices (List[Dict[str, float]]):
                Vertices in the CPD surrounded by compounds.
            comp_facets (Dict[Composition, List[int]] ):
               List of vertex indices comprising each Composition.
            target_composition (Composition):
                Target Composition.
            target_comp_chempot (Dict[str, Dict[str, float]]):
                Chemical potentials at the vertices of the target Composition.
        """
        self.elements = elements
        self.el_refs = el_refs
        self.dim = dim
        self.vertices = vertices
        self.comp_facets = comp_facets
        self.target_composition = target_composition
        self.target_comp_chempot = target_comp_chempot

    @classmethod
    def from_phase_diagram(cls,
                           pd: PhaseDiagram,
                           target_composition: str = None
                           ) -> "ChemPotDiag":
        """

        Args:
            pd (PhaseDiagram):
                Pmg PhaseDiagram object.
            target_composition (Composition):
                Target composition.
        """
        elements = sorted(pd.elements)
        vertices = []
        adjoin_comps = []

        for f in pd.facets:
            comp_list = [pd.qhull_entries[i].composition for i in f]
            energy_list = \
                [pd.get_form_energy_per_atom(pd.qhull_entries[i]) for i in f]
            atomic_frac = [[c.get_atomic_fraction(e) for e in elements]
                           for c in comp_list]
            chempot = dict(zip(elements,
                               np.linalg.solve(atomic_frac, energy_list)))

            vertices.append({str(e): round(chempot[e], 5) for e in elements})
            adjoin_comps.append([pd.qhull_entries[i].composition for i in f])

        comp_facets = {}
        for e in pd.qhull_entries:
            comp_facets[str(e.composition.reduced_composition)] = \
                [i for i, cs in enumerate(adjoin_comps) if e.composition in cs]

        target_composition = str(Composition(target_composition))
        target_comp_chempot = {}
        target_comp_abs_chempot = {}

        if target_composition:
            alphabet_list = list(string.ascii_uppercase)
            for i, f in enumerate(comp_facets[target_composition]):
                target_comp_chempot[alphabet_list[i]] = vertices[f]

        return cls(elements=[str(e) for e in elements],
                   el_refs={str(e): v for e, v in pd.el_refs.items()},
                   dim=len(elements),
                   vertices=vertices,
                   comp_facets=comp_facets,
                   target_composition=target_composition,
                   target_comp_chempot=target_comp_chempot)

    @property
    def target_comp_abs_chempot(self) -> Dict[str, Dict[str, float]]:
        """Same as target_comp_chempot but in absolute scale. """
        abs_chempot = {}
        for k, v in self.target_comp_chempot.items():
            abs_chempot[k] = {e: self.el_refs[e].energy_per_atom + v[e]
                              for e in self.elements}

        return abs_chempot

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
            raise NotImplementedError(f"Drawing {self.dim}-dim is impossible.")

        if title:
            ax.set_title(title)
        else:
            ele = "-".join([e for e in elements])
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
            comp = Composition(comp)
            coords = np.array([[self.vertices[i][e]
                                for i in vertex_indices] for e in elements])
            if len(vertex_indices) == 1:
                if str(comp.elements[0]) == elements[0]:
                    coords = np.append(coords, [[0.0], [draw_range]], axis=1)
                else:
                    coords = np.append(coords, [[draw_range], [0.0]], axis=1)

            color = "black"
            if self.target_composition and str(comp) == self.target_composition:
                color = "blue"

            plt.plot(coords[0], coords[1], zorder=1, color=color, linewidth=3)
            label = latexify(comp.get_reduced_formula_and_factor()[0])
            ax.text(np.average(coords[0]), np.average(coords[1]), label,
                    size='smaller', zorder=2, color=color)

        for label, cp in self.target_comp_chempot.items():
            ax.text(cp[elements[0]], cp[elements[1]], label,
                    zorder=200, color="blue", weight="bold")

        return ax

    def _plot_3d(self, draw_range, elements):
        ax = plt.figure().add_subplot(111, projection='3d')

        for i, (comp, vertex_indices) in enumerate(self.comp_facets.items()):
            comp = Composition(comp)
            coords = np.array([[self.vertices[i][e]
                                for i in vertex_indices] for e in elements])

            comp_elements = [str(e) for e in comp.elements]
            if len(comp_elements) == 1:
                inds = []
                for ind, e in enumerate(elements):
                    if e not in comp_elements:
                        inds.append(ind)

                c1 = coords.copy()
                c2 = coords.copy()
                c1[inds[0], :] = draw_range
                c2[inds[1], :] = draw_range
                c3 = np.array([[0.0], [0.0], [0.0]])
                c3[inds[0], 0] = draw_range
                c3[inds[1], 0] = draw_range
                coords = np.concatenate((coords, c1, c2, c3), axis=1)

            elif len(comp_elements) == 2:
                for ind, e in enumerate(elements):
                    if e not in comp_elements:
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
            if self.target_composition and str(comp) == self.target_composition:
                color = "blue"

            ax.text(center[0], center[1], center[2], label, size='smaller',
                    zorder=100, ha='center', va='center', color=color)

        for label, cp in self.target_comp_chempot.items():
            ax.text(cp[elements[0]], cp[elements[1]], cp[elements[2]], label,
                    zorder=200, color="blue", weight="bold")

        return ax

    @classmethod
    def load_json(cls, filename: str = "cpd.json"):
        return loadfn(filename)

    def to_json_file(self, filename: str = "cpd.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)


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
