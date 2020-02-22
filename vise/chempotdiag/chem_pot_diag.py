# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import json
import string
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from typing import Optional, List, Dict

from monty.json import MSONable, MontyEncoder
from monty.serialization import loadfn

from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.util.string import latexify

from vise.util.logger import get_logger

logger = get_logger(__name__)


class ChemPotDiag(MSONable):
    """chemical potential diagrams"""
    def __init__(self,
                 elements: List[Element],
                 el_ref_list: List[PDEntry],
                 dim: int,
                 vertices: List[List[float]],
                 qhull_entries: List[PDEntry],
                 comp_facets: List[List[int]],
                 target_comp: str,
                 target_comp_chempot: Dict[str, List[float]],
                 unstable_energy: Optional[float] = None):
        """

        Note: Chemical potential is in eV/*atom* even for the complex diagrams.

        Args:
            elements (List[Element]):
                List of pmg Element involved.
            el_ref_list (List[PDEntry]):
                List of PDEntry for each Element.
            dim (int):
                Dimension of chemical potential diagram, or number of elements.
            vertices (List[Dict[Element, float]]):
                Vertices in the CPD surrounded by compounds.
            qhull_entries (List[PDEntry]):
                Entries comprising convex hull.
            comp_facets (List[List[int]]):
               List of vertex indices comprising each composition.
               Indices are the same as qhull_entries, therefore composition
               and name are obtained via qhull_entries.
            target_comp (Composition):
                Target Composition.
            target_comp_chempot (Dict[str, Dict[Element, float]]):
                Chemical potentials at the vertices of the target composition.
                The keys are alphabets.
        """
        self.elements = elements
        self.el_ref_list = el_ref_list
        self.dim = dim
        self.vertices = vertices
        self.qhull_entries = qhull_entries
        self.comp_facets = comp_facets
        self.target_comp = target_comp
        self.target_comp_chempot = target_comp_chempot
        self.unstable_energy = unstable_energy

    @classmethod
    def from_phase_diagram(
            cls,
            pd: PhaseDiagram,
            target_comp: Optional[str] = None,
            allow_unstable_target_chempot: bool = False,
            ) -> "ChemPotDiag":
        """Create class instance object from pymatgen PhaseDiagram object

        Args:
            pd (PhaseDiagram):
                Pmg PhaseDiagram object.
            target_comp (None/str):
                Target composition.
            allow_unstable_target_chempot:
                Whether to allow unstable target chemical potential.
                The vertex is defined at the point where the surrounding
                competing phases in the phase diagram equilibrates.
        """
        elements = sorted(pd.elements)
        vertices = []
        adjoin_comps = []  # List of adjoining compositions to each vertex.

        for fc in pd.facets:
            comp_list = [pd.qhull_entries[i].composition for i in fc]
            e_list = [pd.get_form_energy_per_atom(pd.qhull_entries[i])
                      for i in fc]
            atomic_frac = [[c.get_atomic_fraction(e) for e in elements]
                           for c in comp_list]
            chempot = dict(zip(elements, np.linalg.solve(atomic_frac, e_list)))

            vertices.append([round(chempot[e], 10) for e in elements])
            adjoin_comps.append([pd.qhull_entries[i].composition for i in fc])

        comp_facets = []  # List of vertex indices comprising each composition.
        for e in pd.qhull_entries:
            comp_facets.append(
                [i for i, cs in enumerate(adjoin_comps) if e.composition in cs])

        target_comp_chempot = {}
        unstable_energy = None
        if target_comp:
            target_inds = [i for i, e in enumerate(pd.qhull_entries)
                           if Composition(e.name) == Composition(target_comp)]

            alphabet_list = list(string.ascii_uppercase)
            if not target_inds:
                if not allow_unstable_target_chempot:
                    raise ValueError(f"Target composition {target_comp} "
                                     f"is invalid. Choose from "
                                     f"{[e.name for e in pd.qhull_entries]}")
                else:
                    target_inds = \
                        [e for e in pd.unstable_entries
                         if Composition(e.name) == Composition(target_comp)]
                    if not target_inds:
                        raise ValueError(
                            f"Target composition {target_comp} does not "
                            f"exist even in unstable entries. Choose from "
                            f"{[e.name for e in pd.unstable_entries]}")
                    comp = target_inds[0].composition
                    unstable_energy = pd.get_hull_energy(comp)
                    # Find the facet in phase diagram where comp belongs to.
                    c = pd.pd_coords(comp)
                    facets = []
                    for f, s in zip(pd.facets, pd.simplexes):
                        if s.in_simplex(c, PhaseDiagram.numerical_tol / 10):
                            facets.append(f)

                    # Find the index of vertices in chemical potential diagram.
                    for i, fc in enumerate(facets):
                        ind = [v.tolist() for v in pd.facets].index(fc.tolist())
                        target_comp_chempot[alphabet_list[i]] = vertices[ind]
                        logger.warning(
                            f"Unstable compound {comp} is evaluated in "
                            f"ChemPotDiag. Unstable energy: {unstable_energy}.")
            else:
                if len(comp_facets) > 26:
                    raise ValueError(f"Too many vertices: {len(comp_facets)}.")

                for i, fc in enumerate(comp_facets[target_inds[0]]):
                    target_comp_chempot[alphabet_list[i]] = vertices[fc]

        return cls(elements=elements,
                   el_ref_list=[pd.el_refs[e] for e in elements],
                   dim=len(elements),
                   vertices=vertices,
                   qhull_entries=pd.qhull_entries,
                   comp_facets=comp_facets,
                   target_comp=target_comp,
                   target_comp_chempot=target_comp_chempot,
                   unstable_energy=unstable_energy)

    @property
    def target_comp_abs_chempot(self) -> Dict[str, List[float]]:
        """Same as target_comp_chempot but in absolute scale.

        Mostly used for calculating chemical potential-related properties such
        as point defects.
        """
        abs_chempot = {}
        for point, chempot in self.target_comp_chempot.items():
            abs_chempot[point] = [round(e.energy_per_atom + c, 10)
                                  for e, c in zip(self.el_ref_list, chempot)]
        return abs_chempot

    def draw_diagram(self,
                     title: str = None,
                     filename: Optional[str] = None) -> plt:
        """Draw chemical potential diagram.

        Args:
            title (str):
                Title of the diagram.
            filename (None/str):
                Specify this when saving the diagram as an image file.

        Returns:
            Pyplot
        """
        draw_range = min(-1.0, min(sum([i for i in self.vertices], [])) * 1.1)

        if self.dim == 2 or self.dim == 3:
            ax = self._plot(draw_range, title)
        else:
            raise NotImplementedError(f"Drawing {self.dim} is impossible.")

        if self.dim == 2:
            ax.set_xlim(draw_range, 0)
            ax.set_ylim(draw_range, 0)
            plt.gca().set_aspect('equal', adjustable='box')
        else:
            ax.set_xlim3d(draw_range, 0)
            ax.set_ylim3d(0, draw_range)
            ax.set_zlim3d(draw_range, 0)

        ax.grid(color='b', alpha=0.2, linestyle='dashed', linewidth=0.5)
        plt.tight_layout()

        if filename:
            plt.savefig(filename)
        else:
            plt.show()

    def _plot(self, draw_range: float, title: str):

        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['ps.fonttype'] = 42
        plt.rcParams['svg.fonttype'] = "none"

        if self.dim == 2:
            ax = plt.figure().add_subplot(111)
        else:
            ax = plt.figure().add_subplot(111, projection='3d')

        text_kwargs = {}
        species = ["NA"] * self.dim
        for i, vertex_inds in enumerate(self.comp_facets):
            name = self.qhull_entries[i].name
            label = latexify(name)
            coords = [self.vertices[i] for i in vertex_inds]
            if len(coords) == 0:
                continue
            coords = np.array(coords).transpose()
            comp = self.qhull_entries[i].composition

            comp_color = "black"
            if self.target_comp and name == self.target_comp:
                comp_color = "blue"

            if self.dim == 2:
                if len(comp.elements) == 1:
                    ind = 0 if comp.elements[0] == self.elements[0] else 1
                    # [0, draw_range] & [draw_range, 0] for 1st & 2nd elements.
                    coord = [[draw_range * ind], [draw_range * (1 - ind)]]
                    coords = np.append(coords, coord, axis=1)
                    set_label = f"set_{['x', 'y'][ind]}label"
                    getattr(ax, set_label)(f"Chemical potential of {label}")
                    species[ind] = label

                plt.plot(*coords, zorder=1, color=comp_color, linewidth=3)

            else:  # 3-dim
                if len(comp.elements) == 1:
                    ex_inds = []
                    for ind, e in enumerate(self.elements):
                        if e in comp.elements:
                            label_ind = ind
                            set_label = f"set_{['x', 'y', 'z'][ind]}label"
                        else:
                            ex_inds.append(ind)

                    # Project vertices to the two edges and add one corner
                    c1 = coords.copy()
                    c2 = coords.copy()
                    c1[ex_inds[0], :] = c2[ex_inds[1], :] = draw_range
                    c3 = np.array([[0.0], [0.0], [0.0]])
                    c3[ex_inds[0], 0] = c3[ex_inds[1], 0] = draw_range
                    coords = np.concatenate((coords, c1, c2, c3), axis=1)
                    getattr(ax, set_label)(f"Chemical potential of {label}")
                    species[label_ind] = label

                elif len(comp.elements) == 2:
                    for ind, e in enumerate(self.elements):
                        if not any([str(e) == str(e2) for e2 in comp.elements]):
                            break
                    c = coords.copy()
                    c[ind, :] = draw_range
                    coords = np.append(coords, c, axis=1)

                sorted_t_coords = sort_coords(np.transpose(coords))
                face = Poly3DCollection([sorted_t_coords])

                color = [(comp.get_atomic_fraction(e) + 3) / 4
                         for e in self.elements]
                face.set_color(color)
                face.set_edgecolor("black")
                ax.add_collection3d(face)
                text_kwargs.update({"ha": "center", "va": "center"})

            center = np.average(coords, axis=1)
            ax.text(*center, label, size='smaller', zorder=100,
                    color=comp_color, **text_kwargs)

        title = title or f"Chemical potential diagram of {'-'.join(species)}"
        ax.set_title(title)

        for label, cp in self.target_comp_chempot.items():
            ax.text(*cp, label, zorder=200, color="blue", weight="bold")

        return ax

    @classmethod
    def load_json(cls, filename: str = "cpd.json"):
        return loadfn(filename)

    def to_json_file(self, filename: str = "cpd.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)


def sort_coords(coords: np.ndarray) -> np.ndarray:
    """Sort coordinates based on the angle with first coord from the center.

    Args:
        coords (np.ndarray):
            Coordinates to be sorted. The format of coords is as follows.
            np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]

    Returns:
        np.ndarray for sorted coordinates.
    """
    if len(coords[0]) != 3:
        raise ValueError("Only valid for 3D vector")

    center = np.average(coords, axis=0)
    rel_coords = coords - center
    n = np.cross(rel_coords[0], rel_coords[1])
    if abs(np.linalg.norm(n)) < 1e-8:  # Skip parallel vectors.
        n = np.cross(rel_coords[0], rel_coords[2])
    n = n / np.linalg.norm(n)

    def angle_between_v0(index: int) -> float:
        """
        Args:
            index (int): index of coords.

        Returns (float):
            Angle between rays from the center to rel_coords[0] and
            rel_coords[int].
        """
        v0 = rel_coords[0] / np.linalg.norm(rel_coords[0])
        v = rel_coords[index] / np.linalg.norm(rel_coords[index])
        det = np.linalg.det(np.concatenate(([v0], [v], [n]), axis=0))
        angle = np.arctan2(np.clip(np.dot(v0, v), -1.0, 1.0), det)
        return angle

    indices = [i for i in range(len(coords))]
    indices.sort(key=angle_between_v0)

    return coords[indices]
