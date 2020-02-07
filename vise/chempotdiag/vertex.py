#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.
from typing import Optional, List, Dict, Sequence, Union, Set
from collections import OrderedDict, Iterable
import string

import numpy as np

from pymatgen.core.periodic_table import Element

from chempotdiag.compound import Compound, CompoundsList, ElemOrderType


class Vertex:
    """
        Object for vertex in chemical potential diagram.
    """

    def __init__(self,
                 label: Optional[str],
                 coords: Dict[Union[str, Element], float]):
        """
        Create a Vertex object.
        Args:
            label (None/str): Label of vertex. This is used to draw a diagram.
            coords ({Element or str: float}):
                Coordinates of vertex in chemical potential diagram.
        """
        # for 1 unary system, coords may be float
        if not isinstance(coords, Iterable):
            coords = [coords]
        self.label = label
        self._coords = {Element(k): v for k, v in coords.items()}

    @property
    def elements(self) -> Set[Element]:
        """
            (list of Elements) Elemental elements of the composition vector.
        """
        return set(self.coords.keys())

    @property
    def coords(self) -> Dict[Element, float]:
        """
            ({Element or str: float) Coordinates of vertex.
        """
        return self._coords

    def coords_vector(self, elements: ElemOrderType) -> np.ndarray:
        return np.array([self.coords[Element(e)] for e in elements])

    def as_dict(self) -> Dict:
        d = {"label": self.label,
             "coords": self.coords}
        return d

    @classmethod
    def from_dict(cls, d: Dict) -> "Vertex":
        label = d["label"]
        coords = d["coords"]
        return cls(label, coords)

    def __repr__(self):
        return (f"Vertex("
                f"Label: {self.label} , "
                f"Coordinates: {self.coords})")

    def __eq__(self, other: "Vertex"):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self == other.
        """
        if not isinstance(other, Vertex):
            raise TypeError("Vertex class can not be"
                            " compared with other class.")
        if any([e1 != e2 for e1, e2 in zip(self.elements, other.elements)]):
            return False
        elif any([self.coords[e] != other.coords[e]
                  for e in self.elements]):
            return False
        return True

    def __ne__(self, other: "Vertex"):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self != other.
        """
        return not self == other

    def __lt__(self, other: "Vertex"):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self < other.
        """
        if not isinstance(other, Vertex):
            raise TypeError("Vertex class can not be"
                            " compared with other class.")
        if self.elements != other.elements:
            raise ValueError("Can't be compared because elements differ")
        for e in sorted(self.elements):
            c1 = self.coords[e]
            c2 = other.coords[e]
            if not np.isclose(c1, c2, atol=0, rtol=1e-12):
                return c1 < c2
        else:  # all elements are strictly the same.
            return False

    def __gt__(self, other: "Vertex"):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self > other.
        """
        return (not self < other) and self != other

    def __ge__(self, other: "Vertex"):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self >= other.
        """
        return not self < other

    def __le__(self, other: "Vertex"):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
        Returns (bool): If self > other.
        """
        return not self > other

    def almost_equal(self, other: "Vertex", tol: float = 1e-5):
        """
        Check if two equilibrium_points is almost same.
        Information is label is ignored.
        Args:
            other (Vertex): Compared vertex.
            tol (float): Tolerance for numeric error of coordinates.
        Returns (bool): If self almost equals to other.
        """
        if not isinstance(other, Vertex):
            raise TypeError(f"Vertex class can not be "
                            f"compared with {other} type.")
        if self.elements != other.elements:
            return False
        for e in sorted(self.elements):
            c1 = self.coords[e]
            c2 = other.coords[e]
            if abs(c1 - c2) > tol:
                return False
        return True


class VertexOnBoundary(Vertex):

    def __init__(self,
                 label: Optional[str],
                 coords: Dict[Union[str, Element], float],
                 boundary_range: float,
                 boundary_range_limit: float):
        super().__init__(label, coords)
        self._boundary_range = boundary_range
        self._boundary_range_limit = boundary_range_limit
        self._element_of_boundary_coords =\
            {e: np.abs(c - boundary_range) < 1e-6 for e, c in coords.items()}

    def set_boundary_range(self, boundary_range: float):
        if boundary_range > self._boundary_range_limit:
            raise ValueError("Boundary range is shallower than deepest vertex.")
        for element in self._element_of_boundary_coords:
            if self._element_of_boundary_coords[element]:
                self.coords[Element(element)] = boundary_range

    @property
    def boundary_range_limit(self):
        return self._boundary_range_limit


class VerticesList(list):

    _TYPE_ERROR_MESSAGE = "VerticesArray must be contains only vertex."
    _DIF_ELEM_ERROR_MESSAGE = "Vertices must have the same elements."

    def __init__(self, *args, **kw):
        list.__init__(self, *args, **kw)
        for v in self:
            if not isinstance(v, Vertex):
                raise TypeError(self._TYPE_ERROR_MESSAGE)
            if len(self) != 0 and self[0].elements != v.elements:
                raise ValueError(self._DIF_ELEM_ERROR_MESSAGE)

    def append(self, vertex: Vertex):
        if not isinstance(vertex, Vertex):
            raise TypeError(self._TYPE_ERROR_MESSAGE)
        if len(self) != 0 and vertex.elements != self.elements:
            raise ValueError(self._DIF_ELEM_ERROR_MESSAGE)
        list.append(self, vertex)

    def extend(self, vertices: "VerticesList"):
        if not isinstance(vertices, list) or \
                any([not isinstance(v, Vertex) for v in vertices]):
            raise TypeError(self._TYPE_ERROR_MESSAGE)
        elif len(self) != 0 and len(vertices) != 0 and\
                any(v.elements != self.elements for v in vertices):
            raise ValueError(self._DIF_ELEM_ERROR_MESSAGE)
        else:
            list.extend(self, vertices)

    def __add__(self, vertices: "VerticesList"):
        if not isinstance(vertices, list) or \
                any([not isinstance(v, Vertex) for v in vertices]):
            raise TypeError(self._TYPE_ERROR_MESSAGE)
        elif len(self) != 0 and len(vertices) != 0 and \
                any(v.elements != self.elements for v in vertices):
            raise ValueError(self._DIF_ELEM_ERROR_MESSAGE)
        else:
            return VerticesList(list.__add__(self, vertices))

    def __setitem__(self, key: int, vertex: Vertex):
        if not isinstance(vertex, Vertex):
            raise TypeError(self._TYPE_ERROR_MESSAGE)
        elif len(self) != 0 and self.elements != vertex.elements:
            raise ValueError(self._DIF_ELEM_ERROR_MESSAGE)
        list.__setitem__(self, key, vertex)

    @property
    def elements(self) -> List[Element]:
        return self[0].elements

    def coords_vector(self, elements: ElemOrderType) -> np.ndarray:
        return np.array([v.coords_vector(elements) for v in self])

    def sorted_to_loop_in_3d(self, elements: Optional[ElemOrderType]):
        """
        Args:
            sort indices of vertex to loop (for drawing diagram)
        """
        elements = elements if elements \
            else sorted([str(e) for e in self.elements])
        if any([len(v.coords) != 3 for v in self]):
            raise ValueError("This function can be applied to 3D vector,"
                             "But input vectors contain non-3D vector.")
        vertices_coords = [v.coords_vector(elements) for v in self]
        mean = np.zeros(len(self[0].coords))
        for v in vertices_coords:
            mean += v / len(self)
        from_mean = [v - mean for v in vertices_coords]
        n = np.cross(from_mean[0], from_mean[1])
        if abs(np.linalg.norm(n)) < 1e-8:  # unfortunately, get parallel vectors
            n = np.cross(from_mean[0], from_mean[2])
        n = n / np.linalg.norm(n)

        def angle_between_v0(index: int) -> float:
            """
            Args:
                index (int): index of vertices_coords.
            Returns (float):
                angle between from_mean[index] and from_mean[0]
            """
            v0 = from_mean[0]
            v = from_mean[index]
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

    def clear_label(self):
        """
            Clear labels of all equilibrium_points.
        """
        for i, _ in enumerate(self):
            self[i].label = None

    @property
    def is_labeled(self) -> bool:
        """
            (bool) Is equilibrium_points are labeled.
        """
        return any([v.label for v in self])

    def set_alphabetical_label(self):
        """
            Label equilibrium_points alphabetically.
        """
        name_list = list(string.ascii_uppercase)
        count = 0
        for i in range(len(self)):
            if isinstance(self[i], VertexOnBoundary):
                continue
            self[i].label = name_list[count]
            count += 1

    def get_indices_and_vertices(self, label: Optional["VerticesList"]):
        """
        Find object of Compound from self by name(str) of compound.
        Args:
            label (str):
        Returns (None/VerticesList):
            Matched vertex data from label. If no compounds match, return None.
        """
        matched_list = [(i, v) for i, v in enumerate(self)
                        if v.name == label]
        if len(matched_list) == 0:
            return None
        return VerticesList(matched_list)

    def set_boundary_range(self, boundary_range: float):
        """
            Change coordinates of equilibrium_points on boundary and range of drawing.
        """
        for v in self:
            if isinstance(v, VertexOnBoundary):
                v.set_boundary_range(boundary_range)

    @property
    def boundary_range_limit(self) -> float:
        limit = None
        for v in self:
            if isinstance(v, VertexOnBoundary):
                limit = v.boundary_range_limit
        if limit is None:
            raise ValueError("All equilibrium_points are not on boundary."
                             "It is unexpected error."
                             "Please contact author.")
        return limit

    # TODO: document, unittest
    def almost_equal(self, other: "VerticesList", tol: float = 1e-5):
        """
        Check if this object and other object almost equal.
        If only order of elements like ([Ca, O] and [O, Ca]),
        judges they are same.
        Args:
            other (list of Compound or CompoundsList) : Compared object.
            tol (float) : Numerical tolerance. (Of energy, composition)
        Returns (bool) : If this object and other object almost equal.

        """
        if len(self) != len(other):
            return False
        else:
            for v1, v2 in zip(sorted(self), sorted(other)):
                if not v1.almost_equal(v2, tol=tol):
                    return False
        return True


