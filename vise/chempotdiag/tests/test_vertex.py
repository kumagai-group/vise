#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.
import unittest

import numpy as np

from pymatgen.core.periodic_table import Element

from chempotdiag.vertex import Vertex, VerticesList, VertexOnBoundary


class TestVertex(unittest.TestCase):

    def test_vertex_dict(self):
        v = Vertex("Label", {"Mg": 0.00000, "O": -7.77345})
        d = v.as_dict()
        v2 = Vertex.from_dict(d)
        self.assertEqual(v, v2)

    def test_vertex_list(self):
        v1 = Vertex("Label", {"Mg": 0.00000})
        v2 = Vertex("Label", {"O": 0.00000})
        with self.assertRaises(TypeError):
            VerticesList([v1, 2])
        with self.assertRaises(ValueError):
            VerticesList([v1, v2])

    def test_vertex_on_boundary(self):
        with self.assertRaises(ValueError):
            # boundary_range and chempot differ
            VertexOnBoundary("Label", {"Mg": -1.0}, -2.0, -2.0)
        v = VertexOnBoundary("Label", {"Mg": -2.0}, -2.0, -2.0)
        self.assertAlmostEqual(v.boundary_range_limit, -2.0, 5)
        v.set_boundary_range(-2.1)  # OK
        self.assertAlmostEqual(v.coords[Element("Mg")], -2.1, 5)
        with self.assertRaises(ValueError):
            v.set_boundary_range(-1.9)  # shallower than limit

    def test_3d_loop(self):
        """
        before:
        
                 points[2]
        points[1]          points[0]
                 points[3]
        ---------------------------------
        after:
                  points[2]
        points[3]           points[1]
                  points[0]
        """

        points = np.array([[1, 0, 0],
                           [-1, 0, 0],
                           [0, 1, 0],
                           [0, -1, 0]
                           ])
        elements = ["H", "He", "Li"]
        coords = points
        vl = VerticesList([Vertex(None, {e: c for e, c in zip(elements, vect)})
                           for vect in coords])
        sorted_expected = np.array([[0, -1, 0],
                                    [1, 0, 0],
                                    [0, 1, 0],
                                    [-1, 0, 0]
                                    ])
        after_sorted = vl.sorted_to_loop_in_3d(elements)
        for vect, expected in \
                zip(after_sorted, sorted_expected):
            vect: Vertex
            for actual, exp in zip(vect.coords_vector(elements), expected):
                self.assertAlmostEqual(actual, exp, 5)
        after_sorted.set_alphabetical_label()
        for v, expected in zip(after_sorted, ["A", "B", "C", "D"]):
            self.assertEqual(v.label, expected)


if __name__ == "__main__":
    unittest.main()
