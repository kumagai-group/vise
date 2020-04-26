# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np

from pymatgen.analysis.structure_matcher import StructureMatcher

from vise.util.structure_symmetrizer import (
    StructureSymmetrizer, cell_to_structure)


from vise.util.testing import ViseTest


class CellToStructureTest(ViseTest):

    def test(self):
        cell = (np.array([[1., 0., 0.],
                          [0., 1., 0.],
                          [0., 0., 1.]], dtype='double', order='C'),
                np.array([[0., 0., 0.],
                          [0.5, 0.5, 0.5]], dtype='double', order='C'),
                np.array([1, 2], dtype='intc'))
        print(cell_to_structure(cell))


class StructureSymmetryTest(ViseTest):

    def setUp(self):
        structure = self.get_structure_by_name("conventional_MgO")
        self.ss = StructureSymmetrizer(structure)

    def test(self):
#        print(self.ss.spglib_sym_data)
#        print(self.ss.conventional)
#        print(self.ss.primitive)
        print(self.ss.find_seekpath_data())
        print(self.ss.seekpath_data)
        print(self.ss.band_primitive)


# class StructureToSpglibCellTest(ViseTest):

    # def setUp(self):
    #     structure = self.get_structure_by_name("conventional_MgO")
    #     self.cell = structure_to_spglib_cell(structure)

    # def test(self):
    #     expected_lattice = \
    #         [[4.2468930196252002, 0.0000000000000000, 0.0000000000000000],
    #          [0.0000000000000000, 4.2468930196252002, 0.0000000000000000],
    #          [0.0000000000000000, 0.0000000000000000, 4.2468930196252002]]
    #     expected_1st_position = [0, 0, 0]
    #     expected_atomic_numbers = [12, 12, 12, 12, 8, 8, 8, 8]

        # self.assertTrue(expected_lattice, self.cell[0])
        # self.assertTrue(expected_1st_position, self.cell[1][0])
        # self.assertTrue(expected_atomic_numbers, self.cell[2])


# class FindSpglibStandardConventionalTest(ViseTest):
#     def setUp(self):
#         self.structure = self.get_structure_by_name("MgO")
#         self.cell = structure_to_spglib_cell(self.structure)


# class FindHPKOTPrimitiveTest(ViseTest):
#     def test(self):
#         expected = self.get_structure_by_name("MgO")
#         actual = find_hpkot_primitive(
#             self.get_structure_by_name("conventional_MgO"))

        # sm = StructureMatcher(ltol=0.0001, stol=0.0001, angle_tol=0.001,
        #                       primitive_cell=False, scale=False)
        # self.assertTrue(sm.fit(expected, actual))


# class StructureToSeekpathTest(ViseTest):

    # def setUp(self):
    #     # self.conventional_cell_structure = \
    #     #     Structure.from_file("PPOSCAR-MgO-strange_lattice")
    #     self.primitive_cell_structure = self.get_structure_by_name("MgO")

#     # def test_raise_error(self):
#     #     with self.assertRaises(InvalidStructureError):
#     #         structure_to_seekpath(self.conventional_cell_structure)

    # def test_return_res(self):
    #     res = structure_to_seekpath(self.primitive_cell_structure)
    #     expected = np.array([0.5, 0, 0.5])
    #     np.testing.assert_equal(expected, res["explicit_kpoints_rel"][-1])
    #     self.assertEqual("X", res["explicit_kpoints_labels"][-1])



