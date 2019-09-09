# -*- coding: utf-8 -*-

import numpy as np
import os

from pymatgen.core.structure import Structure
from pymatgen.core.sites import PeriodicSite

from vise.util.structure_handler import (
    get_symmetry_dataset, structure_to_spglib_cell, spglib_cell_to_structure,
    find_hpkot_primitive, structure_to_seekpath )
from vise.util.testing import ViseTest
from vise.core.error_classes import InvalidStructureError

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class GetSymmetryDatasetTest(ViseTest):

    def setUp(self):
        self.structure = Structure.from_file("POSCAR-MgSe300atoms")

    def test(self):
        dataset = get_symmetry_dataset(self.structure, symprec=0.1)
        print(dataset)


class StructureToSpglibCellTest(ViseTest):

    def setUp(self):
        structure = Structure.from_file("BPOSCAR-MgO")
        self.cell = structure_to_spglib_cell(structure)

    def test(self):
        expected_lattice = \
            [[4.2468930196252002, 0.0000000000000000, 0.0000000000000000],
             [0.0000000000000000, 4.2468930196252002, 0.0000000000000000],
             [0.0000000000000000, 0.0000000000000000, 4.2468930196252002]]
        expected_1st_position = [0, 0, 0]
        expected_atomic_numbers = [12, 12, 12, 12, 8, 8, 8, 8]

        self.assertTrue(expected_lattice, self.cell[0])
        self.assertTrue(expected_1st_position, self.cell[1][0])
        self.assertTrue(expected_atomic_numbers, self.cell[2])


class SpglibCellToStructureTest(ViseTest):
    def setUp(self):
        self.structure = Structure.from_file("BPOSCAR-MgO")
        self.cell = structure_to_spglib_cell(self.structure)
        self.returned_structure = spglib_cell_to_structure(self.cell)

    def test(self):
        self.assertTrue(self.structure == self.returned_structure)


class FindSpglibStandardConventionalTest(ViseTest):
    def setUp(self):
        self.structure = Structure.from_file("PPOSCAR-MgO")
        self.cell = structure_to_spglib_cell(self.structure)


class FindHPKOTPrimitiveTest(ViseTest):
    def test(self):
        expected = Structure.from_file("PPOSCAR-MgO")
        actual = find_hpkot_primitive(Structure.from_file("BPOSCAR-MgO"))

        from pymatgen.analysis.structure_matcher import StructureMatcher
        sm = StructureMatcher(ltol=0.0001, stol=0.0001, angle_tol=0.001,
                              primitive_cell=False, scale=False)
        self.assertTrue(sm.fit(expected, actual))


class StructureToSeekpathTest(ViseTest):

    def setUp(self):
        self.conventional_cell_structure = \
            Structure.from_file("PPOSCAR-MgO-strange_lattice")
        self.primitive_cell_structure = Structure.from_file("PPOSCAR-MgO")

    def test_raise_error(self):
        with self.assertRaises(InvalidStructureError):
            structure_to_seekpath(self.conventional_cell_structure)

    def test_return_res(self):
        res = structure_to_seekpath(self.primitive_cell_structure)
        expected = np.array([0.5, 0, 0.5])
        np.testing.assert_equal(expected, res["explicit_kpoints_rel"][-1])
        self.assertEqual("X", res["explicit_kpoints_labels"][-1])


class StructureToSeekpathTest(ViseTest):

    def setUp(self):
        self.conventional_cell_structure = \
            Structure.from_file("PPOSCAR-MgO-strange_lattice")
        self.primitive_cell_structure = Structure.from_file("PPOSCAR-MgO")

    def test_raise_error(self):
        with self.assertRaises(InvalidStructureError):
            structure_to_seekpath(self.conventional_cell_structure)

    def test_return_res(self):
        res = structure_to_seekpath(self.primitive_cell_structure)
        expected = np.array([0.5, 0, 0.5])
        np.testing.assert_equal(expected, res["explicit_kpoints_rel"][-1])
        self.assertEqual("X", res["explicit_kpoints_labels"][-1])


