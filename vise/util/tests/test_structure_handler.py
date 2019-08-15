# -*- coding: utf-8 -*-

import numpy as np
import os
import unittest

from pymatgen.core.structure import Structure
from pymatgen.core.sites import PeriodicSite

from obadb.util.structure_handler \
    import get_symmetry_dataset, get_point_group_from_dataset, \
    structure_to_spglib_cell, spglib_cell_to_structure, find_hpkot_primitive, \
    structure_to_seekpath, get_coordination_environment, \
    get_coordination_distances, NotAppropriatePrimitiveError

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "vasp")


class GetSymmetryDatasetTest(unittest.TestCase):

    def setUp(self):
        self.structure = Structure.from_file("POSCAR-MgSe300atoms")

    def test(self):
        dataset = get_symmetry_dataset(self.structure, symprec=0.1)
        print(dataset)


class GetPointGroupFromDatasetTest(unittest.TestCase):

    def setUp(self):
        structure = Structure.from_file("POSCAR-MgSe300atoms")
        self.sym_dataset = get_symmetry_dataset(structure, symprec=0.1)
        self.lattice = structure.lattice.matrix

    def test(self):
        coords = [0.133333, 0.066667, 0.166252]

        print(get_point_group_from_dataset(self.sym_dataset, coords,
                                           self.lattice, symprec=0.1))


class StructureToSpglibCellTest(unittest.TestCase):

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


class SpglibCellToStructureTest(unittest.TestCase):
    def setUp(self):
        self.structure = Structure.from_file("BPOSCAR-MgO")
        self.cell = structure_to_spglib_cell(self.structure)
        self.returned_structure = spglib_cell_to_structure(self.cell)

    def test(self):
        self.assertTrue(self.structure == self.returned_structure)


class FindSpglibStandardConventionalTest(unittest.TestCase):
    def setUp(self):
        self.structure = Structure.from_file("PPOSCAR-MgO")
        self.cell = structure_to_spglib_cell(self.structure)


class FindHPKOTPrimitiveTest(unittest.TestCase):
    def test(self):
        expected = Structure.from_file("PPOSCAR-MgO")
        actual = find_hpkot_primitive(Structure.from_file("BPOSCAR-MgO"))

        from pymatgen.analysis.structure_matcher import StructureMatcher
        sm = StructureMatcher(ltol=0.0001, stol=0.0001, angle_tol=0.001,
                              primitive_cell=False, scale=False)
        self.assertTrue(sm.fit(expected, actual))


class StructureToSeekpathTest(unittest.TestCase):

    def setUp(self):
        self.conventional_cell_structure = \
            Structure.from_file("PPOSCAR-MgO-strange_lattice")
        self.primitive_cell_structure = Structure.from_file("PPOSCAR-MgO")

    def test_raise_error(self):
        with self.assertRaises(NotAppropriatePrimitiveError) as cm:
            structure_to_seekpath(self.conventional_cell_structure)

    def test_return_res(self):
        res = structure_to_seekpath(self.primitive_cell_structure)
        expected = np.array([0.5, 0, 0.5])
        np.testing.assert_equal(expected, res["explicit_kpoints_rel"][-1])
        self.assertEqual("X", res["explicit_kpoints_labels"][-1])


class StructureToSeekpathTest(unittest.TestCase):

    def setUp(self):
        self.conventional_cell_structure = \
            Structure.from_file("PPOSCAR-MgO-strange_lattice")
        self.primitive_cell_structure = Structure.from_file("PPOSCAR-MgO")

    def test_raise_error(self):
        with self.assertRaises(NotAppropriatePrimitiveError) as cm:
            structure_to_seekpath(self.conventional_cell_structure)

    def test_return_res(self):
        res = structure_to_seekpath(self.primitive_cell_structure)
        expected = np.array([0.5, 0, 0.5])
        np.testing.assert_equal(expected, res["explicit_kpoints_rel"][-1])
        self.assertEqual("X", res["explicit_kpoints_labels"][-1])


class GetCoordinationDistanceTest(unittest.TestCase):
    def setUp(self):
        self.structure_s2n2o7 = Structure.from_file("POSCAR-Sn2Nb2O7")
        self.structure_zns = Structure.from_file("POSCAR-ZnS")

    def test_get_coordination_environment(self):
        s2n2o7_6th_expected = \
            (PeriodicSite("O", [5.2888, 8.5916, 5.2888],
                          self.structure_s2n2o7.lattice,
                          coords_are_cartesian=True),
             2.723862367496558)
        env1 = get_coordination_environment(self.structure_s2n2o7, 0)

        self.assertEqual(env1[5], s2n2o7_6th_expected)

        zns_16th_expected = \
            (PeriodicSite("S", [-1.3619, -1.3619, 1.3619],
                          self.structure_zns.lattice,
                          coords_are_cartesian=True),
             3.8519712615815305)

        env2 = get_coordination_environment(self.structure_zns, 32)

        self.assertEqual(env2[15], zns_16th_expected)

    def test_get_neighbors(self):
        a = get_coordination_distances(self.structure_zns, 0)
        for k, v in a.items():
            print(k + ": " + " ".join([str(round(i, 2)) for i in v]))


