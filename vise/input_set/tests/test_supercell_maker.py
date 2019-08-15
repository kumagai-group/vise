# -*- coding: utf-8 -*-
import os
import unittest

from pymatgen.core.structure import Structure

from vise.input_maker.supercell_maker import Supercell

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "vasp")


class SupercellTest(unittest.TestCase):
    def setUp(self):
        self.structure = Structure.from_file("PPOSCAR-MgO")
#       from_file(os.path.join(test_dir, "PPOSCAR-MgO"))

    def test_init1(self):
        multi = [2, 1, 1]
        s1 = Supercell(structure=self.structure, multi=multi, comment='')
        s1.supercell_structure.to(filename="PPOSCAR-MgO-2x1x1")

    def test_init2(self):
        multi = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
        s2 = Supercell(structure=self.structure, multi=multi, comment='')
        s2.supercell_structure.to(filename="PPOSCAR-MgO-conv")

