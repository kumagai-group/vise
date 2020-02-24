# -*- coding: utf-8 -*-

from vise.input_set.input_set import ViseInputSet
from vise.input_set.xc import Xc
from vise.util.testing import ViseTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InputSetTest(ViseTest):
    def setUp(self) -> None:
        mgo = self.get_structure_by_name("MgO")
        self.input_set = ViseInputSet.make_input(structure=mgo, xc=Xc.hse)

    def test_write(self):
        self.input_set.write_input(output_dir=".")

    def test_dict(self):
        expected = self.input_set.as_dict()
        actual = ViseInputSet.from_dict(expected).as_dict()
        self.assertEqual(expected, actual)
