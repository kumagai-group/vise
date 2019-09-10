# -*- coding: utf-8 -*-

from vise.input_set.new_input_set import InputSet
from vise.util.testing import ViseTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InputSetTest(ViseTest):
    def setUp(self) -> None:
        mgo = self.get_structure_by_name("YMnO3")
        self.input_set = InputSet.make_input(structure=mgo,
                                             xc="hse")

    def test_write(self):
        self.input_set.write_input(output_dir=".")

    def test_dict(self):
        expected = self.input_set.as_dict()
        actual = InputSet.from_dict(expected).as_dict()
        self.assertEqual(expected, actual)
