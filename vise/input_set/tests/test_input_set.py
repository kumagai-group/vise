# -*- coding: utf-8 -*-

from vise.input_set.new_input_set import InputSet
from vise.util.testing import ViseTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InputSetTest(ViseTest):
    def test(self):
        mgo = self.get_structure_by_name("MgO")
        input_set = InputSet.make_input(structure=mgo)
        input_set.write_input(output_dir=".")