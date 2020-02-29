# -*- coding: utf-8 -*-

from vise.util.testing import ViseTest
from vise.input_set.incar import incar_flags

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class IncarFlagsTest(ViseTest):

    def test(self) -> None:
        expected = ["SMASS", "NBLOCK"]
        actual = incar_flags["md"]
        self.assertEqual(expected, actual)


