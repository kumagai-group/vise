# -*- coding: utf-8 -*-

from vise.util.testing import ViseTest
from vise.input_set.incar import incar_flags


class IncarFlagsTest(ViseTest):

    def test(self) -> None:
        expected = {"NBLOCK": None, "SMASS": None}
        actual = incar_flags["md"]
        self.assertEqual(expected, actual)


